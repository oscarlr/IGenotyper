"""Small public synthetic fixtures for assembly, restart and WhatsHap failures."""
import os
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pysam
import pytest
import pyBigWig

from IGenotyper.assembly.scripts import assembly_contigs, create_assemble_script, region_assembled
from IGenotyper.commands.assembly import combine_sequence, combine_assembly_sequences
from IGenotyper.command_lines.assembly import Assembly
from IGenotyper.command_lines.alignments import Align
from IGenotyper.command_lines.clt import CommandLine
from IGenotyper.command_lines.reads import ReadManip
from IGenotyper.command_lines.snps import Snps
from IGenotyper.common.validation import require_usable_reads, validate_vcf

CPU = SimpleNamespace(threads=1, cluster=False)


def write_bam(path, platform='REVIO', reads=(), chroms=('chr1', 'igh')):
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
              'SQ': [{'SN': c, 'LN': 400} for c in chroms],
              'RG': [{'ID': 'rg', 'SM': 'sample', 'PM': platform}]}
    with pysam.AlignmentFile(str(path), 'wb', header=header) as stream:
        for tid, name, alt, flag in reads:
            read = pysam.AlignedSegment()
            read.query_name = name
            bases = list('ACGT' * 50)
            if alt:
                bases[50], bases[100] = 'A', 'C'
            read.query_sequence = ''.join(bases)
            read.flag = flag
            read.reference_id = tid
            read.reference_start = 0
            read.mapping_quality = 60
            read.cigarstring = '200M'
            read.query_qualities = pysam.qualitystring_to_array('I' * 200)
            read.set_tag('RG', 'rg')
            stream.write(read)
    pysam.index(str(path))
    return str(path)


def files_for(tmp_path, platform='REVIO'):
    bam = write_bam(tmp_path / 'input.bam', platform)
    return SimpleNamespace(tmp=str(tmp_path), input_bam=bam, ccs_bam=bam,
        ccs_fastq=str(tmp_path / 'reads.fasta'), ccs_to_ref_phased=bam,
        subreads_to_ref_phased=bam, scripts=str(Path('IGenotyper/scripts').resolve()),
        assembly_script=str(Path('IGenotyper/data/assembly.sh').resolve()),
        assembly_fasta=str(tmp_path / 'assembly.fasta'),
        igh_assembly_fasta=str(tmp_path / 'igh.fasta'))


@pytest.mark.parametrize('platform,relative,mode', [
    ('SEQUEL', 'contigs.fasta', '-pacbio'),
    ('SEQUELII', 'canu/canu.contigs.fasta', '-pacbio-hifi'),
    ('REVIO', 'canu/canu.contigs.fasta', '-pacbio-hifi')])
def test_platform_execution_and_collection(tmp_path, platform, relative, mode):
    files = files_for(tmp_path, platform)
    region = tmp_path / 'assembly/chr1/0_200/1'
    region.mkdir(parents=True)
    (region / 'reads.fasta').write_text('>read\nACGT\n')
    # Stub external tools only; execute the actual generated assembly script.
    bindir = tmp_path / 'bin'; bindir.mkdir()
    canu = bindir / 'canu'
    canu.write_text('#!/bin/bash\nset -e\n' +
        'printf "%s\\n" "$@" > ' + str(tmp_path / 'canu.args') + '\n' +
        'mkdir -p ' + str(region / 'canu') + '\n' +
        "printf '>canu\\nACGT\\n' > " + str(region / 'canu/canu.contigs.fasta') + '\n')
    canu.chmod(0o755)
    # Existing valid polishing output is preserved for SEQUEL.
    if platform == 'SEQUEL':
        (region / relative).write_text('>polished\nTGCA\n')
    script = create_assemble_script(files, CPU, str(region), 'chr1', 0, 200, '1')
    env = dict(os.environ, PATH=str(bindir) + os.pathsep + str(Path(sys.executable).parent) + os.pathsep + os.environ['PATH'])
    subprocess.run(['bash', script], check=True, env=env)
    assert mode in (tmp_path / 'canu.args').read_text().splitlines()
    assert assembly_contigs(str(region), platform) == str(region / relative)
    assert region_assembled(str(region), platform)
    combine_assembly_sequences(files, [('chr1', 0, 200, '1')])
    assert ('TGCA' if platform == 'SEQUEL' else 'ACGT') in Path(files.assembly_fasta).read_text()
    assert not Path(files.igh_assembly_fasta).exists()


@pytest.mark.parametrize('content', [None, '', '>empty\n'])
def test_missing_assembly_preserves_previous_output(tmp_path, content):
    files = files_for(tmp_path)
    region = tmp_path / 'assembly/igh/0_200/1'; (region / 'canu').mkdir(parents=True)
    (region / 'done').touch()
    if content is not None:
        (region / 'canu/canu.contigs.fasta').write_text(content)
    Path(files.assembly_fasta).write_text('>previous\nACTG\n')
    assert not region_assembled(str(region), 'REVIO')
    with pytest.raises(RuntimeError, match='assembly FASTA'):
        combine_sequence(files, [('igh', 0, 200, '1')], files.assembly_fasta, 'fasta')
    assert Path(files.assembly_fasta).read_text() == '>previous\nACTG\n'
    with pytest.raises(RuntimeError, match='No assembly regions'):
        combine_sequence(files, [], files.assembly_fasta, 'fasta')


@pytest.mark.parametrize("exit_code", [0, 7])
def test_failed_assembly_never_writes_done(tmp_path, exit_code):
    files = files_for(tmp_path)
    region = tmp_path / 'region'; region.mkdir()
    (region / 'reads.fasta').write_text('>r\nACGT\n')
    (region / 'done').touch()
    bindir = tmp_path / 'bin'; bindir.mkdir()
    canu = bindir / 'canu'; canu.write_text('#!/bin/bash\nexit %s\n' % exit_code); canu.chmod(0o755)
    script = create_assemble_script(files, CPU, str(region), 'igh', 0, 200, '1')
    with patch.dict(os.environ, PATH=str(bindir) + os.pathsep + os.environ['PATH']):
        with pytest.raises(subprocess.CalledProcessError):
            Assembly(files, CPU, 'sample').run_assembly_scripts([script])
    assert not (region / 'done').exists()


def test_failed_and_legacy_outputs_are_not_reused(tmp_path):
    runner = CommandLine(None, None, None)
    output = tmp_path / 'partial.vcf'
    command = "printf 'partial' > '%s'; exit 1" % output
    for _ in range(2):
        with pytest.raises(subprocess.CalledProcessError):
            runner.run_command(command, output)
        assert not Path(str(output) + '.success.json').exists()
    assert not output.exists()
    runner.run_command("printf 'complete' > '%s'" % output, output)
    assert output.read_text() == 'complete'
    with patch('subprocess.check_call', side_effect=AssertionError('should reuse')):
        runner.run_command("printf 'complete' > '%s'" % output, output)


def test_multiple_outputs_and_pipeline_failure(tmp_path):
    runner = CommandLine(None, None, None)
    bam, index = tmp_path / 'x.data', tmp_path / 'x.data.index'
    command = "printf bam > '%s'; printf index > '%s'" % (bam, index)
    runner.run_command(command, [bam, index])
    bam.unlink()
    runner.run_command(command, [bam, index])
    assert bam.read_text() == 'bam'
    with pytest.raises(subprocess.CalledProcessError):
        runner.run_command("false | printf broken > '%s'" % bam, [bam, index])
    assert bam.read_text() == 'bam' and index.read_text() == 'index'
    with pytest.raises(RuntimeError, match='Missing or empty'):
        runner.run_command("printf bam > '%s'" % bam, [bam, index])


@pytest.mark.parametrize('flag', [None, 256, 2048])
def test_zero_usable_bam_fails_before_conversion(tmp_path, flag):
    bam = write_bam(tmp_path / 'reads.bam', reads=[] if flag is None else [(0, 'r', False, flag)])
    files = SimpleNamespace(ccs_bam=bam, ccs_fastq=str(tmp_path / 'reads.fa'))
    with patch('subprocess.check_call', side_effect=AssertionError('must fail before conversion')):
        with pytest.raises(ValueError, match='zero usable'):
            ReadManip(files, CPU, 'sample').turn_ccs_reads_to_fastq()


def test_empty_track_is_real_zero_bigwig(tmp_path):
    bam = write_bam(tmp_path / 'empty.bam')
    track = tmp_path / 'empty.bw'
    with patch('subprocess.check_call', side_effect=AssertionError('must not run bamCoverage')):
        Align(None, CPU, 'sample').bam_to_bigwig(bam, str(track))
    with pyBigWig.open(str(track)) as bw:
        assert bw.chroms() == {'chr1': 400, 'igh': 400}
        assert bw.values('igh', 0, 400) == [0.0] * 400


def test_nonempty_coverage_failure_propagates(tmp_path):
    bam = write_bam(tmp_path / 'reads.bam', reads=[(0, 'r', False, 0)])
    with patch('subprocess.check_call', side_effect=subprocess.CalledProcessError(1, 'bamCoverage')):
        with pytest.raises(subprocess.CalledProcessError):
            Align(None, CPU, 'sample').bam_to_bigwig(bam, str(tmp_path / 'track.bw'))


def variant_fixture(tmp_path, collision_chrom=None):
    ref = tmp_path / 'ref.fa'
    ref.write_text('>chr1\n' + 'ACGT' * 100 + '\n>igh\n' + 'ACGT' * 100 + '\n')
    pysam.faidx(str(ref))
    reads = []
    for tid, chrom in enumerate(['chr1', 'igh']):
        for i in range(8):
            name = 'collision' if collision_chrom == chrom and i < 2 else 'read%s' % i
            reads.append((tid, name, i % 2, 0))
    bam = write_bam(tmp_path / 'reads.bam', reads=reads)
    return SimpleNamespace(ref=str(ref)), bam


@pytest.mark.parametrize('collision_chrom', ['chr1', 'igh'])
def test_whats_hap_28_exact_crash_and_safe_restart(tmp_path, collision_chrom):
    whatshap = pytest.importorskip('whatshap')
    if whatshap.__version__ != '2.8':
        pytest.skip('Upstream failure reproduction is specific to supported WhatsHap 2.8')
    files, bam = variant_fixture(tmp_path, collision_chrom)
    candidates, output = str(tmp_path / 'candidates.vcf'), str(tmp_path / 'genotypes.vcf')
    caller = Snps(files, CPU, 'sample')
    caller.snp_candidates(bam, candidates)
    assert len(validate_vcf(candidates, 'sample')) == 4
    # Capture the real upstream traceback, then exercise our runner twice.
    failed = subprocess.run(['whatshap', 'genotype', '--sample', 'sample', '--ignore-read-groups',
        '--reference', files.ref, '-o', output, candidates, bam], capture_output=True, text=True)
    assert failed.returncode != 0
    assert 'compute_genotypes' in failed.stderr and 'RuntimeError: No variants present' in failed.stderr
    assert Path(output).stat().st_size > 0
    # The repository adapter repairs only zero-observation internal reads.
    caller.snp_genotypes(bam, candidates, output)
    assert len(validate_vcf(output, 'sample')) == 4
    assert Path(output + '.success.json').exists()


def test_whats_hap_unique_names_succeeds(tmp_path):
    pytest.importorskip('whatshap')
    files, bam = variant_fixture(tmp_path)
    caller = Snps(files, CPU, 'sample')
    candidates, output = str(tmp_path / 'candidates.vcf'), str(tmp_path / 'genotypes.vcf')
    caller.snp_candidates(bam, candidates)
    caller.snp_genotypes(bam, candidates, output)
    assert len(validate_vcf(output, 'sample')) == 4
    with patch('subprocess.check_call', side_effect=AssertionError('must reuse valid VCF')):
        caller.snp_genotypes(bam, candidates, output)
    phased = str(tmp_path / 'phased.vcf')
    caller.phase_snps(phased, output, [bam])
    assert len(validate_vcf(phased, 'sample')) == 4


def test_conversion_preserves_distinct_strand_names(tmp_path):
    names = ['movie/1/ccs/fwd', 'movie/1/ccs/rev']
    bam = write_bam(tmp_path / 'reads.bam', reads=[(0, name, i, 0) for i, name in enumerate(names)])
    files = SimpleNamespace(ccs_bam=bam, ccs_fastq=str(tmp_path / 'reads.fa'))
    # Provide samtools via pysam's bundled samtools for a real BAM-to-FASTA conversion.
    def execute(command, **kwargs):
        assert 'sed' not in command
        import shlex
        destination = shlex.split(command)[-1]
        pysam.fasta("-0", destination, bam, catch_stdout=False)
    with patch('subprocess.check_call', side_effect=execute):
        ReadManip(files, CPU, 'sample').turn_ccs_reads_to_fastq()
    headers = [line[1:] for line in Path(files.ccs_fastq).read_text().splitlines() if line.startswith('>')]
    assert headers == names
    old_headers = [name.replace('/ccs', '/0_8').replace('/fwd', '').replace('/rev', '') for name in names]
    assert len(set(old_headers)) == 1


def test_vcf_validation_rejects_truncation_and_order(tmp_path):
    header = '##fileformat=VCFv4.2\n##contig=<ID=igh,length=400>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n'
    row = 'igh\t%s\t.\tA\tC\t.\t.\t.\tGT\t0/1\n'
    vcf = tmp_path / 'test.vcf'
    vcf.write_text(header + row % 10 + row % 20)
    sites = validate_vcf(str(vcf), 'sample')
    vcf.write_text(header + row % 10)
    with pytest.raises(ValueError, match='retain all'):
        validate_vcf(str(vcf), 'sample', expected=sites)
    for positions in [(20, 10), (10, 10)]:
        vcf.write_text(header + ''.join(row % p for p in positions))
        with pytest.raises(ValueError, match='Unsorted or duplicate'):
            validate_vcf(str(vcf), 'sample')


def test_guard_removes_only_empty_grouped_read_and_retains_positions():
    from whatshap.core import Read, ReadSet, compute_genotypes
    from whatshap.variants import ReadSetReader, AlignedRead
    from IGenotyper.command_lines.whatshap_genotype import compute_genotypes_with_evidence

    left, right = Read('same', 60), Read('same', 60)
    left.add_variant(50, 0, 30)
    right.add_variant(50, 1, 30)
    grouped = ReadSetReader.create_read_from_group([
        AlignedRead(left, False, False, 0, 200),
        AlignedRead(right, False, False, 0, 200)], 100000)
    assert len(grouped) == 0
    empty = ReadSet(); empty.add(grouped)
    with pytest.raises(RuntimeError, match='No variants present'):
        compute_genotypes(empty, [50])
    with pytest.raises(ValueError, match='No informative variant observations'):
        compute_genotypes_with_evidence(empty, [50])
    with pytest.raises(ValueError, match='No informative variant observations'):
        compute_genotypes_with_evidence(ReadSet(), [50])

    informative = Read('informative', 60); informative.add_variant(50, 1, 30)
    evidence = ReadSet(); evidence.add(informative)
    mixed = ReadSet(); mixed.add(grouped); mixed.add(informative)
    expected_gt, expected_gl = compute_genotypes(evidence, [50, 100])
    actual_gt, actual_gl = compute_genotypes_with_evidence(mixed, [50, 100])
    assert actual_gt == expected_gt
    assert [list(gl) for gl in actual_gl] == [list(gl) for gl in expected_gl]
    assert len(actual_gt) == 2  # Uncovered candidate retained with uniform likelihood.
    assert list(actual_gl[1]) == pytest.approx([1/3, 1/3, 1/3])


def test_adapter_does_not_hide_other_errors_or_accept_unknown_version():
    from whatshap.core import Read, ReadSet
    from IGenotyper.command_lines import whatshap_genotype as adapter
    read = Read('r', 60); read.add_variant(5, 0, 30)
    reads = ReadSet(); reads.add(read)
    with patch('whatshap.core.compute_genotypes', side_effect=RuntimeError('unrelated')):
        with pytest.raises(RuntimeError, match='unrelated'):
            adapter.compute_genotypes_with_evidence(reads, [5])
    with patch.object(adapter, 'version', return_value='9.9'):
        with pytest.raises(RuntimeError, match='supports WhatsHap 2.8'):
            adapter.main()


def test_all_conflicting_evidence_fails_without_publishing(tmp_path):
    files, _ = variant_fixture(tmp_path)
    reads = [(0, 'pair%s' % i, alt, 0) for i in range(3) for alt in [False, True]]
    bam = write_bam(tmp_path / 'conflicting.bam', reads=reads)
    caller = Snps(files, CPU, 'sample')
    candidates, output = str(tmp_path / 'candidates.vcf'), str(tmp_path / 'genotypes.vcf')
    caller.snp_candidates(bam, candidates)
    assert len(validate_vcf(candidates, 'sample')) == 2
    for _ in range(2):
        with pytest.raises(subprocess.CalledProcessError):
            caller.snp_genotypes(bam, candidates, output)
        assert not Path(output).exists()
        assert not Path(output + '.success.json').exists()


def test_parseable_partial_vcf_cannot_be_published(tmp_path):
    files, bam = variant_fixture(tmp_path)
    caller = Snps(files, CPU, 'sample')
    candidates = str(tmp_path / 'candidates.vcf')
    output = tmp_path / 'out.vcf'
    caller.snp_candidates(bam, candidates)
    sites = validate_vcf(candidates, 'sample')
    def partial(paths):
        text = Path(candidates).read_text()
        Path(paths[0]).write_text('\n'.join(text.splitlines()[:-1]) + '\n')
    with pytest.raises(ValueError, match='retain all'):
        caller.run_command('truncated-success', output,
            action=partial, validator=lambda paths: validate_vcf(paths[0], 'sample', expected=sites))
    assert not output.exists()
    assert not Path(str(output) + '.success.json').exists()


def test_stale_receipt_and_changed_input_force_rerun(tmp_path):
    runner = CommandLine(None, None, None)
    source, output = tmp_path / 'source', tmp_path / 'output'
    source.write_text('original')
    command = "cat '%s' > '%s'" % (source, output)
    runner.run_command(command, output, inputs=[source])
    output.write_text('corrupted')
    runner.run_command(command, output, inputs=[source])
    assert output.read_text() == 'original'
    source.write_text('new naming')
    runner.run_command(command, output, inputs=[source])
    assert output.read_text() == 'new naming'


def test_staging_handles_spaces_and_preserves_prior_output_on_failure(tmp_path):
    directory = tmp_path / 'space dir'; directory.mkdir()
    output = directory / 'out file.vcf'; output.write_text('previous')
    runner = CommandLine(None, None, None)
    with pytest.raises(subprocess.CalledProcessError):
        runner.run_command("printf partial > '%s'; exit 1" % output, output)
    assert output.read_text() == 'previous'
    runner.run_command("printf validated > '%s'" % output, output)
    assert output.read_text() == 'validated'


def test_phasing_preserves_original_names_and_invalidates_old_bam(tmp_path):
    from IGenotyper.phasing.reads import phase_alignments
    files, _ = variant_fixture(tmp_path)
    names = ['movie/1/ccs/fwd', 'movie/1/ccs/rev'] + ['read%s' % i for i in range(6)]
    bam = write_bam(tmp_path / 'names.bam', reads=[(0, name, i % 2, 0) for i, name in enumerate(names)])
    caller = Snps(files, CPU, 'sample')
    candidates, genotype, phased = [str(tmp_path / name) for name in ['candidates.vcf', 'genotypes.vcf', 'phased.vcf']]
    caller.snp_candidates(bam, candidates)
    caller.snp_genotypes(bam, candidates, genotype)
    caller.phase_snps(phased, genotype, [bam])
    output = write_bam(tmp_path / 'phased.bam', reads=[(0, 'old/collapsed/0_8', False, 0)])
    phase_alignments(phased, bam, 'sample', output)
    with pysam.AlignmentFile(output, 'rb') as reads:
        assert [r.query_name for r in reads] == names
    with patch('IGenotyper.phasing.reads.read_in_phased_vcf', side_effect=AssertionError('must reuse')):
        phase_alignments(phased, bam, 'sample', output)


def test_stale_region_is_archived_after_phasing_changes(tmp_path):
    from IGenotyper.assembly.scripts import get_assembly_scripts, assembly_provenance, record_assembly_success
    import json
    files = files_for(tmp_path)
    directory = tmp_path / 'assembly/igh/0_200/1'
    (directory / 'canu').mkdir(parents=True)
    (directory / 'canu/canu.contigs.fasta').write_text('>old\nACGT\n')
    record_assembly_success(str(directory), "REVIO", assembly_provenance(files, "REVIO"))
    assert get_assembly_scripts(files, CPU, [('igh', 0, 200, '1')]) == []
    write_bam(Path(files.ccs_to_ref_phased), reads=[(0, 'new/fwd', False, 0)])
    assert len(get_assembly_scripts(files, CPU, [('igh', 0, 200, '1')])) == 1
    assert not (directory / 'canu/canu.contigs.fasta').exists()
    assert len(list(directory.parent.glob('1.previous-*/canu/canu.contigs.fasta'))) == 1


def test_unmapped_sequence_is_usable_but_missing_sequence_is_not(tmp_path):
    bam = tmp_path / 'unmapped.bam'
    for sequence in [None, 'ACGT']:
        with pysam.AlignmentFile(str(bam), 'wb', header={'HD': {'VN': '1.6'}}) as stream:
            read = pysam.AlignedSegment()
            read.query_name = 'raw/ccs/fwd'
            read.flag = 4
            read.query_sequence = sequence
            stream.write(read)
        if sequence is None:
            with pytest.raises(ValueError, match='zero usable'):
                require_usable_reads(str(bam))
        else:
            require_usable_reads(str(bam))


def test_real_bam_and_index_are_republished_as_one_completed_step(tmp_path):
    source = write_bam(tmp_path / 'source.bam', reads=[(0, 'movie/1/ccs/fwd', False, 0)])
    output = str(tmp_path / 'output.bam')
    runner = CommandLine(None, None, None)
    calls = []
    def produce(paths):
        calls.append(True)
        shutil.copyfile(source, paths[0])
        pysam.index(paths[0])
    runner.run_command('test-bam', [output, output + '.bai'], inputs=[source], action=produce)
    runner.run_command('test-bam', [output, output + '.bai'], inputs=[source], action=produce)
    assert len(calls) == 1
    Path(output).unlink()
    runner.run_command('test-bam', [output, output + '.bai'], inputs=[source], action=produce)
    assert len(calls) == 2
    with pysam.AlignmentFile(output, 'rb') as reads:
        assert reads.check_index()
        assert next(reads).query_name == 'movie/1/ccs/fwd'
