"""Real BAM extraction, generated shell scripts and collection; stub assemblers."""
import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

from Bio import SeqIO
import pysam
import pytest

from IGenotyper.assembly.scripts import (
    assembly_contigs, assembly_provenance, get_assembly_scripts,
    region_assembled, region_status, region_provenance,
)
from IGenotyper.assembly.workflow import select_assembly_workflow
from IGenotyper.command_lines.assembly import Assembly
from IGenotyper.commands.assembly import combine_assembly_sequences, combine_sequence

CPU = SimpleNamespace(threads=1, cluster=False)
COVERED = ('chr1', 0, 500, '0')
UNCOVERED = ('chr1', 5000, 5500, '0')


def write_reads(path, platform, read_type, qualities, phased=False, ds=None):
    rg = {'ID': 'rg', 'PM': platform, 'SM': 'sample'}
    if ds is not None:
        rg['DS'] = ds
    elif read_type is not None:
        rg['DS'] = 'READTYPE=%s;BINDINGKIT=synthetic' % read_type
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
              'SQ': [{'SN': chrom, 'LN': 10000} for chrom in ['chr1', 'igh']],
              'RG': [dict(rg, ID=hap) for hap in ['0', '1', '2']] if phased else [rg]}
    with pysam.AlignmentFile(str(path), 'wb', header=header) as out:
        for i, quality in enumerate(qualities):
            read = pysam.AlignedSegment()
            read.query_name = 'movie/%s/%s' % (i, '0_200' if read_type == 'SUBREAD' else 'ccs')
            read.query_sequence = 'ACGT' * 50
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 100
            read.mapping_quality = 60
            read.cigarstring = '200M'
            read.query_qualities = pysam.qualitystring_to_array('I' * 200)
            read.set_tag('RG', '0' if phased else 'rg')
            if quality is not None:
                read.set_tag('rq', quality)
            out.write(read)
    pysam.index(str(path))
    return str(path)


def make_files(tmp_path, platform='SEQUEL', read_type='CCS', qualities=(None,), ds=None):
    original = write_reads(tmp_path / 'input.bam', platform, read_type, qualities, ds=ds)
    mapped = write_reads(tmp_path / 'mapped.bam', platform, read_type, qualities, phased=True, ds=ds)
    return SimpleNamespace(tmp=str(tmp_path), input_bam=original,
        ccs_to_ref_phased=mapped if read_type != 'SUBREAD' else str(tmp_path / 'absent-ccs.bam'),
        subreads_to_ref_phased=mapped if read_type == 'SUBREAD' else str(tmp_path / 'absent-subreads.bam'),
        scripts=str(Path('IGenotyper/scripts').resolve()),
        assembly_script=str(Path('IGenotyper/data/assembly.sh').resolve()),
        assembly_fasta=str(tmp_path / 'assembly.fasta'),
        igh_assembly_fasta=str(tmp_path / 'igh.fasta'))


@pytest.fixture
def stub_tools(tmp_path, monkeypatch):
    bindir = tmp_path / 'bin'; bindir.mkdir()
    log = tmp_path / 'tools.jsonl'
    program = '''
import json, os, shutil, sys
from pathlib import Path
import pysam
from Bio import SeqIO
name, args = Path(sys.argv[0]).name, sys.argv[1:]
with open(os.environ['IG_TEST_LOG'], 'a') as stream:
    stream.write(json.dumps([name, args]) + '\\n')
if name == 'canu':
    output = Path(args[args.index('-d') + 1]); output.mkdir(parents=True)
    assert list(SeqIO.parse(args[-1], 'fasta'))
    mode = os.environ.get('IG_TEST_CANU', 'success')
    if mode in ('missing', 'missing_failure'):
        pass
    elif mode == 'zero':
        (output / 'canu.contigs.fasta').touch()
    elif mode == 'empty':
        (output / 'canu.contigs.fasta').write_text('>empty\\n')
    else:
        (output / 'canu.contigs.fasta').write_text('>synthetic\\nACGT\\n')
    if mode in ('failure', 'missing_failure'):
        sys.exit(19)
elif name == 'pbindex':
    Path(args[0] + '.pbi').write_text('stub index')
elif name == 'pbmm2':
    shutil.copyfile(args[-2], args[-1])
elif name == 'gcpp':
    Path(args[args.index('-o') + 1]).write_text('>polished\\nTGCA\\n')
    if os.environ.get('IG_TEST_POLISH_FAILURE'):
        sys.exit(23)
elif name == 'samtools':
    assert args[0] == 'faidx'
    pysam.faidx(args[1])
'''
    for tool in ['canu', 'pbindex', 'pbmm2', 'gcpp', 'samtools']:
        executable = bindir / tool
        executable.write_text('#!' + sys.executable + '\n' + program)
        executable.chmod(0o755)
    monkeypatch.setenv('PATH', str(bindir) + os.pathsep + os.environ['PATH'])
    monkeypatch.setenv('IG_TEST_LOG', str(log))
    return lambda: [json.loads(line) for line in log.read_text().splitlines()] if log.exists() else []


def directory(files, block):
    chrom, start, end, hap = block
    return Path(files.tmp) / 'assembly' / chrom / ('%s_%s' % (start, end)) / hap


def execute(files, blocks):
    scripts = get_assembly_scripts(files, CPU, blocks)
    Assembly(files, CPU, 'sample').run_assembly_scripts(scripts)
    return scripts


@pytest.mark.parametrize('platform,read_type,qualities,expected,flag', [
    ('SEQUEL', 'CCS', (None,), 'ccs', '-pacbio'),
    ('SEQUEL', 'CCS', (0.95,), 'ccs', '-pacbio'),
    ('SEQUEL', 'CCS', (0.999,), 'hifi', '-pacbio-hifi'),
    ('SEQUELII', 'CCS', (0.999,), 'hifi', '-pacbio-hifi'),
    ('REVIO', 'CCS', (0.999,), 'hifi', '-pacbio-hifi'),
    ('REVIO', 'CCS', (0.999, 0.95), 'ccs', '-pacbio'),
    ('SEQUELII', 'CCS', (0.999, None), 'ccs', '-pacbio'),
    ('SEQUEL', 'SUBREAD', (0.85,), 'subreads', '-pacbio'),
    ('SEQUELII', 'SUBREAD', (0.85,), 'subreads', '-pacbio'),
])
def test_dispatch_execution_collection(tmp_path, stub_tools, platform, read_type, qualities, expected, flag):
    files = make_files(tmp_path, platform, read_type, qualities)
    workflow = select_assembly_workflow(files.input_bam)
    assert workflow.name == expected
    provenance = assembly_provenance(files, workflow)
    absent = files.ccs_to_ref_phased if read_type == 'SUBREAD' else files.subreads_to_ref_phased
    assert absent not in provenance['inputs']
    assert not Path(absent).exists()
    execute(files, [COVERED])
    calls = stub_tools()
    canu = next(args for name, args in calls if name == 'canu')
    assert flag in canu
    assert ('gcpp' in [name for name, _ in calls]) == (expected == 'subreads')
    if workflow.polish:
        align_args = next(args for name, args in calls if name == 'pbmm2')
        assert align_args[align_args.index('--preset') + 1] == 'SUBREAD'
    contigs = assembly_contigs(str(directory(files, COVERED)), workflow)
    assert contigs.endswith('/contigs.fasta' if workflow.polish else '/canu/canu.contigs.fasta')
    combine_assembly_sequences(files, [COVERED])
    assert str(next(SeqIO.parse(files.assembly_fasta, 'fasta')).seq) == ('TGCA' if workflow.polish else 'ACGT')
    assert region_assembled(str(directory(files, COVERED)), workflow, region_provenance(provenance, *COVERED))
    # Valid results retain their receipts and do not rerun the assembler.
    assert execute(files, [COVERED]) == []
    assert stub_tools() == calls


def test_mixed_covered_uncovered_and_stale_skip_outputs(tmp_path, stub_tools):
    files = make_files(tmp_path)
    empty = directory(files, UNCOVERED)
    (empty / 'canu').mkdir(parents=True)
    (empty / 'canu/canu.contigs.fasta').write_text('>old\nTTTT\n')
    (empty / 'reads.fasta').write_text('>old\nTTTT\n')
    (empty / 'done').write_text('old untrusted result')
    execute(files, [UNCOVERED, COVERED])
    workflow = select_assembly_workflow(files.input_bam)
    provenance = assembly_provenance(files, workflow)
    assert region_status(str(empty), workflow, region_provenance(provenance, *UNCOVERED)) == 'skipped_no_coverage'
    assert not region_assembled(str(empty), workflow, region_provenance(provenance, *UNCOVERED))
    assert not (empty / 'done').exists()
    assert len(list(empty.parent.glob('0.previous-*/canu/canu.contigs.fasta'))) == 1
    # Even manually restored stale contigs and an old done marker must not win.
    (empty / 'canu').mkdir()
    (empty / 'canu/canu.contigs.fasta').write_text('>stale\nTTTT\n')
    (empty / 'done').write_text((directory(files, COVERED) / 'done').read_text())
    combine_assembly_sequences(files, [UNCOVERED, COVERED])
    records = list(SeqIO.parse(files.assembly_fasta, 'fasta'))
    assert len(records) == 1 and str(records[0].seq) == 'ACGT'
    assert len([name for name, _ in stub_tools() if name == 'canu']) == 1
    assert execute(files, [UNCOVERED, COVERED]) == []


def test_all_uncovered_finishes_without_empty_assembly(tmp_path, stub_tools):
    files = make_files(tmp_path)
    blocks = [UNCOVERED, ('igh', 0, 500, '0')]
    execute(files, blocks)
    assert stub_tools() == []
    assert combine_assembly_sequences(files, blocks) == 0
    assert not Path(files.assembly_fasta).exists()
    Path(files.assembly_fasta).write_text('>previous\nACGT\n')
    assert combine_assembly_sequences(files, blocks) == 0
    assert not Path(files.assembly_fasta).exists()
    assert next(tmp_path.glob('assembly.fasta.previous-*')).read_text() == '>previous\nACGT\n'


def test_empty_igh_subset_does_not_abort_other_assembly(tmp_path, stub_tools):
    files = make_files(tmp_path)
    blocks = [COVERED, ('igh', 0, 500, '0')]
    Path(files.igh_assembly_fasta).write_text('>stale\nTTTT\n')
    execute(files, blocks)
    combine_assembly_sequences(files, blocks)
    assert Path(files.assembly_fasta).is_file()
    assert not Path(files.igh_assembly_fasta).exists()
    assert len(list(tmp_path.glob('igh.fasta.previous-*'))) == 1


@pytest.mark.parametrize('mode', ['failure', 'missing_failure', 'empty'])
def test_actual_tool_failure_is_not_no_coverage(tmp_path, stub_tools, monkeypatch, mode):
    files = make_files(tmp_path)
    monkeypatch.setenv('IG_TEST_CANU', mode)
    for _ in range(2):
        with pytest.raises(subprocess.CalledProcessError):
            execute(files, [COVERED])
        root = directory(files, COVERED)
        assert not (root / 'done').exists()
        assert not (root / 'skipped.json').exists()
        assert not (root / 'canu/canu.contigs.fasta').exists()
    assert len([name for name, _ in stub_tools() if name == 'canu']) == 2
    with pytest.raises(RuntimeError):
        combine_assembly_sequences(files, [COVERED])


def test_polishing_failure_is_fatal(tmp_path, stub_tools, monkeypatch):
    files = make_files(tmp_path, read_type='SUBREAD', qualities=(0.85,))
    monkeypatch.setenv('IG_TEST_POLISH_FAILURE', 'yes')
    with pytest.raises(subprocess.CalledProcessError):
        execute(files, [COVERED])
    root = directory(files, COVERED)
    assert not (root / 'done').exists()
    assert not (root / 'contigs.fasta').exists()
    assert not (root / 'skipped.json').exists()


def test_raw_requires_actionable_missing_subreads_error(tmp_path):
    files = make_files(tmp_path, read_type='SUBREAD')
    Path(files.subreads_to_ref_phased).unlink()
    with pytest.raises(ValueError, match='requires phased subreads BAM.*mapping/phasing'):
        get_assembly_scripts(files, CPU, [COVERED])


@pytest.mark.parametrize('block', [('unknown', 0, 500, '0'), ('chr1', 9000, 11000, '0'), ('chr1', 0, 500, '9')])
def test_invalid_region_is_not_skipped(tmp_path, stub_tools, block):
    files = make_files(tmp_path)
    with pytest.raises(subprocess.CalledProcessError):
        execute(files, [block])
    assert not (directory(files, block) / 'skipped.json').exists()
    assert stub_tools() == []


def test_missing_index_is_not_no_coverage(tmp_path, stub_tools):
    files = make_files(tmp_path)
    Path(files.ccs_to_ref_phased + '.bai').unlink()
    with pytest.raises(ValueError, match='index'):
        execute(files, [UNCOVERED])
    assert stub_tools() == []


def test_skip_receipt_invalidates_when_reads_change(tmp_path, stub_tools):
    files = make_files(tmp_path)
    execute(files, [UNCOVERED])
    # Write a new phased alignment at the previously uncovered locus.
    bam = files.ccs_to_ref_phased
    with pysam.AlignmentFile(bam, 'rb') as source:
        header, reads = source.header, list(source)
    reads[0].reference_start = 5000
    with pysam.AlignmentFile(bam, 'wb', header=header) as out:
        out.write(reads[0])
    pysam.index(bam)
    assert len(execute(files, [UNCOVERED])) == 1
    combine_assembly_sequences(files, [UNCOVERED])
    assert not (directory(files, UNCOVERED) / 'skipped.json').exists()
    assert len([name for name, _ in stub_tools() if name == 'canu']) == 1


@pytest.mark.parametrize('read_type,ds', [(None, None), ('SEGMENT', 'READTYPE=SEGMENT;SOURCE=CCS')])
def test_supported_ccs_metadata_and_name_recovery(tmp_path, stub_tools, read_type, ds):
    files = make_files(tmp_path, read_type=read_type, ds=ds, qualities=(0.999,))
    assert select_assembly_workflow(files.input_bam).name == 'hifi'
    execute(files, [COVERED])
    combine_assembly_sequences(files, [COVERED])


@pytest.mark.parametrize('ds', ['READTYPE=UNKNOWN', 'READTYPE=CCS;READTYPE=SUBREAD'])
def test_ambiguous_metadata_is_actionable(tmp_path, ds):
    files = make_files(tmp_path, ds=ds)
    with pytest.raises(ValueError, match='READTYPE'):
        select_assembly_workflow(files.input_bam)


def test_skip_receipt_cannot_be_copied_to_a_different_region(tmp_path, stub_tools):
    files = make_files(tmp_path)
    execute(files, [UNCOVERED])
    target = directory(files, COVERED)
    target.mkdir(parents=True)
    (target / 'skipped.json').write_text((directory(files, UNCOVERED) / 'skipped.json').read_text())
    assert len(execute(files, [COVERED])) == 1
    combine_assembly_sequences(files, [COVERED, UNCOVERED])
    assert len(list(SeqIO.parse(files.assembly_fasta, 'fasta'))) == 1


def test_malformed_bam_is_not_no_coverage(tmp_path, stub_tools):
    files = make_files(tmp_path)
    Path(files.ccs_to_ref_phased).write_bytes(b'not a BAM')
    with pytest.raises(ValueError):
        execute(files, [UNCOVERED])
    assert stub_tools() == []
    assert not (directory(files, UNCOVERED) / 'skipped.json').exists()


def test_raw_polishing_tools_are_required(tmp_path, stub_tools):
    files = make_files(tmp_path, read_type='SUBREAD', qualities=(0.85,))
    (tmp_path / 'bin/gcpp').unlink()
    scripts = get_assembly_scripts(files, CPU, [COVERED])
    result = subprocess.run(['/bin/bash', scripts[0]], capture_output=True, text=True,
                            env=dict(os.environ, PATH=str(tmp_path / 'bin') + ':/bin'))
    assert result.returncode != 0
    assert 'requires gcpp for polishing' in result.stderr
    assert not (directory(files, COVERED) / 'done').exists()
    assert not (directory(files, COVERED) / 'skipped.json').exists()


def test_dispatch_changes_invalidate_completed_output(tmp_path, stub_tools):
    files = make_files(tmp_path, qualities=(0.999,))
    execute(files, [COVERED])
    write_reads(Path(files.input_bam), 'SEQUEL', 'CCS', (0.9,))
    execute(files, [COVERED])
    modes = [('-pacbio-hifi' if '-pacbio-hifi' in args else '-pacbio')
             for name, args in stub_tools() if name == 'canu']
    assert modes == ['-pacbio-hifi', '-pacbio']
    combine_assembly_sequences(files, [COVERED])


def rewrite_input_groups(files, header_ids, read_id):
    path = files.input_bam
    with pysam.AlignmentFile(path, 'rb') as source:
        header, reads = source.header.to_dict(), list(source)
    header['RG'] = [dict(header['RG'][0], ID=value) for value in header_ids]
    temporary = path + '.tmp'
    with pysam.AlignmentFile(temporary, 'wb', header=header) as out:
        for read in reads:
            read.set_tag('RG', read_id)
            out.write(read)
    os.replace(temporary, path)
    pysam.index(path)


def test_barcode_parent_metadata_dispatch_and_collection(tmp_path, stub_tools, caplog):
    files = make_files(tmp_path, qualities=(0.999, 0.999))
    rewrite_input_groups(files, ['b5706d50'], 'b5706d50/0--0')
    assert select_assembly_workflow(files.input_bam).name == 'hifi'
    assert sum('using metadata from barcode parent' in r.message for r in caplog.records) == 1
    execute(files, [COVERED])
    combine_assembly_sequences(files, [COVERED])
    assert len(list(SeqIO.parse(files.assembly_fasta, 'fasta'))) == 1
    with pysam.AlignmentFile(files.input_bam, 'rb') as bam:
        assert all(read.get_tag('RG') == 'b5706d50/0--0' for read in bam)


def test_exact_barcode_group_takes_precedence():
    from IGenotyper.assembly.workflow import resolve_read_group
    groups = {'b5706d50': None, 'b5706d50/0--0': 'CCS'}
    assert resolve_read_group('b5706d50/0--0', groups) == 'b5706d50/0--0'


@pytest.mark.parametrize('read_id', ['deadbeef/0--0', 'b5706d50/not-a-barcode', 'b5706d50/0--0/extra'])
def test_invalid_barcode_groups_remain_errors(tmp_path, read_id):
    files = make_files(tmp_path)
    rewrite_input_groups(files, ['b5706d50'], read_id)
    with pytest.raises(ValueError, match='unknown RG'):
        select_assembly_workflow(files.input_bam)


def test_duplicate_parent_groups_are_rejected(tmp_path):
    files = make_files(tmp_path)
    rewrite_input_groups(files, ['b5706d50', 'b5706d50'], 'b5706d50/0--0')
    with pytest.raises(ValueError, match='Duplicate @RG'):
        select_assembly_workflow(files.input_bam)


def test_one_6447_base_read_is_not_treated_as_no_coverage(tmp_path, stub_tools, monkeypatch):
    files = make_files(tmp_path, qualities=(0.999,))
    path = files.ccs_to_ref_phased
    with pysam.AlignmentFile(path, 'rb') as source:
        header, read = source.header, next(source)
    read.query_sequence = 'A' * 6447
    read.cigarstring = '6447M'
    with pysam.AlignmentFile(path, 'wb', header=header) as out:
        out.write(read)
    pysam.index(path)
    monkeypatch.setenv('IG_TEST_CANU', 'failure')
    with pytest.raises(subprocess.CalledProcessError):
        execute(files, [COVERED])
    root = directory(files, COVERED)
    assert len(next(SeqIO.parse(root / 'reads.fasta', 'fasta'))) == 6447
    assert not (root / 'skipped.json').exists()
    assert not (root / 'done').exists()


def coverage_files(tmp_path, depth):
    files = make_files(tmp_path, qualities=(0.999,))
    write_reads(Path(files.ccs_to_ref_phased), 'SEQUEL', 'CCS', (0.999,) * depth, phased=True)
    files.target_regions = str(tmp_path / 'targets.bed')
    Path(files.target_regions).write_text('chr1\t100\t300\nchr1\t5000\t6000\n')
    annotations = tmp_path / 'annotations'; annotations.mkdir()
    files.reference_annotations = str(annotations)
    (annotations / 'IG_loci.bed').write_text('chr1\t100\t300\tigh\nchr1\t5000\t6000\ttrb\n')
    files.input_args = str(tmp_path / 'args.json')
    Path(files.input_args).write_text(json.dumps({'sample': 'synthetic'}))
    files.phased_blocks = 'unused'
    return files


def run_coverage_sample(files, monkeypatch, events, coverage_bed=None):
    from IGenotyper.commands import assembly
    monkeypatch.setattr(assembly, 'FileManager', lambda *args, **kwargs: files)
    def plan(*args):
        events.append('plan')
        return [COVERED]
    monkeypatch.setattr(assembly, 'get_phased_blocks', plan)
    monkeypatch.setattr(assembly.Align, 'map_assembly', lambda *args: events.append('map'))
    monkeypatch.setattr(assembly, 'phase_assembly', lambda *args: events.append('phase_assembly'))
    return assembly.run_assembly(False, 1, 8, False, 'unused', 2, files.tmp, None, coverage_bed)


@pytest.mark.parametrize('depth', [0, 1, 19, 20, 21])
def test_20x_sample_gate_before_assembly(tmp_path, stub_tools, monkeypatch, depth):
    from IGenotyper.command_lines.clt import signatures
    files = coverage_files(tmp_path, depth)
    vcf = tmp_path / 'valid-phasing.vcf'; vcf.write_text('preserve phasing')
    protected = [files.ccs_to_ref_phased, files.ccs_to_ref_phased + '.bai', str(vcf)]
    before = signatures(protected)
    events = []
    for _ in range(2):
        result = run_coverage_sample(files, monkeypatch, events)
        assert result['coverage']['mean_depth'] == depth
        assert result['coverage']['target_bases'] == 200
        assert result['retryable'] is False
        assert signatures(protected) == before
        if depth < 20:
            assert result['status'] == 'insufficient_coverage'
            assert result['assembly_completed'] is False
            assert not Path(files.assembly_fasta).exists()
            assert not (tmp_path / 'assembly').exists()
            assert events == []
            assert stub_tools() == []
        else:
            assert result['status'] == 'completed'
            assert result['assembly_completed'] is True
            assert Path(files.assembly_fasta).exists()
            assert 'map' in events and 'phase_assembly' in events
    assert json.loads((tmp_path / 'assembly_status.json').read_text()) == result


def test_coverage_counts_uncovered_bases_and_unions_overlaps(tmp_path):
    from IGenotyper.assembly.coverage import measure_ig_coverage
    files = coverage_files(tmp_path, 20)
    override = tmp_path / 'ig.bed'
    override.write_text('chr1\t100\t300\nchr1\t150\t250\nchr1\t400\t600\n')
    result = measure_ig_coverage(files, files.ccs_to_ref_phased, str(override))
    assert result['target_bases'] == 400
    assert result['aligned_bases'] == 4000
    assert result['mean_depth'] == 10
    assert result['below_threshold'] is True
    assert [region['mean_depth'] for region in result['intervals']] == [20, 0]


def test_coverage_excludes_secondary_supplementary_duplicate_and_qcfail(tmp_path):
    from IGenotyper.assembly.coverage import measure_ig_coverage
    files = coverage_files(tmp_path, 5)
    bam = files.ccs_to_ref_phased
    with pysam.AlignmentFile(bam, 'rb') as source:
        header, reads = source.header, list(source)
    with pysam.AlignmentFile(bam, 'wb', header=header) as out:
        for read, flag in zip(reads, [0, 256, 2048, 1024, 512]):
            read.flag = flag
            out.write(read)
    pysam.index(bam)
    assert measure_ig_coverage(files, bam)['mean_depth'] == 1


def test_coverage_counts_aligned_bases_not_deletions_or_clips(tmp_path):
    from IGenotyper.assembly.coverage import measure_ig_coverage
    files = coverage_files(tmp_path, 1)
    bam = files.ccs_to_ref_phased
    with pysam.AlignmentFile(bam, 'rb') as source:
        header, read = source.header, next(source)
    read.query_sequence = 'A' * 200
    read.cigarstring = '50S50M50D50M50S'
    with pysam.AlignmentFile(bam, 'wb', header=header) as out:
        out.write(read)
    pysam.index(bam)
    result = measure_ig_coverage(files, bam)
    assert result['aligned_bases'] == 100
    assert result['mean_depth'] == .5


def test_adequate_coverage_does_not_hide_canu_failure(tmp_path, stub_tools, monkeypatch):
    files = coverage_files(tmp_path, 20)
    monkeypatch.setenv('IG_TEST_CANU', 'failure')
    with pytest.raises(subprocess.CalledProcessError):
        run_coverage_sample(files, monkeypatch, [])
    result = json.loads((tmp_path / 'assembly_status.json').read_text())
    assert result['status'] == 'failed'
    assert result['retryable'] is True
    assert result['assembly_completed'] is False
    assert result['coverage']['mean_depth'] == 20


def test_improved_coverage_reconsiders_skipped_sample(tmp_path, stub_tools, monkeypatch):
    files = coverage_files(tmp_path, 1)
    assert run_coverage_sample(files, monkeypatch, [])['status'] == 'insufficient_coverage'
    write_reads(Path(files.ccs_to_ref_phased), 'SEQUEL', 'CCS', (0.999,) * 20, phased=True)
    assert run_coverage_sample(files, monkeypatch, [])['status'] == 'completed'


def test_invalid_coverage_coordinates_are_errors(tmp_path, stub_tools, monkeypatch):
    files = coverage_files(tmp_path, 1)
    override = tmp_path / 'wrong.bed'; override.write_text('missing\t0\t100\n')
    with pytest.raises(ValueError, match='outside the BAM reference'):
        run_coverage_sample(files, monkeypatch, [], str(override))
    assert json.loads((tmp_path / 'assembly_status.json').read_text())['status'] == 'failed'
    assert stub_tools() == []


def test_coverage_rhesus_uses_ig_only_target_bed(tmp_path):
    from IGenotyper.assembly.coverage import measure_ig_coverage
    files = coverage_files(tmp_path, 20)
    files.rhesus = True
    Path(files.target_regions).write_text('chr1\t100\t300\n')
    assert measure_ig_coverage(files, files.ccs_to_ref_phased)['mean_depth'] == 20


def test_coverage_threshold_does_not_round_up(tmp_path, stub_tools, monkeypatch):
    files = coverage_files(tmp_path, 20)
    bam = files.ccs_to_ref_phased
    with pysam.AlignmentFile(bam, 'rb') as source:
        header, reads = source.header, list(source)
    reads[-1].query_sequence = 'A' * 199
    reads[-1].cigarstring = '199M'
    with pysam.AlignmentFile(bam, 'wb', header=header) as out:
        for read in reads:
            out.write(read)
    pysam.index(bam)
    result = run_coverage_sample(files, monkeypatch, [])
    assert result['coverage']['mean_depth'] == 19.995
    assert result['status'] == 'insufficient_coverage'
    assert stub_tools() == []


@pytest.mark.parametrize('mode', ['missing', 'zero'])
@pytest.mark.parametrize('read_type', ['CCS', 'SUBREAD'])
def test_no_canu_contigs_is_reusable_skip(tmp_path, stub_tools, monkeypatch, mode, read_type):
    files = make_files(tmp_path, read_type=read_type)
    monkeypatch.setenv('IG_TEST_CANU', mode)
    execute(files, [COVERED])
    root = directory(files, COVERED)
    workflow = select_assembly_workflow(files.input_bam)
    expected = region_provenance(assembly_provenance(files, workflow), *COVERED)
    assert region_status(str(root), workflow, expected) == 'skipped_no_contigs'
    assert not (root / 'done').exists()
    assert not region_assembled(str(root), workflow, expected)
    # A valid skip must take precedence over stale contigs, even if restored later.
    stale = Path(assembly_contigs(str(root), workflow))
    stale.parent.mkdir(parents=True, exist_ok=True)
    stale.write_text('>stale\nTTTT\n')
    assert combine_assembly_sequences(files, [COVERED]) == 0
    assert not Path(files.assembly_fasta).exists()
    assert execute(files, [COVERED]) == []
    assert [name for name, _ in stub_tools()] == ['canu']
    # Changed source evidence invalidates a no-contigs receipt.
    write_reads(Path(files.input_bam), 'SEQUEL', read_type, (None, None))
    assert len(get_assembly_scripts(files, CPU, [COVERED])) == 1


def test_mixed_no_contigs_and_assembled_regions(tmp_path, stub_tools, monkeypatch):
    files = make_files(tmp_path)
    second = ('chr1', 100, 600, '0')
    monkeypatch.setenv('IG_TEST_CANU', 'missing')
    execute(files, [COVERED])
    monkeypatch.setenv('IG_TEST_CANU', 'success')
    execute(files, [second])
    assert combine_assembly_sequences(files, [COVERED, second]) == 1
    records = list(SeqIO.parse(files.assembly_fasta, 'fasta'))
    assert len(records) == 1 and 'c=chr1:100-600' in records[0].id


@pytest.mark.parametrize('mode', ['missing', 'zero'])
def test_sample_no_contigs_exits_cleanly_without_mapping(tmp_path, stub_tools, monkeypatch, mode):
    from IGenotyper.command_lines.clt import signatures
    files = coverage_files(tmp_path, 20)
    before = signatures([files.ccs_to_ref_phased, files.ccs_to_ref_phased + '.bai'])
    monkeypatch.setenv('IG_TEST_CANU', mode)
    events = []
    Path(files.assembly_fasta).write_text('>old\nACGT\n')
    for _ in range(2):
        result = run_coverage_sample(files, monkeypatch, events)
        assert result['status'] == 'no_contigs'
        assert result['assembly_completed'] is False
        assert result['retryable'] is False
        assert not Path(files.assembly_fasta).exists()
        assert 'map' not in events and 'phase_assembly' not in events
        assert signatures(before) == before
    assert json.loads((tmp_path / 'assembly_status.json').read_text()) == result
    assert [name for name, _ in stub_tools()] == ['canu']
