"""Synthetic final BAMs/VCFs and stub pipeline actions for phase restart checks."""
import json
from pathlib import Path
from types import SimpleNamespace

import pysam
import pytest

from IGenotyper.command_lines.clt import signatures
from IGenotyper.phasing.completion import (
    final_outputs, provenance, record_completion, phasing_complete, receipt_path,
)


def fixture_files(tmp_path):
    files = SimpleNamespace(log=str(tmp_path), rhesus=False, data_directory=str(tmp_path), tmp=str(tmp_path))
    for attr in ['input_bam', 'ref', 'target_regions', 'gene_coords', 'sv_coords',
                 'ccs_to_ref_phased', 'phased_snps_vcf', 'phased_blocks', 'report',
                 'stats_json', 'gene_cov', 'plot_phasing', 'plot_gene_cov',
                 'plot_sv_gene_cov', 'input_args']:
        setattr(files, attr, str(tmp_path / attr))
        Path(getattr(files, attr)).write_text('synthetic\n')
    files.ccs_bam = files.input_bam
    Path(files.ref + '.fai').write_text('chr1\t1000\t6\t1000\t1001\n')
    header = {'HD': {'SO': 'coordinate'}, 'SQ': [{'SN': 'chr1', 'LN': 1000}],
              'RG': [{'ID': '0'}, {'ID': '1'}, {'ID': '2'}]}
    for path in [files.input_bam, files.ccs_to_ref_phased]:
        with pysam.AlignmentFile(path, 'wb', header=header) as bam:
            read = pysam.AlignedSegment()
            read.query_name = 'movie/1/ccs/fwd'
            read.query_sequence = 'ACGT'
            read.reference_id = 0
            read.reference_start = 10
            read.cigarstring = '4M'
            read.set_tag('RG', '0')
            bam.write(read)
        pysam.index(path)
    Path(files.phased_snps_vcf).write_text(
        '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n'
        'chr1\t11\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\n')
    Path(files.input_args).write_text(json.dumps(dict(bam=files.input_bam, sample='sample', input_vcf=None)))
    return files


def invoke(files, monkeypatch, threads=1):
    from IGenotyper.commands import phase
    monkeypatch.setattr(phase, 'FileManager', lambda *args: files)
    return phase.run_phasing(files.input_bam, files.log, 'sample', threads, 20,
                             False, 'unused', 1, None, files.tmp, False, files.data_directory)


def stub_pipeline(monkeypatch, events, fail=False):
    from IGenotyper.commands import phase
    monkeypatch.setattr(phase.ReadManip, 'turn_ccs_reads_to_fastq', lambda self: events.append('convert'))
    monkeypatch.setattr(phase.Align, 'map_ccs_reads', lambda self: events.append('map'))
    monkeypatch.setattr(phase, 'generate_phased_snps', lambda *args: events.append('vcf'))
    monkeypatch.setattr(phase, 'phase_ccs', lambda *args: events.append('phase'))
    monkeypatch.setattr(phase.Snps, 'phased_blocks_from_ccs_snps', lambda self: events.append('blocks'))
    def stats(*args):
        events.append('stats')
        if fail:
            raise RuntimeError('plot failure')
    monkeypatch.setattr(phase, 'phasing_stats', stats)


def test_finished_pipeline_skips_all_steps_even_with_new_threads(tmp_path, monkeypatch):
    files = fixture_files(tmp_path)
    events = []
    stub_pipeline(monkeypatch, events)
    invoke(files, monkeypatch)
    assert events == ['convert', 'map', 'vcf', 'phase', 'blocks', 'stats']
    assert receipt_path(files).exists()
    before = signatures(final_outputs(files))
    events.clear()
    invoke(files, monkeypatch, threads=16)
    assert events == []
    assert signatures(final_outputs(files)) == before


@pytest.mark.parametrize('attr', ['ccs_to_ref_phased', 'phased_snps_vcf', 'phased_blocks',
                                  'report', 'stats_json', 'gene_cov', 'plot_phasing',
                                  'plot_gene_cov', 'plot_sv_gene_cov', 'input_args'])
def test_missing_or_changed_final_output_invalidates_completion(tmp_path, attr):
    files = fixture_files(tmp_path)
    expected = provenance(files, 'sample', None)
    record_completion(files, expected)
    Path(getattr(files, attr)).write_text('')
    assert not phasing_complete(files, expected)


def test_index_and_inputs_and_sample_invalidate_completion(tmp_path):
    files = fixture_files(tmp_path)
    expected = provenance(files, 'sample', None)
    record_completion(files, expected)
    assert not phasing_complete(files, provenance(files, 'other', None))
    Path(files.target_regions).write_text('changed annotation')
    assert not phasing_complete(files, provenance(files, 'sample', None))
    Path(files.ccs_to_ref_phased + '.bai').unlink()
    assert not phasing_complete(files, expected)
    with pytest.raises(ValueError, match='Missing or empty'):
        record_completion(files, provenance(files, 'sample', None))


def test_failed_final_step_never_marks_pipeline_complete(tmp_path, monkeypatch):
    files = fixture_files(tmp_path)
    events = []
    stub_pipeline(monkeypatch, events, fail=True)
    with pytest.raises(RuntimeError, match='plot failure'):
        invoke(files, monkeypatch)
    assert not receipt_path(files).exists()
    assert not phasing_complete(files, provenance(files, 'sample', None))


def write_legacy_receipts(files):
    # Real schema-2 receipt structures; subprocess execution is stubbed here.
    for output in [files.ccs_to_ref_phased, files.phased_snps_vcf, files.phased_blocks,
                   files.plot_phasing, files.plot_gene_cov]:
        outputs = [output]
        if output == files.ccs_to_ref_phased:
            outputs.append(output + '.bai')
        if output == files.plot_gene_cov:
            outputs.append(files.plot_sv_gene_cov)
        state = dict(schema=2, command='synthetic successful command',
                     inputs=signatures([files.input_bam, files.ref]), outputs=signatures(outputs))
        Path(output + '.success.json').write_text(json.dumps(state))
    Path(files.input_args).touch()  # The prior pipeline writes args last.


def test_adopt_existing_successful_run_without_executing_pipeline(tmp_path, monkeypatch):
    files = fixture_files(tmp_path)
    write_legacy_receipts(files)
    events = []
    stub_pipeline(monkeypatch, events)
    invoke(files, monkeypatch)
    assert events == [] and receipt_path(files).exists()
    # Future checks need only the final receipt, inputs, and final outputs.
    for path in tmp_path.glob('*.success.json'):
        if path != receipt_path(files):
            path.unlink()
    invoke(files, monkeypatch)
    assert events == []


def test_nonempty_partial_vcf_without_success_receipt_is_not_adopted(tmp_path):
    files = fixture_files(tmp_path)
    write_legacy_receipts(files)
    Path(files.phased_snps_vcf + '.success.json').unlink()
    assert not phasing_complete(files, provenance(files, 'sample', None))


def test_stale_multifile_receipt_is_not_adopted(tmp_path):
    files = fixture_files(tmp_path)
    write_legacy_receipts(files)
    Path(files.plot_sv_gene_cov).write_text('changed')
    assert not phasing_complete(files, provenance(files, 'sample', None))


def test_corrupt_bam_cannot_be_marked_complete(tmp_path):
    files = fixture_files(tmp_path)
    Path(files.ccs_to_ref_phased).write_bytes(b'not a BAM')
    with pytest.raises(Exception):
        record_completion(files, provenance(files, 'sample', None))
    assert not receipt_path(files).exists()


def test_file_manager_respects_explicit_new_bam(tmp_path, monkeypatch):
    from IGenotyper.files import FileManager
    logs = tmp_path / 'logs'; logs.mkdir()
    (logs / 'args.json').write_text(json.dumps(dict(bam='old.bam', tmp='old-tmp')))
    monkeypatch.setattr(FileManager, 'file_structure', lambda *args: None)
    files = FileManager(str(tmp_path), bam='new.bam')
    assert files.input_bam == 'new.bam'
    assert files.tmp == str(tmp_path / 'tmp')
    assert FileManager(str(tmp_path)).input_bam == 'old.bam'
