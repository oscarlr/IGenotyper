"""Whole-pipeline completion independent of disposable phasing intermediates."""
import json
import os
from pathlib import Path
import tempfile

import pysam

from IGenotyper.command_lines.clt import signatures, non_emptyfile
from IGenotyper.common.validation import validate_vcf


def final_outputs(files):
    return [files.ccs_to_ref_phased, files.ccs_to_ref_phased + '.bai',
            files.phased_snps_vcf, files.phased_blocks, files.report,
            files.stats_json, files.gene_cov, files.plot_phasing,
            files.plot_gene_cov, files.plot_sv_gene_cov, files.input_args]


def provenance(files, sample, input_vcf):
    inputs = [files.input_bam, files.ref, files.ref + '.fai', files.target_regions,
              files.gene_coords, files.sv_coords]
    if input_vcf:
        inputs.append(input_vcf)
    # Bump this version when phasing semantics change. Assembly-only edits and
    # resource settings must not invalidate successful phasing.
    return dict(schema=1, pipeline='phasing-completion:v1', sample=sample,
                rhesus=files.rhesus, input_vcf=os.path.abspath(input_vcf) if input_vcf else None,
                inputs=signatures(list(dict.fromkeys(inputs))))


def validate_final_outputs(files, sample):
    for path in final_outputs(files):
        if not non_emptyfile(path):
            raise ValueError('Missing or empty final phasing output: %s' % path)
    pysam.quickcheck(files.ccs_to_ref_phased)
    with pysam.AlignmentFile(files.ccs_to_ref_phased, 'rb') as bam:
        bam.check_index()
        # Read all records once at completion/adoption, catching corrupt blocks
        # and checking that this is a phased BAM (haplotype 0 is valid).
        count = 0
        for read in bam.fetch():
            if not read.has_tag('RG') or read.get_tag('RG') not in ('0', '1', '2'):
                raise ValueError('Invalid haplotype group in phased BAM')
            count += 1
        if not count:
            raise ValueError('Phased BAM has no mapped reads')
    validate_vcf(files.phased_snps_vcf, sample)


def receipt_path(files):
    return Path(files.log, 'phasing.success.json')


def record_completion(files, expected):
    if signatures(expected['inputs']) != expected['inputs']:
        raise RuntimeError('Phasing inputs changed during execution')
    validate_final_outputs(files, expected['sample'])
    result = dict(expected, outputs=signatures(final_outputs(files)))
    fd, temporary = tempfile.mkstemp(prefix='.phasing-success-', dir=files.log)
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(result, stream, indent=2)
        os.replace(temporary, receipt_path(files))
    finally:
        Path(temporary).unlink(missing_ok=True)


def valid_command_receipt(path, visiting=None):
    """Validate an existing receipt chain before adopting a pre-marker run."""
    visiting = set() if visiting is None else visiting
    if path in visiting:
        raise ValueError('Cyclic completion receipts')
    visiting.add(path)
    try:
        state = json.loads(Path(path + '.success.json').read_text())
        if (state['schema'] != 2 or path not in state['outputs'] or
                signatures(state['inputs']) != state['inputs'] or
                signatures(state['outputs']) != state['outputs']):
            raise ValueError('Stale command receipt: %s' % path)
        sources = set(state['inputs'])
        for source in state['inputs']:
            if Path(source + '.success.json').exists():
                sources.update(valid_command_receipt(source, visiting))
        return sources
    finally:
        visiting.remove(path)


def phasing_complete(files, expected):
    marker = receipt_path(files)
    try:
        if marker.exists():
            state = json.loads(marker.read_text())
            # Signatures cover the validated BAM/index and all final outputs.
            # Do not scan millions of reads again when every file is unchanged.
            return state == dict(expected, outputs=signatures(final_outputs(files)))
        args = json.loads(Path(files.input_args).read_text())
        if (args['sample'] != expected['sample'] or
                os.path.abspath(args['bam']) != os.path.abspath(files.input_bam) or
                (os.path.abspath(args['input_vcf']) if args.get('input_vcf') else None) != expected['input_vcf']):
            return False
        # Adoption requires proof of successful commands, not merely a parseable
        # partial VCF. These receipts also lead back to mapping and conversion.
        sources = set()
        for path in [files.ccs_to_ref_phased, files.phased_snps_vcf,
                     files.phased_blocks, files.plot_phasing, files.plot_gene_cov]:
            sources.update(valid_command_receipt(path))
        if not {files.input_bam, files.ref}.issubset(sources):
            return False
        # args.json was written last by the old pipeline. Changed final files or
        # newer annotations cannot be adopted using an older successful run.
        completed_at = Path(files.input_args).stat().st_mtime_ns
        if any(Path(path).stat().st_mtime_ns > completed_at
               for path in set(expected['inputs']) | set(final_outputs(files))):
            return False
        record_completion(files, expected)
        return True
    except (OSError, ValueError, RuntimeError, KeyError, TypeError):
        return False
