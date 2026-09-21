"""Whole-pipeline completion independent of disposable phasing intermediates."""
import json
import os
from pathlib import Path
import tempfile
import subprocess

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


def valid_command_receipt(path, visiting=None, recover=None):
    """Validate an existing receipt chain before adopting a pre-marker run."""
    visiting = set() if visiting is None else visiting
    if path in visiting:
        raise ValueError('Cyclic completion receipts')
    visiting.add(path)
    try:
        state = json.loads(Path(path + '.success.json').read_text())
        if (state['schema'] != 2 or path not in state['outputs'] or
                signatures(state['outputs']) != state['outputs']):
            raise ValueError('Stale command receipt: %s' % path)
        mismatches = set()
        for source, signature in state['inputs'].items():
            try:
                current = signatures([source])[source]
            except OSError:
                current = None
            if current != signature:
                mismatches.add(source)
        if mismatches and not (recover and recover(path, state, mismatches)):
            raise ValueError('Stale command inputs: %s' % path)
        sources = set(state['inputs'])
        for source in state['inputs']:
            if Path(source + '.success.json').exists():
                sources.update(valid_command_receipt(source, visiting, recover))
        return sources
    finally:
        visiting.remove(path)


def recover_legacy_blocks(files, sample, path, state, mismatches):
    """Only the old shared lengths dependency may be replaced by revalidation."""
    from IGenotyper.command_lines.snps import phased_blocks_command
    legacy = getattr(files, 'legacy_chr_lengths', None)
    if (path != files.phased_blocks or legacy is None or mismatches != {legacy}
            or set(state['inputs']) != {legacy, files.phased_snps_vcf}
            or set(state['outputs']) != {files.phased_blocks}
            or state['command'] != phased_blocks_command(sample, path, legacy, files.phased_snps_vcf)):
        return False
    return validate_existing_blocks(files, sample)


def validate_existing_blocks(files, sample):
    path = files.phased_blocks
    from IGenotyper.command_lines.snps import Snps, write_chromosome_lengths
    validate_vcf(files.phased_snps_vcf, sample)
    before = signatures([files.ref, files.ref + '.fai', files.phased_snps_vcf, path])
    with tempfile.TemporaryDirectory(prefix='.check-legacy-blocks-', dir=files.log) as work:
        lengths = str(Path(work, 'chr_lengths.txt'))
        regenerated = str(Path(work, 'phased_blocks.txt'))
        write_chromosome_lengths(files.ref, lengths)
        Snps(files, None, sample).phased_blocks(regenerated, lengths, files.phased_snps_vcf)
        if Path(regenerated).read_bytes() != Path(path).read_bytes():
            return False
    if signatures(before) != before:
        return False
    print('Validated legacy phase blocks against reference-derived lengths; preserving existing outputs.')
    return True


class LegacyPhasingError(RuntimeError):
    """Unverifiable old results must not be automatically overwritten."""


def adopt_untracked_phasing(files, expected):
    """Validate old completed runs that predate command success receipts."""
    try:
        validate_final_outputs(files, expected['sample'])
        completed_at = Path(files.input_args).stat().st_mtime_ns
        if any(Path(path).stat().st_mtime_ns > completed_at for path in expected['inputs']):
            raise ValueError('An input or reference is newer than the saved completed run')
        # Compare the phased VCF with retained source variants when available.
        # A parseable but truncated phased VCF is not adequate evidence.
        source_vcf = expected['input_vcf'] or getattr(files, 'snps_vcf', None)
        if not source_vcf or not non_emptyfile(source_vcf):
            raise ValueError('Missing source VCF needed to verify legacy phased variant completeness')
        if Path(source_vcf + '.success.json').exists():
            valid_command_receipt(source_vcf)
        sites = validate_vcf(source_vcf, expected['sample'])
        validate_vcf(files.phased_snps_vcf, expected['sample'], expected=sites)
        with pysam.AlignmentFile(files.ccs_to_ref_phased, 'rb') as bam:
            with open(files.ref + '.fai') as stream:
                reference = {fields[0]: int(fields[1]) for fields in
                             (line.rstrip().split('\t') for line in stream)}
            if dict(zip(bam.references, bam.lengths)) != reference:
                raise ValueError('Phased BAM reference dictionary does not match the reference index')
        if not validate_existing_blocks(files, expected['sample']):
            raise ValueError('Legacy phase blocks do not match the phased VCF')
        record_completion(files, expected)
        print('Adopted completed legacy phasing; existing BAM, variants and reports are unchanged.')
        return True
    except (OSError, ValueError, RuntimeError, KeyError, TypeError, subprocess.CalledProcessError, pysam.SamtoolsError) as error:
        raise LegacyPhasingError(
            'Existing legacy phasing could not be validated: %s. '
            'No rephasing was started and existing results were preserved. '
            'Restore the missing completion evidence or use a new output directory to rephase.' % error
        ) from error


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
        terminal_paths = [files.ccs_to_ref_phased, files.phased_snps_vcf,
                          files.phased_blocks, files.plot_phasing, files.plot_gene_cov]
        if not any(Path(path + '.success.json').exists() for path in terminal_paths):
            return adopt_untracked_phasing(files, expected)
        # Adoption requires proof of successful commands, not merely a parseable
        # partial VCF. These receipts also lead back to mapping and conversion.
        sources = set()
        for path in [files.ccs_to_ref_phased, files.phased_snps_vcf,
                     files.phased_blocks, files.plot_phasing, files.plot_gene_cov]:
            sources.update(valid_command_receipt(path, recover=lambda path, state, mismatches:
                recover_legacy_blocks(files, expected['sample'], path, state, mismatches)))
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
    except LegacyPhasingError:
        raise
    except (OSError, ValueError, RuntimeError, KeyError, TypeError, subprocess.CalledProcessError):
        return False
