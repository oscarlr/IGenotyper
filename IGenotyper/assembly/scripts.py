#!/usr/bin/env python3
import os
import json
import hashlib
import uuid
import tempfile
import sys
from dataclasses import asdict
from pathlib import Path
from shlex import quote

import pysam

from IGenotyper.common.helper import create_directory
from IGenotyper.common.validation import fasta_records
from IGenotyper.command_lines.clt import signatures
from IGenotyper.assembly.workflow import select_assembly_workflow


def assembly_contigs(directory, workflow, filetype='fasta'):
    relative = 'contigs' if workflow.polish else 'canu/canu.contigs'
    return '%s/%s.%s' % (directory, relative, filetype)


def assembly_bam(files, workflow):
    return files.subreads_to_ref_phased if workflow.polish else files.ccs_to_ref_phased


def assembly_provenance(files, workflow):
    source = assembly_bam(files, workflow)
    if not os.path.isfile(source):
        requirement = 'phased subreads BAM (run subread mapping/phasing first)' if workflow.polish else 'phased CCS BAM (run phase first)'
        raise ValueError('%s workflow requires %s: %s' % (workflow.name, requirement, source))
    indexes = [source + '.bai', str(Path(source).with_suffix('.bai')), source + '.csi', str(Path(source).with_suffix('.csi'))]
    indexes = list(dict.fromkeys(path for path in indexes if os.path.isfile(path)))
    if not indexes:
        raise ValueError('Missing index for assembly BAM; run samtools index %s' % source)
    with pysam.AlignmentFile(source, 'rb') as bam:
        bam.check_index()
    inputs = list(dict.fromkeys([files.input_bam, source] + indexes))
    code = [Path(files.assembly_script), Path(__file__), Path(__file__).with_name('regions.py'), Path(__file__).with_name('workflow.py')]
    return {'schema': 2, 'workflow': json.loads(json.dumps(asdict(workflow))),
            'inputs': signatures(inputs),
            'script': hashlib.sha256(b''.join(path.read_bytes() for path in code)).hexdigest()}


def region_provenance(provenance, chrom, start, end, hap):
    return dict(provenance, region=[chrom, int(start), int(end), str(hap)], flank=1000)


def region_failure(directory, provenance=None):
    """Read a current Canu failure receipt; failure is never assembly success."""
    try:
        result = json.loads(Path(directory, 'failed.json').read_text())
        if result['status'] != 'failed_canu' or not isinstance(result['exit_code'], int) or result['exit_code'] <= 0:
            return None
        if provenance is not None and result['provenance'] != provenance:
            return None
        if signatures([result['log']]) != result['logs']:
            return None
        return result
    except (OSError, ValueError, KeyError, TypeError):
        return None


def record_canu_failure(directory, provenance, exit_code, workdir):
    validate_assembly_inputs(provenance)
    log = str(Path(workdir, 'canu.log').resolve())
    result = dict(status='failed_canu', provenance=provenance,
                  region=provenance['region'], exit_code=int(exit_code),
                  workdir=str(Path(workdir).resolve()), log=log, logs=signatures([log]))
    fd, temporary = tempfile.mkstemp(prefix='.failure-', dir=directory)
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(result, stream, indent=2)
        for marker in ('done', 'skipped.json'):
            Path(directory, marker).unlink(missing_ok=True)
        os.replace(temporary, Path(directory, 'failed.json'))
    finally:
        Path(temporary).unlink(missing_ok=True)
    print('Canu failed for %s (exit %s); log: %s' % (provenance['region'], exit_code, log), file=sys.stderr)


def region_status(directory, workflow, provenance=None):
    """Only validated assembled or skipped outcomes are reusable."""
    if Path(directory, 'failed.json').exists():
        return 'failed_canu' if region_failure(directory, provenance) is not None else None
    try:
        # A skip takes precedence over any old contigs/done file in this directory.
        skip = Path(directory, 'skipped.json')
        state = json.loads((skip if skip.exists() else Path(directory, 'done')).read_text())
        status = state.get('status')
        allowed = ('skipped_no_coverage', 'skipped_no_contigs') if skip.exists() else ('assembled',)
        if status not in allowed:
            return None
        outputs = {} if skip.exists() else signatures([assembly_contigs(directory, workflow)])
        if state.get('outputs') != outputs:
            return None
        if provenance is not None and state != dict(provenance, status=status, outputs=outputs):
            return None
        if not skip.exists():
            fasta_records(assembly_contigs(directory, workflow))
        return status
    except (OSError, ValueError, RuntimeError, AttributeError):
        return None


def region_assembled(directory, workflow, provenance=None):
    return region_status(directory, workflow, provenance) == 'assembled'


def validate_assembly_inputs(provenance):
    if signatures(provenance['inputs']) != provenance['inputs']:
        raise RuntimeError('Assembly inputs changed during execution; rerun assembly')


def record_region_result(directory, provenance, status, contigs=None):
    validate_assembly_inputs(provenance)
    if status == 'assembled':
        fasta_records(contigs)
        outputs, marker, other = signatures([contigs]), 'done', 'skipped.json'
    elif status == 'skipped_no_coverage':
        reads = Path(directory, 'reads.fasta')
        if not reads.is_file() or reads.stat().st_size != 0:
            raise RuntimeError('Cannot record no coverage without a validated empty extraction')
        outputs, marker, other = {}, 'skipped.json', 'done'
    elif status == 'skipped_no_contigs':
        if contigs is None or (os.path.exists(contigs) and (not os.path.isfile(contigs) or os.path.getsize(contigs) != 0)):
            raise RuntimeError('Cannot record no contigs with a nonempty or invalid output')
        outputs, marker, other = {}, 'skipped.json', 'done'
    else:
        raise ValueError('Unknown assembly status: %s' % status)
    fd, temporary = tempfile.mkstemp(prefix='.status-', dir=directory)
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(dict(provenance, status=status, outputs=outputs), stream)
        Path(directory, other).unlink(missing_ok=True)
        os.replace(temporary, os.path.join(directory, marker))
    finally:
        Path(temporary).unlink(missing_ok=True)


def record_assembly_success(directory, workflow, provenance):
    record_region_result(directory, provenance, 'assembled', assembly_contigs(directory, workflow))


def create_assemble_script(files, cpu, directory, chrom, start, end, hap, workflow=None, provenance=None):
    workflow = workflow or select_assembly_workflow(files.input_bam)
    provenance = region_provenance(provenance or assembly_provenance(files, workflow), chrom, start, end, hap)
    params = {
        'hap': hap, 'assembly_bam': assembly_bam(files, workflow),
        'chrom': chrom, 'start': int(start), 'end': int(end),
        'region': '%s:%s-%s' % (chrom, max(1, int(start) - 1000 + 1), int(end) + 1000),
        'output': directory, 'threads': int(cpu.threads),
        'size': int(end) - int(start) + 2000,
        'subreads': files.input_bam, 'python_scripts': files.scripts,
        'data_setting': workflow.canu_flag, 'polish': int(workflow.polish),
        'contigs': assembly_contigs(directory, workflow),
        'completion': json.dumps(provenance), 'python': sys.executable,
    }
    bashfile = '%s/assemble.sh' % directory
    with open(bashfile, 'w') as stream:
        stream.write('#!/bin/bash\nset -euo pipefail\n')
        for name, value in params.items():
            stream.write('%s=%s\n' % (name, quote(str(value))))
        stream.write(Path(files.assembly_script).read_text())
    return bashfile


def get_assembly_scripts(files, cpu, phased_blocks, workflow=None):
    assembly_scripts = []
    workflow = workflow or select_assembly_workflow(files.input_bam)
    provenance = assembly_provenance(files, workflow)
    print('Assembly workflow: %s (%s, %s)' % (workflow.name, workflow.read_type, workflow.accuracy))
    for chrom, start, end, hap in phased_blocks:
        directory = '%s/assembly/%s/%s_%s/%s' % (files.tmp, chrom, start, end, hap)
        expected = region_provenance(provenance, chrom, start, end, hap)
        if region_status(directory, workflow, expected) in ('assembled', 'skipped_no_coverage', 'skipped_no_contigs'):
            continue
        if os.path.isdir(directory):
            os.replace(directory, directory + '.previous-' + uuid.uuid4().hex)
        create_directory(directory)
        assembly_scripts.append(create_assemble_script(files, cpu, directory, chrom, start, end, hap, workflow, provenance))
    return assembly_scripts
