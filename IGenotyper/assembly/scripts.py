#!/usr/bin/env python3
import os
import json
import hashlib
import uuid
import tempfile
from pathlib import Path
from shlex import quote

from IGenotyper.common.helper import create_directory, run_type
from IGenotyper.common.validation import fasta_records
from IGenotyper.command_lines.clt import signatures


def assembly_contigs(directory, platform, filetype='fasta'):
    if platform in ('SEQUELII', 'REVIO'):
        return '%s/canu/canu.contigs.%s' % (directory, filetype)
    if platform == 'SEQUEL':
        return '%s/contigs.%s' % (directory, filetype)
    raise ValueError('Unsupported PacBio platform: %s' % platform)


def assembly_provenance(files, platform):
    inputs = [files.ccs_to_ref_phased]
    if platform == 'SEQUEL':
        inputs += [files.subreads_to_ref_phased, files.input_bam]
    return {'schema': 1, 'platform': platform, 'inputs': signatures(inputs),
            'script': hashlib.sha256(Path(files.assembly_script).read_bytes()).hexdigest()}


def region_assembled(directory, platform, provenance=None):
    if not os.path.isfile('%s/done' % directory):
        return False
    try:
        fasta_records(assembly_contigs(directory, platform))
    except RuntimeError:
        return False
    try:
        state = json.loads(Path(directory, "done").read_text())
        if state.get('outputs') != signatures([assembly_contigs(directory, platform)]):
            return False
        if provenance is not None and state != dict(provenance, outputs=state['outputs']):
            return False
    except (OSError, ValueError, AttributeError):
        return False
    return True


def record_assembly_success(directory, platform, provenance):
    contigs = assembly_contigs(directory, platform)
    fasta_records(contigs)
    fd, temporary = tempfile.mkstemp(prefix='.done-', dir=directory)
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(dict(provenance, outputs=signatures([contigs])), stream)
        os.replace(temporary, os.path.join(directory, 'done'))
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def create_assemble_script(files, cpu, directory, chrom, start, end, hap):
    platform = run_type(files.input_bam)
    contigs = assembly_contigs(directory, platform)
    flank = 1000
    params = {
        'hap': hap,
        'ccs_to_ref': files.ccs_to_ref_phased,
        'region': '%s:%s-%s' % (chrom, max(1, int(start) - flank + 1), int(end) + flank),
        'output': directory,
        'threads': int(cpu.threads),
        'size': int(end) - int(start) + 2 * flank,
        'subreads': files.input_bam,
        'subreads_to_ref': files.subreads_to_ref_phased,
        'python_scripts': files.scripts,
        'pacbio_machine': platform,
        'contigs': contigs,
        'completion': json.dumps(assembly_provenance(files, platform)),
    }
    bashfile = '%s/assemble.sh' % directory
    with open(bashfile, 'w') as stream:
        stream.write('#!/bin/bash\nset -euo pipefail\n')
        for name, value in params.items():
            stream.write('%s=%s\n' % (name, quote(str(value))))
        stream.write(Path(files.assembly_script).read_text())
    return bashfile


def get_assembly_scripts(files, cpu, phased_blocks):
    assembly_scripts = []
    platform = run_type(files.input_bam)
    provenance = assembly_provenance(files, platform)
    for chrom, start, end, hap in phased_blocks:
        directory = '%s/assembly/%s/%s_%s/%s' % (files.tmp, chrom, start, end, hap)
        if region_assembled(directory, platform, provenance):
            continue
        if os.path.isdir(directory):
            # Preserve untrusted/stale intermediate results for recovery.
            os.replace(directory, directory + ".previous-" + uuid.uuid4().hex)
        create_directory(directory)
        assembly_scripts.append(create_assemble_script(files, cpu, directory, chrom, start, end, hap))
    return assembly_scripts
