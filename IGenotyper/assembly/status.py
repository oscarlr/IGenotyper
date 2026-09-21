"""Atomic sample outcomes for humans and automatic retry schedulers."""
import json
import os
from pathlib import Path
import tempfile


def write_sample_status(files, sample, status, provenance, coverage, error=None):
    destination = Path(files.assembly_fasta).parent / 'assembly_status.json'
    result = {'schema': 1, 'sample': sample, 'status': status,
              'assembly_completed': status == 'completed',
              'retryable': status in ('running', 'failed'),
              'provenance': provenance, 'coverage': coverage}
    if status == 'insufficient_coverage':
        result['message'] = ('Assembly skipped: mean IG target coverage %.3fx is below %sx. '
                             'Do not automatically retry unchanged input. Phasing outputs are preserved.'
                             % (coverage['mean_depth'], coverage['minimum_mean_depth']))
    if status == 'no_contigs':
        result['message'] = ('Assembly finished without contigs; mapping and assembly phasing skipped. '
                             'Phasing outputs are preserved. Do not automatically retry unchanged input.')
    if error is not None:
        result['error'] = '%s: %s' % (type(error).__name__, error)
    fd, temporary = tempfile.mkstemp(prefix='.assembly-status-', dir=destination.parent)
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(result, stream, indent=2)
        os.replace(temporary, destination)
    finally:
        Path(temporary).unlink(missing_ok=True)
    print('Assembly status: %s (%s)' % (status, destination))
    if 'message' in result:
        print(result['message'])
    return result
