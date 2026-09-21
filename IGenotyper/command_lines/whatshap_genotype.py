"""WhatsHap 2.8 genotype adapter; no installed dependency files are modified.

Only zero-observation reads are removed, after WhatsHap groups alignments and
before its initial genotype computation. Candidate positions and all informative
observations remain unchanged. Other failures propagate normally.
"""
import logging
from importlib.metadata import version

ADAPTER_VERSION = 'empty-read-guard-v1'
logger = logging.getLogger(__name__)


def compute_genotypes_with_evidence(readset, positions=None):
    from whatshap.core import compute_genotypes

    keep = [i for i, read in enumerate(readset) if len(read) > 0]
    removed = len(readset) - len(keep)
    if removed:
        logger.warning('Excluded %d grouped reads with zero variant observations; retained %d informative reads', removed, len(keep))
    if not keep:
        # Do not invent calls or silently skip a chromosome with candidate sites.
        raise ValueError('No informative variant observations remain after read grouping; '
                         'genotyping cannot proceed. Candidate positions were not discarded.')
    return compute_genotypes(readset.subset(keep), positions)


def main():
    installed = version('whatshap')
    if installed != '2.8':
        raise RuntimeError('IGenotyper genotype adapter supports WhatsHap 2.8; found %s. '
                           'Use environment.yml or validate a new version before updating the adapter.' % installed)
    from whatshap.cli import genotype
    from whatshap.__main__ import main as whatshap_main

    original = genotype.compute_genotypes
    genotype.compute_genotypes = compute_genotypes_with_evidence
    try:
        whatshap_main()
    finally:
        genotype.compute_genotypes = original


if __name__ == '__main__':
    main()
