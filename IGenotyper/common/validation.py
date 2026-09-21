"""Validation shared by command execution and synthetic regression tests."""
import pysam
from Bio import SeqIO


def require_usable_reads(path):
    # Match samtools fasta's default exclusion of secondary/supplementary reads.
    # Unmapped primary reads are valid input for mapping; SEQ must be present.
    with pysam.AlignmentFile(path, 'rb', check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            if not (read.flag & 0x900) and read.query_sequence:
                return
    raise ValueError('BAM contains zero usable primary reads with sequence: %s' % path)


def fasta_records(path):
    try:
        records = list(SeqIO.parse(path, 'fasta'))
    except (OSError, ValueError) as error:
        raise RuntimeError('Missing or invalid assembly FASTA: %s' % path) from error
    if not records or any(not record.seq for record in records):
        raise RuntimeError('Empty assembly FASTA: %s' % path)
    return records


def validate_vcf(path, sample=None, expected=None, reference=None, bam=None):
    """Reject malformed, unordered or duplicate sites; never repair/drop records.

    Return site identities to verify that genotype/phase retained every candidate.
    Empty candidate sets are allowed: successful exit, not size, proves completion.
    """
    sites = []
    seen_chroms = set()
    last_chrom, last_pos = None, 0
    with pysam.VariantFile(path) as vcf:
        if sample is not None and sample not in vcf.header.samples:
            raise ValueError('VCF %s lacks sample %s' % (path, sample))
        for record in vcf:
            if record.chrom != last_chrom:
                if record.chrom in seen_chroms:
                    raise ValueError('VCF chromosomes are not contiguous: %s' % path)
                seen_chroms.add(record.chrom)
                last_chrom, last_pos = record.chrom, 0
            if record.pos <= last_pos:
                raise ValueError('Unsorted or duplicate VCF position: %s:%s in %s' % (record.chrom, record.pos, path))
            if not record.alts or record.pos < 1:
                raise ValueError('Invalid variant in %s' % path)
            last_pos = record.pos
            sites.append((record.chrom, record.pos, record.ref, record.alts))
    if expected is not None and sites != expected:
        raise ValueError('VCF output does not retain all input sites: %s' % path)
    if reference is not None:
        with pysam.FastaFile(reference) as ref, pysam.AlignmentFile(bam, 'rb') as reads:
            reads.check_index()
            for chrom, pos, base, alts in sites:
                if chrom not in reads.references or chrom not in ref.references:
                    raise ValueError('Candidate contig %s is absent from BAM/reference' % chrom)
                if ref.fetch(chrom, pos - 1, pos - 1 + len(base)).upper() != base.upper():
                    raise ValueError('Candidate REF mismatch at %s:%s' % (chrom, pos))
    return sites
