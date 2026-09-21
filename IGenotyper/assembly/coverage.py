"""Pre-assembly mean depth over IG targets, including uncovered target bases."""
from collections import defaultdict
from pathlib import Path

import pysam

from IGenotyper.command_lines.clt import signatures

MIN_ASSEMBLY_COVERAGE = 20
IG_LOCI = {'igh', 'ighc', 'igk', 'igl'}


def read_intervals(path, ig_only=False):
    intervals = []
    with open(path) as stream:
        for number, line in enumerate(stream, 1):
            if not line.strip() or line.startswith(('#', 'track ', 'browser ')):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError('Invalid BED row %s in %s' % (number, path))
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError('Invalid BED interval on row %s in %s' % (number, path))
            if ig_only and (len(fields) < 4 or fields[3].lower() not in IG_LOCI):
                continue
            intervals.append((fields[0], start, end))
    return intervals


def merge_intervals(intervals):
    grouped = defaultdict(list)
    for chrom, start, end in intervals:
        grouped[chrom].append((start, end))
    merged = []
    for chrom in sorted(grouped):
        spans = []
        for start, end in sorted(grouped[chrom]):
            if spans and start <= spans[-1][1]:
                spans[-1][1] = max(end, spans[-1][1])
            else:
                spans.append([start, end])
        merged.extend((chrom, start, end) for start, end in spans)
    return merged


def coverage_beds(files, override=None):
    if override:
        return [str(override)]
    if getattr(files, 'rhesus', False):
        # This bundled BED contains only the three rhesus IG loci.
        return [files.target_regions]
    return [files.target_regions, str(Path(files.reference_annotations) / 'IG_loci.bed')]


def ig_target_intervals(beds):
    targets = read_intervals(beds[0])
    if len(beds) == 2:
        loci = read_intervals(beds[1], ig_only=True)
        targets = [(chrom, max(start, lo), min(end, hi))
                   for chrom, start, end in targets
                   for locus_chrom, lo, hi in loci
                   if chrom == locus_chrom and max(start, lo) < min(end, hi)]
    result = merge_intervals(targets)
    if not result:
        raise ValueError('No IG target intervals for coverage; supply a reference-matched --coverage-bed')
    return result


def measure_ig_coverage(files, bam_path, override=None):
    beds = coverage_beds(files, override)
    before = signatures([bam_path] + beds)
    intervals = ig_target_intervals(beds)
    measured = []
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        bam.check_index()
        for chrom, start, end in intervals:
            if chrom not in bam.references or end > bam.get_reference_length(chrom):
                raise ValueError('IG coverage interval %s:%s-%s is outside the BAM reference; check --coverage-bed' % (chrom, start, end))
            aligned_bases = 0
            for read in bam.fetch(chrom, start, end):
                # Same usable alignment flags as assembly extraction; all haplotypes.
                if read.flag & 3844:
                    continue
                if not read.query_sequence:
                    raise ValueError('Missing sequence for coverage read %s' % read.query_name)
                # Aligned query bases only: M/= /X contribute; D/N and clips do not.
                aligned_bases += sum(max(0, min(end, hi) - max(start, lo))
                                     for lo, hi in read.get_blocks())
            measured.append({'chrom': chrom, 'start': start, 'end': end,
                             'target_bases': end - start, 'aligned_bases': aligned_bases,
                             'mean_depth': aligned_bases / (end - start)})
    if signatures([bam_path] + beds) != before:
        raise RuntimeError('Coverage inputs changed during calculation; rerun assembly')
    target_bases = sum(region['target_bases'] for region in measured)
    aligned_bases = sum(region['aligned_bases'] for region in measured)
    return {'mean_depth': aligned_bases / target_bases,
            'minimum_mean_depth': MIN_ASSEMBLY_COVERAGE,
            'below_threshold': aligned_bases < MIN_ASSEMBLY_COVERAGE * target_bases,
            'target_bases': target_bases, 'aligned_bases': aligned_bases,
            'method': 'aligned_query_bases / union_IG_target_bases; all haplotypes; exclude flags 3844',
            'inputs': before, 'intervals': measured}
