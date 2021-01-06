#!/bin/env python
import pysam

from pybedtools import BedTool

def total_reads(bam):
    sam = pysam.AlignmentFile(bam)
    assert int(sam.mapped) == 0
    reads = int(sam.unmapped)
    return reads    

def num_target_reads(bam):
    sam = pysam.AlignmentFile(bam)
    reads = int(sam.mapped)
    return reads

def target_region_coverage(files,min_cov=10):
    regions = []
    sam = pysam.AlignmentFile(files.ccs_to_ref)
    capture_regions = load_bed_regions(files.target_regions)
    for chrom,start,end in capture_regions:
        bases = 0
        for pileupcolumn in samfile.pileup(chrom,start,end):
            if pileupcolumn.n >= min_cov:
                bases += 1
        regions.append([chrom,start,end,bases])
    return regions

def input_stats(files):
    num_subreads = total_subreads(files.input_bam)
    num_ccs = total_ccs_reads(files.ccs_bam)
    on_target_count = num_target_reads(files.ccs_to_ref)
    ref_bases_ccs_coverage = target_region_coverage(files)

def phased_snps(files):
    num_phased_snps,num_unphased_snps = {}
    snps = snps_from_reads(files)
    for chrom in snps:
        if chrom not in num_phased_snps:
            num_phased_snps[chrom] = 0
        if chrom not in num_unphased_snps:
            num_unphased_snps[chrom] = 0
        for pos in snps[chrom]:
            if "/" in snps[chrom][pos]:
                num_unphased_snps[chrom] += 1
            if "|" in snps[chrom][pos]:
                num_phased_snps[chrom] += 1
    return (num_phased_snps,num_unphased_snps)

def phased_bases_per_chrom(files):
    #capture_regions = load_bed_regions(files.target_regions)
    phased_regions = BedTool(get_phased_blocks(files,files.phased_blocks))
    regions_with_coverage = BedTool(target_region_coverage(files))
    phased_regions_with_cov = phased_regions.intersect(regions_with_coverage)
    for region in phased_regions_with_cov:
        print region
    

def phased_stats(files):
    num_phased_snps,num_unphased_snps = phased_snps(files)
    phased_bases = phased_bases_per_chrom(files)


def phasing_stats(files):
    input_stats(files)
    phased_stats(files)
