#!/bin/env python
import pysam
from pybedtools import BedTool

from IGenotyper.common.plot import plot_histogram
from IGenotyper.common.helper import load_bed_regions,snps_from_reads,get_phased_blocks,show_value,skip_read,get_phased_regions,write_to_bashfile

def total_reads(bam):
    sam = pysam.AlignmentFile(bam,check_sq=False)
    assert int(sam.mapped) == 0
    reads = int(sam.unmapped)
    return reads    

def num_target_reads(bam):
    sam = pysam.AlignmentFile(bam)
    count = 0
    for read in sam:
        if skip_read(read):
            continue
        count += 1
    return count

def target_region_coverage(files,min_cov=10):
    regions = []
    sam = pysam.AlignmentFile(files.ccs_to_ref)
    capture_regions = load_bed_regions(files.target_regions)
    for chrom,start,end in capture_regions:
        bases = 0
        for pileupcolumn in sam.pileup(chrom,start,end):
            if pileupcolumn.n >= min_cov:
                bases += 1
        regions.append([chrom,start,end,bases])
    return regions

def read_lengths(bam_file):
    lengths = []
    samfile = pysam.AlignmentFile(bam_file, "rb",check_sq=False)
    for read in samfile:
        lengths.append(read.query_length)
    return lengths

def read_quality(bam_file):
    quality = []
    samfile = pysam.AlignmentFile(bam_file, "rb",check_sq=False)
    for read in samfile:
        score = sum(map(float,read.query_qualities))/read.query_length
        quality.append(score)
    return quality
        
def input_stats(files):
    ccs_read_lengths = read_lengths(files.ccs_bam)
    ccs_read_length_avg = sum(ccs_read_lengths)/len(ccs_read_lengths)    
    plot_histogram(ccs_read_lengths,"CCS read lengths",files.plot_ccs_read_lengths)
    ccs_quals = read_quality(files.ccs_bam)
    ccs_quals_avg = sum(ccs_quals)/len(ccs_quals)
    plot_histogram(ccs_quals,"CCS quality scores",files.plot_ccs_read_quals)
    regions_covered = target_region_coverage(files)
    igh_bases = None
    for chrom,start,end,num_bases in regions_covered:
        if chrom == "igh":
            igh_bases = num_bases
    stats = {
        "num_subreads": total_reads(files.input_bam),
        "num_ccs": total_reads(files.ccs_bam),
        "on_target_count": num_target_reads(files.ccs_to_ref),
        "ref_bases_ccs_coverage": igh_bases,
        "ccs_read_length_avg": ccs_read_length_avg,
        "ccs_quals_avg": round(ccs_quals_avg,2)
    }
    return stats
    
def phased_snps(files):
    num_phased_snps = {}
    num_unphased_snps = {}
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

def min_coverage(feature,cov = 10):
    return int(feature.name) >= 10
    
def covered_regions(bam):
    a = BedTool(bam)
    b = a.genomecov(bg=True)
    min_cov = b.filter(min_coverage)
    return min_cov
    
def phased_bases_per_chrom(files):
    phased_bases = {}
    phased_regions = BedTool(get_phased_blocks(files,files.phased_blocks))
    regions_with_coverage = covered_regions(files.ccs_to_ref)
    phased_regions_with_cov = phased_regions.intersect(regions_with_coverage)
    for region in phased_regions_with_cov:
        hap = show_value(region.name)
        if hap == "2":
            continue
        if hap == "0":
            continue
        chrom = show_value(region.chrom)
        start = int(show_value(region.start))
        end = int(show_value(region.end))
        if chrom not in phased_bases:
            phased_bases[chrom] = 0
        phased_bases[chrom] += (end - start)
    return phased_bases

def phased_stats(files,stats):
    num_phased_snps,num_unphased_snps = phased_snps(files)
    stats["num_phased_snps"] = num_phased_snps["igh"]
    stats["num_unphased_snps"] = num_unphased_snps["igh"]
    stats["phased_bases"] = phased_bases_per_chrom(files)["igh"]
    return stats

def phasing_plot(files,plot_tools_command_line,align_command_line):
    hapall_bw = "%s/hapall.bw" % files.tmp
    align_command_line.bam_to_bigwig(files.ccs_to_ref,hapall_bw)

    unphased_snps_bed = "%s/unphased_snps.bed" % files.tmp
    phased_snps_bed = "%s/phased_snps.bed" % files.tmp
    unphased_snps_bed_fh = open(unphased_snps_bed,'w')
    phased_snps_bed_fh = open(phased_snps_bed,'w')    
    snps = snps_from_reads(files)
    for chrom in snps:
        for pos in snps[chrom]:
            gt = snps[chrom][pos]
            out = [chrom,pos,pos+1]
            if gt in ["1/1"]:
                unphased_snps_bed_fh.write("%s\n" % "\t".join(map(str,out)))
            if gt in ["0|1","1|0"]:
                phased_snps_bed_fh.write("%s\n" % "\t".join(map(str,out)))
    unphased_snps_bed_fh.close()
    phased_snps_bed_fh.close()

    phased_blocks_bed = "%s/phased_blocks_bed" % files.tmp
    #phased_regions = get_phased_regions(files.phased_blocks,min_length=0,min_variants=0)
    with open(phased_blocks_bed,'w') as fh:        
        with open(files.phased_blocks,'r') as fofh:
            for line in fofh:
                line = line.rstrip().split('\t')
                if "#" in line[0]:
                    continue
                region = [line[1]] + line[3:]
                fh.write("%s\n" % "\t".join(map(str,region)))                   
    
    for i in ["0","1","2"]:
        hap_bw = "%s/hap%s.bw" % (files.tmp,i)
        align_command_line.hap_bam_to_bigwig(files.ccs_to_ref,i,hap_bw)

    params = {
        "hapall_bw": hapall_bw,
        "unphased_snps_bed": unphased_snps_bed,
        "phased_snps_bed": phased_snps_bed,
        "phased_blocks_bed": phased_blocks_bed,
        "hap0_bw": "%s/hap0.bw" % files.tmp,
        "hap1_bw": "%s/hap1.bw" % files.tmp,
        "hap2_bw": "%s/hap2.bw" % files.tmp,
        "genes": files.gene_coords,
        "sv": files.sv_coords
    }

    template = "%s/templates/pygenometracks.ini" % files.package_directory
    pygenometracks_config = "%s/pygenometracks.ini" % files.tmp
    write_to_bashfile(template,pygenometracks_config,params)            

    plot_tools_command_line.run_pygenometracks(pygenometracks_config,files.plot_phasing)
    
def phasing_stats(sample_name,files,plot_tools_command_line,align_command_line):
    stats = input_stats(files)
    stats = phased_stats(files,stats)
    phasing_plot(files,plot_tools_command_line,align_command_line)

    template = "%s/templates/report.html" % files.package_directory
    stats["sample"] = sample_name
    write_to_bashfile(template,files.report,stats)            
    
