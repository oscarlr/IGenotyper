#!/bin/env python
import pysam
import numpy as np
from pybedtools import BedTool

from IGenotyper.common.plot import plot_histogram,plot_barplot
from IGenotyper.common.helper import load_bed_regions,snps_from_reads,get_phased_blocks,show_value,skip_read,get_phased_regions,write_to_bashfile,get_igh_region,non_emptyfile

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

def target_region_coverage(files,primary_alignment_bam,min_cov=10):
    regions = []
    sam = pysam.AlignmentFile(primary_alignment_bam)
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
        
def coverage_stats(files,bamfn):
    bam = BedTool(bamfn)
    bedgraph = bam.genomecov(bg=True)
    targets = BedTool(files.target_regions)
    bedgraph = bedgraph.intersect(targets)
    coverage = []    
    for row in bedgraph:
        if row[0] == "igh":
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            cov = int(row[3])
            for _ in range(start,end):
                coverage.append(cov)
    igh_chrom,igh_start,igh_end = get_igh_region(files.target_regions)
    igh = BedTool([(igh_chrom,igh_start,igh_end)])    
    igh = igh.subtract(bedgraph)
    for row in igh:
        for _ in range(row.start,row.end):
            coverage.append(cov)
    mean = round(np.mean(coverage),2)                  
    std = round(np.std(coverage),2)
    cv = round(std/mean,2)
    return (mean,std,cv)

def input_stats(files,primary_alignment_bam):
    ccs_read_lengths = read_lengths(files.ccs_bam)
    ccs_read_length_avg = sum(ccs_read_lengths)/len(ccs_read_lengths)    
    plot_histogram(ccs_read_lengths,"CCS read lengths",files.plot_ccs_read_lengths)
    ccs_quals = read_quality(files.ccs_bam)
    ccs_quals_avg = sum(ccs_quals)/len(ccs_quals)
    plot_histogram(ccs_quals,"CCS quality scores",files.plot_ccs_read_quals)
    regions_covered = target_region_coverage(files,primary_alignment_bam)
    igh_bases = None
    for chrom,start,end,num_bases in regions_covered:
        if chrom == "igh":
            igh_bases = num_bases
    target_read_count = num_target_reads(primary_alignment_bam)
    all_read_count = total_reads(files.ccs_bam)
    mean_cov, std_cov, cv_cov = coverage_stats(files,primary_alignment_bam)
    stats = {
        "num_subreads": total_reads(files.input_bam),
        "num_ccs": all_read_count,
        "on_target_count": target_read_count,
        "percent_on_target": round((float(target_read_count)/all_read_count),2) * 100,
        "ref_bases_ccs_coverage": igh_bases,
        "ccs_read_length_avg": ccs_read_length_avg,
        "ccs_quals_avg": round(ccs_quals_avg,2),
        "mean_coverage": mean_cov,
        "std_coverage": std_cov,
        "coeff_of_var": cv_cov
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
    
def phased_bases_per_chrom(files,primary_alignment_bam):
    phased_bases = {}
    phased_regions = BedTool(get_phased_blocks(files,files.phased_blocks))
    regions_with_coverage = covered_regions(primary_alignment_bam)
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

def feat_coverage(features,bam):
    sam = pysam.AlignmentFile(bam)
    features = BedTool(features)
    coverage = []
    for row in features:
        chrom = str(row[0])
        start = int(row[1])
        end = int(row[2])
        feature = str(row[3])
        row_coverage = [chrom,start,end,feature,0,0,0,0]
        for read in sam.fetch(reference=chrom,start=start,end=end):
            if skip_read(read):
                continue
            if read.reference_start > start:
                continue
            if read.reference_end < end:
                continue
            read_hap = int(read.get_tag("RG",True)[0])
            row_coverage[4] += 1
            row_coverage[5 + read_hap] += 1
        coverage.append(row_coverage)
    return coverage

# def feat_coverage(feat_fn,primary_alignment_bam):
#     ccs_alignment = BedTool(primary_alignment_bam)
#     ccs_bedgraph = ccs_alignment.genomecov(bg=True)
#     feat = BedTool(feat_fn)
#     feat_cov = feat.intersect(ccs_bedgraph,wao=True)
#     cov = {}
#     for gc in feat_cov:
#         chrom = str(gc[0])
#         start = int(gc[1])
#         end = int(gc[2])
#         name = str(gc[3])
#         coverage = gc[7]
#         if coverage != ".":
#             coverage = int(gc[7])
#         bases = int(gc[8])
#         if (chrom,start,end,name) not in cov:
#             cov[(chrom,start,end,name)] = 0
#         if coverage != ".":
#             cov[(chrom,start,end,name)] += (coverage*bases)
#     feat_coverage = []
#     for feat in cov:
#         chrom = feat[0]
#         start = feat[1]
#         end = feat[2]
#         name = feat[3]
#         coverage = float(cov[feat])/(feat[2] - feat[1])
#         feat_coverage.append([chrom,start,end,name,coverage])
#     return feat_coverage
    
def plot_gene_cov(files,gene_cov,plot_tools_command_line):
    with open(files.gene_cov,'w') as fh:
        for gene in gene_cov:
            fh.write("%s\n" % "\t".join(map(str,gene)))
    plot_tools_command_line.rplot_gene_cov()
    
    
def phased_stats(files,stats,plot_tools_command_line,primary_alignment_bam):
    num_phased_snps,num_unphased_snps = phased_snps(files)
    stats["num_phased_snps"] = num_phased_snps["igh"]
    stats["num_unphased_snps"] = num_unphased_snps["igh"]
    phased_bases = phased_bases_per_chrom(files,primary_alignment_bam)
    if "igh" in phased_bases:
        stats["phased_bases"] = phased_bases["igh"]
    else:
        stats["phased_bases"] = 0
    gene_cov = feat_coverage(files.gene_coords,primary_alignment_bam)
    gene_cov.sort(key = lambda x: x[-1])
    name = [i[3] for i in gene_cov]
    cov = [i[-1] for i in gene_cov]
    # plot_barplot(xvals,yvals,xlab,ylab,plotfn):
    # plot_barplot(name,cov,"Gene coverage","",files.plot_gene_cov)
    plot_gene_cov(files,gene_cov,plot_tools_command_line)
    stats["num_genes_no_cov"] = len([i for i in gene_cov if i[-1] == 0])
    return stats

def phasing_plot(files,plot_tools_command_line,align_command_line,primary_alignment_bam):
    hapall_bw = "%s/hapall.bw" % files.tmp
    align_command_line.bam_to_bigwig(primary_alignment_bam,hapall_bw)
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
                if line[3] == line[4]:
                    continue
                region = [line[1]] + line[3:]
                fh.write("%s\n" % "\t".join(map(str,region)))                   
    
    for i in ["0","1","2"]:
        hap_bw = "%s/hap%s.bw" % (files.tmp,i)
        align_command_line.hap_bam_to_bigwig(primary_alignment_bam,i,hap_bw)

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
    primary_alignment_bam = "%s/primary.bam" % files.tmp
    align_command_line.primary_alignments(files.ccs_to_ref_phased,primary_alignment_bam)

    stats = input_stats(files,primary_alignment_bam)
    stats = phased_stats(files,stats,plot_tools_command_line,primary_alignment_bam)
    if not non_emptyfile(files.plot_phasing):
        phasing_plot(files,plot_tools_command_line,align_command_line,primary_alignment_bam)

    template = "%s/templates/report.html" % files.package_directory
    stats["sample"] = sample_name
    write_to_bashfile(template,files.report,stats)            
    
    template = "%s/templates/stats.json" % files.package_directory
    write_to_bashfile(template,files.stats_json,stats)

