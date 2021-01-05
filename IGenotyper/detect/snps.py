#!/bin/env python
import json
import pysam
from pybedtools import BedTool

from IGenotyper.common.helper import load_bed_regions,vcf_header,phased_blocks_merged_seq,intervals_overlapping


# from IGenotyper.helper import assembly_location,intervals_overlapping,get_haplotype,load_bed_regions,vcf_header,get_phased_blocks,coords_not_overlapping,assembly_coords,create_directory,get_ref_seq,skip_read #non_overlapping,interval_intersection,contig_coords,hap_coords,non_overlapping_hap_coords,skip_read,coords_not_overlapping

#from IGenotyper.commands.msa.msa_to_variants_without_hmm import path_to_variants

def snps_from_reads(files):
    snps = {}
    with open(files.phased_snps_vcf,'r') as vcf_fh:
        for line in vcf_fh:
            if line.startswith("#"):
                continue
            line = line.rstrip().split("\t")
            chrom = line[0]
            position = line[1]
            genotype = line[9].split(":")[0]
            if chrom not in snps:
                snps[chrom] = {}                
            snps[chrom][int(position) - 1] = genotype
    return snps

def snp_read_genotype(ccs_snps,chrom,pos):
    genotype = False
    if chrom in ccs_snps:
        if pos in ccs_snps[chrom]:
            genotype = ccs_snps[chrom][pos]
    return genotype
    
def snp_read_support(ccs_snps,chrom,pos):
    support = False
    if chrom in ccs_snps:
        if pos in ccs_snps[chrom]:
            support = True
    return support
    
def snp_feat_overlap(feat_coords,chrom,pos):
    feat = False
    for feat_coord in feat_coords:
        if chrom != feat_coord[0]:
            return feat
        feat_pos = feat_coord[:3] #[feat_coord[1],feat_coord[2]]
        snp_pos = [chrom,pos,pos+1]
        if intervals_overlapping(feat_pos,snp_pos):
            feat = feat_coord[3]
    return feat

def alternate_allele(reads,ref_base,in_phased_region):
    aa = set()
    for read_name,read_base,read_qual,read_hap in reads:
        if read_base != ref_base:
            if in_phased_region:
                if read_hap != "0":
                    aa.add(read_base)
            else:
                if read_hap == "0":
                    aa.add(read_base)
    if len(aa) == 0:
        aa = set(".")
    return ",".join(list(aa))

def quality(reads):
    qual = 0.0
    count = 0.0
    for read_name,read_base,read_qual,read_hap in reads:
        qual += read_qual
        count += 1
    return int(round(qual/count,0))

def unphased_genotype(reads,ref_base):
    genotype = None
    bases = set()
    for read_name,read_base,read_qual,read_hap in reads:
        if read_base == ref_base:
	    genotype = "0/0"
        elif read_base != ref_base:
            genotype = "1/1"
        bases.add(read_base)
    assert len(bases) == 1
    return genotype

def phased_hap_bases(reads):
    hap1_base = None
    hap2_base = None
    for read_name,read_base,read_qual,read_hap in reads:
        if read_hap == "1":
            hap1_base = read_base
        if read_hap == "2":
            hap2_base = read_base    
    return (hap1_base,hap2_base)

def missing_hap_genotype(hap1_base,hap2_base,ref_base):
    genotype = None
    if hap1_base == None and hap2_base == None:
        genotype = "./."
    elif hap1_base == None:
        if hap2_base != ref_base:
            genotype = "./1"
        else:
            genotype = "./0"
    elif hap2_base == None:
        if hap1_base != ref_base:
            genotype = "./1"
        else:
            genotype = "./0"
    return genotype

def hap_genotype(hap1_base,hap2_base,ref_base):
    genotype = None
    #assert not ((ref_base == hap1_base) and (ref_base == hap2_base))
    if ref_base == hap1_base:
        if ref_base == hap2_base:
            genotype = "0/0"
        else:
            genotype = "0/1"
    elif ref_base == hap2_base:
        genotype = "0/1"
    elif hap1_base == hap2_base:
        genotype = "1/1"
    elif hap1_base != hap2_base:
        genotype = "1/2"
    return genotype
        
def genotype(reads,in_phased_region,ref_base):
    if in_phased_region == False:
        gt = unphased_genotype(reads,ref_base)
    else:
        hap1_base,hap2_base = phased_hap_bases(reads)
        if None in [hap1_base,hap2_base]:        
            gt = missing_hap_genotype(hap1_base,hap2_base,ref_base)
        else:
            gt = hap_genotype(hap1_base,hap2_base,ref_base)
    assert gt != None
    return gt

def snp_in_feature(feature_fn,chrom,pos):
    bed = load_bed_regions(feature_fn,True)
    return snp_feat_overlap(bed,chrom,pos)

def reads_in_pos(pileupcolumn):
    reads = []
    for pileupread in pileupcolumn.pileups:
        if pileupread.is_del:
            continue
        if pileupread.is_refskip:
            continue
        if pileupread.alignment.query_qualities == None:
            qual = 60
        else:
            qual = pileupread.alignment.query_qualities[pileupread.query_position]
        reads.append([pileupread.alignment.query_name,
                      pileupread.alignment.query_sequence[pileupread.query_position],
                      qual,
                      pileupread.alignment.get_tag("RG",True)[0]])
    return reads

def snp_info(bases,chrom,pos,**kwargs):
    read_support = snp_read_support(kwargs["ccs_snps"],chrom,pos)
    contig = ",".join([i[0] for i in bases])
    sv = snp_feat_overlap(kwargs["sv"],chrom,pos)
    intronic = snp_feat_overlap(kwargs["intronic"],chrom,pos)
    lp1 = snp_feat_overlap(kwargs["lp1"],chrom,pos)
    rss = snp_feat_overlap(kwargs["rss"],chrom,pos)
    gene = snp_feat_overlap(kwargs["gene"],chrom,pos)
    vdj = snp_feat_overlap(kwargs["vdj"],chrom,pos)
    read_genotype = snp_read_genotype(kwargs["ccs_snps"],chrom,pos)
    info = ("read_support=%s;"
            "contig=%s;"
            "SV=%s;"
            "intronic=%s;"
            "LP1=%s;"
            "RSS=%s;"
            "gene=%s;"
            "VDJ=%s;"
            "read_genotype=%s" %
            (read_support,
             contig,
             sv,
             intronic,                
             lp1,
             rss,
             gene,
             vdj,
             read_genotype))
    return info
    
def vcf_snp(bases,chrom,pos,**kwargs):    
    in_phased_region = snp_feat_overlap(kwargs["phased_regions"],chrom,pos)
    ref_base = kwargs["ref"].fetch(chrom,pos,pos + 1).upper()
    snp = [
        chrom,
        pos + 1,
        ".",
        ref_base,
        alternate_allele(bases,ref_base,in_phased_region),
        quality(bases),
        "PASS",
        snp_info(bases,chrom,pos,**kwargs),
        "GT",
        genotype(bases,in_phased_region,ref_base)
        ]
    return "\t".join(map(str,snp))

def not_valid_pos(bases,chrom,pos,**kwargs):
    not_valid = False
    in_phased_region = snp_feat_overlap(kwargs["phased_regions"],chrom,pos)
    hap_bases = {}
    b = set()
    for name,base,qual,hap in bases:
        b.add(base)
        if in_phased_region and hap == "0":
            continue
        if (not in_phased_region) and (hap in ["1","2"]):
            continue        
        if hap not in hap_bases:
            hap_bases[hap] = set()
        hap_bases[hap].add(base)
    for hap in hap_bases:
        if len(hap_bases[hap]) > 1:
            not_valid = True
    if not in_phased_region:
        if len(b) != 1:
             not_valid = True
    return not_valid
    
def detect_snp(chrom,pileupcolumn,**kwargs):    
    bases = reads_in_pos(pileupcolumn)
    if len(bases) == 0:
        return None
    if not_valid_pos(bases,chrom,pileupcolumn.pos,**kwargs):
        return None
    snp = vcf_snp(bases,chrom,pileupcolumn.pos,**kwargs)
    return snp

def just_phased_regions(files):
    regions = []
    phased_regions = phased_blocks_merged_seq(files)
    for region in phased_regions:
        r = region[0:3]
        h = region[3]
        if h == "0":
            r.append(False)
        else:
            r.append(True)
        if r not in regions:
            regions.append(r)
    return regions
    
def detect_snps(files,sample):
    samfile = pysam.AlignmentFile(files.merged_assembly_to_ref_phased)
    capture_regions = load_bed_regions(files.target_regions)
    header = vcf_header(sample)
    vcf_fh = open(files.snps_assembly_vcf, 'w')
    vcf_fh.write("%s\n" % "\n".join(header))

    ref = pysam.FastaFile(files.ref)
    phased_regions = just_phased_regions(files) #phased_blocks_merged_seq(files)
    ccs_snps = snps_from_reads(files)
    bedfh = {
        "sv": load_bed_regions(files.sv_coords,True),
        "intronic": load_bed_regions(files.introns,True),
        "lp1": load_bed_regions(files.lpart1,True),
        "rss": load_bed_regions(files.rss,True),
        "gene": load_bed_regions(files.gene_coords,True),
        "vdj": load_bed_regions(files.vdj_coords,True)
        }
    
    for chrom,start,end in capture_regions:
        for pileupcolumn in samfile.pileup(chrom,start,end):
            snp = detect_snp(chrom,pileupcolumn,ref=ref,phased_regions=phased_regions,ccs_snps=ccs_snps,**bedfh)
            if snp != None:                    
                vcf_fh.write("%s\n" % snp)
    vcf_fh.close()
