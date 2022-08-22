#!/bin/env python
import json
import pysam
from pybedtools import BedTool

from IGenotyper.common.helper import load_bed_regions,vcf_header,intervals_overlapping,snps_from_reads,get_phased_blocks,get_haplotype


# from IGenotyper.helper import assembly_location,intervals_overlapping,get_haplotype,load_bed_regions,vcf_header,get_phased_blocks,coords_not_overlapping,assembly_coords,create_directory,get_ref_seq,skip_read #non_overlapping,interval_intersection,contig_coords,hap_coords,non_overlapping_hap_coords,skip_read,coords_not_overlapping

#from IGenotyper.commands.msa.msa_to_variants_without_hmm import path_to_variants

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
            continue
        feat_pos = feat_coord[:3] #[feat_coord[1],feat_coord[2]]
        snp_pos = [chrom,pos,pos+1]
        if intervals_overlapping(feat_pos,snp_pos):
            feat = feat_coord[3]
    return feat

def alternate_allele(reads,ref_base,in_phased_region):
    aa = set()
    for read_name,read_base,read_qual,read_hap in reads:
        if in_phased_region and	read_hap == "0":
            continue
        if read_base != ref_base:
            aa.add(read_base)
    if len(aa) == 0:        
        aa = set(".")
    aa = list(aa)
    if len(aa) > 1:
        if "DEL" in aa:
            aa.remove("DEL")
    return ",".join(aa)

def quality(reads):
    qual = 0.0
    count = 0.0
    for read_name,read_base,read_qual,read_hap in reads:
        qual += read_qual
        count += 1
    return int(round(qual/count,0))

def unphased_genotype(reads,ref_base,sv):
    genotype = None
    bases = set()
    for read_name,read_base,read_qual,read_hap in reads:
        if read_base == ref_base:
            if sv == False:
                genotype = "0/0"
            else:
                genotype = "0/0" #"./0"
        elif read_base != ref_base:
            if read_base == "DEL":
                alt = "."
            else:
                alt = "1"
            if sv == False:
                genotype = "%s/%s" % (alt,alt)
            else:
                genotype = "%s/%s" % (alt,alt)
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
    # if hap1_base == None and hap2_base == None:
    #     genotype = "./."
    if hap1_base != None:
        if hap1_base != ref_base:
            genotype = "./1"
        else:
            genotype = "./0"
    elif hap2_base != None:
        if hap2_base != ref_base:
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
            if hap2_base == "DEL":
                genotype = "./0"
            else:
                genotype = "0/1"
    elif ref_base == hap2_base:
        if hap1_base == "DEL":
            genotype = "./0"
        else:
            genotype = "0/1"
    elif hap1_base == hap2_base:
        if hap1_base == "DEL":
            genotype = "./."
        else:
            genotype = "1/1"
    elif hap1_base != hap2_base:
        if "DEL" not in [hap1_base,hap2_base]:
            genotype = "1/2"
        else:
            genotype = "./1"
    return genotype
        
def genotype(reads,in_phased_region,ref_base,sv):
    if in_phased_region == False:
        gt = unphased_genotype(reads,ref_base,sv)
    else:
        hap1_base,hap2_base = phased_hap_bases(reads)
        assert None not in [hap1_base,hap2_base]
        # if None in [hap1_base,hap2_base]:        
        #     gt = missing_hap_genotype(hap1_base,hap2_base,ref_base)
        # else:
        gt = hap_genotype(hap1_base,hap2_base,ref_base)
    assert gt != None
    return gt

def snp_in_feature(feature_fn,chrom,pos):
    bed = load_bed_regions(feature_fn,True)
    return snp_feat_overlap(bed,chrom,pos)

def reads_in_pos(pileupcolumn,in_phased_region):
    reads = []
    # Iteration over every ref pos
    for pileupread in pileupcolumn.pileups:
        #hap = pileupread.alignment.get_tag("RG",True)[0]
        hap = get_haplotype(pileupread.alignment.query_name)
        if in_phased_region:
            if hap == "0":
                continue
        else:
            if hap in ["1","2"]:
                continue
        if pileupread.indel > 0:
            base = pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position + pileupread.indel]
        elif pileupread.is_del == 1:
            base = "DEL"
        else:
            base = pileupread.alignment.query_sequence[pileupread.query_position]
        if pileupread.alignment.query_qualities == None:
            qual = 60
        else:
            qual = pileupread.alignment.query_qualities[pileupread.query_position]
        reads.append([pileupread.alignment.query_name,
                      base,
                      qual,
                      hap])
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
    aa = alternate_allele(bases,ref_base,in_phased_region)
    sv = snp_feat_overlap(kwargs["sv"],chrom,pos)
    gt = genotype(bases,in_phased_region,ref_base,sv)
    if gt == "0/0":
        assert aa == "."
    else:
        assert aa != "."
    if aa == "DEL":
        assert "." in gt
        aa = "."
    snp = [
        chrom,
        pos + 1,
        ".",
        ref_base,
        aa,
        quality(bases),
        "PASS",
        snp_info(bases,chrom,pos,**kwargs),
        "GT",
        gt
        ]
    return "\t".join(map(str,snp))

def not_valid_pos(bases,chrom,pos,**kwargs):
    # Not allowed
    # 1.
    # 2.
    # 3. Positions with only unphased contigs in phased regions
    not_valid = False
    in_phased_region = snp_feat_overlap(kwargs["phased_regions"],chrom,pos)
    hap_bases = {}
    #b = set()
    for name,base,qual,hap in bases:
        #b.add(base)
        if in_phased_region and hap == "0":
            continue
        if (not in_phased_region) and (hap in ["1","2"]):
            continue        
        if hap not in hap_bases:
            hap_bases[hap] = set()
        hap_bases[hap].add(base)
    for hap in hap_bases:
        if len(hap_bases[hap]) != 1:
            not_valid = True
    # if not in_phased_region:
    #     #if len(b) != 1:
    #     if len(hap_bases["0"]) != 1:
    #          not_valid = True
    if in_phased_region:
        if "1" not in hap_bases:
            not_valid = True
        if "2" not in hap_bases:
            not_valid = True
    else:
        if "0" not in hap_bases:
            not_valid = True
    return not_valid
    
def detect_snp(chrom,pileupcolumn,**kwargs):    
    in_phased_region = snp_feat_overlap(kwargs["phased_regions"],chrom,pileupcolumn.pos)
    bases = reads_in_pos(pileupcolumn,in_phased_region)
    if len(bases) == 0:
        return None
    if not_valid_pos(bases,chrom,pileupcolumn.pos,**kwargs):
        return None
    snp = vcf_snp(bases,chrom,pileupcolumn.pos,**kwargs)
    return snp

def just_phased_regions(files):
    regions = []
    phased_regions = get_phased_blocks(files,files.phased_blocks) #phased_blocks_merged_seq(files)
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
    #samfile = pysam.AlignmentFile(files.merged_assembly_to_ref_phased)
    capture_regions = load_bed_regions(files.target_regions)

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

    header = vcf_header(sample)
    for type_ in ["assembly","igh_assembly"]:
        if type_ == "assembly":            
            vcf_fh = open(files.snps_assembly_vcf, 'w')
            vcf_fh.write("%s\n" % "\n".join(header))    
            samfile = pysam.AlignmentFile(files.assembly_to_ref_phased)
        if type_ == "igh_assembly":
            vcf_fh = open(files.snps_igh_assembly_vcf, 'w')
            vcf_fh.write("%s\n" % "\n".join(header))    
            samfile = pysam.AlignmentFile(files.igh_assembly_to_ref_subs_phased)
        for chrom,start,end in capture_regions:
            if type_ == "igh_assembly":
                if chrom != "igh":
                    continue
            for pileupcolumn in samfile.pileup(chrom,start,end):
                snp = detect_snp(chrom,pileupcolumn,ref=ref,phased_regions=phased_regions,ccs_snps=ccs_snps,**bedfh)
                if snp != None:                    
                    vcf_fh.write("%s\n" % snp)
        vcf_fh.close()
