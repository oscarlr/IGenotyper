#!/bin/env python
import json
import pysam
from pybedtools import BedTool

from IGenotyper.files import FileManager
from IGenotyper.cpu import CpuManager
from IGenotyper.clt import CommandLine

# from IGenotyper.helper import assembly_location,intervals_overlapping,get_haplotype,load_bed_regions,vcf_header,get_phased_blocks,coords_not_overlapping,assembly_coords,create_directory,get_ref_seq,skip_read #non_overlapping,interval_intersection,contig_coords,hap_coords,non_overlapping_hap_coords,skip_read,coords_not_overlapping

#from IGenotyper.commands.msa.msa_to_variants_without_hmm import path_to_variants

def snps_from_reads(files):
    snps = {}
    with open(files.phased_snvs_vcf,'r') as vcf_fh:
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

# def snp_info(files,snps):
#     snp_info_fields = {}
#     ccs_snps = snps_from_reads(files)
#     sv_coords = load_bed_regions(files.sv_coords,True)
#     introns = load_bed_regions(files.introns,True)
#     lpart1 = load_bed_regions(files.lpart1,True)
#     rss = load_bed_regions(files.rss,True)
#     gene_coords = load_bed_regions(files.gene_coords,True)
#     vdj_coords = load_bed_regions(files.vdj_coords,True)
#     for chrom in snps:
#         if chrom not in snp_info_fields:
#             snp_info_fields[chrom] = {}
#         for pos in snps[chrom]:
#             snp_info_fields[chrom][pos] = {
#                 "read_support": snp_read_support(ccs_snps,chrom,pos),
#                 "contig": snps[chrom][pos][4],
#                 "SV": snp_feat_overlap(sv_coords,chrom,pos),
#                 "intronic": snp_feat_overlap(introns,chrom,pos),
#                 "LP1": snp_feat_overlap(lpart1,chrom,pos),
#                 "RSS": snp_feat_overlap(rss,chrom,pos),
#                 "gene": snp_feat_overlap(gene_coords,chrom,pos),
#                 "VDJ": snp_feat_overlap(vdj_coords,chrom,pos),
#                 "read_genotype": snp_read_genotype(ccs_snps,chrom,pos)
#             }
#     return snp_info_fields

# def parse_info(info_field,g):
#     if g in ["./0","0/0"]:
#         read_support = "."
#         read_genotype = "."
#     else:
#         read_support = info_field["read_support"]
#         read_genotype = info_field["read_genotype"]
#     info = ("read_support=%s;"
#             "contig=%s;SV=%s;"
#             "intronic=%s;"
#             "LP1=%s;"
#             "RSS=%s;"
#             "gene=%s;"
#             "VDJ=%s;"
#             "read_genotype=%s" % (
#                 read_support,
#                 info_field["contig"],
#                 info_field["SV"],
#                 info_field["intronic"],
#                 info_field["LP1"],
#                 info_field["RSS"],
#                 info_field["gene"],
#                 info_field["VDJ"],
#                 read_genotype
#             ))
#     return info
    
# def write_snps(files,snps,snps_info_field,sample):
#     output_lines = vcf_header(sample)
#     with open(files.snvs_assembly_vcf, 'w') as vcf_fh:
#         vcf_fh.write("%s\n" % "\n".join(output_lines))
#         for chrom in snps:
#             snps_positions = snps[chrom].keys()
#             snps_positions.sort()
#             for pos in snps_positions:
#                 r,g,a,q,n = snps[chrom][pos]
#                 info = parse_info(snps_info_field[chrom][pos],g)
#                 rs_id = "."
#                 out = [chrom,pos + 1,rs_id,r,a,q,"PASS",info,"GT",g]
#                 vcf_fh.write("%s\n" % "\t".join(map(str,out)))

def get_alternate_allele(reads,ref_base):
    # What happens in one hap has two alleles?
    aa = set()
    for read_base,read_qual,read_hap in reads:
        if read_base != ref_base:
            aa.add(read_base)
    if len(aa) == 0:
        aa = set(".")
    return ",".join(list(aa))

def get_quality(reads):
    qual = 0.0
    count = 0.0
    for read_base,read_qual,read_hap in reads:
        qual += read_qual
        count += 1
    return int(round(qual/count,0))

def unphased_genotype(alternate_allele):
    alternate_alleles = alternate_allele.split(",")
    if len(alternate_alleles) == 1:
        return "1/1"
    return "./."

def phased_hap_bases(reads):
    hap1_base = None
    hap2_base = None
    for read_base,read_qual,read_hap in reads:        
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
    assert genotype != None
    return genotype

def hap_genotype(hap1_base,hap2_base,ref_base):
    genotype = None
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
        
def get_genotype(reads,in_phased_region,alternate_allele,ref_base):
    if alternate_allele == ".":
        return "0/0"
    if in_phased_region == False:
        return unphased_genotype(alternate_allele)
    hap1_base,hap2_base = phased_hap_bases(reads)
    if None in [hap1_base,hap2_base]:        
        genotype = missing_hap_genotype(hap1_base,hap2_base,ref_base)
    else:
        gentype = hap_genotype(hap1_base,hap2_base,ref_base)
    assert genotype != None
    return genotype            

def snp_in_feature(feature_fn,chrom,pos):
    bed = load_bed_regions(feature_fn,True)
    return snp_feat_overlap(bed,chrom,pos)

def read_bases(files,pileupcolumn):
    reads = {}
    for pileupread in pileupcolumn.pileups:
        if pileupread.is_del:
            continue
        if pileupread.is_refskip:
            continue                
        read_base = pileupread.alignment.query_sequence[pileupread.query_position]
        read_qual = pileupread.alignment.query_qualities[pileupread.query_position]
        read_hap = read.get_tag("RG",True)[0]
        reads[read.query_name] = [read_base,read_qual,read_hap]
        contigs.append(read.query_name)
    return reads

def snp_info(reads,files,chrom,pos):
    ref = pysam.FastaFile(files.ref)
    ref_base = ref.fetch(chrom,pos,pos + 1).upper()        
    ccs_snps = snps_from_reads(files)
    phased_regions = phased_blocks_merged_seq()
    in_phased_region = snp_feat_overlap(phased_regions,chrom,pos)
    alternate_allele = get_alternate_allele(reads,ref_base)
    quality = get_quality(reads,ref_base)
    genotype = get_genotype(reads,in_phased_region,alternate_allele,ref_base)
    read_support = snp_read_support(ccs_snps,chrom,pos)
    contig = ",".join(reads.keys())
    SV = snp_in_feature(files.sv_coords,pos)
    intronic = snp_in_feature(files.introns,pos)
    LP1 = snp_in_feature(files.lpart1,pos)
    RSS = snp_in_feature(files.rss,pos)
    gene = snp_in_feature(files.gene_coords,pos)
    VDJ = snp_in_feature(files.vdj_coords,pos)

def detect_snp(files,chrom,pileupcolumn):    
    reads = read_bases(files,pileupcolumn)
    if len(reads) == 0:
        return None
    info = snp_info(reads,files,chrom,pileupcolumn.pos)
    snp = []
    return snp
    
def detect_snps(files,sample):
    samfile = pysam.AlignmentFile(files.merged_assembly_to_ref_phased)
    capture_regions = load_bed_regions(files.target_regions,True)
    header = vcf_header(sample)
    vcf_fh = open(files.snvs_assembly_vcf, 'w')
    vcf_fh.write("%s\n" % "\n".join(header))
    for chrom,start,end in capture_regions:
        for pileupcolumn in samfile.pileup(chrom,start,end):
            snp = detect_snp(files,chrom,pileupcolumn)
            if snp != None:                    
                vcf_fh.write("%s\n" % snp)
    vcf_fh.close()
    
# def unphased_coords(chrom,hap0_coords,hap_coords):
#     phased_intervals = []
#     unphased_intervals = []
#     new_lists = [phased_intervals,unphased_intervals]
#     old_lists = [hap_coords,hap0_coords]
#     for new_list,list_ in zip(new_lists,old_lists):
#         for start,end in list_:
#             new_list.append((chrom,start,end))
#     phased_intervals = pybedtools.BedTool(phased_intervals)
#     unphased_intervals = pybedtools.BedTool(unphased_intervals)
#     no_phased_coords = unphased_intervals.subtract(phased_intervals)    
#     intervals = []
#     for feat in no_phased_coords:
#         intervals.append([feat[0],feat[1],feat[2]])
#     return intervals
    
# def overlapping_contig_coords(files):
#     coords = non_overlapping_hap_coords(files)
#     hap_1_chroms = coords["hap1"].keys()
#     hap_2_chroms = coords["hap2"].keys()
#     chroms = set(hap_1_chroms).intersect(set(hap_2_chroms))
#     chrom_intersect_coords = {"phase": {},"unphased":{}}
#     for chrom in chroms:
#         hap_1 = coords["hap1"][chrom]
#         hap_2 = coords["hap2"][chrom]
#         intersect = interval_intersection(hap_1,hap_2)
#         chrom_intersect_coords["phase"][chrom] = intersect
#         if chrom in coords["hap0"]:
#             chrom_intersect_coords["unphased"][chrom] = unphased_coords(coords["hap0"][chrom],intersect)
#     return chrom_intersect_coords

# def get_haps_coords(coords):
#     hap_coords = {}
#     for coord in coords:
#         chrom = coord[0]
#         start = coord[1]
#         end = coord[2]
#         hap = coord[3]
#         if hap not in hap_coords:
#             hap_coords[hap] = []
#         hap_coords[hap].append([chrom,start,end])
#     return hap_coords

# def is_coverage(feature,cov = 1):
#     return int(feature.name) == cov
    
# def msa_coords(files):
#     phased_regions = get_haps_coords(get_phased_blocks(files))
#     contig_coords = get_haps_coords(assembly_coords(files))
#     filtered_coords = {}
#     for hap in ["0","1","2"]:
#         contig_coords_bed = BedTool(contig_coords[hap])
#         contig_cov = contig_coords_bed.genomecov(bg=True,g="%s.fai" % files.ref)
#         contig_cov_filtered = contig_cov.filter(is_coverage)
#         phased_contig_cov_filtered = contig_cov_filtered.intersect(phased_regions[hap])
#         filtered_coords[hap] = phased_contig_cov_filtered
#     coords = {
#         "phased": filtered_coords["1"].intersect(filtered_coords["2"]),
#         "unphased": filtered_coords["0"]
#     }
#     return coords

# def extract_assembly_sequence(files,chrom,start,end,hap):
#     hap_sequence = None
#     samfile = pysam.AlignmentFile(files.assembly_to_ref,'rb')
#     for contig in samfile.fetch(chrom,start,end):
#         if skip_read(contig):
#             continue
#         if contig.reference_start > start:
#             continue
#         if contig.reference_end < end:
#             continue
#         if get_haplotype(contig.query_name) != hap:
#             continue
#         aligned_pairs = contig.get_aligned_pairs()
#         query_start = None
#         query_end = None
#         for query_pos, ref_pos in aligned_pairs:
#             if query_pos == None:
#                 continue
#             if ref_pos == None:
#                 continue
#             if int(ref_pos) <= int(start):
#                 query_start = query_pos
#             query_end = query_pos
#             if int(ref_pos) > int(end):
#                 break
#         assert query_start != None
#         assert query_end != None
#         hap_sequence = contig.query_sequence[query_start:query_end]
#         #sequence_qual = contig.query_qualities[query_start:query_end]
#         #assert len(sequence_qual) == len(hap_sequence)
#     return hap_sequence

# def extract_sequence(files,chrom,start,end,phase,haps):
#     outfasta = "%s/%s/fasta/%s_%s_%s.fasta" % (files.msa,phase,chrom,start,end)
#     with open(outfasta,'w') as outfasta_fh:
#         outfasta_fh.write(">ref\n%s\n" % get_ref_seq(files,chrom,start,end))
#         for h in haps:
#             seq = extract_assembly_sequence(files,chrom,start,end,h)
#             if h != "0":
#                 outfasta_fh.write(">hap%s\n%s\n" % (h,seq))
#             else:
#                 outfasta_fh.write(">hap1\n%s\n" % seq)
#                 outfasta_fh.write(">hap2\n%s\n" % seq)
#     return outfasta

# def call_variants(files,fastafile,chrom,start,end,phase,command_line_tools,variants):
#     msa_fn = "%s/%s/msa/%s_%s_%s.clu" % (files.msa,phase,chrom,start,end)
#     variants_fn = "%s/%s/variants/%s_%s_%s.bed" % (files.msa,phase,chrom,start,end)
#     command_line_tools.run_kalign(fastafile,msa_fn)
#     path_to_variants(msa_fn,chrom,start,end,variants_fn)    
#     with open(variants_fn,'r') as fh:
#         for line in fh:
#             line = line.strip().split('\t')
#             variants.append(line)
    
# def detect_non_snp_variants(files,sample,command_line_tools):
#     variants = []
#     coords = msa_coords(files)
#     for phase in coords:
#         if phase == "phased":
#             haps = ["1","2"]
#         else:
#             haps = ["0"]
#         for coord in coords[phase]:
#             chrom = str(coord[0])
#             start = int(coord[1])
#             end = int(coord[2])            
#             create_directory("%s/%s/fasta" % (files.msa,phase))
#             fastafile = extract_sequence(files,chrom,start,end,phase,haps)
#             create_directory("%s/%s/msa" % (files.msa,phase))
#             create_directory("%s/%s/variants" % (files.msa,phase))
#             call_variants(files,fastafile,chrom,start,end,phase,command_line_tools,variants)
#     variants.sort(key=lambda x: int(x[2]))
#     indels_fh = open(files.indels_assembly_bed,'w')
#     sv_fh = open(files.sv_assembly_bed,'w')
#     for variant in variants:
#         variant_size = int(variant[5])
#         if variant_size > 50:
#             sv_fh.write("%s\n" % "\t".join(variant))
#         else:
#             indels_fh.write("%s\n" % "\t".join(variant))
#     indels_fh.close()
#     sv_fh.close()
            
# def run_detect(outdir,hom):
#     files = FileManager(outdir)
#     cpu = CpuManager()
#     command_line_tools = CommandLine(files,cpu)
#     with open(files.input_args,'r') as fh:
#         phasing_args = json.load(fh)
#     sample = phasing_args["sample"]
#     #detect_snps(files,sample)
#     detect_non_snp_variants(files,sample,command_line_tools)
    
# def main(args):
#     run_detect(**vars(args))
    

