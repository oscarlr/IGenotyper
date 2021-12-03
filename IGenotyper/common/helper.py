#!/bin/env python
import os
import sys
import heapq
import pysam
import shutil
import datetime
import pybedtools
from string import Template
from collections import namedtuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def create_folders(folders):
    for folder in folders:
        if not os.path.exists(folder):
            os.makedirs(folder)

def non_emptyfile(checkfile):
    return os.path.isfile(checkfile) and os.path.getsize(checkfile) > 0

def remove_files(file_lists):
    for f in file_lists:
        os.remove(f)

def remove_vcfs(files):
    vcfs = [files.snp_candidates,files.snps_vcf,files.phased_snps_vcf]
    remove_files(vcfs)

def show_value(s):
    if sys.version_info.major == 2:
        if isinstance(s, unicode):
            return str(s)
    return s

def create_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def write_to_bashfile(template_bash,bashfile,params):
    filein = open(template_bash)
    src = Template(filein.read())
    output_lines = src.safe_substitute(params)
    bashfh = open(bashfile,'w')
    bashfh.write(output_lines)
    filein.close()
    bashfh.close()

def read_is_unphased(read):
    haplotype = read.get_tag("RG",True)[0]
    if haplotype == "0":
        return True
    return False

def clean_up(files):
    if os.path.isdir(files.tmp):
        shutil.rmtree(files.tmp)

def assembly_location(read_name):
    read_origin = read_name.split("_")[0].split('=')[1]
    chrom = read_origin.split(":")[0]
    start = int(read_origin.split(":")[1].split("-")[0])
    end = int(read_origin.split(":")[1].split("-")[1])
    return [chrom,start,end]

def get_haplotype(entry):
    #return read_name.split("_")[1].split('=')[1]
    hap_start = entry.index("h=") + len("h=")
    return entry[hap_start:hap_start + 1]

def intervals_overlapping(a, b):
    if a[0] != b[0]:
        return False
    overlapping = False
    num_overlapping = max(0, min(a[2], b[2]) - max(a[1], b[1]))
    if num_overlapping > 0:
        overlapping = True
    return overlapping

def load_bed_regions(bedfile,add_fourth=False):
    bed_regions = []
    with open(bedfile,'r') as bedfh:
        for line in bedfh:
            line = line.rstrip().split('\t')
            chrom = line[0]
            start = int(line[1])
            end = int(line[2])
            if add_fourth:
                annotation = True
                if len(line) == 4:
                    annotation = line[3]
                bed_regions.append([chrom,start,end,annotation])
            else:
                bed_regions.append([chrom,start,end])
    return bed_regions

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

def vcf_header(sample_name="sample"):
    i = datetime.datetime.now()
    line = [ "##fileformat=VCFv4.2",
             "##fileDate=%s%s%s" % (i.year,i.month,i.day),
             "##source=IGenotyper",
             "##INFO=<ID=SV,Number=1,Type=String,Description=\"Type of structural variant\">",
             "##INFO=<ID=contig,Number=2,Type=String,Description=\"Contig containing SNP\">",
             "##INFO=<ID=VDJ,Number=1,Type=String,Description=\"Type of region\">",
             "##INFO=<ID=read_support,Number=1,Type=String,Description=\"Support from CCS reads\">",
             "##INFO=<ID=intronic,Number=1,Type=String,Description=\"SNP in intron of gene\">",
             "##INFO=<ID=LP1,Number=1,Type=String,Description=\"SNP in leader part 1 sequence of gene\">",
             "##INFO=<ID=RSS,Number=1,Type=String,Description=\"SNP in recombination signal sequence of gene\">",
             "##INFO=<ID=gene,Number=1,Type=String,Description=\"SNP in gene\">",
             "##INFO=<ID=igh_region,Number=1,Type=String,Description=\"SNP in IGHV, IGHD or IGHJ\">",
             "##INFO=<ID=read_genotype,Number=1,Type=String,Description=\"Phased genotype only in phase in specified haplotype block\">",
             "##INFO=<ID=haplotype_block,Number=1,Type=String,Description=\"Haplotype block containing phased SNPs\">",
             "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",             
             "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">",
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s" % sample_name]
    return line


def non_overlapping(data):
    # from https://stackoverflow.com/questions/62102830/get-non-overlapping-distinct-intervals-from-a-set-of-intervals-python
    out = []
    starts = sorted([(i[0], 1) for i in data])  # start of interval adds a layer of overlap
    ends = sorted([(i[1], -1) for i in data])   # end removes one
    layers = 0
    current = []
    for value, event in heapq.merge(starts, ends):    # sorted by value, then ends (-1) before starts (1)
        layers += event
        if layers ==1:  # start of a new non-overlapping interval
            current.append(value)
        elif current:  # we either got out of an interval, or started an overlap
            current.append(value)
            out.append(current)
            current = []
    return out

# def interval_intersection(arr1,arr2):
#     i = 0
#     j = 0
    
#     n = len(arr1) 
#     m = len(arr2) 

#     intersect = []
#     # Loop through all intervals unless one  
#     # of the interval gets exhausted 
#     while i < n and j < m: 
          
#         # Left bound for intersecting segment 
#         l = max(arr1[i][0], arr2[j][0]) 
          
#         # Right bound for intersecting segment 
#         r = min(arr1[i][1], arr2[j][1]) 
          
#         # If segment is valid print it 
#         if l <= r:
#             intersect.append([l,r])
  
#         # If i-th interval's right bound is  
#         # smaller increment i else increment j 
#         if arr1[i][1] < arr2[j][1]: 
#             i += 1
#         else: 
#             j += 1

#     return intersect


def assembly_coords(assembly_to_ref):
    coords = []
    samfile = pysam.AlignmentFile(assembly_to_ref)
    for read in samfile:
        if skip_read(read):
            continue
        chrom = samfile.get_reference_name(read.reference_id)
        hap = read.get_tag("RG",True)[0]        
        # hap = get_haplotype(read.query_name)
        # if coords_not_overlapping(read,chrom):
        #     continue
        coords.append([chrom,int(read.reference_start),int(read.reference_end),hap,read.query_name])
    return coords

# def hap_coords(files):
#     coords = contig_coords(files)
#     hap_coords = {}
#     for contig_name in coords:
#         hap = get_haplotype(contig_name)
#         if hap not in hap_coords:
#             hap_coords[hap] = {}
#         hap_coords[hap][contig_name] = coords[contig_name]
#     return hap_coords

# def non_overlapping_hap_coords(files):
#     coords = hap_coords(files)
#     non_overlapping_coords = {}
#     for hap in coords:
#         chrom_coords = {}
#         for contig_name in coords[hap]:
# 	    for chrom,start,end in coords[hap][contig_name]:
#                 if chrom not in chrom_coords:
#                     chrom_coords[chrom] = []
#                 chrom_coords[chrom].append([start,end])
# 	non_overlapping_coords[hap] = {}
#         for chrom in chrom_coords:
#             non_overlapping_coords[hap][chrom] = non_overlapping(chrom_coords[chrom])
#     return non_overlapping_coords

def skip_read(read):
    skip = False
    if read.is_unmapped:
        skip = True
    if read.is_supplementary:
        skip = True
    if read.is_secondary:
        skip = True
    return skip

def coords_not_overlapping(read,mapped_chrom):
    assembled_region = assembly_location(read.query_name)
    mapped_start = int(read.reference_start)
    mapped_end = int(read.reference_end)
    mapped_region = [mapped_chrom,mapped_start,mapped_end]
    not_overlapping = True
    if intervals_overlapping(assembled_region,mapped_region):
        not_overlapping = False
    return not_overlapping

def get_phased_regions(phased_blocks,min_length=500,min_variants=2):
    blocks = []
    Block = namedtuple('Block', ['sample','chrom','start_1','start','end','num_variants'])
    with open(phased_blocks, 'r') as fh:
        header = fh.readline()
        for line in fh:
            line = line.rstrip().split('\t')
            block = Block._make(line)
            if int(block.num_variants) < min_variants:
                continue
            if (int(block.end) - int(block.start)) < min_length:
                continue
            blocks.append([block.chrom, int(block.start), int(block.end)])
    return sorted(blocks, key=lambda x: x[1])

def add_haplotype_to_blocks(phased_blocks,regions,haplotype):
    for region in regions:
        block = [
            show_value(region.chrom),
            show_value(region.start),
            show_value(region.end),
            haplotype
            ]
        phased_blocks.append(block)
    return phased_blocks

def get_phased_blocks(files,phased_blocks_source):
    phased_blocks = []
    target_regions = pybedtools.BedTool(files.target_regions)
    phased_regions = target_regions.intersect(pybedtools.BedTool(get_phased_regions(phased_blocks_source)))
    unphased_regions = target_regions.subtract(phased_regions)
    for haplotype in ["1","2"]:
        phased_blocks = add_haplotype_to_blocks(phased_blocks,phased_regions,haplotype)
    phased_blocks = add_haplotype_to_blocks(phased_blocks,unphased_regions,"0")
    return phased_blocks

# def phased_blocks_merged_seq(files):
#     return get_phased_blocks(files,files.phased_blocks_merged_seq)
    
def get_ref_seq(files,chrom,start,end):
    fasta = pysam.FastaFile(files.ref)
    return fasta.fetch(reference=chrom,start=start,end=end)

def get_mapping_pos(bam):
    pos = {}
    samfile = pysam.AlignmentFile(bam)
    for read in samfile:
        if skip_read(read):
            continue
        chrom = samfile.get_reference_name(read.reference_id)
        pos[read.query_name] = [chrom,read.reference_start,read.reference_end]
    return pos

def get_igh_region(target_regions_bed):
    igh_chrom = "igh"
    igh_start = None
    igh_end = None
    target_regions = pybedtools.BedTool(target_regions_bed)
    for target_region in target_regions:
        if str(target_region.chrom) == igh_chrom:
            igh_start = int(target_region.start)
            igh_end = int(target_region.end)
    assert igh_start != None
    assert igh_end != None
    return (igh_chrom,igh_start,igh_end)

def extract_sequence_from(read,chrom,start,end):
    read_start = None
    read_end = None
    aligned_pairs = read.get_aligned_pairs()
    for query_pos, ref_pos in aligned_pairs:
        if query_pos == None:
            continue
        if ref_pos == None:
            continue
        if ref_pos <= start:
            read_start = query_pos
        if ref_pos > end:
            break        
        read_end = query_pos
    if read_start == None:
        return ""
    if read_end == None:
        return ""
    return read.query_sequence[read_start:read_end].upper()

def extract_sequence(phased_bam,bed,fasta):
    records = []
    coords = load_bed_regions(bed,True)
    samfile = pysam.AlignmentFile(phased_bam)
    for chrom,start,end,feat in coords:
        i = 0
        for read in samfile.fetch(chrom,start,end):
            if skip_read(read):
                continue
            if "h=" not in read.query_name:
                hap = read.get_tag("RG",True)[0]
            else:
                hap = get_haplotype(read.query_name)
            name = "feat=%s_hap=%s_pos=%s:%s-%s_i=%s" % (feat,hap,chrom,start,end,i)
            seq = extract_sequence_from(read,chrom,start,end)
            if len(seq) == 0:
                continue
            record = SeqRecord(Seq(seq),id=name,name=name,description="")
            records.append(record)
            i += 1
    SeqIO.write(records,fasta,"fasta") 

def run_type(bam):
    sam = pysam.AlignmentFile(bam,check_sq=False)
    return dict(sam.header)['RG'][0]['PM']

