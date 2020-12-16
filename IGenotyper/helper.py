#!/bin/env python
import os
import sys
import shutil
import datetime
from string import Template

def create_folders(folders):
    for folder in folders:
        if not os.path.exists(folder):
            os.makedirs(folder)

def non_emptyfile(checkfile):
    return os.path.isfile(checkfile) and os.path.getsize(checkfile) > 0

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

def get_haplotype(read_name):
    return read_name.split("_")[1].split('=')[1]

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
                annotation = line[3]
                bed_regions.append([chrom,start,end,annotation])
            else:
                bed_regions.append([chrom,start,end])
    return bed_regions

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
