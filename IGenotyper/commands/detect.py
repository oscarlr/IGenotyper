#!/bin/env python
import json
import pysam
import os
from pybedtools import BedTool

from IGenotyper.files import FileManager
from IGenotyper.common.cpu import CpuManager

from IGenotyper.command_lines.variants import VariantTools

from IGenotyper.detect.snps import detect_snps
from IGenotyper.detect.msa_indels_svs import detect_msa_variants
from IGenotyper.detect.alleles import detect_alleles


#from IGenotyper.clt import CommandLine
#from IGenotyper.helper import assembly_location,intervals_overlapping,get_haplotype,load_bed_regions,vcf_header,get_phased_blocks,coords_not_overlapping,assembly_coords,create_directory,get_ref_seq,skip_read #non_overlapping,interval_intersection,contig_coords,hap_coords,non_overlapping_hap_coords,skip_read,coords_not_overlapping

#from IGenotyper.commands.msa.msa_to_variants_without_hmm import path_to_variants

def add_arguments(subparser):
    subparser.add_argument('--rhesus', default=False, action='store_true')
    subparser.add_argument('--hom', metavar='HOM', help='Add homozygous reference genotype')
    subparser.add_argument('outdir', metavar='OUTDIR', help='Directory for output')

def run_detect(outdir, hom, rhesus):
    files = FileManager(outdir, rhesus=rhesus)

    with open(files.input_args, 'r') as fh:
        phasing_args = json.load(fh)
    sample = str(phasing_args["sample"])  # Convert sample to string

    cpu = CpuManager()
    variants_command_line = VariantTools(files, cpu, sample)

    detect_snps(files, sample)
    #detect_alleles(files)
    #detect_msa_variants(files,sample,variants_command_line)

    os.system("python %s/AIM_pipeline.py %s %s" % (files.ancestry, sample, outdir))

def main(args):
    run_detect(**vars(args))
