#!/bin/env python
import json
import os
import pysam

import warnings
from Bio import BiopythonWarning

# Filter out Biopython warnings
warnings.filterwarnings("ignore", category=BiopythonWarning)

from pybedtools import BedTool
from IGenotyper.files import FileManager
from IGenotyper.common.cpu import CpuManager
from IGenotyper.command_lines.variants import VariantTools
from IGenotyper.detect.snps import detect_snps
from IGenotyper.detect.alleles import detect_alleles

def add_arguments(subparser):
    subparser.add_argument('--rhesus', default=False, action='store_true')
    subparser.add_argument('--hom', metavar='HOM', help='Add homozygous reference genotype')
    subparser.add_argument('outdir', metavar='OUTDIR', help='Directory for output')

def run_detect(outdir, hom, rhesus):
    print("Initializing FileManager...")
    files = FileManager(outdir, rhesus=rhesus)
    print("FileManager initialized")

    with open(files.input_args, 'r') as fh:
        phasing_args = json.load(fh)
    sample = str(phasing_args["sample"])  # Convert sample to string explicitly

    print("Sample:", sample)

    print("Initializing CPU Manager...")
    cpu = CpuManager()
    print("CPU Manager initialized")

    print("Initializing VariantTools...")
    variants_command_line = VariantTools(files, cpu, sample)
    print("VariantTools initialized")

    print("Detecting SNPs...")
    detect_snps(files, sample)
    print("SNPs detected")

    print("Detecting Alleles...")
    detect_alleles(files)
    print("Alleles detected")

    print("Executing AIM_pipeline.py...")
    os.system("python %s/AIM_pipeline.py %s %s" % (files.ancestry, sample, outdir))
    print("AIM_pipeline.py executed")

def main(args):
    run_detect(**vars(args))
    print("Execution complete")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Description of your program")
    add_arguments(parser)
    args = parser.parse_args()
    main(args)
