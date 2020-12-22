#!/bin/env python
from IGenotyper.cpu import CpuManager
from IGenotyper.files import FileManager

from IGenotyper.helper import non_emptyfile,clean_up,remove_files
from IGenotyper.command_lines.reads import ReadManip
from IGenotyper.command_lines.alignments import Align

from IGenotyper.phasing.snps import generate_phased_snps
from IGenotyper.phasing.reads import phase_subreads,phase_ccs

import os
import json
from shutil import copyfile

def add_arguments(subparser):
    subparser.add_argument('--sample',metavar='SAMPLE',default="sample",help='Name of sample')
    subparser.add_argument('--threads',metavar='THREADS',default=1,help='Number of threads')
    subparser.add_argument('--mem',metavar='MEM',default=20,help='Memory for cluster')
    subparser.add_argument('--cluster',default=False, action='store_true', help='Use cluster')
    subparser.add_argument('--queue', metavar='QUEUE',default="premium",help='Queue for cluster')
    subparser.add_argument('--walltime', metavar='WALLTIME',default=24,help='Walltime for cluster')
    subparser.add_argument('--tmp', metavar='TMP', default="tmp", help='Temporary folder')
    subparser.add_argument('--input_vcf', metavar='VCF', help='Phased VCF file to phase reads')
    subparser.add_argument('bam', metavar='BAM', help='PacBio bam file')
    subparser.add_argument('outdir',metavar='OUTDIR',help='Directory for output')

def save_parameters(files,sample,input_vcf):
    paramaters = {
        "bam": files.input_bam,
        "sample": sample,
        "input_vcf": input_vcf,
        "tmp": files.tmp
    }
    with open(files.input_args,'w') as fh:
        json.dump(paramaters,fh,sort_keys=True, indent=4)

def run_phasing(
        bam,
        outdir,
        sample,
        threads,
        mem,
        cluster,
        queue,
        walltime,
        input_vcf,
        tmp
):
    files = FileManager(outdir,bam,tmp)
    if non_emptyfile(files.input_args):
        sys.exit(0)
        
    cpu = CpuManager(threads,mem,cluster,queue,walltime)
    reads_command_line = ReadManip(files,cpu,sample)
    align_command_line = Align(files,cpu,sample)
    
    reads_command_line.generate_ccs_reads()
    
    if not non_emptyfile(files.phased_snvs_vcf):
        reads_command_line.turn_ccs_reads_to_fastq()
        align_command_line.map_ccs_reads()
        align_command_line.map_subreads()
        
        if input_vcf is None:
            generate_phased_snps(files,cpu,sample)
        else:
            copyfile(input_vcf,files.phased_snvs_vcf)

    phase_ccs(files,sample)
    phase_subreads(files,sample)
    
    iterations = 2
    for iteration in range(0,iterations):
        remove_vcfs(files)
        fix_ccs_alignment(reads_command_line,files,iteration)
        fix_subread_alignment(reads_command_line,files,iteration)
        phase_ccs(files,sample)
        phase_subreads(files,sample)

    save_parameters(files,sample,input_vcf)
    clean_up(files)
    
def main(args):
    run_phasing(**vars(args))
