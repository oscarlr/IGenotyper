#!/bin/env python
from IGenotyper.cpu import CpuManager
from IGenotyper.files import FileManager
from IGenotyper.clt import CommandLine
#from IGenotyper.vcffn import read_in_phased_vcf
#from IGenotyper.bam import create_phased_bam_header
from IGenotyper.helper import non_emptyfile,clean_up
#from rephase import fix_phased_alignments

import os
import json
import pysam
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

# def create_tag(hap):
#     haptag = ("RG", str(hap), "Z")
#     return haptag

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
    command_line_tools = CommandLine(files,cpu)

    command_line_tools.generate_ccs_reads()
    
    if not non_emptyfile(files.phased_snvs_vcf):
        command_line_tools.turn_ccs_reads_to_fastq()
        command_line_tools.map_ccs_reads()
        command_line_tools.map_subreads()
        
        if input_vcf is None:
            detect_snps_from_reads(files.ccs_to_ref,files.ref,sample,files.snp_candidates,files.snvs_vcf)
            phase_snps_from_reads(sample,files.ref,files.hased_snvs_vcf,files.snvs_vcf,files.ccs_to_ref)
        else:
            copyfile(input_vcf,files.phased_snvs_vcf)

    phase_alignments(files.phased_snvs_vcf,files.ccs_to_ref,sample,files.ccs_to_ref_phased)
    phase_alignments(files.phased_snvs_vcf,files.subreads_to_ref,sample,files.subreads_to_ref_phased)
    
    iterations = 2
    for iteration in range(0,iterations):
        vcfs = [files.snp_candidates,files.snvs_vcf,files.phased_snvs_vcf]
        remove_files(vcfs)
        for phased_bam,unphased_bam in zip([files.ccs_to_ref_phased,files.subreads_to_ref_phased],
                                          [files.ccs_to_ref,files.subreads_to_ref]):
            samfile = fix_alignments(files.tmp,phased_bam,iteration)
            os.remove(unphased_bam)
            os.remove("%s.bai" % unphased_bam)
            command_line_tools.sam_to_sorted_bam(samfile,unphased_bam)
        detect_snps_from_reads(files.ccs_to_ref,files.ref,sample,files.snp_candidates,files.snvs_vcf)
        phase_snps_from_reads(sample,files.ref,files.phased_snvs_vcf,files.snvs_vcf,files.ccs_to_ref)
        phase_alignments(files.phased_snvs_vcf,files.ccs_to_ref,sample,files.ccs_to_ref_phased)
        phase_alignments(files.phased_snvs_vcf,files.subreads_to_ref,sample,files.subreads_to_ref_phased)

    save_parameters(files,sample,input_vcf)
    clean_up(files)
    
def main(args):
    run_phasing(**vars(args))
