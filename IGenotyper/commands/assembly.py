#!/bin/env python
import os
import json
import pybedtools
from Bio import SeqIO

from IGenotyper.files import FileManager

from IGenotyper.common.cpu import CpuManager
from IGenotyper.common.helper import non_emptyfile,get_phased_blocks

from IGenotyper.command_lines.assembly import Assembly
from IGenotyper.command_lines.alignments import Align
from IGenotyper.command_lines.snps import Snps

from IGenotyper.phasing.reads import phase_merged_seqs

from IGenotyper.assembly.scripts import get_assembly_scripts
from IGenotyper.assembly.merge_assembly import merge_assembly

def add_arguments(subparser):
    subparser.add_argument('--threads', metavar='THREADS', default=1, help='Number of threads')
    subparser.add_argument('--mem', metavar='MEM', default=8, help='Memory for cluster')
    subparser.add_argument('--cluster', default=False, action='store_true', help='Use cluster')
    subparser.add_argument('--queue', metavar='QUEUE', default="premium", help='Queue for cluster')
    subparser.add_argument('--walltime', metavar='WALLTIME', default=2, help='Walltime for cluster')
    subparser.add_argument('outdir',metavar='OUTDIR',help='Directory for output')

def combine_sequence(files,phased_blocks,outfile,type_):
    seqs = []
    for chrom, start, end, hap in phased_blocks:
        dir = "%s/assembly/%s/%s_%s/%s" % (files.tmp, chrom, start, end, hap)
        contig = "%s/contigs.%s" % (dir,type_)
        if os.path.isfile(contig):
            contigs = list(SeqIO.parse(contig,type_))
            total_contigs = len(contigs)
            for i,record in enumerate(contigs):
                record.id = "c=%s:%s-%s_h=%s_i=%s_t=%s_/0/0_0" % (chrom,start,end,hap,i,total_contigs)
                record.description = ""
                seqs.append(record)
    SeqIO.write(seqs,outfile,type_)

def combine_assembly_sequences(files,phased_blocks):
    combine_sequence(files,phased_blocks,files.assembly_fasta,"fasta")
    combine_sequence(files,phased_blocks,files.assembly_fastq,"fastq")

def run_assembly(
        threads,
        mem,
        cluster,
        queue,
        walltime,
        outdir
):
    files = FileManager(outdir)

    with open(files.input_args,'r') as fh:
        phasing_args = json.load(fh)
    sample = phasing_args["sample"]

    cpu = CpuManager(threads, mem, cluster, queue, walltime)
    snps = Snps(files,cpu,sample)
    assembly_command_line = Assembly(files,cpu,sample)
    align_command_line = Align(files,cpu,sample)
    #snps_command_line = Snps(files,cpu,sample)
    
    if not non_emptyfile(files.assembly_fastq):
        phased_blocks = get_phased_blocks(files)
        assembly_scripts = get_assembly_scripts(files,cpu,phased_blocks)
        assembly_command_line.run_assembly_scripts(assembly_scripts)
        combine_assembly_sequences(files,phased_blocks)

    align_command_line.map_assembly()

    merge_assembly(files,align_command_line,sample)
    snps.phase_snvs_with_merged_seq()
    snps.phased_blocks_from_merged_seq()
    phase_merged_seqs(files,sample)

    
def main(args):
    run_assembly(**vars(args))
