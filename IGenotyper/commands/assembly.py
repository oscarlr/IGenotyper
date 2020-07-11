#!/bin/env python
from IGenotyper.cpu import CpuManager
from IGenotyper.files import FileManager
from IGenotyper.clt import CommandLine
from IGenotyper.helper import show_value,create_directory,write_to_bashfile

import os
import pybedtools
from collections import namedtuple

def add_arguments(subparser):
    subparser.add_argument('--threads', metavar='THREADS', default=1, help='Number of threads')
    subparser.add_argument('--mem', metavar='MEM', default=8, help='Memory for cluster')
    subparser.add_argument('--cluster', default=False, action='store_true', help='Use cluster')
    subparser.add_argument('--queue', metavar='QUEUE', default="premium", help='Queue for cluster')
    subparser.add_argument('--walltime', metavar='WALLTIME', default=2, help='Walltime for cluster')
    subparser.add_argument('--sample', metavar='SAMPLE', default="sample", help='Name of sample')
    subparser.add_argument('bam', metavar='BAM', help='PacBio bam file')
    subparser.add_argument('outdir',metavar='OUTDIR',help='Directory for output')

def get_phased_regions(files,min_length=500,min_variants=2):
    blocks = []
    Block = namedtuple('Block', ['sample','chrom','start_1','start','end','num_variants'])
    with open(files.phased_blocks, 'r') as fh:
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

def get_phased_blocks(files):
    phased_blocks = []
    target_regions = pybedtools.BedTool(files.target_regions)
    phased_regions = pybedtools.BedTool(get_phased_regions(files))
    unphased_regions = target_regions.subtract(phased_regions)
    for haplotype in ["1","2"]:
        phased_blocks = add_haplotype_to_blocks(phased_blocks,phased_regions,haplotype)
    phased_blocks = add_haplotype_to_blocks(phased_blocks,unphased_regions,"0")
    return phased_blocks

def region_assembled(directory):
    assembled = False
    if os.path.isfile("%s/done" % directory):
        assembled = True
    return assembled

def create_assemble_script(files,cpu,dir,chrom,start,end,hap):
    flank = 1000
    samtools_hap = "-r %s" % hap
    bashfile = "%s/assemble.sh" % dir
    params = {
        "hap": samtools_hap,
        "ccs_to_ref": files.ccs_to_ref_phased,
        "chrom": chrom,
        "start": max(0, int(start) - flank),
        "end": int(end) + flank,
        "output": dir,
        "threads": cpu.threads,
        "size": (int(end) - int(start)) + (flank * 2),
        "subreads": files.input_bam,
        "subreads_to_ref": files.subreads_to_ref_phased,
        "python_scripts": files.scripts,
        "ref": files.ref
    }
    write_to_bashfile(files.assembly_script,bashfile,params)
    return bashfile

def get_assembly_scripts(files,cpu,phased_blocks):
    assembly_scripts = []
    for chrom,start,end,hap in phased_blocks:
        dir = "%s/assembly/%s/%s_%s/%s" % (files.tmp,chrom,start,end,hap)
        if region_assembled(dir):
            continue
        create_directory(dir)
        assembly_script = create_assemble_script(files,cpu,dir,chrom,start,end,hap)
        assembly_scripts.append(assembly_script)
    return assembly_scripts

def run_assembly(
        threads,
        mem,
        cluster,
        queue,
        walltime,
        sample,
        bam,
        outdir
):
    files = FileManager(outdir,bam=bam)
    cpu = CpuManager(threads, mem, cluster, queue, walltime)
    command_line_tools = CommandLine(files,cpu)
    command_line_tools.phase_blocks(sample)
    phased_blocks = get_phased_blocks(files)
    assembly_scripts = get_assembly_scripts(files,cpu,phased_blocks)
    command_line_tools.run_assembly_scripts(assembly_scripts)

def main(args):
    run_assembly(**vars(args))