#!/bin/env python
import os

from IGenotyper.common.helper import create_directory,write_to_bashfile

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
