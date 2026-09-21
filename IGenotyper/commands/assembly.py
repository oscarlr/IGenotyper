#!/bin/env python
import os
import json
import tempfile
import filecmp
import pybedtools
from Bio import SeqIO

from IGenotyper.files import FileManager

from IGenotyper.common.cpu import CpuManager
from IGenotyper.common.helper import non_emptyfile,get_phased_blocks,run_type

from IGenotyper.command_lines.assembly import Assembly
from IGenotyper.command_lines.alignments import Align
#from IGenotyper.command_lines.snps import Snps

from IGenotyper.phasing.reads import phase_assembly

from IGenotyper.assembly.scripts import get_assembly_scripts, assembly_contigs, assembly_provenance, region_assembled
from IGenotyper.common.validation import fasta_records
#from IGenotyper.assembly.merge_assembly import merge_assembly

def add_arguments(subparser):
    subparser.add_argument('--rhesus',default=False, action='store_true')
    subparser.add_argument('--threads', metavar='THREADS', default=1, help='Number of threads')
    subparser.add_argument('--mem', metavar='MEM', default=8, help='Memory for cluster')
    subparser.add_argument('--cluster', default=False, action='store_true', help='Use cluster')
    subparser.add_argument('--queue', metavar='QUEUE', default="premium", help='Queue for cluster')
    subparser.add_argument('--walltime', metavar='WALLTIME', default=2, help='Walltime for cluster')
    subparser.add_argument('--data-dir', help='Directory containing reference.fasta')
    subparser.add_argument('outdir',metavar='OUTDIR',help='Directory for output')

def combine_sequence(files,phased_blocks,outfile,type_,chrom_select=None):
    seqs = []
    platform = run_type(files.input_bam)
    selected = [block for block in phased_blocks if chrom_select is None or block[0] == chrom_select]
    if not selected:
        raise RuntimeError("No assembly regions selected for %s" % outfile)
    for chrom, start, end, hap in selected:
        directory = "%s/assembly/%s/%s_%s/%s" % (files.tmp, chrom, start, end, hap)
        contig = assembly_contigs(directory, platform, type_)
        contigs = fasta_records(contig)
        if not region_assembled(directory, platform, assembly_provenance(files, platform)):
            raise RuntimeError("Assembly region has no valid completion record: %s" % directory)
        for i, record in enumerate(contigs):
            record.id = "c=%s:%s-%s_h=%s_i=%s_t=%s_/0/0_0" % (chrom, start, end, hap, i, len(contigs))
            record.description = ""
            seqs.append(record)
    fd, temporary = tempfile.mkstemp(dir=os.path.dirname(outfile) or ".", suffix=".fasta")
    try:
        with os.fdopen(fd, "w") as stream:
            SeqIO.write(seqs, stream, type_)
        if not os.path.isfile(outfile) or not filecmp.cmp(temporary, outfile, shallow=False):
            os.replace(temporary, outfile)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)

def combine_assembly_sequences(files,phased_blocks):
    combine_sequence(files,phased_blocks,files.assembly_fasta,"fasta")
    if any(block[0] == "igh" for block in phased_blocks):
        combine_sequence(files,phased_blocks,files.igh_assembly_fasta,"fasta","igh")
    #combine_sequence(files,phased_blocks,files.assembly_fastq,"fastq")

def run_assembly(
        rhesus,
        threads,
        mem,
        cluster,
        queue,
        walltime,
        outdir,
        data_dir
):
    files = FileManager(outdir,rhesus=rhesus,data_dir=data_dir)

    with open(files.input_args,'r') as fh:
        phasing_args = json.load(fh)
    sample = phasing_args["sample"]

    cpu = CpuManager(threads, mem, cluster, queue, walltime)
    #snps = Snps(files,cpu,sample)
    assembly_command_line = Assembly(files,cpu,sample)
    align_command_line = Align(files,cpu,sample)
    #snps_command_line = Snps(files,cpu,sample)

    phased_blocks = get_phased_blocks(files,files.phased_blocks)    
    
    #if not non_emptyfile(files.assembly_fastq): # CHANGED to GET CONSTANT
    assembly_scripts = get_assembly_scripts(files,cpu,phased_blocks)
    assembly_command_line.run_assembly_scripts(assembly_scripts)
    
    combine_assembly_sequences(files,phased_blocks)
    align_command_line.map_assembly()
    phase_assembly(files,sample)

    # merge_assembly(files,align_command_line,sample)
    # snps.phase_snvs_with_merged_seq()
    # snps.phased_blocks_from_merged_seq()
    # phase_merged_seqs(files,sample)

    
def main(args):
    run_assembly(**vars(args))
