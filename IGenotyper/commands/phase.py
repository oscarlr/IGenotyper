#!/bin/env python
from IGenotyper.files import FileManager
from IGenotyper.common.cpu import CpuManager


from IGenotyper.command_lines.snps import Snps
from IGenotyper.command_lines.reads import ReadManip
from IGenotyper.command_lines.alignments import Align
from IGenotyper.command_lines.plot import PlotTools

from IGenotyper.phasing.snps import generate_phased_snps
from IGenotyper.phasing.reads import phase_ccs
from IGenotyper.phasing.stats import phasing_stats

import os
import json
from shutil import copyfile
from IGenotyper.common.validation import require_usable_reads, validate_vcf
from IGenotyper.phasing.completion import provenance, phasing_complete, record_completion, receipt_path

def add_arguments(subparser):
    subparser.add_argument('--rhesus',default=False, action='store_true')
    subparser.add_argument('--sample',metavar='SAMPLE',default="sample",help='Name of sample')
    subparser.add_argument('--threads',metavar='THREADS',default=1,help='Number of threads')
    subparser.add_argument('--mem',metavar='MEM',default=20,help='Memory for cluster')
    subparser.add_argument('--cluster',default=False, action='store_true', help='Use cluster')
    subparser.add_argument('--queue', metavar='QUEUE',default="premium",help='Queue for cluster')
    subparser.add_argument('--walltime', metavar='WALLTIME',default=24,help='Walltime for cluster')
    subparser.add_argument('--tmp', metavar='TMP', default="tmp", help='Temporary folder')
    subparser.add_argument('--data-dir', help='Directory containing reference.fasta')
    subparser.add_argument('--input_vcf', metavar='VCF', help='Phased VCF file to phase reads')
    subparser.add_argument('bam', metavar='BAM', help='PacBio bam file')
    subparser.add_argument('outdir',metavar='OUTDIR',help='Directory for output')

def save_parameters(files,sample,input_vcf):
    paramaters = {
        "bam": files.input_bam,
        "sample": sample,
        "input_vcf": input_vcf,
        "tmp": files.tmp,
        "data_dir": files.data_directory
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
        tmp,
        rhesus,
        data_dir
):    
    files = FileManager(outdir,bam,tmp,rhesus,data_dir)

    expected = provenance(files, sample, input_vcf)
    if phasing_complete(files, expected):
        print('Phasing already completed and validated; reusing final outputs.')
        return
    receipt_path(files).unlink(missing_ok=True)

    cpu = CpuManager(threads,mem,cluster,queue,walltime)
    reads_command_line = ReadManip(files,cpu,sample)
    align_command_line = Align(files,cpu,sample)
    plot_command_line = PlotTools(files,cpu,sample)
    snps_command_line = Snps(files,cpu,sample)
    
    require_usable_reads(files.ccs_bam)
    reads_command_line.turn_ccs_reads_to_fastq()
    align_command_line.map_ccs_reads()
    if input_vcf is None:
        # Each command checks its success receipt; file size cannot prove completion.
        generate_phased_snps(files, cpu, sample)
    else:
        validate_vcf(input_vcf, sample)
        if os.path.abspath(input_vcf) != os.path.abspath(files.phased_snps_vcf):
            snps_command_line.run_command(
                "import_phased_vcf:v1 sample=%s" % sample, files.phased_snps_vcf,
                inputs=[input_vcf],
                validator=lambda paths: validate_vcf(paths[0], sample),
                action=lambda paths: copyfile(input_vcf, paths[0]))
    phase_ccs(files, sample)

    snps_command_line.phased_blocks_from_ccs_snps()
    phasing_stats(sample,files,plot_command_line,align_command_line)

    save_parameters(files,sample,input_vcf)
    record_completion(files, expected)
    #clean_up(files)
    
def main(args):
    run_phasing(**vars(args))
