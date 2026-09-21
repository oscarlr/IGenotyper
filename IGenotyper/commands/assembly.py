#!/bin/env python
import os
import json
import tempfile
import filecmp
import uuid
from Bio import SeqIO

from IGenotyper.files import FileManager

from IGenotyper.common.cpu import CpuManager
from IGenotyper.common.helper import get_phased_blocks

from IGenotyper.command_lines.assembly import Assembly
from IGenotyper.command_lines.alignments import Align
#from IGenotyper.command_lines.snps import Snps

from IGenotyper.phasing.reads import phase_assembly

from IGenotyper.assembly.scripts import get_assembly_scripts, assembly_contigs, assembly_provenance, region_status, region_provenance, region_failure
from IGenotyper.assembly.workflow import select_assembly_workflow
from IGenotyper.common.validation import fasta_records
from IGenotyper.assembly.scripts import assembly_bam, validate_assembly_inputs
from IGenotyper.assembly.coverage import measure_ig_coverage
from IGenotyper.assembly.status import write_sample_status
#from IGenotyper.assembly.merge_assembly import merge_assembly

def add_arguments(subparser):
    subparser.add_argument('--rhesus',default=False, action='store_true')
    subparser.add_argument('--threads', metavar='THREADS', default=1, help='Number of threads')
    subparser.add_argument('--mem', metavar='MEM', default=8, help='Memory for cluster')
    subparser.add_argument('--cluster', default=False, action='store_true', help='Use cluster')
    subparser.add_argument('--queue', metavar='QUEUE', default="premium", help='Queue for cluster')
    subparser.add_argument('--walltime', metavar='WALLTIME', default=2, help='Walltime for cluster')
    subparser.add_argument('--data-dir', help='Directory containing reference.fasta')
    subparser.add_argument('--coverage-bed', help='Reference-matched IG loci BED for the 20x mean-coverage assembly gate')
    subparser.add_argument('outdir',metavar='OUTDIR',help='Directory for output')

def combine_sequence(files,phased_blocks,outfile,type_,chrom_select=None,workflow=None,allow_empty=False):
    seqs = []
    workflow = workflow or select_assembly_workflow(files.input_bam)
    provenance = assembly_provenance(files, workflow)
    selected = [block for block in phased_blocks if chrom_select is None or block[0] == chrom_select]
    if not selected:
        raise RuntimeError("No assembly regions selected for %s" % outfile)
    for chrom, start, end, hap in selected:
        directory = "%s/assembly/%s/%s_%s/%s" % (files.tmp, chrom, start, end, hap)
        status = region_status(directory, workflow, region_provenance(provenance, chrom, start, end, hap))
        if status in ("skipped_no_coverage", "skipped_no_contigs", "failed_canu"):
            continue
        contig = assembly_contigs(directory, workflow, type_)
        contigs = fasta_records(contig)
        if status != "assembled":
            raise RuntimeError("Assembly region has no valid completion record: %s" % directory)
        for i, record in enumerate(contigs):
            record.id = "c=%s:%s-%s_h=%s_i=%s_t=%s_/0/0_0" % (chrom, start, end, hap, i, len(contigs))
            record.description = ""
            seqs.append(record)
    if not seqs:
        if allow_empty:
            # A checked empty result must not leave a stale FASTA at its path.
            if os.path.exists(outfile):
                os.replace(outfile, outfile + '.previous-' + uuid.uuid4().hex)
            return 0
        raise RuntimeError("No valid contigs overall: all selected assembly regions were skipped (%s)" % outfile)
    fd, temporary = tempfile.mkstemp(dir=os.path.dirname(outfile) or ".", suffix=".fasta")
    try:
        with os.fdopen(fd, "w") as stream:
            SeqIO.write(seqs, stream, type_)
        if not os.path.isfile(outfile) or not filecmp.cmp(temporary, outfile, shallow=False):
            os.replace(temporary, outfile)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return len(seqs)

def combine_assembly_sequences(files,phased_blocks,workflow=None):
    workflow = workflow or select_assembly_workflow(files.input_bam)
    count = combine_sequence(files,phased_blocks,files.assembly_fasta,"fasta",workflow=workflow,allow_empty=True)
    if any(block[0] == "igh" for block in phased_blocks):
        combine_sequence(files,phased_blocks,files.igh_assembly_fasta,"fasta","igh",workflow=workflow,allow_empty=True)
    return count

def run_assembly(
        rhesus,
        threads,
        mem,
        cluster,
        queue,
        walltime,
        outdir,
        data_dir,
        coverage_bed=None
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

    # Gate before region planning or any assembly command. Phasing is read-only here.
    workflow = select_assembly_workflow(files.input_bam)
    provenance = assembly_provenance(files, workflow)
    coverage = None
    failed_regions = []
    write_sample_status(files, sample, 'running', provenance, coverage)
    try:
        coverage = measure_ig_coverage(files, assembly_bam(files, workflow), coverage_bed)
        validate_assembly_inputs(provenance)
        if coverage['below_threshold']:
            return write_sample_status(files, sample, 'insufficient_coverage', provenance, coverage)
        phased_blocks = get_phased_blocks(files,files.phased_blocks)
        assembly_scripts = get_assembly_scripts(files,cpu,phased_blocks,workflow)
        assembly_command_line.run_assembly_scripts(assembly_scripts)
        for chrom, start, end, hap in phased_blocks:
            directory = "%s/assembly/%s/%s_%s/%s" % (files.tmp, chrom, start, end, hap)
            failure = region_failure(directory, region_provenance(provenance, chrom, start, end, hap))
            if failure is not None:
                failed_regions.append(failure)
        if combine_assembly_sequences(files,phased_blocks,workflow) == 0:
            if failed_regions:
                raise RuntimeError('No valid contigs recovered; Canu failed in %s region(s). See assembly_status.json for logs.' % len(failed_regions))
            return write_sample_status(files, sample, 'no_contigs', provenance, coverage)
        align_command_line.map_assembly()
        phase_assembly(files,sample)
    except Exception as error:
        write_sample_status(files, sample, 'failed', provenance, coverage, error, failed_regions)
        raise
    status = 'completed_with_failures' if failed_regions else 'completed'
    return write_sample_status(files, sample, status, provenance, coverage, failed_regions=failed_regions)

    # merge_assembly(files,align_command_line,sample)
    # snps.phase_snvs_with_merged_seq()
    # snps.phased_blocks_from_merged_seq()
    # phase_merged_seqs(files,sample)

    
def main(args):
    run_assembly(**vars(args))
