#!/bin/env python
from IGenotyper.cpu import CpuManager
from IGenotyper.files import FileManager
from IGenotyper.clt import CommandLine
from IGenotyper.vcffn import read_in_phased_vcf
from IGenotyper.bam import create_phased_bam_header

import pysam

def add_arguments(subparser):
    subparser.add_argument('--sample',metavar='SAMPLE',default="sample",help='Name of sample')
    subparser.add_argument('--threads',metavar='THREADS',default=1,help='Number of threads')
    subparser.add_argument('--mem',metavar='MEM',default=20,help='Memory for cluster')
    subparser.add_argument('--cluster',default=False, action='store_true', help='Use cluster')
    subparser.add_argument('--queue', metavar='QUEUE',default="premium",help='Queue for cluster')
    subparser.add_argument('--walltime', metavar='WALLTIME',default=24,help='Walltime for cluster')
    subparser.add_argument('--tmp', metavar='TMP', default="tmp", help='Temporary folder')
    subparser.add_argument('bam', metavar='BAM', help='PacBio bam file')
    subparser.add_argument('outdir',metavar='OUTDIR',help='Directory for output')

def pre_phase_processing(command_line_tools,sample):
    command_line_tools.generate_ccs_reads()
    command_line_tools.turn_ccs_reads_to_fastq()
    command_line_tools.map_subreads()
    command_line_tools.map_ccs_reads()
    command_line_tools.genotype_snvs_from_ccs(sample)
    command_line_tools.phase_genotype_snvs(sample)

def create_tag(hap):
    haptag = ("RG", str(hap), "Z")
    return haptag

def calculate_tag(basetups,hap_dict):
    '''
    Returns a 0, 1 or -1 as the value of a read
    '''
    unphased_tag = 0
    haplotype1_tag = 1
    haplotype2_tag = 2
    b1_bases = 0.0
    b0_bases = 0.0
    thres = .9
    for basetup in basetups:
        rpos, b, qual = basetup
        b0, b1 = hap_dict[rpos]
        if b == b1 or b == b0: # check to make sure the base is called
            if b == b0:
                b0_bases += 1.0
            else:
                b1_bases += 1.0
    if (b1_bases + b0_bases) == 0:
        return create_tag(unphased_tag)
    b0_frac = b0_bases/(b1_bases + b0_bases)
    b1_frac = b1_bases/(b1_bases + b0_bases)
    if b0_frac >= thres:
        return create_tag(haplotype1_tag)
    elif b1_frac >= thres:
        return create_tag(haplotype2_tag)
    else:
        return create_tag(unphased_tag)

def phase_read(read,snps,chrom):
    unphased_tag = 0
    haplotype1_tag = 1
    haplotype2_tag = 2
    base_tuples = []
    if chrom in snps:
        for read_pos, ref_pos in read.get_aligned_pairs():
            if ref_pos is None:
                continue
            if ref_pos in snps[chrom]:
                if read_pos is not None:
                    qual = 8
                    if read.query_qualities != None:
                        qual = read.query_qualities[read_pos]
                    base = read.query_sequence[read_pos]
                    base_tuples.append((ref_pos, base, qual))
    if len(base_tuples) > 0:
        read_group_tag = calculate_tag(base_tuples, snps[chrom])
    else:
        read_group_tag = create_tag(unphased_tag)
    read_tags = read.get_tags()
    tags_to_add = []
    for tag in read_tags:
        if tag[0] != "RG":
            tags_to_add.append(tag)
    tags_to_add.append(read_group_tag)
    read.set_tags(tags_to_add)
    return read

def phase_reads(bam,outbam,vcf):
    phase_snps = vcf.phased_variants()
    unphased_bam = pysam.AlignmentFile(bam, 'rb')
    phased_bam_header = create_phased_bam_header(bam)
    phased_bam = pysam.AlignmentFile(outbam,'wb',header=phased_bam_header)
    for read in unphased_bam.fetch():
        if read.is_unmapped:
            continue
        chrom = unphased_bam.get_reference_name(read.reference_id)
        tagged_read = phase_read(read,phase_snps,chrom)
        phased_bam.write(tagged_read)
    unphased_bam.close()
    phased_bam.close()

def phase_sequencing_data(files,sample):
    vcf = read_in_phased_vcf(files.phased_snvs_vcf,sample)
    phase_reads(files.ccs_to_ref,files.ccs_to_ref_phased,vcf)
    phase_reads(files.subreads_to_ref,files.subreads_to_ref_phased,vcf)

def run_phasing(
        bam,
        outdir,
        sample,
        threads,
        mem,
        cluster,
        queue,
        walltime,
        tmp
):
    files = FileManager(outdir,bam,tmp)
    cpu = CpuManager(threads,mem,cluster,queue,walltime)
    command_line_tools = CommandLine(files,cpu)
    pre_phase_processing(command_line_tools,sample)
    phase_sequencing_data(files,sample)

def main(args):
    run_phasing(**vars(args))
