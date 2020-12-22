#!/bin/env

from IGenotyper.common.command_line import wh_find_snv_candidates,wh_genotype,wh_phase

def detect_snps_from_reads(bam,ref,sample,vcf_no_genotype,vcf_genotype):
    args = [sample,ref,bam,vcf_no_genotype]
    wh_find_snv_candidates(args)

    args = [sample,ref,vcf_genotype,vcf_no_genotype,bam]
    wh_genotype(args)

def phase_snps_from_reads(sample,ref,phased_vcf,vcf,bam):
    args = [sample,ref,phased_vcf,vcf,bam]
    wh_phase(args)
