#!/bin/env python

from IGenotyper.command_line.alignments import Snps

def generate_phased_snps(files,cpu,sample):
    snps_command_line = Snps(files,cpu,sample)
    snps_command_line.snp_candidates_from_ccs()
    snps_command_line.snp_genotypes_from_ccs()
    snps_command_line.phase_ccs_snvs()
