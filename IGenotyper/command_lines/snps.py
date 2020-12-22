#!/bin/env python
import os
from lsf.lsf import Lsf

from IGenotyper.command_lines.clt import CommandLine

class Snps(CommandLine):
    def __init__(self,files,cpu,sample):
        CommandLine.__init__(self,files,cpu,sample)

    def snp_candidates(self,bam,snp_candidates):
        args = [self.sample,
                self.files.ref,
                bam,
                snp_candidates]
        command = ("CONDA_BASE=$(conda info --base) \n"
                   "source ${CONDA_BASE}/etc/profile.d/conda.sh \n"
                   "conda activate whatshap-latest \n"
                   "whatshap find_snv_candidates "
                   "--sample %s "
                   "%s "
                   "%s "
                   "--pacbio "
                   "-o %s #> /dev/null 2>&1 "
                   "conda deactivate " % tuple(args))
        self.run_command(command,snp_candidates)
        
    def snp_candidates_from_ccs(self):
        self.snp_candidates(self.files.ccs_to_ref,self.files.snp_candidates)

    def snp_genotypes(self,bam,snp_candidates,snps_vcf):
        args = [self.sample,
                self.files.ref,
                snvs_vcf,
                snp_candidates,
                bam]
        command = ("CONDA_BASE=$(conda info --base) \n"
                   "source ${CONDA_BASE}/etc/profile.d/conda.sh \n"
                   "conda activate whatshap-latest \n"
                   "whatshap genotype "
                   "--sample %s "
                   "--ignore-read-groups "
                   "--reference %s "
                   "-o %s "
                   "%s "
                   "%s #> /dev/null 2>&1 " 
                   "conda deactivate " % tuple(args))
        self.run_command(command,snvs_vcf)

    def snp_genotypes_from_ccs(self):
        self.snp_genotypes(self,self.files.ccs_to_ref,self.files.snp_candidates,self.files.snps_vcf)
        
    def phase_snvs(self,phased_snvs_vcf,snvs_vcf,bams):
        bams = " ".join(bams)
        args = [self.sample,
                self.files.ref,
                phased_snvs_vcf,
                snvs_vcf,
                bams]
        command = ("CONDA_BASE=$(conda info --base) \n"
                   "source ${CONDA_BASE}/etc/profile.d/conda.sh \n"
                   "conda activate whatshap-latest \n"
                   "whatshap phase "
                   "--sample %s "
                   "--reference %s "
                   "--ignore-read-groups "
                   "--distrust-genotypes "
                   "-o %s "
                   "%s "
                   "%s #> /dev/null 2>&1\n"
                   "conda deactivate" % tuple(args))
        self.run_command(command,phased_snvs_vcf)

    def phase_ccs_snvs(self):
        phase_snvs(self,self.files.phased_snvs_vcf,self.files.snvs_vcf,[self.files.ccs_to_ref])
        
    def phased_blocks(self,phased_blocks,phased_snvs_vcf):
        args = [self.sample,
                phased_blocks,
                phased_snvs_vcf]
        command = ("CONDA_BASE=$(conda info --base) \n"
                   "source ${CONDA_BASE}/etc/profile.d/conda.sh \n"
                   "conda activate whatshap-latest \n"
                   "whatshap stats "
                   "--sample %s "
                    "--block-list %s "
                   "%s #> /dev/null 2>&1\n"
                   "conda deactivate" % tuple(args))
        self.run_command(command,phased_blocks)

    def phased_blocks_from_ccs_snps(self):
        self.phased_blocks(self.files.phased_blocks,self.files.phased_snvs_vcf)    
