#!/bin/env python
import os
from lsf.lsf import Lsf
from helper import non_emptyfile

class CommandLine:
    def __init__(self,files,cpu):
        self.files = files
        self.cpu = cpu

    def run_command(self,command,output_file):
        if not non_emptyfile(output_file):
            os.system(command)

    def generate_ccs_reads(self):
        print "Generating CCS reads..."
        min_passes = 2
        args = [self.cpu.threads,
                min_passes,
                self.files.input_bam,
                self.files.ccs_bam]
        command = ("ccs "
                   "--num-threads %s "
                   "--min-passes %s "               
                   "%s "
                   "%s > /dev/null 2>&1" % tuple(args))
        output_file = "%s.pbi" % self.files.ccs_bam
        self.run_command(command,output_file)

    def turn_ccs_reads_to_fastq(self):
        args = [self.files.ccs_fastq_unedited,
                self.files.ccs_bam,
                self.files.ccs_fastq_unedited,
                self.files.ccs_fastq]
        command = ("bam2fastq "
                   "-o %s %s\n"
                   "zcat %s.fastq.gz | sed 's/ccs/0_8/g' > %s\n" % tuple(args))
        self.run_command(command,self.files.ccs_fastq)

    def map_reads_with_blasr(self,reads,prefix,ref,opts=""):
        args = [reads,
                ref,
                ref,
                prefix,
                opts,
                self.cpu.threads]
        command = ("blasr "
                   "%s "
                   "%s "
                   "--sa %s.sa "
                   "--out %s.sam "
                   "--sam "
                   "%s "
                   "--nproc %s " % tuple(args))
        output_file = "%s.sam" % prefix
        self.run_command(command,output_file)

    def sam_to_sorted_bam(self,prefix,sorted_bam):
        sam = "%s.sam" % prefix
        bam = "%s.bam" % prefix
        args = [sam,bam,bam,sorted_bam,sorted_bam]
        command = ("samtools view -Sbh %s > %s \n"
                   "samtools sort %s -o %s > /dev/null 2>&1 \n"
                   "samtools index %s" % tuple(args))
        sorted_bam_bai = "%s.bai" % sorted_bam
        self.run_command(command,sorted_bam_bai)

    def map_subreads(self):
        print "Mapping subreads..."
        prefix = "%s/subreads_to_ref" % self.files.tmp
        sorted_bam_tmp = "%s.sorted.bam" % prefix
        self.map_reads_with_blasr(self.files.input_bam,prefix,self.files.ref)
        self.sam_to_sorted_bam(prefix,sorted_bam_tmp)
        self.select_target_reads(sorted_bam_tmp,self.files.subreads_to_ref)

    def map_assembly(self):
        print "Mapping assembly..."
        prefix = "%s/assembly_to_ref" % self.files.tmp
        self.map_reads_with_blasr(self.files.assembly_fastq,prefix,self.files.ref)
        self.sam_to_sorted_bam(prefix,self.files.assembly_to_ref)

    def select_target_reads(self,bam_file,igh_bam_file):
        ## add aim regions
        args = [bam_file,self.files.target_regions,igh_bam_file,
                igh_bam_file]
        command = ("samtools view -Sbh %s -L %s > %s \n"
                   "samtools index %s" % tuple(args))
        self.run_command(command,"%s.bai" % igh_bam_file)

    def map_ccs_reads(self):
        print "Mapping CCS reads..."
        prefix = "%s/ccs_to_ref" % self.files.tmp
        sorted_bam_tmp = "%s.sorted.bam" % prefix
        self.map_reads_with_blasr(self.files.ccs_fastq,prefix,self.files.ref)
        self.sam_to_sorted_bam(prefix,sorted_bam_tmp)
        self.select_target_reads(sorted_bam_tmp,self.files.ccs_to_ref)

    def genotype_snvs_from_ccs(self,sample_name):
        args = [sample_name,
                self.files.ref,
                self.files.ccs_to_ref,
                self.files.snp_candidates,
                sample_name,
                self.files.ref,
                self.files.snvs_vcf,
                self.files.snp_candidates,
                self.files.ccs_to_ref]
        command = ("CONDA_BASE=$(conda info --base) \n"
                   "source ${CONDA_BASE}/etc/profile.d/conda.sh \n"
                   "conda activate whatshap-latest \n"
                   "whatshap find_snv_candidates "
                   "--sample %s "
                   "%s "
                   "%s "
                   "--pacbio "
                   "-o %s #> /dev/null 2>&1 \n"
                   "whatshap genotype "
                   "--sample %s "
                   "--ignore-read-groups "
                   "--reference %s "
                   "-o %s "
                   "%s "
                   "%s #> /dev/null 2>&1 " 
                   "conda deactivate " % tuple(args))
        self.run_command(command, self.files.snvs_vcf)

    def phase_genotype_snvs(self,sample_name):
        args = [sample_name,
                self.files.ref,
                self.files.phased_snvs_vcf,
                self.files.snvs_vcf,
                self.files.ccs_to_ref]
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
        self.run_command(command, self.files.phased_snvs_vcf)

    def phase_blocks(self,sample_name):
        args = [sample_name,
                self.files.phased_blocks,
                self.files.phased_snvs_vcf]
        command = ("CONDA_BASE=$(conda info --base) \n"
                   "source ${CONDA_BASE}/etc/profile.d/conda.sh \n"
                   "conda activate whatshap-latest \n"
                   "whatshap stats "
                   "--sample %s "
                    "--block-list %s "
                   "%s #> /dev/null 2>&1\n"
                   "conda deactivate" % tuple(args))
        self.run_command(command, self.files.phased_blocks)

    def run_assembly_scripts(self,assembly_scripts):
        if not self.cpu.cluster:
            for script in assembly_scripts:
                command = "sh %s" % script
                os.system(command)
        else:
            hpc = Lsf()
            for job in assembly_scripts:
                hpc.config(cpu=self.cpu.threads,walltime=self.cpu.walltime,
                           memory=self.cpu.mem,queue=self.cpu.queue)
                hpc.submit("%s" % job)
            hpc.wait()