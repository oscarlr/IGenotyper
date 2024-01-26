#!/bin/env python
import os
from lsf.lsf import Lsf
from IGenotyper.common.helper import non_emptyfile


class CommandLine:
    def __init__(self,files,cpu,sample):
        self.files = files
        self.cpu = cpu
        self.sample = sample
    
    def run_command(self,command,output_file):
        print "-----------------"
        print "Checking %s" % output_file        
        if not non_emptyfile(output_file):    
            print "\tRunning command... \n%s" % command            
            os.system(command)
        print "-----------------"
            
# def generate_ccs_reads(cpu,input_bam,ccs_bam,min_passes=2):
#     print "Generating CCS reads..."
#     args = [cpu.threads,min_passes,input_bam,ccs_bam]
#     command = ("ccs "
#                "--num-threads %s "
#                "--min-passes %s "               
#                "%s "
#                "%s > /dev/null 2>&1" % tuple(args))
#     output_file = "%s.pbi" % ccs_bam
#     run_command(command,output_file)

# def turn_ccs_reads_to_fastq(ccs_fastq_unedited,ccs_bam,ccs_fastq):
#     args = [ccs_fastq_unedited,ccs_bam,ccs_fastq_unedited,ccs_fastq]
#     command = ("bam2fastq "
#                "-o %s %s\n"
#                "zcat %s.fastq.gz | sed 's/ccs/0_8/g' > %s\n" % tuple(args))
#     run_command(command,files.ccs_fastq)

#     def map_reads_with_blasr(self,reads,prefix,ref,opts=""):
#         args = [reads,
#                 ref,
#                 ref,
#                 prefix,
#                 opts,
#                 self.cpu.threads]
#         command = ("blasr "
#                    "%s "
#                    "%s "
#                    "--sa %s.sa "
#                    "--out %s.sam "
#                    "--sam "
#                    "%s "
#                    "--nproc %s " % tuple(args))
#         output_file = "%s.sam" % prefix
#         self.run_command(command,output_file)

#     def sam_to_sorted_bam(self,prefix,sorted_bam):
#         sam = "%s.sam" % prefix
#         bam = "%s.bam" % prefix
#         args = [sam,bam,bam,sorted_bam,sorted_bam]
#         command = ("samtools view -Sbh %s > %s \n"
#                    "samtools sort %s -o %s > /dev/null 2>&1 \n"
#                    "samtools index %s" % tuple(args))
#         sorted_bam_bai = "%s.bai" % sorted_bam
#         self.run_command(command,sorted_bam_bai)

#     def map_subreads(self):
#         print "Mapping subreads..."
#         prefix = "%s/subreads_to_ref" % self.files.tmp
#         sorted_bam_tmp = "%s.sorted.bam" % prefix
#         if not non_emptyfile(self.files.subreads_to_ref):
#             if not non_emptyfile("%s.bai" % sorted_bam_tmp):
#                 self.map_reads_with_blasr(self.files.input_bam,prefix,self.files.ref)
#                 self.sam_to_sorted_bam(prefix,sorted_bam_tmp)
#             self.select_target_reads(sorted_bam_tmp,self.files.subreads_to_ref)

#     def map_assembly(self):
#         print "Mapping assembly..."
#         prefix = "%s/assembly_to_ref" % self.files.tmp
#         self.map_reads_with_blasr(self.files.assembly_fastq,prefix,self.files.ref)
#         self.sam_to_sorted_bam(prefix,self.files.assembly_to_ref)

#     def select_target_reads(self,bam_file,igh_bam_file):
#         ## add aim regions
#         args = [bam_file,self.files.target_regions,igh_bam_file,
#                 igh_bam_file]
#         command = ("samtools view -Sbh %s -L %s > %s \n"
#                    "samtools index %s" % tuple(args))
#         self.run_command(command,"%s.bai" % igh_bam_file)

def map_ccs_reads(self):
    print "Mapping CCS reads..."
    prefix = "%s/ccs_to_ref" % self.files.tmp
    sorted_bam_tmp = "%s.sorted.bam" % prefix
    if not non_emptyfile(self.files.ccs_to_ref):
        if not non_emptyfile("%s.bai" % sorted_bam_tmp):
            self.map_reads_with_blasr(self.files.ccs_fastq,prefix,self.files.ref)
            self.sam_to_sorted_bam(prefix,sorted_bam_tmp)
        self.select_target_reads(sorted_bam_tmp,self.files.ccs_to_ref)

# # wh_find_snv_candidates,wh_genotype,wh_phase

#     def wh_find_snv_candidates(self,sample_name):
#         args = [sample_name,
#                 self.files.ref,
#                 self.files.ccs_to_ref,
#                 self.files.snp_candidates,
#                 sample_name,
#                 self.files.ref,
#                 self.files.snvs_vcf,
#                 self.files.snp_candidates,
#                 self.files.ccs_to_ref]
#         command = ("CONDA_BASE=$(conda info --base) \n"
#                    "source ${CONDA_BASE}/etc/profile.d/conda.sh \n"
#                    "conda activate whatshap-latest \n"
#                    "whatshap find_snv_candidates "
#                    "--sample %s "
#                    "%s "
#                    "%s "
#                    "--pacbio "
#                    "-o %s #> /dev/null 2>&1 \n"
#                    "whatshap genotype "
#                    "--sample %s "
#                    "--ignore-read-groups "
#                    "--reference %s "
#                    "-o %s "
#                    "%s "
#                    "%s #> /dev/null 2>&1 " 
#                    "conda deactivate " % tuple(args))
#         self.run_command(command, self.files.snvs_vcf)

#     def phase_genotype_ccs_snvs(self,sample_name):
#         args = [sample_name,
#                 self.files.ref,
#                 self.files.phased_snvs_vcf,
#                 self.files.snvs_vcf,
#                 "%s" % self.files.ccs_to_ref]
#         self.phase_genotype_snvs(sample_name,args)
        
#     def phase_genotype_snvs(self,sample_name,args):
#         command = ("CONDA_BASE=$(conda info --base) \n"
#                    "source ${CONDA_BASE}/etc/profile.d/conda.sh \n"
#                    "conda activate whatshap-latest \n"
#                    "whatshap phase "
#                    "--sample %s "
#                    "--reference %s "
#                    "--ignore-read-groups "
#                    "--distrust-genotypes "
#                    "-o %s "
#                    "%s "
#                    "%s #> /dev/null 2>&1\n"
#                    "conda deactivate" % tuple(args))
#         self.run_command(command, args[2])

#     def phase_ccs_blocks(self,sample_name):
#         args = [sample_name,
#                 self.files.phased_blocks,
#                 self.files.phased_snvs_vcf]
#         self.phase_blocks(sample_name,args)
    
#     def phase_blocks(self,sample_name,args):
#         command = ("CONDA_BASE=$(conda info --base) \n"
#                    "source ${CONDA_BASE}/etc/profile.d/conda.sh \n"
#                    "conda activate whatshap-latest \n"
#                    "whatshap stats "
#                    "--sample %s "
#                     "--block-list %s "
#                    "%s #> /dev/null 2>&1\n"
#                    "conda deactivate" % tuple(args))
#         self.run_command(command,args[1])

#     def run_assembly_scripts(self,assembly_scripts):
#         if not self.cpu.cluster:
#             for script in assembly_scripts:
#                 command = "sh %s" % script
#                 os.system(command)
#         else:
#             hpc = Lsf()
#             for job in assembly_scripts:
#                 hpc.config(cpu=self.cpu.threads,walltime=self.cpu.walltime,
#                            memory=self.cpu.mem,queue=self.cpu.queue)
#                 hpc.submit("%s" % job)
#             hpc.wait()

#     def run_kalign(self,fastafn,clufn):
#         #| sed 's/Kalign/CLUSTAL/g' > ${dir}/msa.clu
#         # s=100
#         # e=0.85
#         # t=0.45
#         # m=0
#         # args = [s,e,t,m,fastafn,clufn]
#         #args = [fastafn,clufn]
#         command = "kalign -i %s -f clu -o %s" % (fastafn,clufn)
#         self.run_command(command,clufn)

    def blast_seq(self,fastafn,blast_out):
        args = [fastafn,fastafn,blast_out]
        command = ("blastn -query %s -subject %s "
                   "-outfmt \"6 length pident nident mismatch gapopen gaps qseqid qstart qend qlen sseqid sstart send slen sstrand\" "
                   "> %s " % tuple(args))
        self.run_command(command,blast_out)
