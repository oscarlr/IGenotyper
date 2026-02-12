#!/bin/env python
import os
from lsf.lsf import Lsf
from IGenotyper.common.helper import non_emptyfile

from IGenotyper.command_lines.clt import CommandLine

class ReadManip(CommandLine):
    def __init__(self,files,cpu,sample):
        CommandLine.__init__(self,files,cpu,sample)

    def generate_ccs_reads(self):
        #print "Generating CCS reads..."
        min_passes = 2
        args = [self.cpu.threads,
                min_passes,
                self.files.input_bam,
                self.files.ccs_bam]
        command = ("ccs "
                   "--num-threads %s "
                   "--min-passes %s "               
                   "%s "
                   "%s #> /dev/null 2>&1" % tuple(args))
        output_file = "%s.pbi" % self.files.ccs_bam
        self.run_command(command,output_file)
        command = "samtools index %s" % self.files.ccs_bam
        output_file = "%s.bai" % self.files.ccs_bam
        self.run_command(command,output_file)

    def turn_ccs_reads_to_fastq(self):
        args = [self.files.ccs_bam,
                self.files.ccs_fastq]
        command = (
            "samtools fastq %s | "
            "sed 's/ccs/0_8/g' | "
            "sed 's/\\/fwd//g' | "
            "sed 's/\\/rev//g' "
            "> %s\n" % tuple(args)
        )
        self.run_command(command, self.files.ccs_fastq)

    # def turn_ccs_reads_to_fastq(self):
    #     args = [self.files.ccs_fastq_unedited,
    #             self.files.ccs_bam,
    #             self.files.ccs_fastq_unedited,
    #             self.files.ccs_fastq]
    #     command = ("bam2fasta "
    #                "-o %s %s\n"
    #                "zcat %s.fasta.gz | sed 's/ccs/0_8/g' | sed 's/\/fwd//g' | sed 's/\/rev//g' > %s\n" % tuple(args))
    #     self.run_command(command,self.files.ccs_fastq)
