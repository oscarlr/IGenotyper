#!/bin/env python
import os
from lsf.lsf import Lsf
from IGenotyper.helper import non_emptyfile

from IGenotyper.command_lines.clt import CommandLine

class Align(CommandLine):
    def __init__(self,files,cpu,sample):
        CommandLine.__init__(self,files,cpu,sample)

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
        if not non_emptyfile(self.files.subreads_to_ref):
            if not non_emptyfile("%s.bai" % sorted_bam_tmp):
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
        if not non_emptyfile(self.files.ccs_to_ref):
            if not non_emptyfile("%s.bai" % sorted_bam_tmp):
                self.map_reads_with_blasr(self.files.ccs_fastq,prefix,self.files.ref)
                self.sam_to_sorted_bam(prefix,sorted_bam_tmp)
            self.select_target_reads(sorted_bam_tmp,self.files.ccs_to_ref)

