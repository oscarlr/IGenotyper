#!/bin/env python
import os
from lsf.lsf import Lsf
from IGenotyper.common.helper import non_emptyfile

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
                   "--noSplitSubreads "
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

    def create_igh_ref(self):
        print "Creating igh reference..."
        igh_ref = "%s/igh_ref.fasta" % self.files.tmp        
        args = [self.files.ref,igh_ref,
                igh_ref,
                igh_ref]
        command = ("samtools faidx %s igh > %s \n"
                   "samtools faidx %s \n"
                   "sawriter %s \n" % tuple(args))
        self.run_command(command,"%s.sa" % igh_ref)
        return igh_ref

    def map_igh_assembly(self):
        print "Mapping igh assembly..."
        igh_ref = self.create_igh_ref()
        prefix = "%s/igh_assembly_to_ref" % self.files.tmp
        self.map_reads_with_blasr(self.files.igh_assembly_fasta,prefix,igh_ref,"--insertion 0 --deletion 0 --minMatch 15 --maxMatch 30")
        self.sam_to_sorted_bam(prefix,self.files.igh_assembly_to_ref)

        prefix = "%s/igh_assembly_to_ref_subs" % self.files.tmp
        self.map_reads_with_blasr(self.files.igh_assembly_fasta,prefix,igh_ref,"--insertion 16 --deletion 20 --minMatch 15 --maxMatch 30")
        self.sam_to_sorted_bam(prefix,self.files.igh_assembly_to_ref_subs)
        
    def map_assembly(self):
        print "Mapping assembly..."
        prefix = "%s/assembly_to_ref" % self.files.tmp
        self.map_reads_with_blasr(self.files.assembly_fasta,prefix,self.files.ref,"--insertion 16 --deletion 20 --minMatch 15 --maxMatch 30")
        #self.map_reads_with_blasr(self.files.assembly_fasta,prefix,self.files.ref,"--insertion 0 --deletion 0 --minMatch 15 --maxMatch 30")
        self.sam_to_sorted_bam(prefix,self.files.assembly_to_ref)
        #self.map_igh_assembly()

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
        #sorted_bam_tmp = "%s.sorted.bam" % prefix
        if not non_emptyfile(self.files.ccs_to_ref):
            if not non_emptyfile("%s.bai" % self.files.ccs_to_ref):
                #self.map_reads_with_blasr(self.files.ccs_fastq,prefix,self.files.ref)
                self.map_reads_with_minimap2(self.files.ccs_fastq,prefix,self.files.ref)
                self.sam_to_sorted_bam(prefix,self.files.ccs_to_ref)
                #self.sam_to_sorted_bam(prefix,sorted_bam_tmp)
            #self.select_target_reads(sorted_bam_tmp,self.files.ccs_to_ref)

    def blast_seq(self,fastafn,blast_out):
        args = [fastafn,fastafn,blast_out]
        command = ("blastn -query %s -subject %s "
                   "-outfmt \"6 length pident nident mismatch gapopen gaps qseqid qstart qend qlen sseqid sstart send slen sstrand\" "
                   "> %s " % tuple(args))
        self.run_command(command,blast_out)

    def map_merged_assembly(self):
        print "Mapping merged assembly..."
        prefix = "%s/merged_assembly_to_ref" % self.files.tmp
        #opt="--insertion 0 --deletion 0 --minMatch 35 --maxMatch 50 --scoreMatrix \"-100 50 50 50 50 50 -100 50 50 50 50 50 -100 50 50 50 50 50 -100 50 50 50 50 50 -100\""
        opt="--minMatch 35 --maxMatch 50"
        self.map_reads_with_blasr(self.files.merged_assembly,prefix,self.files.ref,opt)
        self.sam_to_sorted_bam(prefix,self.files.merged_assembly_to_ref)

    def bam_to_bigwig(self,bam,bigwig):
        args =  [bam,bigwig]
        command = ("CONDA_BASE=$(conda info --base) \n"
                   "source ${CONDA_BASE}/etc/profile.d/conda.sh \n"
                   "conda activate pygenometracks \n"
                   "bamCoverage -b %s -o %s " % tuple(args))
        self.run_command(command,bigwig)

    def select_hap_sequence(self,bam,hap,outbam):
        args = [bam,hap,outbam,
                outbam]
        command= ("samtools view -Sbh -F 3884 %s -r %s > %s \n"
                  "samtools index %s " % tuple(args))
        self.run_command(command,"%s.bai" % outbam)
        
    def hap_bam_to_bigwig(self,bam,hap,bigwig):
        outbam = "%s/%s.bam" % (self.files.tmp,hap)
        self.select_hap_sequence(bam,hap,outbam)
        self.bam_to_bigwig(outbam,bigwig)
        
    def primary_alignments(self,inbam,outbam):
        args = [inbam,outbam,
                outbam]
        command = ("samtools view -Sbh -F 3884 %s > %s \n"
                   "samtools index %s " % tuple(args))
        self.run_command(command,"%s.bai" % outbam)

    def map_reads_with_minimap2(self, reads, prefix, ref, opts=""):
        """
        Map reads/contigs to a reference using minimap2.
        Default preset: HiFi reads (-x map-hifi)
        """

        preset = opts.strip()
        if preset == "":
            preset = "-x map-hifi"

        args = [self.cpu.threads,
                preset,
                ref,
                reads,
                prefix]

        command = (
            "minimap2 "
            "-t %s "
            "-a "
            "%s "
            "%s "
            "%s "
            "> %s.sam " % tuple(args)
        )

        output_file = "%s.sam" % prefix
        self.run_command(command, output_file)
