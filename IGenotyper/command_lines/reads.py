#!/usr/bin/env python3
from shlex import quote
from IGenotyper.common.validation import require_usable_reads

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
        self.run_command(command, [self.files.ccs_bam, output_file])
        command = "samtools index %s %s" % (quote(self.files.ccs_bam), quote(self.files.ccs_bam + ".bai"))
        output_file = "%s.bai" % self.files.ccs_bam
        self.run_command(command, output_file, inputs=[self.files.ccs_bam])

    def turn_ccs_reads_to_fastq(self):
        require_usable_reads(self.files.ccs_bam)
        command = (
            "samtools fasta -n -@ %s %s > %s"
            % (
                int(self.cpu.threads),
                quote(self.files.ccs_bam),
                quote(self.files.ccs_fastq),
            )
        )
        self.run_command(command, self.files.ccs_fastq, inputs=[self.files.ccs_bam])
