#!/usr/bin/env python3
from shlex import quote

from IGenotyper.command_lines.clt import CommandLine, non_emptyfile


class Align(CommandLine):
    def map_reads_with_minimap2(self, reads, sorted_bam, ref, preset="map-hifi"):
        """Map reads and stream directly into a sorted, indexed BAM."""
        if reads.lower().endswith((".bam", ".cram")):
            query = "-"
            input_command = "samtools fasta -@ %s %s | " % (
                int(self.cpu.threads),
                quote(reads),
            )
        else:
            query = quote(reads)
            input_command = ""
        minimap2_ref = ref
        if ref == self.files.ref:
            minimap2_ref = getattr(self.files, "minimap2_ref", ref)
        command = (
            "set -o pipefail; %sminimap2 -t %s -a -x %s %s %s | "
            "samtools sort -@ %s -o %s - && samtools index %s"
        ) % (
            input_command,
            int(self.cpu.threads),
            quote(preset),
            quote(minimap2_ref),
            query,
            int(self.cpu.threads),
            quote(sorted_bam),
            quote(sorted_bam),
        )
        self.run_command(command, "%s.bai" % sorted_bam)

    def sam_to_sorted_bam(self, prefix, sorted_bam):
        sam = "%s.sam" % prefix
        command = "samtools sort -@ %s -o %s %s && samtools index %s" % (
            int(self.cpu.threads),
            quote(sorted_bam),
            quote(sam),
            quote(sorted_bam),
        )
        self.run_command(command, "%s.bai" % sorted_bam)

    def map_subreads(self):
        print("Mapping subreads...")
        prefix = "%s/subreads_to_ref" % self.files.tmp
        sorted_bam_tmp = "%s.sorted.bam" % prefix
        if not non_emptyfile(self.files.subreads_to_ref):
            if not non_emptyfile("%s.bai" % sorted_bam_tmp):
                self.map_reads_with_minimap2(
                    self.files.input_bam, sorted_bam_tmp, self.files.ref, "map-pb"
                )
            self.select_target_reads(sorted_bam_tmp, self.files.subreads_to_ref)

    def create_igh_ref(self):
        print("Creating IGH reference...")
        igh_ref = "%s/igh_ref.fasta" % self.files.tmp
        command = "samtools faidx %s igh > %s && samtools faidx %s" % (
            quote(self.files.ref),
            quote(igh_ref),
            quote(igh_ref),
        )
        self.run_command(command, "%s.fai" % igh_ref)
        return igh_ref

    def map_igh_assembly(self):
        print("Mapping IGH assembly...")
        igh_ref = self.create_igh_ref()
        self.map_reads_with_minimap2(
            self.files.igh_assembly_fasta,
            self.files.igh_assembly_to_ref_subs,
            igh_ref,
            "asm20",
        )

    def map_assembly(self):
        print("Mapping assembly...")
        self.map_reads_with_minimap2(
            self.files.assembly_fasta,
            self.files.assembly_to_ref,
            self.files.ref,
            "asm20",
        )

    def select_target_reads(self, bam_file, target_bam_file):
        command = "samtools view -bh %s -L %s -o %s && samtools index %s" % (
            quote(bam_file),
            quote(self.files.target_regions),
            quote(target_bam_file),
            quote(target_bam_file),
        )
        self.run_command(command, "%s.bai" % target_bam_file)

    def map_ccs_reads(self):
        print("Mapping CCS reads...")
        prefix = "%s/ccs_to_ref" % self.files.tmp
        if not non_emptyfile(self.files.ccs_to_ref):
            if not non_emptyfile("%s.bai" % self.files.ccs_to_ref):
                self.map_reads_with_minimap2(
                    self.files.ccs_fastq,
                    self.files.ccs_to_ref,
                    self.files.ref,
                    "map-hifi",
                )

    def blast_seq(self, fastafn, blast_out):
        command = (
            "blastn -query %s -subject %s "
            "-outfmt '6 length pident nident mismatch gapopen gaps qseqid "
            "qstart qend qlen sseqid sstart send slen sstrand' > %s"
            % (quote(fastafn), quote(fastafn), quote(blast_out))
        )
        self.run_command(command, blast_out)

    def map_merged_assembly(self):
        print("Mapping merged assembly...")
        self.map_reads_with_minimap2(
            self.files.merged_assembly,
            self.files.merged_assembly_to_ref,
            self.files.ref,
            "asm20",
        )

    def bam_to_bigwig(self, bam, bigwig):
        command = "bamCoverage -b %s -o %s" % (quote(bam), quote(bigwig))
        self.run_command(command, bigwig)

    def select_hap_sequence(self, bam, hap, outbam):
        command = "samtools view -bh -F 3884 -r %s -o %s %s && samtools index %s" % (
            quote(hap),
            quote(outbam),
            quote(bam),
            quote(outbam),
        )
        self.run_command(command, "%s.bai" % outbam)

    def hap_bam_to_bigwig(self, bam, hap, bigwig):
        outbam = "%s/%s.bam" % (self.files.tmp, hap)
        self.select_hap_sequence(bam, hap, outbam)
        self.bam_to_bigwig(outbam, bigwig)

    def primary_alignments(self, inbam, outbam):
        command = "samtools view -bh -F 3884 -o %s %s && samtools index %s" % (
            quote(outbam),
            quote(inbam),
            quote(outbam),
        )
        self.run_command(command, "%s.bai" % outbam)
