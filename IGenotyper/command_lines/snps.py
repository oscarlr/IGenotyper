#!/usr/bin/env python3
from shlex import quote

from IGenotyper.command_lines.clt import CommandLine


class Snps(CommandLine):
    """WhatsHap 2.x variant calling, genotyping, and phasing commands."""

    def snp_candidates(self, bam, snp_candidates):
        command = (
            "whatshap find_snv_candidates --sample %s --pacbio -o %s %s %s"
            % (
                quote(self.sample),
                quote(snp_candidates),
                quote(self.files.ref),
                quote(bam),
            )
        )
        self.run_command(command, snp_candidates)

    def snp_candidates_from_ccs(self):
        self.snp_candidates(self.files.ccs_to_ref, self.files.snp_candidates)

    def snp_genotypes(self, bam, snp_candidates, snps_vcf):
        command = (
            "whatshap genotype --sample %s --ignore-read-groups "
            "--reference %s -o %s %s %s"
            % (
                quote(self.sample),
                quote(self.files.ref),
                quote(snps_vcf),
                quote(snp_candidates),
                quote(bam),
            )
        )
        self.run_command(command, snps_vcf)

    def snp_genotypes_from_ccs(self):
        self.snp_genotypes(
            self.files.ccs_to_ref, self.files.snp_candidates, self.files.snps_vcf
        )

    def phase_snps(self, phased_snps_vcf, snps_vcf, bams):
        command = (
            "whatshap phase --sample %s --reference %s --ignore-read-groups "
            "--distrust-genotypes -o %s %s %s"
            % (
                quote(self.sample),
                quote(self.files.ref),
                quote(phased_snps_vcf),
                quote(snps_vcf),
                " ".join(quote(bam) for bam in bams),
            )
        )
        self.run_command(command, phased_snps_vcf)

    def phase_ccs_snvs(self):
        self.phase_snps(
            self.files.phased_snps_vcf,
            self.files.snps_vcf,
            [self.files.ccs_to_ref],
        )

    def phased_blocks(self, phased_blocks, chr_lengths, phased_snps_vcf):
        command = (
            "whatshap stats --sample %s --block-list %s --chr-lengths %s %s"
            % (
                quote(self.sample),
                quote(phased_blocks),
                quote(chr_lengths),
                quote(phased_snps_vcf),
            )
        )
        self.run_command(command, phased_blocks)

    def phased_blocks_from_ccs_snps(self):
        with open(self.files.chr_lengths, "w") as outfh:
            with open("%s.fai" % self.files.ref, "r") as infh:
                for line in infh:
                    fields = line.rstrip().split("\t")
                    outfh.write("%s\t%s\n" % (fields[0], fields[1]))
        self.phased_blocks(
            self.files.phased_blocks,
            self.files.chr_lengths,
            self.files.phased_snps_vcf,
        )

    def phase_snvs_with_merged_seq(self):
        bams = [self.files.ccs_to_ref, self.files.merged_assembly_to_ref]
        self.phase_snps(
            self.files.merged_phased_snps_vcf,
            self.files.phased_snps_vcf,
            bams,
        )

    def phased_blocks_from_merged_seq(self):
        self.phased_blocks(
            self.files.phased_blocks_merged_seq,
            self.files.chr_lengths,
            self.files.merged_phased_snps_vcf,
        )
