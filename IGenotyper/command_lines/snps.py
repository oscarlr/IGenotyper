#!/usr/bin/env python3
from shlex import quote
import sys
from importlib.metadata import version
from IGenotyper.command_lines.whatshap_genotype import ADAPTER_VERSION

from IGenotyper.command_lines.clt import CommandLine
from IGenotyper.common.validation import validate_vcf


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
        self.run_command(command, snp_candidates,
                         validator=lambda outputs: validate_vcf(outputs[0], self.sample),
                         inputs=[bam, self.files.ref])

    def snp_candidates_from_ccs(self):
        self.snp_candidates(self.files.ccs_to_ref, self.files.snp_candidates)

    def snp_genotypes(self, bam, snp_candidates, snps_vcf):
        sites = validate_vcf(snp_candidates, self.sample, reference=self.files.ref, bam=bam)
        command = (
            "%s -m IGenotyper.command_lines.whatshap_genotype genotype --sample %s --ignore-read-groups "
            "--reference %s -o %s %s %s"
            % (
                quote(sys.executable),
                quote(self.sample),
                quote(self.files.ref),
                quote(snps_vcf),
                quote(snp_candidates),
                quote(bam),
            )
        )
        command += " # adapter=%s whatshap=%s" % (ADAPTER_VERSION, version("whatshap"))
        self.run_command(command, snps_vcf,
                         validator=lambda outputs: validate_vcf(outputs[0], self.sample, expected=sites),
                         inputs=[bam, snp_candidates, self.files.ref])

    def snp_genotypes_from_ccs(self):
        self.snp_genotypes(
            self.files.ccs_to_ref, self.files.snp_candidates, self.files.snps_vcf
        )

    def phase_snps(self, phased_snps_vcf, snps_vcf, bams):
        sites = validate_vcf(snps_vcf, self.sample)
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
        self.run_command(command, phased_snps_vcf,
                         validator=lambda outputs: validate_vcf(outputs[0], self.sample, expected=sites),
                         inputs=[snps_vcf, self.files.ref] + bams)

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
        self.run_command(command, phased_blocks, inputs=[phased_snps_vcf, chr_lengths])

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
