#!/bin/env python
import os
import json
from IGenotyper.common.helper import create_folders,non_emptyfile,run_type
from IGenotyper.common.reference import validate_reference

class FileManager():
    def __init__(self,outdir,bam=None,tmp = "tmp",rhesus=False,data_dir=None):
        self.outdir = outdir
        self.input_bam = bam
        self.tmp = tmp
        self.rhesus = rhesus

        if tmp == "tmp":
            self.tmp = "%s/%s" % (self.outdir,tmp)

        self.folder_structure()
        
        if non_emptyfile("%s/args.json" % self.log):
            with open("%s/args.json" % self.log, 'r') as fh:
                phasing_args = json.load(fh)
            if bam is None:
                self.input_bam = phasing_args["bam"]
                self.tmp = phasing_args["tmp"]
            if data_dir is None:
                data_dir = phasing_args.get("data_dir")

        self.file_structure(data_dir)

    def folder_structure(self):
        self.preprocess = "%s/preprocessed" % self.outdir
        self.variants = "%s/variants" % self.outdir
        self.alignments = "%s/alignments" % self.outdir
        self.assembly = "%s/assembly" % self.outdir
        self.assembly_merge = "%s/assembly/merge" % self.outdir
        self.alleles = "%s/alleles" % self.outdir
        self.msa = "%s/variants/msa" % self.outdir
        self.msa_fasta = "%s/variants/msa/fasta" % self.outdir
        self.msa_msa = "%s/variants/msa/msa" % self.outdir
        self.msa_variants = "%s/variants/msa/variants" % self.outdir
        self.plots = "%s/plots" % self.outdir
        self.log = "%s/logs" % self.outdir
        self.package_directory = os.path.dirname(os.path.abspath(__file__))
        self.scripts = "%s/scripts" % self.package_directory        
        self.templates = "%s/templates" % self.package_directory
        self.ancestry = "%s/ancestry" % self.package_directory
        
        folders = [
            self.outdir,
            self.preprocess,
            self.alignments,
            self.variants,
            self.alleles,
            self.assembly,
            self.tmp,
            self.log,
            self.msa,
            self.assembly_merge,
            self.msa_fasta,
            self.msa_msa,
            self.msa_variants,
            self.plots
        ]

        create_folders(folders)

    def file_structure(self, data_dir=None):
        #pacbio_machine = run_type(self.input_bam)
        pacbio_machine = "SEQUELII" 
        self.ccs_bam = self.input_bam #"%s/ccs.bam" % self.preprocess
        self.ccs_pbi = "%s.pbi" % self.input_bam #"%s/ccs.bam.pbi" % self.preprocess

        self.ccs_fastq = "%s/ccs.fasta" % self.tmp
        self.ccs_fastq_unedited = "%s/ccs.fasta" % self.tmp

        package_data = "%s/data" % self.package_directory
        self.reference_annotations = "%s/immune_receptor_genomics/251106" % package_data

        if self.rhesus:
            reference_data = "%s/rhesus" % package_data
            compatibility_data = reference_data
        else:
            requested_data = data_dir or os.environ.get("IGENOTYPER_DATA_DIR")
            default_data = os.path.join(
                os.environ.get("XDG_DATA_HOME", os.path.expanduser("~/.local/share")),
                "igenotyper",
            )
            reference_data = os.path.abspath(
                os.path.expanduser(requested_data or default_data)
            )
            # Keep existing development checkouts working, while installed
            # packages use the external user data directory by default.
            legacy_reference = os.path.join(package_data, "reference.fasta")
            if not requested_data and not os.path.isfile(
                os.path.join(reference_data, "reference.fasta")
            ) and os.path.isfile(legacy_reference):
                reference_data = package_data
            compatibility_data = package_data

        self.data_directory = reference_data
        self.ref = "%s/reference.fasta" % reference_data
        self.minimap2_ref = (
            "%s.mmi" % self.ref if non_emptyfile("%s.mmi" % self.ref) else self.ref
        )
        self.target_regions = "%s/target_regions.bed" % compatibility_data
        self.sv_coords = "%s/sv_coords.bed" % compatibility_data
        self.lpart1 = "%s/lpart1.bed" % compatibility_data
        self.rss = "%s/rss.bed" % compatibility_data
        self.vdj_coords = "%s/vdj_coords.bed" % compatibility_data
        self.allele_db = "%s/alleles.fasta" % compatibility_data

        # These tracks are versioned with the current reference FASTA.
        self.introns = "%s/intron.bed" % self.reference_annotations
        self.gene_coords = "%s/gene.bed" % self.reference_annotations
        self.constant_gene_coords = "%s/constant.bed" % self.reference_annotations

        if not self.rhesus:
            validate_reference(
                self.ref,
                "%s/reference.fasta.fai" % self.reference_annotations,
                [
                    self.target_regions,
                    self.sv_coords,
                    self.introns,
                    self.lpart1,
                    self.rss,
                    self.gene_coords,
                    self.constant_gene_coords,
                    self.vdj_coords,
                ],
            )

        self.subreads_to_ref = "%s/subreads_to_ref.sorted.bam" % self.preprocess
        self.subreads_to_ref_phased = "%s/subreads_to_ref_phased.sorted.bam" % self.alignments
        self.ccs_to_ref = "%s/ccs_to_ref.sorted.bam" % self.preprocess
        self.ccs_to_ref_phased = "%s/ccs_to_ref_phased.sorted.bam" % self.alignments

        self.snp_candidates = "%s/snvs_candidates.vcf" % self.tmp
        self.snps_vcf = "%s/snvs_from_ccs.vcf" % self.variants
        self.phased_snps_vcf = "%s/snvs_phased_from_ccs.vcf" % self.variants
        self.snps_assembly_vcf = "%s/snvs_assembly.vcf" % self.variants
        self.snps_igh_assembly_vcf = "%s/snvs_igh_assembly_subs.vcf" % self.variants
        self.indels_assembly_bed = "%s/indel_assembly.bed" % self.variants
        self.sv_assembly_bed = "%s/sv_assembly.bed" % self.variants
        
        self.alleles_matches_in_ccs = "%s/alleles_matches_in_ccs.txt" % self.alleles

        self.assembly_fasta = "%s/contigs.fasta" % self.assembly
        self.assembly_fastq = "%s/contigs.fastq" % self.assembly
        self.assembly_to_ref = "%s/contigs_to_ref.sorted.bam" % self.alignments
        self.assembly_to_ref_phased = "%s/contigs_to_ref_phased.sorted.bam" % self.alignments

        self.igh_assembly_fasta = "%s/igh_contigs.fasta" % self.assembly
        self.igh_assembly_fastq = "%s/igh_contigs.fastq" % self.assembly
        self.igh_assembly_to_ref = "%s/igh_contigs_to_ref.sorted.bam" % self.alignments
        self.igh_assembly_to_ref_phased = "%s/igh_contigs_to_ref_phased.sorted.bam" % self.alignments

        self.igh_assembly_to_ref_subs = "%s/igh_contigs_to_ref_subs.sorted.bam" % self.alignments
        self.igh_assembly_to_ref_subs_phased = "%s/igh_contigs_to_ref_subs_phased.sorted.bam" % self.alignments

        self.phased_blocks = "%s/phased_blocks.txt" % self.variants
        self.legacy_chr_lengths = "%s/chr_lengths.txt" % compatibility_data
        self.chr_lengths = "%s/chr_lengths.txt" % self.log
        self.input_args = "%s/args.json" % self.log
        self.assembly_script = "%s/data/assembly.sh" % self.package_directory

        # self.merged_assembly = "%s/merged_assembly.fasta" % self.assembly_merge
        # self.contigs_blasted = "%s/blast.txt" % self.assembly_merge
        # self.merged_contigs = "%s/merged_contigs.txt" % self.assembly_merge
        # self.merged_assembly_to_ref = "%s/merged_assembly_to_ref.sorted.bam" % self.alignments
        # self.merged_assembly_to_ref_phased = "%s/merged_assembly_to_ref_phased.sorted.bam" % self.alignments
        # self.merged_phased_snps_vcf = "%s/snvs_phased_from_merged_seq.vcf" % self.variants
        # self.phased_blocks_merged_seq = "%s/phased_blocks_merged_seq.txt" % self.variants

        self.plot_ccs_read_lengths = "%s/ccs_read_lengths.png" % self.plots
        self.plot_ccs_read_quals = "%s/ccs_read_quals.png" % self.plots
        self.plot_phasing = "%s/phasing.png" % self.plots
        self.plot_gene_cov = "%s/gene_cov.png" % self.plots
        self.plot_sv_gene_cov = "%s/sv_gene_cov.png" % self.plots

        self.report = "%s/report.html" % self.outdir
        self.stats_json = "%s/stats.json" % self.log
        self.gene_cov = "%s/gene_cov.txt" % self.log
        
        self.assembly_genes_fasta = "%s/assembly_genes.fasta" % self.alleles
        self.assembly_genes_alleles = "%s/assembly_alleles.bed" % self.alleles
        self.ccs_genes_fasta = "%s/ccs_genes.fasta" % self.alleles
        self.ccs_genes_alleles = "%s/ccs_alleles.bed" % self.alleles

        self.constant_assembly_genes_fasta = "%s/constant_assembly_genes.fasta" % self.alleles
        self.constant_assembly_genes_alleles = "%s/constant_assembly_alleles.bed" % self.alleles
        self.constant_ccs_genes_fasta = "%s/constant_ccs_genes.fasta" % self.alleles
        self.constant_ccs_genes_alleles = "%s/constant_ccs_alleles.bed" % self.alleles
