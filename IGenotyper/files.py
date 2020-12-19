#!/bin/env python
import os
import json
from helper import create_folders,non_emptyfile

class FileManager():
    def __init__(self,outdir,bam=None,tmp = "tmp"):
        self.outdir = outdir
        self.input_bam = bam
        self.tmp = tmp

        if tmp == "tmp":
            self.tmp = "%s/%s" % (self.outdir,tmp)

        self.folder_structure()
        self.file_structure()

        if non_emptyfile(self.input_args):
            with open(self.input_args, 'r') as fh:
                phasing_args = json.load(fh)
            self.input_bam = phasing_args["bam"]
            self.tmp = phasing_args["tmp"]

    def folder_structure(self):
        self.preprocess = "%s/preprocessed" % self.outdir
        self.variants = "%s/variants" % self.outdir
        self.alignments = "%s/alignments" % self.outdir
        self.assembly = "%s/assembly" % self.outdir
        self.alleles = "%s/alleles" % self.outdir
        self.msa = "%s/variants/msa" % self.outdir
        self.log = "%s/logs" % self.outdir
        self.package_directory = os.path.dirname(os.path.abspath(__file__))
        self.scripts = "%s/scripts" % self.package_directory

        folders = [
            self.outdir,
            self.preprocess,
            self.alignments,
            self.variants,
            self.alleles,
            self.assembly,
            self.tmp,
            self.log,
            self.msa
        ]

        create_folders(folders)

    def file_structure(self):
        self.ccs_bam = "%s/ccs.bam" % self.preprocess
        self.ccs_pbi = "%s/ccs.bam.pbi" % self.preprocess
        self.ccs_fastq = "%s/ccs.fastq" % self.tmp
        self.ccs_fastq_unedited = "%s/ccs.fastq" % self.tmp

        self.ref = "%s/data/reference.fasta" % self.package_directory
        self.target_regions = "%s/data/target_regions.bed" % self.package_directory
        self.sv_coords = "%s/data/sv_coords.bed" % self.package_directory
        self.introns = "%s/data/introns.bed" % self.package_directory
        self.lpart1 = "%s/data/lpart1.bed" % self.package_directory
        self.rss = "%s/data/rss.bed" % self.package_directory
        self.gene_coords = "%s/data/gene_coords.bed" % self.package_directory
        self.vdj_coords = "%s/data/vdj_coords.bed" % self.package_directory
        
        self.subreads_to_ref = "%s/subreads_to_ref.sorted.bam" % self.preprocess
        self.subreads_to_ref_phased = "%s/subreads_to_ref_phased.sorted.bam" % self.alignments
        self.ccs_to_ref = "%s/ccs_to_ref.sorted.bam" % self.preprocess
        self.ccs_to_ref_phased = "%s/ccs_to_ref_phased.sorted.bam" % self.alignments

        self.snp_candidates = "%s/snvs_candidates.vcf" % self.tmp
        self.snvs_vcf = "%s/snvs_from_ccs.vcf" % self.variants
        self.phased_snvs_vcf = "%s/snvs_phased_from_ccs.vcf" % self.variants
        self.snvs_assembly_vcf = "%s/snvs_assembly.vcf" % self.variants
        self.indels_assembly_bed = "%s/indel_assembly.bed" % self.variants
        self.sv_assembly_bed = "%s/sv_assembly.bed" % self.variants
        
        self.alleles_matches_in_ccs = "%s/alleles_matches_in_ccs.txt" % self.alleles

        self.assembly_fasta = "%s/contigs.fasta" % self.assembly
        self.assembly_fastq = "%s/contigs.fastq" % self.assembly
        self.assembly_to_ref = "%s/contigs_to_ref.sorted.bam" % self.alignments

        self.phased_blocks = "%s/phased_blocks.txt" % self.variants
        self.input_args = "%s/args.json" % self.log
        self.assembly_script = "%s/data/assembly.sh" % self.package_directory

