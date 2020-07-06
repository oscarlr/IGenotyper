#!/bin/env python
import vcf

class Variant:
    def __init__(self):
        self.chrom = None
        self.pos = None
        self.phased = False
        self.allele_bases = []

    def from_record_allele_bases(self,variant_record,sample_name):
        alleles = variant_record.genotype(sample_name)['GT'].split("|")
        alleles = map(int, alleles)
        allele_bases = [variant_record.REF]
        for alt_base in variant_record.ALT:
            allele_bases.append(alt_base)
        self.allele_bases = map(lambda x: allele_bases[x], alleles)

    def from_record(self,variant_record,sample_name):
        self.chrom = variant_record.CHROM
        self.pos = variant_record.POS
        if "|" in variant_record.genotype(sample_name)['GT']:
            self.phased = True
            self.from_record_allele_bases(variant_record,sample_name)

class Vcf:
    def __init__(self):
        self.variants_by_pos = {}

    def add_variant(self,variant):
        if variant.chrom not in self.variants_by_pos:
            self.variants_by_pos[variant.chrom] = {}
        self.variants_by_pos[variant.chrom][variant.pos] = variant

    def phased_variants(self):
        phased_variants_by_pos = {}
        for chrom in self.variants_by_pos:
            for pos in self.variants_by_pos[chrom]:
                variant = self.variants_by_pos[chrom][pos]
                if not variant.phased:
                    continue
                if chrom not in phased_variants_by_pos:
                    phased_variants_by_pos[chrom] = {}
                phased_variants_by_pos[chrom][pos] = variant
        return phased_variants_by_pos

def read_in_phased_vcf(vcffn,sample_name):
    vcf_file = Vcf()
    vcf_reader = vcf.Reader(open(vcffn, 'r'))
    for record in vcf_reader:
        variant = Variant()
        variant.from_record(record,sample_name)
        vcf_file.add_variant(variant)
    return vcf_file
