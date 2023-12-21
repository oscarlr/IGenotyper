#!/bin/env python

def create_phased_bam_header(unphased_bam):
    unphased_tag = 0
    haplotype1_tag = 1
    haplotype2_tag = 2
    readgroup_unphased = { "ID": unphased_tag }
    readgroup_hap1 = { "ID": haplotype1_tag }
    readgroup_hap2 = { "ID": haplotype2_tag }
    if "RG" in unphased_bam.header:
        for group in unphased_bam.header["RG"]:
            for item in group:
                if item != "ID":
                    readgroup_unphased[item] = group[item]
                    readgroup_hap1[item] = group[item]
                    readgroup_hap2[item] = group[item]
    phased_bam_header = unphased_bam.header.to_dict()
    phased_bam_header["RG"] = [readgroup_unphased,readgroup_hap1,readgroup_hap2]
    return phased_bam_header

def create_tag(hap):
    haptag = ("RG", str(hap), "Z")
    return haptag
