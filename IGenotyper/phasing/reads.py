#!/bin/env python
import pysam

from IGenotyper.common.helper import non_emptyfile
from IGenotyper.common.vcffn import read_in_phased_vcf
from IGenotyper.common.bamfn import create_phased_bam_header,create_tag

def calculate_tag(basetups,hap_dict):
    '''
    Returns a 0, 1 or -1 as the value of a read
    '''
    unphased_tag = 0
    haplotype1_tag = 1
    haplotype2_tag = 2
    b1_bases = 0.0
    b0_bases = 0.0
    thres = .9
    for basetup in basetups:
        rpos, b, qual = basetup
        b0, b1 = hap_dict[rpos].allele_bases
        if b == b1 or b == b0: # check to make sure the base is called
            if b == b0:
                b0_bases += 1.0
            else:
                b1_bases += 1.0
    if (b1_bases + b0_bases) == 0:
        return create_tag(unphased_tag)
    b0_frac = b0_bases/(b1_bases + b0_bases)
    b1_frac = b1_bases/(b1_bases + b0_bases)
    if b0_frac >= thres:
        return create_tag(haplotype1_tag)
    elif b1_frac >= thres:
        return create_tag(haplotype2_tag)
    else:
        return create_tag(unphased_tag)

def phase_read(read,snps,chrom):
    unphased_tag = 0
    haplotype1_tag = 1
    haplotype2_tag = 2
    base_tuples = []
    if chrom in snps:
        for read_pos, ref_pos in read.get_aligned_pairs():
            if ref_pos is None:
                continue
            if ref_pos in snps[chrom]:
                if read_pos is not None:
                    qual = 8
                    if read.query_qualities != None:
                        qual = read.query_qualities[read_pos]
                    base = read.query_sequence[read_pos].upper()
                    base_tuples.append((ref_pos, base, qual))
    if len(base_tuples) > 0:
        read_group_tag = calculate_tag(base_tuples, snps[chrom])
    else:
        read_group_tag = create_tag(unphased_tag)
    read_tags = read.get_tags()
    tags_to_add = []
    for tag in read_tags:
        if tag[0] != "RG":
            tags_to_add.append(tag)
    tags_to_add.append(read_group_tag)
    read.set_tags(tags_to_add)
    return read

def phase_alignments(vcffn,bam,sample,phased_bam):
    vcf = read_in_phased_vcf(vcffn,sample)
    phased_bam_index = "%s.bai" % phased_bam
    if non_emptyfile(phased_bam_index):
        return None
    phase_snps = vcf.phased_variants()
    unphased_bam = pysam.AlignmentFile(bam, 'rb')
    phased_bam_header = create_phased_bam_header(unphased_bam)
    phased_bam = pysam.AlignmentFile(outbam,'wb',header=phased_bam_header)
    for read in unphased_bam.fetch():
        if read.is_unmapped:
            continue
        chrom = unphased_bam.get_reference_name(read.reference_id)
        tagged_read = phase_read(read,phase_snps,chrom)
        phased_bam.write(tagged_read)
    unphased_bam.close()
    phased_bam.close()
    pysam.index(outbam)
