#!/bin/env python
import json
import pysam
from pybedtools import BedTool

from IGenotyper.common.helper import phased_blocks_merged_seq,assembly_coords,get_ref_seq,skip_read

from IGenotyper.detect.msa_to_variants_without_hmm import path_to_variants

def get_haps_coords(coords):
    hap_coords = {}
    for coord in coords:
        chrom = coord[0]
        start = coord[1]
        end = coord[2]
        hap = coord[3]
        if hap not in hap_coords:
            hap_coords[hap] = []
        hap_coords[hap].append([chrom,start,end])
    return hap_coords

def is_coverage(feature,cov = 1):
    return int(feature.name) == cov

def filter_coords(files,contig_regions,phased_regions):
    contig_coords_bed = BedTool(contig_regions)
    contig_cov = contig_coords_bed.genomecov(bg=True,g="%s.fai" % files.ref)
    contig_cov_filtered = contig_cov.filter(is_coverage)
    phased_contig_cov_filtered = contig_cov_filtered.intersect(phased_regions)
    return phased_contig_cov_filtered        
    
def msa_coords(files):
    phased_regions = phased_blocks_merged_seq(files)
    phased_regions = get_haps_coords(phased_regions)
    contig_coords = assembly_coords(files.merged_assembly_to_ref_phased)
    contig_coords = get_haps_coords(contig_coords)
    filtered_coords = {}
    for hap in ["0","1","2"]:
        filtered_coords[hap] = filter_coords(files,contig_coords[hap],phased_regions[hap])
    coords = {
        "phased": filtered_coords["1"].intersect(filtered_coords["2"]),
        "unphased": filtered_coords["0"]
    }
    return coords

def extract_assembly_sequence(files,chrom,start,end,hap):
    hap_sequence = None
    samfile = pysam.AlignmentFile(files.merged_assembly_to_ref_phased,'rb')
    for contig in samfile.fetch(chrom,start,end):
        if skip_read(contig):
            continue
        if contig.reference_start > start:
            continue
        if contig.reference_end < end:
            continue
        if contig.get_tag("RG",True)[0] != hap:
            continue
        aligned_pairs = contig.get_aligned_pairs()
        query_start = None
        query_end = None
        for query_pos, ref_pos in aligned_pairs:
            if query_pos == None:
                continue
            if ref_pos == None:
                continue
            if int(ref_pos) <= int(start):
                query_start = query_pos
            query_end = query_pos
            if int(ref_pos) > int(end):
                break
        assert query_start != None
        assert query_end != None
        hap_sequence = contig.query_sequence[query_start:query_end]
    return hap_sequence

def extract_sequence(files,chrom,start,end,haps):
    outfasta = "%s/%s_%s_%s.fasta" % (files.msa_fasta,chrom,start,end)
    with open(outfasta,'w') as outfasta_fh:
        outfasta_fh.write(">ref\n%s\n" % get_ref_seq(files,chrom,start,end))
        for h in haps:
            seq = extract_assembly_sequence(files,chrom,start,end,h)
            if h != "0":
                outfasta_fh.write(">hap%s\n%s\n" % (h,seq))
            else:
                outfasta_fh.write(">hap1\n%s\n" % seq)
                outfasta_fh.write(">hap2\n%s\n" % seq)
    return outfasta

def call_variants(files,fastafile,chrom,start,end,command_line_tools,variants):
    msa_fn = "%s/%s_%s_%s.clu" % (files.msa_msa,chrom,start,end)
    variants_fn = "%s/%s_%s_%s.bed" % (files.msa_variants,chrom,start,end)
    command_line_tools.run_kalign(fastafile,msa_fn)
    path_to_variants(msa_fn,chrom,start,end,variants_fn)    
    with open(variants_fn,'r') as fh:
        for line in fh:
            line = line.strip().split('\t')
            variants.append(line)
    
def detect_msa_variants(files,sample,command_line_tools):
    variants = []
    coords = msa_coords(files)
    for phase in coords:
        if phase == "phased":
            haps = ["1","2"]
        else:
            haps = ["0"]
        for coord in coords[phase]:
            chrom = str(coord[0])
            start = int(coord[1])
            end = int(coord[2])            
            fastafile = extract_sequence(files,chrom,start,end,haps)
            call_variants(files,fastafile,chrom,start,end,command_line_tools,variants)
    variants.sort(key=lambda x: int(x[2]))
    header = ["chrom","start","end","event","genotype","size","ref","hap1","hap2","seq_start","seq_end","msa"]
    indels_fh = open(files.indels_assembly_bed,'w')
    sv_fh = open(files.sv_assembly_bed,'w')
    indels_fh.write("%s\n" % "\t".join(header))
    sv_fh.write("%s\n" % "\t".join(header))
    for variant in variants:
        variant_size = int(variant[5])
        if variant_size >= 50:
            sv_fh.write("%s\n" % "\t".join(variant))
        else:
            indels_fh.write("%s\n" % "\t".join(variant))
    indels_fh.close()
    sv_fh.close()
