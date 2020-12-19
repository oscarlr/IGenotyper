#!/bin/env python
#from MsPAC.python_scripts.hmm import *
#from hmm2 import *
from collections import Counter
#from pomegranate import *
from Bio import AlignIO
#import numpy as np
#import pysam
import sys

def get_msa_sequence(clufn):
    alignment = AlignIO.read(clufn,"clustal")
    sequences = {}
    for i,sequence in enumerate(alignment):        
        if sequence.id not in ["ref","hap1","hap2"]:
            if i == 0:
                sequence.id = "ref"
            if i == 1:
                sequence.id = "hap1"
            if i == 2:
                sequence.id = "hap2"
        sequences[sequence.id] = str(sequence.seq).upper()
    return sequences

def get_observations(sequence):
    obs = []
    if len(sequence) == 3:
        h1_seq = sequence["hap1"]
        h2_seq = sequence["hap2"]
        chrom_seq = sequence["ref"]
        for h1,h2,chrom in zip(h1_seq,h2_seq,chrom_seq):
            obs.append(observations["3"][h1][h2][chrom])
        return obs

# def get_quality_scores(sequence,quality_scores,start_index,end_index,padding):
#     quality_start = start_index - sequence[:start_index].count("-") #- padding
#     quality_end = end_index - sequence[:end_index].count("-") #+ padding
#     #print quality_start,quality_end
#     if quality_end - quality_start < (padding*2):
#         quality_start = max(0,quality_start - padding)
#         quality_end = quality_end + padding
#     scores = quality_scores[quality_start:quality_end]
#     #print scores,start_index,end_index,quality_start,quality_end
#     mean = float(sum(scores))/len(scores)
#     return mean

def get_state(ref,h1,h2):
    state = "NORMAL"
    if ref == "-":
        if h1 != "-" and h2 != "-":
            state = "INS_1|1"
        elif h1 != "-":
            state = "INS_1|0"
        elif h2 != "-":
            state = "INS_0|1"
    else:
        if h1 == "-" and h2 == "-":
            state = "DEL_1|1"
        elif h1 == "-":
            state = "DEL_1|0"
        elif h2 == "-":
            state = "DEL_0|1"
    return state
            
def three_sequence_msa_variants(clufn,chrom,ref_start,ref_end,outfn): #,quality_scores,padding):
    sequence = get_msa_sequence(clufn)
    start = False
    current_sv = None
    current_sv_start_index = None
    current_sv_current_index = None
    current_sv_hap1_seq = []    
    current_sv_hap2_seq = []
    current_sv_ref_seq = []
    msa_states = []
    ref_index = -1
    ref_sv_start = None
    ref_sv_index = None
    contain_indel = False
    outfh = open(outfn,'w')
    for i,(h1,h2,ref) in enumerate(zip(sequence["hap1"],sequence["hap2"],sequence["ref"])):
        if ref != "-":
            ref_index += 1
        if start == False:
            if "-" not in (h1,h2,ref):
                start = True
            continue
        msa_state = get_state(ref,h1,h2)                              
        if msa_state != "NORMAL":
            msa_states.append(msa_state)
            if current_sv == None:
                current_sv = msa_state
                current_sv_start_index = i
                current_sv_current_index = i
                ref_sv_start = ref_index
                ref_sv_index = ref_index
                continue
            current_sv_current_index += 1 
            ref_sv_index = ref_index
            if "-" in [h1,h2,ref]:
                contain_indel = True
            if h1 != "-":
                current_sv_hap1_seq.append(h1)
            if h2 != "-":
                current_sv_hap2_seq.append(h2)
            if ref != "-":
                current_sv_ref_seq.append(ref)
        if i != current_sv_current_index and current_sv != None:
            if i + 1 == len(sequence["hap1"]) and "-" in (h1,h2,ref):
                continue
            sv_type, genotype = Counter(msa_states).most_common(1)[0][0].split("_") #current_sv.split("_")
            sv_len = max(len(current_sv_hap1_seq),len(current_sv_hap2_seq),len(current_sv_ref_seq))
            # if quality_scores != None:
            #     hap1_qual_score = get_quality_scores(sequence["hap1"],quality_scores["hap1"][0],current_sv_start_index,current_sv_current_index,padding)
            #     hap2_qual_score = get_quality_scores(sequence["hap2"],quality_scores["hap2"][0],current_sv_start_index,current_sv_current_index,padding)
            # else:
            #     hap1_qual_score = 60
            #     hap2_qual_score = 60
            if contain_indel == False:
                continue
            output = [chrom,                              
                      ref_sv_start + ref_start, # 0-based/UCSC Genome format
                      ref_sv_index + ref_start + 1, # 0-based/UCSC Genome format    
                      sv_type,
                      genotype,
                      sv_len,
                      "".join(current_sv_ref_seq) if len("".join(current_sv_ref_seq)) > 0 else ".",
                      "".join(current_sv_hap1_seq) if len("".join(current_sv_hap1_seq)) > 0 else ".",
                      "".join(current_sv_hap2_seq) if len("".join(current_sv_hap2_seq)) > 0 else ".",
                      current_sv_start_index,
                      current_sv_current_index,
                      clufn]
            outfh.write("%s\n" % "\t".join(map(str,output)))
            current_sv = None
            contain_indel = False
            current_sv_hap1_seq = []    
            current_sv_hap2_seq = []
            current_sv_ref_seq = []
            msa_states = []
    outfh.close()

def path_to_variants(clufn,chrom,start,end,outfn): #,quality_scores,padding):
    three_sequence_msa_variants(clufn,chrom,start,end,outfn) #,quality_scores,padding)

# def load_qual_scores(qual_scoresfn):
#     qual_scores = {}
#     with open(qual_scoresfn,'r') as fh:
#         for line in fh:
#             if ">" in line:
#                 name = line[1:].rstrip()[:4]
#                 contig = line[1:].rstrip()[5:]
#             else:
#                 qual =  line.rstrip().split(',')
#                 qual_scores[name] = (map(int,qual),contig)
#     return qual_scores

# def load_sv_regions(sv_regions_bed):
#     sv_regions = []
#     with open(sv_regions_bed,'w') as fh:
#         for line in fh:
#             line = line.rstrip().split('\t')
#             sv_regions.append(line)
#     return sv_regions

# def main():
#     clufn = sys.argv[1]
#     chrom = sys.argv[2]
#     start = int(sys.argv[3])
#     end = int(sys.argv[4])
#     #qual_scoresfn = sys.argv[5]
#     #qual_scores = load_qual_scores(qual_scoresfn)
#     #padding = int(sys.argv[6])
#     sequence = get_msa_sequence(clufn)
#     #obs = get_observations(sequence)
#     #log, path = model[str(len(sequence))].viterbi(obs)
#     path_to_variants(sequence,clufn,chrom,start,end) #,qual_scores,padding)

# if __name__ == "__main__":
#     main()
