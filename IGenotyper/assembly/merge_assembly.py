#!/bin/env python
import pysam
import networkx
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import namedtuple

from networkx.algorithms.components.connected import connected_components

from IGenotyper.common.helper import skip_read,assembly_location,get_haplotype,get_mapping_pos

def assembly_seq(files):
    seqfn = "%s/seq.fasta" % files.tmp
    sam = pysam.AlignmentFile(files.assembly_to_ref)
    with open(seqfn,'w') as fh:
        for read in sam:
            if skip_read(read):
                continue
            fh.write(">%s\n%s\n" % (read.query_name,read.query_sequence))
    return seqfn

def not_overlapping(alignment):
    flank = 500
    if not (int(alignment.qstart) < flank and (int(alignment.slen) - int(alignment.send)) < flank):
        if not (int(alignment.sstart) < flank and (int(alignment.qlen) - int(alignment.qend)) < flank):
            return True
    return False

def overlapping_contig(alignment):
    #  ----------------
    # -------
    # or
    # -------
    #  ----------------
    contig = None
    thres = 100
    if int(alignment.qstart) < 100 and (int(alignment.qlen) - int(alignment.qend)) < 100:
        contig = alignment.qseqid
    if int(alignment.sstart) < 100 and (int(alignment.slen) - int(alignment.send)) < 100:
        contig = alignment.sseqid
    return contig

def filter_criteria(alignment,mismatch_max,length_min,contigs_to_ignore):
    filter_ = False
    if int(alignment.length) < length_min:
        filter_ = True
    if alignment.qseqid == alignment.sseqid:
        filter_ = True
    if float(alignment.mismatch) > mismatch_max:
        filter_ = True
    if alignment.sstrand != "plus":
        filter_ = True
    return filter_

def wrong_placement(alignment):
    filter_ = False
    if int(alignment.sstart) > int(alignment.qstart):
        filter_ = True
    if not filter_:
        q_chrom,q_start,q_end = assembly_location(alignment.qseqid)
        s_chrom,s_start,s_end = assembly_location(alignment.sseqid)
        start_node,end_node = get_start_end_nodes(alignment)
        if start_node == alignment.qseqid:
            if q_start > s_start:
                filter_ = True
        if start_node == alignment.sseqid:
            if s_start > q_start:
                filter_ = True
    return filter_
    
def filter_blast(blastout,mismatch_max=3,length_min=1000):
    pairs = []
    contigs_to_ignore = []
    alignments = {}
    columns = ["length","pident","nident","mismatch","gapopen","gaps","qseqid",
               "qstart","qend","qlen","sseqid","sstart","send","slen","sstrand"]
    Alignment = namedtuple('Alignment',columns)
    with open(blastout,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            alignment = Alignment._make(line)
            entry = tuple(sorted([alignment.qseqid,alignment.sseqid]))
            if entry in alignments:
                continue
            if filter_criteria(alignment,mismatch_max,length_min,contigs_to_ignore):
                continue
            contig = overlapping_contig(alignment)
            if contig != None:
                contigs_to_ignore.append(contig)
                continue
            if alignment.qseqid in contigs_to_ignore:
                continue
            if alignment.sseqid in contigs_to_ignore:
                continue
            if not_overlapping(alignment):
                continue
            if wrong_placement(alignment):
                continue
            alignments[entry] = alignment
    return (alignments,contigs_to_ignore)

def get_start_end_nodes(alignment):
    #      ---------------- e
    # ---------- s
    # or
    # ------------- s
    #       ---------------- e
    if int(alignment.qstart) >= int(alignment.sstart):
        assert int(alignment.qlen) - int(alignment.qend) < int(alignment.slen) - int(alignment.send)
        start_node = alignment.qseqid
        end_node = alignment.sseqid
    if int(alignment.sstart) > int(alignment.qstart):
        assert int(alignment.slen) - int(alignment.send) < int(alignment.qlen) - int(alignment.qend)
        start_node = alignment.sseqid
        end_node = alignment.qseqid
    return (start_node,end_node)
    
def create_graph(alignments,type_):
    if type_ == "directed":
        G = networkx.DiGraph()
    elif type_ == "undirected":
        G = networkx.Graph()
    for qseqid,sseqid in alignments:
        if not G.has_node(qseqid):
            G.add_node(qseqid)
        if not G.has_node(sseqid):
            G.add_node(sseqid)
        start_node,end_node = get_start_end_nodes(alignments[(qseqid,sseqid)])
        G.add_edge(start_node,end_node)
    return G

def subset_alignments(alignments,edges):
    subset = {}
    for edge in edges:
        edge = tuple(sorted(edge))
        subset[edge] = alignments[edge]
    return subset

def get_contig_coords(alignments,path):
    if len(path) == 1:
	return [[path[0],1,-1]]
    coords = []
    i = 0
    for start_contig, end_contig in zip(path,path[1:]):
        entry = tuple(sorted([start_contig,end_contig]))
        alignment = alignments[entry]
        if alignment.qseqid != start_contig:
            continue
        if alignment.sseqid != end_contig:
            continue
        merge_coords = [alignment.qseqid,1,int(alignment.qstart) - 1]
        coords.append(merge_coords)
        i += 1
        if i == len(path) - 1:
            merge_coords = [alignment.sseqid,1,alignment.slen]
            coords.append(merge_coords)
    return coords

def path_sequence(files,alignments,contigs,merged_contigs_fh):
    seqfn = "%s/seq.fasta" % files.tmp
    coords = get_contig_coords(alignments,contigs)
    inseq = SeqIO.to_dict(SeqIO.parse(seqfn,"fasta"))
    seqs = []
    for contig_name,start,end in coords:
        merged_contigs_fh.write("\t%s\t%s\t%s\n" % (contig_name,start,end))
        if start == 1 and end == -1:
            seq = str(inseq[contig_name].seq[0:])
        else:
            seq = str(inseq[contig_name].seq[int(start):int(end)])
        seqs.append(seq)
    return "".join(seqs)

def valid_paths(contigs):
    # 0 --> 1 --> 2 --> 0 X --> 2
    # 1 ---> 0 X --> 1
    paths = []
    path = []
    first_hap = None
    visited_unphased = False
    current_hap = None
    current_coord = None
    for contig in contigs:
        hap = get_haplotype(contig)
        coord = assembly_location(contig)
        if hap == "0":
            if first_hap != None:
                visited_unphased = True
        if first_hap == None:
            first_hap = hap
            current_hap = hap
            current_coord = assembly_location(contig)
        else:
            hap_switch = False
            if current_hap == "1" and hap == "2":
                hap_switch = True
            if current_hap == "2" and hap == "1":
                hap_switch = True
            if hap_switch:
                if coord == current_coord:
                    paths.append(path)
                    path = []
                    visited_unphased = False
                    first_hap = None
                    current_hap = None
                    current_coord = None
                    path.append(contig)
                    continue
        if visited_unphased:
            if hap in ["1","2"]:
                paths.append(path)
                path = []
                visited_unphased = False
                first_hap = None
                current_hap = None
                current_coord = None
                path.append(contig)
                continue
        current_hap = hap
        current_coord = coord
        path.append(contig)
    paths.append(path)
    return paths                
        
def merge_contigs(files,alignments,contig_pos):
    seqs = {}
    contigs_used = []
    undirected_graph = create_graph(alignments,"undirected")
    merged_contigs_fh = open(files.merged_contigs,'w')
    for i,group in enumerate(connected_components(undirected_graph)):
        edges =  undirected_graph.subgraph(group).edges()
        alignment_edges = subset_alignments(alignments,edges)
        directed_graph = create_graph(alignment_edges,"directed")
        pos1 = 99999999999
        pos2 = 0
        start_contig = None
        end_contig = None
        for edge in edges:
            for e in edge:
                chrom,start,end = contig_pos[e]
                if start < pos1:
                    pos1 = start
                    start_contig = e
                if end > pos2:
                    pos2 = end
                    end_contig = e
        contigs = networkx.shortest_path(directed_graph,source=start_contig,target=end_contig)
        paths = valid_paths(contigs)
        for j,path in enumerate(paths):
            merged_contigs_fh.write("%s.%s\n" % (i,j))
            contig_seq = path_sequence(files,alignments,path,merged_contigs_fh)
            seqs["%s.%s" % (i,j)] = contig_seq
        contigs_used += contigs
    merged_contigs_fh.close()
    return (seqs,contigs_used)

def write_merge_seqs(files,seqs,contigs_used,contigs_to_ignore):
    seqfn = "%s/seq.fasta" % files.tmp
    inseq = SeqIO.to_dict(SeqIO.parse(seqfn,"fasta"))
    index = max(map(float,seqs.keys())) + 1
    for seq in inseq:
        if seq in contigs_to_ignore:
            continue
        hap = get_haplotype(seq)
        if hap != "0":
            if seq in contigs_used:
                continue
        seqs[index] = str(inseq[seq].seq)
        index += 1
    records = []
    for seq in seqs:
        record = SeqRecord(Seq(seqs[seq],"fasta"),id="%s/0/0_8" % seq,name="",description="")
        records.append(record)
    SeqIO.write(records,files.merged_assembly,"fasta")

def blast_assembly(files,align_command_line):
    forward_strand_seq = assembly_seq(files)
    align_command_line.blast_seq(forward_strand_seq,files.contigs_blasted)    

def merge_sequences(files):
    alignments,contigs_to_ignore = filter_blast(files.contigs_blasted)
    contig_pos = get_mapping_pos(files.assembly_to_ref)
    seqs,contigs_used = merge_contigs(files,alignments,contig_pos)
    write_merge_seqs(files,seqs,contigs_used,contigs_to_ignore)    

def merge_assembly(files,align_command_line,sample):
    blast_assembly(files,align_command_line)
    merge_sequences(files)
    align_command_line.map_merged_assembly()    
    
# def refine_assembly(files,command_tools,sample_name):
#     seqfn = assembly_seq(files)
#     create_directory("%s/refinement" % files.assembly)
#     blastout = "%s/refinement/blast_out.txt" % files.assembly
#     command_tools.blast_seq(seqfn,blastout)
#     map_merge_seq(files,command_tools)
#     outbam = "%s/refinement/merge_seq_to_igh.sorted.bam" % files.assembly
#     phased_outbam = "%s/refinement/merge_seq_to_igh_phased.sorted.bam" % files.assembly
#     outvcf = "%s/refinement/merge_seq_to_igh_phased.vcf" % files.assembly
#     out_phased_blocks = "%s/refinement/merge_seq_to_igh_phased.txt" % files.assembly
#     args = [sample_name,
# 	    files.ref,
#             outvcf,
#             files.snvs_vcf,
#             "%s %s" % (outbam,files.ccs_to_ref_phased)]
#     command_tools.phase_genotype_snvs(sample_name,args)
#     args = [sample_name,
#             out_phased_blocks,
#             outvcf]
#     command_tools.phase_blocks(sample_name,args)
#     vcf = read_in_phased_vcf(outvcf,sample_name)
#     phase_reads(outbam,phased_outbam,vcf)
