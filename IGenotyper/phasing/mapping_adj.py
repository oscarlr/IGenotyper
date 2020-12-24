#!/bin/env python
import os
import pysam

from IGenotyper.helper import read_is_unphased

def get_new_bam_file(tmp,bam,iteration):
    new_bam_file = "%s/%s.sam" % (tmp,iteration)
    old_bam_file = pysam.AlignmentFile(bam)
    new_sam_file = pysam.AlignmentFile(new_bam_file,'w',template=old_bam_file)
    old_bam_file.close()
    return new_sam_file

def get_secondary_alignment_reads(old_bam_file):
    secondary_alignments = {}
    samfile = pysam.AlignmentFile(old_bam_file)
    for read in samfile.fetch("igh"):
        if read.is_secondary:
            if read.query_name not in secondary_alignments:
                secondary_alignments[read.query_name] = []
            secondary_alignments[read.query_name].append(read)
    samfile.close()
    return secondary_alignments

def change_alignment_skip_read(read):
    skip = False
    if read.is_secondary:
        skip = True
    if read.is_supplementary:
        skip = True
    if not read_is_unphased(read):
        skip = True
    return skip

def supplementary_score_diff(read,secondary_reads,thres=500):
    diffs = [] 
    primary_read_score = float(read.get_tag("AS",True)[0])
    for secondary_read in secondary_reads:
        secondary_read_score = float(secondary_read.get_tag("AS",True)[0])
        diffs.append(primary_read_score-secondary_read_score)
    if abs(max(diffs)) > thres:
        return True
    return False

def change_read(primary_read,secondary_reads):
    thres = 2500
    secondary_to_primary_read = None
    secondary_to_primary_flag = None
    secondary_to_primary_mapq = None
    index_to_ignore = None
    secondary_to_primary_score = 0
    for i,secondary_read in enumerate(secondary_reads):
        secondary_read_score = secondary_read.get_tag("AS",True)[0]
        if secondary_read_score < secondary_to_primary_score:
            secondary_to_primary_read = secondary_read
            secondary_to_primary_score = secondary_read_score
            secondary_to_primary_flag = secondary_read.flag
            secondary_to_primary_mapq = secondary_read.mapping_quality
            index_to_ignore = i
    if abs(float(primary_read.get_tag("AS",True)[0]) - float(secondary_to_primary_score)) < thres:
        # Change flag to primary alignment
        #secondary_to_primary_read.flag = primary_read.flag
        secondary_to_primary_read.flag = secondary_to_primary_read.flag - 256
        secondary_to_primary_read.mapping_quality = primary_read.mapping_quality
        # Change flag to secondary alignment
        #primary_read.flag = secondary_to_primary_flag
        primary_read.flag = primary_read.flag + 256
        primary_read.mapping_quality = secondary_to_primary_mapq
        #reads_to_return = [primary_read,secondary_to_primary_read]
        reads_to_return = [secondary_to_primary_read]
        for i,read in enumerate(secondary_reads):
            if i == index_to_ignore:
                continue
            reads_to_return.append(read)
    else:
        reads_to_return = secondary_reads
        reads_to_return.append(primary_read)
    return reads_to_return


def change_primary_alignments(new_sam_file,old_bam_file):
    secondary_alignment_reads = get_secondary_alignment_reads(old_bam_file)
    old_sam_file = pysam.AlignmentFile(old_bam_file)
    changed_reads = set()
    for read in old_sam_file.fetch():#"igh"):
        if change_alignment_skip_read(read):
            continue
        if read.query_name in secondary_alignment_reads:            
            secondary_reads = secondary_alignment_reads[read.query_name]
            if read.mapping_quality > 40:
                if supplementary_score_diff(read,secondary_reads):
                    continue
            output_reads = change_read(read,secondary_reads)            
            for output_read in output_reads:
                changed_reads.add(output_read.query_name)
                new_sam_file.write(output_read)
    old_sam_file.close()
    return changed_reads
            
def fix_alignments(tmp,bam,iteration):
    new_sam_file = get_new_bam_file(tmp,bam,iteration)
    moved_reads = change_primary_alignments(new_sam_file,bam)
    old_sam_file = pysam.AlignmentFile(bam)
    for read in old_sam_file:
        if read.query_name in moved_reads:
            continue
        new_sam_file.write(read)
    old_sam_file.close()
    new_sam_file.close()
    return new_sam_file

def fix_ccs_alignment(files,align_command_line,iteration):
    sam_file = fix_alignments(files.tmp,files.ccs_to_ref_phased,iteration)
    os.remove(files.ccs_to_ref)
    os.remove("%s.bai" % files.ccs_to_ref)
    os.remove(files.ccs_to_ref_phased)
    os.remove("%s.bai" % files.ccs_to_ref_phased)
    align_command_line.sam_to_sorted_bam("%s/%s" % (files.tmp,iteration),files.ccs_to_ref)

def fix_subread_alignment(files,align_command_line,iteration):
    sam_file = fix_alignments(files.tmp,files.subreads_to_ref_phased,iteration)
    os.remove(files.subreads_to_ref)
    os.remove("%s.bai" % files.subreads_to_ref)
    os.remove(files.subreads_to_ref_phased)
    os.remove("%s.bai" % files.subreads_to_ref_phased)
    align_command_line.sam_to_sorted_bam("%s/%s" % (files.tmp,iteration),files.subreads_to_ref)
