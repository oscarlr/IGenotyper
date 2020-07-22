#!/bin/env python
from IGenotyper.files import FileManager

import pysam

def add_arguments(subparser):
    subparser.add_argument('--database', metavar='DB', help='Fasta DB with alleles')
    subparser.add_argument('--num_reads', metavar='NUM_READS',default=5,help='Number of reads to support allele call')
    subparser.add_argument('ccs', metavar='BAM', help='PacBio CCS bam file')
    subparser.add_argument('outdir',metavar='OUTDIR',help='Directory for output')

def create_allele_db_trie(database):
    pass

def create_allele_name_db(database):
    pass

def get_allele_matches(read_seq,allele_trie):
    pass

def run_alleles(
        ccs,
        outdir,
        database,
        num_reads
):
    allele_counts = {}
    allele_trie_db = create_allele_db_trie(database)
    allele_seq_to_names = create_allele_name_db(database)
    samfile = pysam.AlignmentFile(bam)
    for read in samfile.fetch():
        if read.is_unmapped:
            continue
        if read.is_secondary:
            continue
        if read.is_supplementary:
            continue
        read_seq = read.query_sequence
        allele_matches = get_allele_matches(read_seq,allele_trie_db)
        for allele_seq in allele_matches:
            allele_name = allele_seq_to_names[allele_seq]
            if allele_name not in allele_counts:
                allele_counts[allele_name] = 0
            allele_counts[allele_name] += 1
    files = FileManager(outdir,ccs)
    with open(files.alleles_matches_in_ccs,'w') as fh:
        for allele_name in allele_counts:
            fh.write("%s\t%s" % (allele_name,allele_counts[allele_name]))

def main(args):
    run_alleles(**vars(args))


# For every single read
#   Check if there is sequence in the read that aligns to an IMGT databse
