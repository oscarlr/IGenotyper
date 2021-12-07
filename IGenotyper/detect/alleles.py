#!/bin/env python
from Bio import SeqIO
from collections import Counter

from IGenotyper.common.helper import extract_sequence

def get_matches(query_entries,database_entries):
    matches = {}
    for query_entry in query_entries:
        query_seq = query_entries[query_entry].seq.upper()
        matches[(query_entry,query_seq)] = set()
        for database_entry in database_entries:
            database_seq = database_entries[database_entry].seq.upper()
            if query_seq.count(database_seq) > 0:
                matches[(query_entry,query_seq)].add(database_entry)
            if query_seq.count(database_seq.reverse_complement()) > 0:
                matches[(query_entry,query_seq)].add(database_entry)
    return matches

def query_gene_hap(entry):
    # gene_name = entry.split("_")[0].split("=")[1]
    # hap = entry.split("_")[1].split("=")[1]
    # chrom = entry.split("_")[2].split(":")[0]    
    # start = entry.split("_")[2].split(":")[1].split("-")[0]
    # end = entry.split("_")[2].split(":")[1].split("-")[1]
    gene_name = entry[entry.index("feat=") + len("feat="):entry.index("_hap=")]
    hap = entry[entry.index("hap=") + len("hap="):entry.index("_pos=")]
    pos = entry[entry.index("pos=") + len("pos="):entry.index("_i=")]
    chrom = pos.split(":")[0]
    start = pos.split(":")[1].split("-")[0]
    end = pos.split(":")[1].split("-")[1]
    return (gene_name,hap,chrom,start,end)

def hit_gene_allele(entry):
    hit_gene = entry.split("_")[0].split("=")[1]
    hit_allele = int(entry.split("_")[1].split("=")[1])
    return (hit_gene,hit_allele)

def genotype_genes(seq_fa,db):
    query_entries = SeqIO.to_dict(SeqIO.parse(seq_fa,"fasta"))
    db_entries = SeqIO.to_dict(SeqIO.parse(db,"fasta"))
    matches = get_matches(query_entries,db_entries)
    alleles = {}
    for query_gene, query_seq in matches:
        gene_name,hap,chrom,start,end = query_gene_hap(query_gene)
        added = False
        if (gene_name,hap,chrom,start,end) not in alleles:
            alleles[(gene_name,hap,chrom,start,end)] = []
        for hit in matches[(query_gene,query_seq)]:
            hit_gene, hit_allele = hit_gene_allele(hit)
            if hit_gene != gene_name:
                continue
            alleles[(gene_name,hap,chrom,start,end)].append(("MATCH",hit_allele))
            added = True
        if added:
            continue
        alleles[(gene_name,hap,chrom,start,end)].append(("NOVEL",query_seq))
    return alleles

def write_genotypes(alleles,fn):
    header = ["#chrom","start","end","gene","hap","category","allele","obs_count"]
    with open(fn,'w') as fh:
        fh.write("%s\n" % "\t".join(header))
        for gene_name,hap,chrom,start,end in alleles:
            all_alleles = alleles[(gene_name,hap,chrom,start,end)]
            allele_counts = Counter(all_alleles)
            for allele_count in allele_counts:
                cat,allele = allele_count
                output = [chrom,start,end,gene_name,hap,cat,allele,allele_counts[allele_count]]
                fh.write("%s\n" % "\t".join(map(str,output)))

def detect_alleles(files):
    extract_sequence(files.assembly_to_ref_phased,files.gene_coords,files.assembly_genes_fasta)
    extract_sequence(files.ccs_to_ref_phased,files.gene_coords,files.ccs_genes_fasta)
    
    assembly_alleles = genotype_genes(files.assembly_genes_fasta,files.allele_db)
    ccs_alleles = genotype_genes(files.ccs_genes_fasta,files.allele_db)
    
    write_genotypes(assembly_alleles,files.assembly_genes_alleles)
    write_genotypes(ccs_alleles,files.ccs_genes_alleles)
