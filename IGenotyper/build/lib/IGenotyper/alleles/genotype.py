#!/bin/env python
from Bio import SeqIO
from collections import Counter

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
    gene_name = entry.split("_")[0].split("=")[1]
    hap = entry.split("_")[1].split("=")[1]
    return (gene_name,hap)

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
        gene_name,hap = query_gene_hap(query_gene)
        added = False
        alleles[(gene_name,hap)] = []
        for hit in matches[(query_gene,query_seq)]:
            hit_gene, hit_allele = hit_gene_allele(hit)
            if hit_gene != gene_name:
                continue
            alleles[(gene_name,hap)].append(("MATCH",hit_allele))
            added = True
        if added:
            continue
        alleles[(gene_name,hap)].append(("NOVEL",query_seq))
    return alleles
            
def write_genotypes(alleles,fn):
    header = ["gene","hap","category","allele","obs_count"]    
    with open(fn,'w') as fh:
        fh.write("%s\n" % "\t".join(header))
        for gene_name,hap in alleles:
            all_alleles = alleles[(gene_name,hap)]
            allele_counts = Counter(all_alleles)
            for allele_count in allele_counts:                
                cat,allele = allele_count
                output = [gene_name,hap,cat,allele,allele_counts[allele_count]]
                fh.write("%s\n" % "\t".join(map(str,output)))


#alleles = genotype_genes(seq_fa,db)
#write_genotypes(alleles,fn)
