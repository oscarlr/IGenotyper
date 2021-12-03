#!/bin/env python
import json
import pysam

from IGenotyper.files import FileManager
from IGenotyper.alleles.genotype import genotype_genes,write_genotypes

from IGenotyper.common.helper import extract_sequence

def add_arguments(subparser):
    subparser.add_argument('--database', metavar='DB', help='Fasta DB with alleles')
    subparser.add_argument('--num_reads', metavar='NUM_READS',default=5,help='Number of reads to support allele call')
    subparser.add_argument('outdir',metavar='OUTDIR',help='Directory for output')

def run_alleles(
        database,
        num_reads,
        outdir
):
    files = FileManager(outdir)

    with open(files.input_args,'r') as fh:
        phasing_args = json.load(fh)
    sample = phasing_args["sample"]

    extract_sequence(files.assembly_to_ref_phased,files.gene_coords,files.assembly_genes_fasta)    
    assembly_alleles = genotype_genes(files.assembly_genes_fasta,database)
    write_genotypes(assembly_alleles,files.assembly_genes_alleles)

    extract_sequence(files.ccs_to_ref_phased,files.gene_coords,files.ccs_genes_fasta)    
    ccs_alleles = genotype_genes(files.ccs_genes_fasta,database)
    write_genotypes(ccs_alleles,files.ccs_genes_alleles)
    

def main(args):
    run_alleles(**vars(args))
