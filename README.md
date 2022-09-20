# IGenotyper

[Introduction](#introduction)  
[Installation](#installation)  
[Getting IGH specific reference](#getting-igh-specific-reference)<br>
[Testing IGenotyper installation](#testing-igenotyper-installation)<br>
[Usage](#usage)<br>
[Running IGenotyper](#running-igenotyper)<br>
[Explanation of steps](#explanation-of-steps)<br>
[Output directories](#output-directories)<br>
[Output files](#output-files)<br>
[Notes](#notes)

## Installation

```

git clone https://github.com/oscarlr/IGv2.git
cd IGv2

conda env create -f pygenometracks.yml

conda env create -f environment.yml
conda activate IGv2
python setup.py install
conda deactivate

conda env create -n whatshap-latest python=3.7
conda activate whatshap-latest
pip install git+https://github.com/whatshap/whatshap
conda deactivate

cd ..
git clone https://github.com/oscarlr/cluster
cd cluster
python setup.py install
export SJOB_DEFALLOC=NONE
```

## Getting IGH specific reference IGH specific reference

```
wget https://watsonlab.blob.core.windows.net/public/reference/hg38_igh_trab/reference.fasta
wget https://watsonlab.blob.core.windows.net/public/reference/hg38_igh_trab/reference.fasta.fai

sawriter reference.fasta
cp reference.fasta* ~/anaconda3/envs/IGv2/lib/python2.7/site-packages/IGenotyper-1.1-py2.7.egg/IGenotyper/data/

# For rheMac10
sawriter reference.fasta
cp reference.fasta* ~/anaconda3/envs/IGv2/lib/python2.7/site-packages/IGenotyper-1.1-py2.7.egg/IGenotyper/data/rhesus/

```

## Testing IGenotyper installation
```
Todo
```

## Running IGenotyper
```
# <pacbio bam file> must have a pbi and bai index. They can be created using pbindex and samtools index, respectively.

IG phase <pacbio bam file> <output> 
IG assemble <output> 
IG detect <output> 
```

## Usage
```
IG phase --help
usage: IG phase [-h] [--rhesus] [--sample SAMPLE] [--threads THREADS]
                [--mem MEM] [--cluster] [--queue QUEUE] [--walltime WALLTIME]
                [--tmp TMP] [--input_vcf VCF]
                BAM OUTDIR

positional arguments:
  BAM                  PacBio bam file
  OUTDIR               Directory for output

optional arguments:
  -h, --help           show this help message and exit
  --rhesus
  --sample SAMPLE      Name of sample
  --threads THREADS    Number of threads
  --mem MEM            Memory for cluster
  --cluster            Use cluster
  --queue QUEUE        Queue for cluster
  --walltime WALLTIME  Walltime for cluster
  --tmp TMP            Temporary folder
  --input_vcf VCF      Phased VCF file to phase reads
```
```
IG assembly --help
usage: IG assembly [-h] [--threads THREADS] [--mem MEM] [--cluster]
                   [--queue QUEUE] [--walltime WALLTIME]
                   OUTDIR

positional arguments:
  OUTDIR               Directory for output

optional arguments:
  -h, --help           show this help message and exit
  --threads THREADS    Number of threads
  --mem MEM            Memory for cluster
  --cluster            Use cluster
  --queue QUEUE        Queue for cluster
  --walltime WALLTIME  Walltime for cluster
```
```
IG detect --help
usage: IG detect [-h] [--hom HOM] OUTDIR

positional arguments:
  OUTDIR      Directory for output

optional arguments:
  -h, --help  show this help message and exit
  --hom HOM   Add homozygous reference genotype
```

## Explanation of steps
### Phase
In the first step `--phase`, the subreads and CCS reads are phased and aligned to the IGH specific reference. Each read has a read group annotation. A read group annotation of 1 and 2 corresponds to haplotype 1 and 2. The read group annotation of 0 corresponds to unassignable reads. In IGV, you can seperate these reads by left clicking and selecting group by read group.

### Assemble
In the second step `--assembly`, the haplotypes are assembled. During this process folders will be created for each region/haplotype block. Within each folder there is a bash script that runs the assembly process. These can be submitted as a single job into the cluster (this speeds up the process).

### Detect
In the third step `--detect`, SNVs, indels, SVs and gene/alleles are genotyped. A VCF file is created for the SNVs, a BED file for the indels and SVs, a TAB-delimited file for the gene/alleles calls.  

## Output directories
alignments  alleles  assembly  logs  plots  preprocessed  report.html  tmp  variants
| Directories            | Description                                          |
|------------------------|------------------------------------------------------|
| `<output>/alignments`  | Alignments of CCS, subreads and contigs (phased and unphased) |
| `<output>/assembly`    | Assembly of IGH locus                                |
| `<output>/variants`    | SNVs, indels and SVs                                 |
| `<output>/alleles`     | Alleles in sample                                    |
| `<output>/logs`        | Log files with input parameters                   |
| `<output>/tmp`         | Temporary files. Could be deleted.                   |

## Output files
1. `alignments/`
    1. `ccs_to_ref*`: CCS reads aligned to reference
    2. `contigs_to_ref*`: All assembled contigs aligned to refernece
    3. `igh_contigs_to_ref*`: IGH assembled contigs aligned to igh reference
2. `assembly/`
    1. `contigs.fasta`: All assembled contigs
    2. `igh_contigs.fasta`: IGH assembled contigs
3. `alleles/`
    1. `assembly_alleles.bed`: Alleles extracted from the assembly for each gene
    2. `assembly_genes.fasta`: Fasta sequence from the assembly for each gene/allele
    3. `ccs_alleles.bed`: Alleles extracted from the CCS reads for each gene
    4. `ccs_genes.fasta`: Fasta sequence from the CCS reads for each gene/allele
4. `logs/`
    1. `gene_cov.txt`: Haplotype coverage for each gene
5. `variants/`
    1. `snvs_phased_from_ccs.vcf`: Phased SNVs detected from the CCS reads
    2. `snvs_assembly.vcf`: SNVS detected from the assembly
    3. `indel_assembly.bed`: Indels detected from the assembly
    4. `sv_assembly.bed`: SVs detected from the assembly
    5. `phased_blocks.txt`: Phased haplotype blocks
