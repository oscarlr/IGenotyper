# IGenotyper

[Introduction](#introduction)  
[Installation](#installation)  
[Getting IGH specific reference](#getting-igh-specific-reference)<br>
[Testing IGenotyper installation](#testing-igenotyper-installation)<br>
[Running IGenotyper](#running-igenotyper)<br>
[Explanation of steps](#explanation-of-steps)<br>
[Output directories](#output-directories)<br>
[Output files](#output-files)<br>
[Notes](#notes)

## Installation

```
git clone https://github.com/oscarlr/IGv2.git
cd IGv2
conda env create -f whatshap.yml
conda env create -f pygenometracks.yml

conda env create -f environment.yml
conda activate IGv2
python setup.py install

cd ..
git clone https://github.com/oscarlr/cluster
cd cluster
python setup.py install
export SJOB_DEFALLOC=NONE

wget https://watsonlab.blob.core.windows.net/public/reference/hg38_igh_trab/reference.fasta
wget https://watsonlab.blob.core.windows.net/public/reference/hg38_igh_trab/reference.fasta.fai
```

## Getting IGH specific reference or download IGH specific reference

```
sawriter reference.fasta
cp reference.fasta* ~/anaconda3/envs/IGv2_testing/lib/python2.7/site-packages/IGenotyper-1.1-py2.7.egg/IGenotyper/data/

# For rheMac10
sawriter reference.fasta
cp reference.fasta* ~/anaconda3/envs/IGv2_testing/lib/python2.7/site-packages/IGenotyper-1.1-py2.7.egg/IGenotyper/data/rhesus/

```

## Testing IGenotyper installation
```
Todo
```

## Running IGenotyper
```
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
  
IG detect --help
usage: IG detect [-h] [--hom HOM] OUTDIR

positional arguments:
  OUTDIR      Directory for output

optional arguments:
  -h, --help  show this help message and exit
  --hom HOM   Add homozygous reference genotype
```
