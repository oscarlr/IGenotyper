# IGenotyper

[Introduction](#introduction)  
[Requirements](#requirements)  
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
