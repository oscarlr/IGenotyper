# IGv2
## Installation

```
git clone https://github.com/oscarlr/IGv2.git
cd IGv2
conda env create -f whatshap.yml

conda env create -f environment.yml
conda activate IGv2_testing
python setup.py install

cd ..
git clone https://github.com/oscarlr/cluster
cd cluster
python setup.py install
export SJOB_DEFALLOC=NONE

wget https://watsonlab.blob.core.windows.net/public/reference/hg38_igh_trab/reference.fasta
wget https://watsonlab.blob.core.windows.net/public/reference/hg38_igh_trab/reference.fasta.fai

sawriter reference.fasta
cp reference.fasta* 


```
