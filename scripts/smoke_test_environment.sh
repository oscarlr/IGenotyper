#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
work=$(mktemp -d)
trap 'rm -rf -- "${work}"' EXIT
cd "${work}"

python - <<'PY'
import Bio
import matplotlib
import networkx
import numpy
import pandas
import pybedtools
import pysam
import reportlab
import vcf
import whatshap
print("Python imports: OK")
PY

cp "${repo_dir}/test/simulated/reference.fasta" reference.fasta
cp "${repo_dir}/test/simulated/reads.fasta" reads.fasta

samtools faidx reference.fasta
minimap2 -a -x map-hifi reference.fasta reads.fasta > mapped.sam
samtools view -b mapped.sam | samtools sort -o mapped.bam
samtools index mapped.bam
samtools quickcheck mapped.bam

cat > variants.vcf <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=igh,length=2000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample
igh	500	.	A	T	60	PASS	.	GT:PS	0|1:500
EOF
bcftools view variants.vcf -Ov -o viewed.vcf
whatshap --help > whatshap-help.txt

printf 'igh\t0\t10\ta\n' > a.bed
printf 'igh\t5\t15\tb\n' > b.bed
bedtools intersect -a a.bed -b b.bed > overlap.bed
test -s overlap.bed

blastn -query reads.fasta -subject reference.fasta -outfmt 6 > blast.tsv
test -s blast.tsv
cat > alignment.fasta <<'EOF'
>one
ACGTACGT
>two
ACGTTCGT
EOF
kalign -i alignment.fasta -o alignment.aln.fasta
test -s alignment.aln.fasta

bamCoverage -b mapped.bam -o coverage.bw --binSize 10 >/dev/null
test -s coverage.bw

cat > tracks.ini <<'EOF'
[test]
file = a.bed
title = test
height = 1
EOF
pyGenomeTracks --tracks tracks.ini --region igh:1-2000 -o tracks.png >/dev/null
test -s tracks.png

Rscript -e 'stopifnot(requireNamespace("ggplot2", quietly=TRUE)); stopifnot(requireNamespace("reshape2", quietly=TRUE))'
canu -version >/dev/null
IG phase --help >/dev/null

printf '%s\n' \
    "minimap2 $(minimap2 --version)" \
    "$(samtools --version | head -1)" \
    "$(bcftools --version | head -1)" \
    "WhatsHap $(whatshap --version)" \
    "$(canu -version 2>&1 | head -1)" \
    "$(kalign --version 2>&1 | head -1)" \
    "Environment smoke test: PASS"
