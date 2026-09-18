#!/usr/bin/env bash
set -euo pipefail

test_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
work=$(mktemp -d)
trap 'rm -rf -- "${work}"' EXIT

cp "${test_dir}/reference.fasta" "${test_dir}/reads.fasta" "${test_dir}/truth.vcf" "${work}/"
cd "${work}"

samtools faidx reference.fasta
minimap2 -a -x map-hifi -R '@RG\tID:simulated\tSM:simulated' \
    reference.fasta reads.fasta | samtools sort -o reads.bam
samtools index reads.bam

canu -p simulated -d canu genomeSize=2k useGrid=false maxThreads=4 \
    minInputCoverage=0 stopOnLowCoverage=0 -pacbio-hifi reads.fasta >/dev/null 2>&1
test -s canu/simulated.unassembled.fasta

whatshap find_snv_candidates --pacbio --sample simulated \
    -o candidates.vcf reference.fasta reads.bam
whatshap genotype --sample simulated --ignore-read-groups \
    --reference reference.fasta -o genotyped.vcf candidates.vcf reads.bam
whatshap phase --sample simulated --ignore-read-groups \
    --reference reference.fasta -o phased.vcf genotyped.vcf reads.bam

expected=$(bcftools query -f '%POS\n' truth.vcf | sort -n | tr '\n' ',')
observed=$(bcftools query -f '%POS\n' phased.vcf | sort -n | tr '\n' ',')
if [[ "${observed}" != "${expected}" ]]; then
    echo "Expected variant positions ${expected}; observed ${observed}" >&2
    exit 1
fi

phased=$(bcftools query -f '[%GT\n]' phased.vcf | grep -c '|')
test "${phased}" -eq 4
echo "Simulated mapping, genotyping, and phasing test: PASS"
