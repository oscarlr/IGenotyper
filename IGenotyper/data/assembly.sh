#!/bin/bash
set -x #-e

output=$1
ccs_to_ref=$2
hap=$3
chrom=$4
start=$5
end=$6
threads=$7
size=$8
subreads_to_ref=$9
python_scripts=${10}
subreads=${11}

# Get CCS reads
if [ ! -s ${output}/reads.fasta ]
then
    samtools view -F 3844 ${ccs_to_ref} ${hap} ${chrom}:${start}-${end} \
    | awk '{ print ">"$1"\n"$10}' > ${output}/reads.fasta
fi

# Assemble with CCS reads
#  contigFilter="2 0 1.0 0.5 0"
if [ ! -s ${output}/canu/canu.contigs.fasta ]
then
    canu \
      -p canu \
    	-d ${output}/canu \
    	corOutCoverage=200 \
    	minThreads=${threads} \
    	genomeSize=${size} \
    	useGrid=0 \
    	-pacbio-raw ${output}/reads.fasta
fi

# Asseemble with subreads
# contigFilter="2 0 1.0 0.5 0"
if [ ! -s ${output}/canu/canu.contigs.fasta ]
then
    rm -fr ${output}/canu
    samtools view -F 3844 ${subreads_to_ref} ${hap} ${chrom}:${start}-${end} \
    | awk '{ print ">"$1"\n"$10}' > ${output}/reads.fasta
    canu \
    -p canu \
    -d ${output}/canu \
    minThreads=${threads} \
    genomeSize=${size} \
    useGrid=0 \
    -pacbio-raw ${output}/reads.fasta
fi

# Error correct contigs
if [ -s ${output}/canu/canu.contigs.fasta ]
then
    if [ ! -s ${output}/subreads.bam.pbi ]
    then
        samtools view -F 3844 ${subreads_to_ref} ${hap} ${chrom}:${start}-${end} \
        | awk '{ print $1 }' | sort | uniq > ${output}/reads.names
        python ${python_scripts}/extract_reads.py \
        -b ${subreads} \
        -n ${output}/reads.names \
        -o ${output}/subreads.bam
        pbindex ${output}/subreads.bam
    fi
    if [ ! -s ${output}/canu/reads_to_canu_contigs.sorted.bam.pbi ]
    then
        pbmm2 align \
        --sort \
        -J ${threads} \
        -j ${threads} \
        ${output}/canu/canu.contigs.fasta \
        ${output}/subreads.bam \
        ${output}/canu/reads_to_canu_contigs.sorted.bam
        pbindex ${output}/canu/reads_to_canu_contigs.sorted.bam
    fi
    if [ ! -s ${output}/contigs.fastq ]
    then
        samtools faidx ${output}/canu/canu.contigs.fasta
        gcpp \
        --reference ${output}/canu/canu.contigs.fasta \
        -j ${threads} \
        -o ${output}/contigs.fasta \
        ${output}/canu/reads_to_canu_contigs.sorted.bam
        gcpp \
        --reference ${output}/canu/canu.contigs.fasta \
        -j ${threads} \
        -o ${output}/contigs.fastq \
        ${output}/canu/reads_to_canu_contigs.sorted.bam
    fi
fi

echo "" > ${output}/done