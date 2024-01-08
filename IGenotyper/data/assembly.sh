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
#pacbio_machine=${12}	ud_edited
pacbio_machine="SEQUELII"
# Get CCS reads
if [ ! -s ${output}/reads.fasta ]
then
    samtools view -F 3844 ${ccs_to_ref} ${hap} ${chrom}:${start}-${end} \
    | awk '{ print ">"$1"\n"$10}' > ${output}/reads.fasta
fi

if [ ${pacbio_machine} != "SEQUELII" ]
then
    data_setting="-pacbio"
else
    data_setting="-pacbio-hifi"
fi

# Assemble with CCS reads
#  contigFilter="2 0 1.0 0.5 0"
if [ ! -s ${output}/canu/canu.contigs.fasta ]
then
    rm -fr ${output}/canu
    canu \
      -p canu \
    	-d ${output}/canu \
    	corOutCoverage=200 \
    	minThreads=${threads} \
    	genomeSize=${size} \
    	useGrid=0 \
	minInputCoverage=0 \
        stopOnLowCoverage=0 \
    	${data_setting} ${output}/reads.fasta
fi

echo "" > ${output}/done
