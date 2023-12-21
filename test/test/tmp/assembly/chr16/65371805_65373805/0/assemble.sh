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
pacbio_machine=${12}

# Get CCS reads
if [ ! -s test/tmp/assembly/chr16/65371805_65373805/0/reads.fasta ]
then
    samtools view -F 3844 test/alignments/ccs_to_ref_phased.sorted.bam -r 0 chr16:65370805-65374805 \
    | awk '{ print ">"$1"\n"$10}' > test/tmp/assembly/chr16/65371805_65373805/0/reads.fasta
fi

if [ SEQUELII != "SEQUELII" ]
then
    data_setting="-pacbio"
else
    data_setting="-pacbio-hifi"
fi

# Assemble with CCS reads
#  contigFilter="2 0 1.0 0.5 0"
if [ ! -s test/tmp/assembly/chr16/65371805_65373805/0/canu/canu.contigs.fasta ]
then
    rm -fr test/tmp/assembly/chr16/65371805_65373805/0/canu
    canu \
      -p canu \
    	-d test/tmp/assembly/chr16/65371805_65373805/0/canu \
    	corOutCoverage=200 \
    	minThreads=1 \
    	genomeSize=4000 \
    	useGrid=0 \
	minInputCoverage=0 \
        stopOnLowCoverage=0 \
    	${data_setting} test/tmp/assembly/chr16/65371805_65373805/0/reads.fasta
fi

# Asseemble with subreads
# contigFilter="2 0 1.0 0.5 0"
# if [ ! -s test/tmp/assembly/chr16/65371805_65373805/0/canu/canu.contigs.fasta ]
# then
#     rm -fr test/tmp/assembly/chr16/65371805_65373805/0/canu
#     samtools view -F 3844 test/alignments/subreads_to_ref_phased.sorted.bam -r 0 chr16:65370805-65374805 \
#     | awk '{ print ">"$1"\n"$10}' > test/tmp/assembly/chr16/65371805_65373805/0/reads.fasta
#     canu \
#     -p canu \
#     -d test/tmp/assembly/chr16/65371805_65373805/0/canu \
#     minThreads=1 \
#     genomeSize=4000 \
#     useGrid=0 \
#     -pacbio-raw test/tmp/assembly/chr16/65371805_65373805/0/reads.fasta
# fi

if [ SEQUELII != "SEQUELII" ]
then
    # Error correct contigs
    if [ -s test/tmp/assembly/chr16/65371805_65373805/0/canu/canu.contigs.fasta ]
    then
	if [ ! -s test/tmp/assembly/chr16/65371805_65373805/0/subreads.bam.pbi ]
	then
            samtools view -F 3844 test/alignments/subreads_to_ref_phased.sorted.bam -r 0 chr16:65370805-65374805 \
		| awk '{ print $1 }' | sort | uniq > test/tmp/assembly/chr16/65371805_65373805/0/reads.names
            python /home/u0jana01/anaconda3/envs/IGv2/lib/python2.7/site-packages/IGenotyper-1.1-py2.7.egg/IGenotyper/scripts/extract_reads.py \
		-b subset_bam/reads.bam \
		-n test/tmp/assembly/chr16/65371805_65373805/0/reads.names \
		-o test/tmp/assembly/chr16/65371805_65373805/0/subreads.bam
            pbindex test/tmp/assembly/chr16/65371805_65373805/0/subreads.bam
	fi
	if [ ! -s test/tmp/assembly/chr16/65371805_65373805/0/canu/reads_to_canu_contigs.sorted.bam.pbi ]
	then
            pbmm2 align \
		--sort \
		-J 1 \
		-j 1 \
		test/tmp/assembly/chr16/65371805_65373805/0/canu/canu.contigs.fasta \
		test/tmp/assembly/chr16/65371805_65373805/0/subreads.bam \
		test/tmp/assembly/chr16/65371805_65373805/0/canu/reads_to_canu_contigs.sorted.bam
            pbindex test/tmp/assembly/chr16/65371805_65373805/0/canu/reads_to_canu_contigs.sorted.bam
	fi
	if [ ! -s test/tmp/assembly/chr16/65371805_65373805/0/contigs.fastq ]
	then
            samtools faidx test/tmp/assembly/chr16/65371805_65373805/0/canu/canu.contigs.fasta
            gcpp \
		--reference test/tmp/assembly/chr16/65371805_65373805/0/canu/canu.contigs.fasta \
		-j 1 \
		-o test/tmp/assembly/chr16/65371805_65373805/0/contigs.fasta \
		test/tmp/assembly/chr16/65371805_65373805/0/canu/reads_to_canu_contigs.sorted.bam
            gcpp \
		--reference test/tmp/assembly/chr16/65371805_65373805/0/canu/canu.contigs.fasta \
		-j 1 \
		-o test/tmp/assembly/chr16/65371805_65373805/0/contigs.fastq \
		test/tmp/assembly/chr16/65371805_65373805/0/canu/reads_to_canu_contigs.sorted.bam
	    if [ ! -s test/tmp/assembly/chr16/65371805_65373805/0/contigs.fasta ]
	    then
		variantCaller \
		    --referenceFilename test/tmp/assembly/chr16/65371805_65373805/0/canu/canu.contigs.fasta \
		    -j 1 \
		    -o test/tmp/assembly/chr16/65371805_65373805/0/contigs.fastq \
		    -o test/tmp/assembly/chr16/65371805_65373805/0/contigs.fasta \
		    test/tmp/assembly/chr16/65371805_65373805/0/canu/reads_to_canu_contigs.sorted.bam
	    fi
	fi
    fi
fi

echo "" > test/tmp/assembly/chr16/65371805_65373805/0/done
