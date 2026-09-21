# Variables are supplied by assembly/scripts.py. Run the generated script with bash.
# An old done marker must not survive a failed rebuild.
rm -f "${output}/done"
if [ ! -s "${output}/reads.fasta" ]; then
    samtools view -F 3844 "${ccs_to_ref}" -r "${hap}" "${region}" \
        | awk '$10 != "*" { print ">"$1"\n"$10 }' > "${output}/reads.fasta.tmp"
    if [ ! -s "${output}/reads.fasta.tmp" ]; then
        echo "Zero usable reads for assembly: ${region}, haplotype ${hap}" >&2
        exit 1
    fi
    mv "${output}/reads.fasta.tmp" "${output}/reads.fasta"
fi

case "${pacbio_machine}" in
    SEQUEL) data_setting="-pacbio" ;;
    SEQUELII|REVIO) data_setting="-pacbio-hifi" ;;
    *) echo "Unsupported PacBio platform: ${pacbio_machine}" >&2; exit 1 ;;
esac

if [ ! -s "${output}/canu/canu.contigs.fasta" ]; then
    canu -p canu -d "${output}/canu" corOutCoverage=200 \
        minThreads="${threads}" genomeSize="${size}" useGrid=0 \
        minInputCoverage=0 stopOnLowCoverage=0 \
        "${data_setting}" "${output}/reads.fasta"
fi
if [ ! -s "${output}/canu/canu.contigs.fasta" ]; then
    echo "Missing or empty Canu contigs: ${output}/canu/canu.contigs.fasta" >&2
    exit 1
fi

# Older SEQUEL data require subread polishing. Never silently use unpolished
# Canu contigs in place of this platform's expected final output.
if [ "${pacbio_machine}" = SEQUEL ] && [ ! -s "${contigs}" ]; then
    samtools view -F 3844 "${subreads_to_ref}" -r "${hap}" "${region}" \
        | awk '{ print $1 }' | sort -u > "${output}/reads.names"
    python "${python_scripts}/extract_reads.py" -b "${subreads}" \
        -n "${output}/reads.names" -o "${output}/subreads.bam"
    pbindex "${output}/subreads.bam"
    pbmm2 align --sort -J "${threads}" -j "${threads}" \
        "${output}/canu/canu.contigs.fasta" "${output}/subreads.bam" \
        "${output}/canu/reads_to_canu_contigs.sorted.bam"
    pbindex "${output}/canu/reads_to_canu_contigs.sorted.bam"
    samtools faidx "${output}/canu/canu.contigs.fasta"
    gcpp --reference "${output}/canu/canu.contigs.fasta" -j "${threads}" \
        -o "${output}/contigs.tmp.fasta" "${output}/canu/reads_to_canu_contigs.sorted.bam"
    python -c 'import sys; from IGenotyper.common.validation import fasta_records; fasta_records(sys.argv[1])' "${output}/contigs.tmp.fasta"
    mv "${output}/contigs.tmp.fasta" "${contigs}"
fi
python -c 'import sys; from IGenotyper.common.validation import fasta_records; fasta_records(sys.argv[1])' "${contigs}"
python -c 'import json, sys; from IGenotyper.assembly.scripts import record_assembly_success; record_assembly_success(sys.argv[1], sys.argv[2], json.loads(sys.argv[3]))' "${output}" "${pacbio_machine}" "${completion}"
