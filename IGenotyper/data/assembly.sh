# Variables are supplied by assembly/scripts.py. Run the generated script with bash.
# Never reuse an old extraction, partial Canu output, or completion marker.
rm -f "${output}/done" "${output}/skipped.json"
count=$("${python}" -m IGenotyper.assembly.regions "${assembly_bam}" \
    "${chrom}" "${start}" "${end}" "${hap}" "${output}/reads.fasta")
if [ "${count}" -eq 0 ]; then
    "${python}" -c 'import json, sys; from IGenotyper.assembly.scripts import record_region_result; record_region_result(sys.argv[1], json.loads(sys.argv[2]), "skipped_no_coverage")' "${output}" "${completion}"
    echo "Skipped assembly (no coverage): ${region}, haplotype ${hap}"
    exit 0
fi

if [ "${polish}" -eq 1 ]; then
    for tool in pbindex pbmm2 gcpp; do
        if ! command -v "${tool}" >/dev/null; then
            echo "SUBREAD assembly requires ${tool} for polishing; install legacy PacBio polishing tools." >&2
            exit 1
        fi
    done
fi

work=$(mktemp -d "${output}/.assembly-XXXXXX")
trap 'rm -rf "${work}"' EXIT
canu -p canu -d "${work}/canu" corOutCoverage=200 \
    minThreads="${threads}" genomeSize="${size}" useGrid=0 \
    minInputCoverage=0 stopOnLowCoverage=0 \
    "${data_setting}" "${output}/reads.fasta"
"${python}" -c 'import sys; from IGenotyper.common.validation import fasta_records; fasta_records(sys.argv[1])' "${work}/canu/canu.contigs.fasta"

if [ "${polish}" -eq 1 ]; then
    awk 'NR % 2 == 1 {print substr($0, 2)}' "${output}/reads.fasta" > "${work}/reads.names"
    "${python}" "${python_scripts}/extract_reads.py" -b "${subreads}" \
        -n "${work}/reads.names" -o "${work}/subreads.bam"
    "${python}" -c 'import sys; from IGenotyper.common.validation import require_usable_reads; require_usable_reads(sys.argv[1])' "${work}/subreads.bam"
    pbindex "${work}/subreads.bam"
    pbmm2 align --preset SUBREAD --sort -J "${threads}" -j "${threads}" \
        "${work}/canu/canu.contigs.fasta" "${work}/subreads.bam" \
        "${work}/canu/reads_to_canu_contigs.sorted.bam"
    pbindex "${work}/canu/reads_to_canu_contigs.sorted.bam"
    samtools faidx "${work}/canu/canu.contigs.fasta"
    gcpp --reference "${work}/canu/canu.contigs.fasta" -j "${threads}" \
        -o "${work}/contigs.fasta" "${work}/canu/reads_to_canu_contigs.sorted.bam"
    "${python}" -c 'import sys; from IGenotyper.common.validation import fasta_records; fasta_records(sys.argv[1])' "${work}/contigs.fasta"
fi
# Publish only validated tool outputs; the completion receipt is written last.
"${python}" -c 'import json, sys; from IGenotyper.assembly.scripts import validate_assembly_inputs; validate_assembly_inputs(json.loads(sys.argv[1]))' "${completion}"
if [ -d "${output}/canu" ]; then
    mv "${output}/canu" "${work}/previous-canu"
fi
mv "${work}/canu" "${output}/canu"
if [ "${polish}" -eq 1 ]; then
    mv "${work}/contigs.fasta" "${contigs}"
fi
"${python}" -c 'import json, sys; from IGenotyper.assembly.scripts import record_region_result; record_region_result(sys.argv[1], json.loads(sys.argv[2]), "assembled", sys.argv[3])' "${output}" "${completion}" "${contigs}"
