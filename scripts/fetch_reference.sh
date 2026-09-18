#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
default_root="${XDG_DATA_HOME:-${HOME}/.local/share}/igenotyper"
data_dir="${IGENOTYPER_DATA_DIR:-${default_root}}"
if [[ "${1:-}" == "--data-dir" ]]; then
    [[ -n "${2:-}" ]] || { echo "--data-dir requires a path" >&2; exit 2; }
    data_dir="${2}"
elif [[ $# -ne 0 ]]; then
    echo "Usage: $0 [--data-dir PATH]" >&2
    exit 2
fi
source_dir="${repo_dir}/IGenotyper/data/immune_receptor_genomics/251106"
url=$(tr -d '[:space:]' < "${source_dir}/reference_path.txt")
output="${data_dir}/reference.fasta"
partial="${output}.part"

mkdir -p "${data_dir}"

if [[ -s "${output}" ]] && cmp -s "${output}.fai" "${source_dir}/reference.fasta.fai"; then
    if [[ ! -s "${output}.mmi" ]]; then
        echo "Building reusable minimap2 index ${output}.mmi"
        minimap2 -x map-hifi -d "${output}.mmi.part" "${output}"
        mv "${output}.mmi.part" "${output}.mmi"
    fi
    echo "Reference is already installed and validated: ${output}"
    exit 0
fi

echo "Downloading the approximately 3.1 GB reference from ${url}"
curl --fail --location --continue-at - --output "${partial}" "${url}"
samtools faidx "${partial}"

if ! cmp -s "${partial}.fai" "${source_dir}/reference.fasta.fai"; then
    echo "Downloaded FASTA does not match the upstream FASTA index" >&2
    exit 1
fi

mv "${partial}" "${output}"
mv "${partial}.fai" "${output}.fai"
echo "Building reusable minimap2 index ${output}.mmi"
minimap2 -x map-hifi -d "${output}.mmi.part" "${output}"
mv "${output}.mmi.part" "${output}.mmi"
echo "Installed and validated ${output}"
