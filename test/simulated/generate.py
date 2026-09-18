#!/usr/bin/env python3
"""Generate a deterministic diploid HiFi-like IGenotyper test dataset."""

import random
from pathlib import Path


OUTDIR = Path(__file__).resolve().parent
SAMPLE = "simulated"
SNPS = {500: "T", 800: "C", 1100: "A", 1400: "G"}


def wrap(sequence, width=80):
    return "\n".join(sequence[i:i + width] for i in range(0, len(sequence), width))


def main():
    random.seed(20260918)
    reference = "".join(random.choice("ACGT") for _ in range(2000))
    haplotype2 = list(reference)
    for position, alternate in SNPS.items():
        if alternate == reference[position - 1]:
            alternate = next(base for base in "ACGT" if base != alternate)
        haplotype2[position - 1] = alternate
    haplotype2 = "".join(haplotype2)

    (OUTDIR / "reference.fasta").write_text(">igh\n%s\n" % wrap(reference))
    with (OUTDIR / "reads.fasta").open("w") as handle:
        for haplotype, sequence in ((1, reference), (2, haplotype2)):
            for replicate in range(6):
                handle.write(
                    ">sim/%d/%d_2000 haplotype=%d\n%s\n"
                    % (haplotype, replicate, haplotype, wrap(sequence))
                )

    with (OUTDIR / "truth.vcf").open("w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##contig=<ID=igh,length=2000>\n")
        handle.write(
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        )
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % SAMPLE)
        for position, requested_alternate in sorted(SNPS.items()):
            ref = reference[position - 1]
            alt = requested_alternate
            if alt == ref:
                alt = next(base for base in "ACGT" if base != alt)
            handle.write(
                "igh\t%d\t.\t%s\t%s\t60\tPASS\t.\tGT\t0/1\n"
                % (position, ref, alt)
            )


if __name__ == "__main__":
    main()
