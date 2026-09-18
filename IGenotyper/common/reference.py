#!/usr/bin/env python3
import os


def read_fasta_index(path):
    lengths = {}
    with open(path) as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError("Invalid FASTA index line %s in %s" % (line_number, path))
            lengths[fields[0]] = int(fields[1])
    return lengths


def validate_bed(path, contig_lengths):
    errors = []
    with open(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                errors.append("line %s has fewer than three columns" % line_number)
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                errors.append("line %s has non-integer coordinates" % line_number)
                continue
            length = contig_lengths.get(fields[0])
            if length is None:
                errors.append("line %s uses unknown contig %s" % (line_number, fields[0]))
            elif start < 0 or end <= start or end > length:
                errors.append(
                    "line %s has invalid interval %s:%s-%s (length %s)"
                    % (line_number, fields[0], start, end, length)
                )
    if errors:
        raise ValueError("Invalid BED file %s:\n%s" % (path, "\n".join(errors)))


def validate_reference(reference, expected_index, bed_files):
    index = "%s.fai" % reference
    if not os.path.isfile(reference):
        raise FileNotFoundError(
            "Reference FASTA not found at %s. Run scripts/fetch_reference.sh "
            "or pass --data-dir/IGENOTYPER_DATA_DIR." % reference
        )
    if not os.path.isfile(index):
        raise FileNotFoundError("Reference FASTA index not found at %s" % index)

    observed = read_fasta_index(index)
    expected = read_fasta_index(expected_index)
    if observed != expected:
        raise ValueError(
            "Reference index %s does not match the bundled annotation version" % index
        )
    for bed_file in bed_files:
        validate_bed(bed_file, observed)
