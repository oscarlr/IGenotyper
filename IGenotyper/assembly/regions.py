"""Extract a validated assembly window; zero local coverage is not an error."""
import os
from pathlib import Path
import tempfile

import pysam


def extract_region_reads(bam_path, chrom, start, end, hap, output):
    start, end = int(start), int(end)
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        bam.check_index()
        if chrom not in bam.references:
            raise ValueError('Assembly contig %s is absent from %s' % (chrom, bam_path))
        length = bam.get_reference_length(chrom)
        if not 0 <= start < end <= length:
            raise ValueError('Invalid assembly coordinates %s:%s-%s (length %s)' % (chrom, start, end, length))
        groups = {str(group['ID']) for group in bam.header.to_dict().get('RG', [])}
        if str(hap) not in groups:
            raise ValueError('Missing haplotype RG=%s in phased BAM %s' % (hap, bam_path))
        fd, temporary = tempfile.mkstemp(prefix='.reads-', dir=Path(output).parent)
        count = 0
        try:
            with os.fdopen(fd, 'w') as stream:
                for read in bam.fetch(chrom, max(0, start - 1000), min(length, end + 1000)):
                    if read.flag & 3844:
                        continue
                    if not read.has_tag('RG') or read.get_tag('RG') not in groups:
                        raise ValueError('Missing or unknown haplotype RG for read %s' % read.query_name)
                    if read.get_tag('RG') != str(hap):
                        continue
                    if not read.query_sequence:
                        raise ValueError('Missing sequence for assembly read %s' % read.query_name)
                    stream.write('>%s\n%s\n' % (read.query_name, read.query_sequence))
                    count += 1
            os.replace(temporary, output)
        finally:
            Path(temporary).unlink(missing_ok=True)
    return count


if __name__ == '__main__':
    import sys
    print(extract_region_reads(*sys.argv[1:]))
