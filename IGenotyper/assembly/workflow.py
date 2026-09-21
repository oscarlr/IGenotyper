"""Choose assembly from biological read type and accuracy, not instrument model."""
from dataclasses import dataclass
import math
import re

import pysam


@dataclass(frozen=True)
class AssemblyWorkflow:
    name: str
    canu_flag: str
    polish: bool
    platforms: tuple
    read_type: str
    accuracy: str


def select_assembly_workflow(path):
    """Require CCS/SUBREAD metadata; recover stripped CCS headers from QNAMEs.

    HiFi requires rq >= 0.99 on every usable primary CCS record. Missing/lower
    accuracy uses Canu's correction pipeline, without requiring external subreads.
    No read is discarded to achieve the HiFi threshold.
    """
    with pysam.AlignmentFile(path, 'rb', check_sq=False) as bam:
        groups = bam.header.to_dict().get('RG', [])
        types = {}
        for group in groups:
            metadata = {}
            for item in group.get('DS', '').split(';'):
                key, separator, value = item.partition('=')
                if separator:
                    key, value = key.strip().upper(), value.strip().upper()
                    if key in metadata and metadata[key] != value:
                        raise ValueError('Conflicting %s metadata in BAM %s' % (key, path))
                    metadata[key] = value
            read_type = metadata.get('READTYPE')
            if read_type == 'SEGMENT' and metadata.get('SOURCE') == 'CCS':
                read_type = 'CCS'
            if read_type not in (None, 'CCS', 'SUBREAD'):
                raise ValueError('Unsupported READTYPE=%s in %s; supply CCS or SUBREAD input' % (read_type, path))
            types[group['ID']] = read_type
        declared = {value for value in types.values() if value is not None}
        if len(declared) > 1:
            raise ValueError('Mixed CCS/SUBREAD BAM; split read types before assembly: %s' % path)
        observed = set()
        count, all_hifi = 0, True
        for read in bam.fetch(until_eof=True):
            if read.flag & 0x900 or not read.query_sequence:
                continue
            if read.has_tag('RG'):
                group_id = read.get_tag('RG')
                if group_id not in types:
                    raise ValueError('Read refers to unknown RG %s in %s' % (group_id, path))
                read_type = types[group_id]
            else:
                read_type = next(iter(declared)) if len(declared) == 1 else None
            if read_type is None:
                if re.fullmatch(r'[^/]+/\d+/ccs(?:/(?:fwd|rev))?(?:/\d+_\d+)?', read.query_name):
                    read_type = 'CCS'
                else:
                    raise ValueError('Cannot determine READTYPE in %s; restore PacBio @RG DS READTYPE=CCS or READTYPE=SUBREAD metadata' % path)
            observed.add(read_type)
            count += 1
            if not read.has_tag('rq'):
                all_hifi = False
            else:
                quality = read.get_tag('rq')
                if not isinstance(quality, (float, int)) or not math.isfinite(quality) or not 0 <= quality <= 1:
                    raise ValueError('Invalid rq predicted accuracy for read %s in %s' % (read.query_name, path))
                all_hifi &= quality >= 0.99
        if count == 0:
            raise ValueError('BAM contains zero usable primary reads with sequence: %s' % path)
        if len(observed | declared) != 1:
            raise ValueError('Mixed CCS/SUBREAD BAM; split read types before assembly: %s' % path)
        read_type = observed.pop()
        platforms = tuple(sorted({group.get('PM', 'unknown') for group in groups}))
    if read_type == 'SUBREAD':
        return AssemblyWorkflow('subreads', '-pacbio', True, platforms, read_type, 'raw')
    if all_hifi:
        return AssemblyWorkflow('hifi', '-pacbio-hifi', False, platforms, read_type, 'all-rq-at-least-0.99')
    return AssemblyWorkflow('ccs', '-pacbio', False, platforms, read_type, 'lower-or-unknown-rq')
