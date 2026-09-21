# Pipeline failure fixes and recovery

## WhatsHap diagnosis

The supported environment pins WhatsHap **2.8**. Synthetic regression tests run
its real candidate generation, genotyping and phasing code, without private data.

The previous FASTA conversion stripped `/fwd` and `/rev` and changed `/ccs` to
`/0_8`. Distinct reads such as `movie/48234718/ccs/fwd` and
`movie/48234718/ccs/rev` therefore acquired the same QNAME. WhatsHap 2.8 groups
by name, source and sample. Its `create_read_from_group` discards disagreeing
observations. When all observations disagree, the resulting read is empty but
remains in the ReadSet. `compute_genotypes` constructs a ColumnIterator that
calls `Read::firstPosition()` / `lastPosition()`, throwing `No variants present`.

Sources inspected at tag v2.8 (commit f07423a7ba0c7609e7d4d5c73fb5e4e240057830):

- [Read grouping](https://github.com/whatshap/whatshap/blob/v2.8/whatshap/variants.py)
- [Exception](https://github.com/whatshap/whatshap/blob/v2.8/src/read.cpp)
- [Column iteration](https://github.com/whatshap/whatshap/blob/v2.8/src/columniterator.cpp)
- [Genotyping](https://github.com/whatshap/whatshap/blob/v2.8/whatshap/cli/genotype.py)

Confirmed locally: a pair with opposite alleles at its sole observed variant
creates an empty grouped read and reproduces the exception. Nonempty synthetic
BAMs/VCFs also reproduce the CLI crash on both `chr1` and `igh`, even with other
informative reads. The same fixtures complete with the repository adapter.
The user independently reported this mechanism in one real sample (27,068
grouped reads, 8,856 candidates, one empty internal read), and the same exception
in 252 failures. Those private results were not independently examined here;
this is not evidence of a complete rerun of those samples.

Candidate generation uses coordinate-ordered pileups, MAPQ >=20, base quality
>=5, minimum ALT depth 3 and relative ALT depth 0.25 under `--pacbio`. Genotyping
ignores read groups as before. WhatsHap's candidate writer omits `##contig`
declarations; this produces warnings, but the unique-name controls genotype
successfully with that same formatting. It is not the cause of the reproduced
exception. IGenotyper now validates sample columns, position order, duplicate
positions, BAM/reference contig membership and REF bases, and verifies exact
site identity/order after genotyping and phasing. It does not sort away,
discard or silently skip problematic records/chromosomes.

FASTA conversion now preserves original QNAMEs with `samtools fasta -n`.
Downstream inspection found no CCS suffix parsing requirement: mapping keeps
FASTA IDs; phasing copies QNAMEs and changes RG tags; assembly extraction copies
QNAMEs into FASTA; subread extraction matches names exactly. Parsers for generated
assembly/gene IDs operate on separately generated IDs, not CCS `/ccs` suffixes.

`python -m IGenotyper.command_lines.whatshap_genotype genotype ...` wraps only
WhatsHap **2.8** in the current interpreter. Immediately before initial genotyping
it removes zero-observation grouped reads, logs their count, and passes every
original candidate position and every informative read to the original kernel.
Other errors propagate. With no informative reads left, it fails explicitly and
publishes no VCF; it does not invent confident genotypes. Unsupported versions
fail with an explicit version message. No site-packages files are modified, and
no later WhatsHap version is assumed to have fixed this bug.

## Completion and recovery

Run the updated repository in the environment from `environment.yml`. Rerun
`IG phase` with the original BAM, sample, output directory and reference options,
then `IG assembly` with that output directory. Do not bypass conversion/mapping
and invoke genotyping directly against an old, renamed BAM as the recovery path:
the adapter prevents the crash but cannot recover evidence already lost by
colliding names.

- External and Python output-producing commands stage files in a temporary
  directory beside the destination. A successful exit and validation of all
  outputs are required before publication. Pipelines use `errexit` and `pipefail`.
- Outputs are individually replaced atomically. A single `*.success.json`
  receipt is published last for the entire output set. It records the command
  (including genotype adapter/version), input and output sizes, nanosecond mtimes
  and ctimes. BAM/index pairs and the two R plot outputs are tracked together.
- Genotype/phased VCF validation compares every output site with the input set;
  a parseable truncated VCF is insufficient. BAMs are checked for intact EOF and
  associated indexes. Receipts are revalidated on reuse.
- Missing, stale or legacy receipts trigger rebuilding. Prior outputs remain
  untouched if execution or validation fails. An interruption between multiple
  file replacements leaves no valid receipt, so the whole step reruns. Concurrent
  writers to the same output directory are not supported.
- Original nonempty legacy files are **not** automatically trusted. The first
  retry regenerates FASTA, mappings, variants and phased BAMs; dependency receipts
  propagate invalidation. The old `args.json` and nonempty-VCF shortcuts no longer
  bypass this recovery. Subsequent successful outputs can be reused.
- Assembly regions with legacy/stale completion records are moved intact to
  sibling `<hap>.previous-<uuid>` directories before rebuilding, preserving older
  contigs for inspection. Keep these until recovery is verified. A fresh output
  directory is also an option if a complete old/new comparison is needed.

## Assembly and coverage behavior

Assembly dispatch now separates instrument (`@RG PM`) from read type
(`@RG DS:READTYPE`) and per-read predicted accuracy (`rq`). In particular,
`PM=SEQUEL; READTYPE=CCS` does not require a subreads BAM. The selector reads the
original input BAM because FASTA conversion/mapping can discard PacBio metadata.
It scans usable primary records once per assembly invocation; it does not sample
quality or discard lower-quality reads to qualify for HiFi mode.

| Input evidence | Canu flag | Polishing | Collected output |
| --- | --- | --- | --- |
| CCS, all usable primary reads have `rq >= 0.99` | `-pacbio-hifi` | None | `canu/canu.contigs.fasta` |
| CCS, any lower or missing `rq` | `-pacbio` (correction, trimming, assembly) | No external subreads required | `canu/canu.contigs.fasta` |
| Explicit SUBREAD | `-pacbio` | `pbindex`, `pbmm2 --preset SUBREAD`, `gcpp` | `contigs.fasta` |

The quality cutoff follows the [PacBio BAM specification](https://pacbiofileformats.readthedocs.io/en/13.1/BAM.html).
`READTYPE=SEGMENT;SOURCE=CCS` is treated as CCS. If READTYPE was stripped, canonical
`movie/hole/ccs[/fwd|rev][/start_end]` names can recover CCS identity. Missing or
lower `rq` is not evidence for HiFi; the correction route is the conservative
choice, consistent with [Canu's input modes](https://canu.readthedocs.io/en/latest/tutorial.html).
Ambiguous/unsupported or mixed CCS/SUBREAD inputs fail with an actionable message.
For incomplete demultiplexed headers, an undeclared `8-hex-ID/forward--reverse`
read group may use metadata from its explicitly declared unsuffixed parent. Exact
header matches take precedence. This recovery logs a warning once per recovered
ID and changes neither RG tags nor QNAMEs. Unknown parents, malformed suffixes,
and duplicate header IDs remain errors; this does not repair the BAM itself.
Instrument model alone does not choose a workflow. Base qualities and pass counts
are not substituted for `rq`.

Explicit SUBREAD assembly extracts from a supplied, indexed, phased subreads BAM
and uses the original subreads for polishing. It does not depend on a separate
CCS BAM. The current CCS-focused `phase` command does not create that legacy
phased subreads input automatically: supply results from a subread mapping/phasing
workflow. Missing BAMs/indexes and missing polishing executables are errors.
The extra polishing tools are not in the modern HiFi environment; the
[explicit SUBREAD preset](https://github.com/PacificBiosciences/pbmm2) avoids the
newer pbmm2 default of CCS.

Region selection partitions the target BED by phase blocks; it does not imply
coverage in every window/haplotype. Extraction validates the contig, coordinates,
haplotype groups and BAM index before checking local coverage. It uses the same
1-kb flanks and alignment-flag exclusions as before, with flanks clamped to the
reference bounds. A legitimate zero-read extraction creates an atomic
`skipped.json` with status `skipped_no_coverage`, **not** a successful `done` marker.
Malformed inputs, missing sequences and tool failures are not converted to skips.

Both assembled and skipped region receipts bind the read workflow, input BAMs
and index, code, region coordinates, flank and haplotype. A changed input, read
mode or old receipt causes the existing region directory to be archived before
rebuilding. A valid skip always takes precedence over stale contigs: they are
never collected, even if old files are restored into that directory. A skip can
be reused on retry while its provenance still matches; it is reconsidered when
inputs change.

Every generated script extracts reads afresh and runs Canu/polishing in a new
staging directory, preventing nonempty partial files from being reused even if
the script itself is retried. Nonempty FASTA validation and successful tool exits
are required before publication; completion is recorded last. Collection ignores
only validated no-coverage skips. Mixed covered/uncovered samples can finish;
all-empty samples fail with `No valid contigs overall` without publishing an
empty assembly or overwriting an older valid combined FASTA. An optional IGH-only
subset with no contigs is omitted (a stale subset FASTA is archived), provided
the overall assembly has valid contigs.

Tests use real synthetic BAM headers/records/indexes, actual extraction and
collection, and the generated shell scripts. Canu, pbindex, pbmm2 and gcpp are
**stubs**, so these tests verify dispatch, commands, publication and failure
handling, not assembly accuracy or compatibility on real datasets. No real
Canu/polishing run or rerun of the reported samples was possible. The conservative
non-HiFi CCS mode may be slower and still legitimately fail at insufficient
coverage; those failures remain explicit.

A haplotype BAM with zero mapped reads produces a real BigWig with explicit
zero-valued intervals spanning every BAM reference contig. It is neither a
zero-byte placeholder nor a missing-data track. Nonempty BAMs still run
`bamCoverage`, and unrelated errors propagate. Input BAMs with no primary reads
containing sequence fail before FASTA conversion; unmapped primary reads remain
valid mapping inputs.

## Synthetic tests

From the repository root in the supported environment:

```sh
python -m pytest test -q
```

Tests require pysam, Biopython, pybedtools, PyVCF3, pyBigWig, pytest and WhatsHap
2.8. The `whatshap` executable must be on PATH from that same environment.
The local verification used Python 3.14 and WhatsHap 2.8, while the Conda
production environment pins Python 3.11. Full production phasing/assembly and
plot rendering are outside these small fixture tests.
