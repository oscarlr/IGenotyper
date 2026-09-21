# IGenotyper

[Introduction](#introduction)  
[Installation](#installation)  
[Reference data](#reference-data)<br>
[Testing IGenotyper installation](#testing-igenotyper-installation)<br>
[Usage](#usage)<br>
[Running IGenotyper](#running-igenotyper)<br>
[Explanation of steps](#explanation-of-steps)<br>
[Output directories](#output-directories)<br>
[Output files](#output-files)<br>
[Notes](#notes)

## Installation

```bash
git clone https://github.com/oscarlr/IGenotyper.git
cd IGenotyper
CONDA_SAFETY_CHECKS=enabled conda env create -f environment.yml
conda activate igenotyper
python -m pip install .
./scripts/fetch_reference.sh
```

Conda safety checks stop installation if cached package files are damaged.
If Conda reports a `SafetyError`, retry with a fresh package cache rather than
disabling the checks:

```bash
CONDA_PKGS_DIRS="$(mktemp -d "${TMPDIR:-/tmp}/igenotyper-conda.XXXXXX")" \
  CONDA_SAFETY_CHECKS=enabled conda env create -f environment.yml
```

On shared filesystems, creating the environment can take several minutes after
downloads finish. Wait for Conda to finish before installing the Python package
or running the tests below.

### Optional LSF integration

The `--cluster` assembly backend uses LSF and requires access to the separate
`Watson-IG/cluster` repository. Skip this section for local execution or Slurm.

```bash
cd ..
git clone https://github.com/Watson-IG/cluster.git
cd cluster
python -m pip install .
export SJOB_DEFALLOC=NONE
```

### Slurm

On Slurm, request a compute allocation and run IGenotyper without `--cluster`.
For example, using your site's partition name and resource limits:

```bash
conda activate igenotyper
srun --partition=compute --ntasks=1 --cpus-per-task=8 --mem=24G --time=02:00:00 \
  IG phase --sample SAMPLE --threads 8 reads.bam output
```

## Reference data

```bash
./scripts/fetch_reference.sh
```

The script downloads the current reference to
`${XDG_DATA_HOME:-$HOME/.local/share}/igenotyper`, validates it against the
FASTA index published by Watson-IG/immune_receptor_genomics, and builds a
reusable minimap2 index. To use another location, run
`./scripts/fetch_reference.sh --data-dir /path/to/data` and pass the same path
to `IG phase --data-dir /path/to/data` (or set `IGENOTYPER_DATA_DIR`). Upstream
versioned annotations remain bundled with the package so they always match the
selected reference release.

The download is approximately 3.1 GB. The installer also creates
`reference.fasta.fai` and `reference.fasta.mmi`. Existing valid downloads are
reused, and an absent minimap2 index is built without downloading the FASTA
again.

Custom location example:

```bash
./scripts/fetch_reference.sh --data-dir /data/igenotyper
export IGENOTYPER_DATA_DIR=/data/igenotyper
```

The command-line `--data-dir` option takes precedence over the environment
variable. IGenotyper validates the FASTA index and all selected BED coordinates
before starting an analysis.

## Testing IGenotyper installation

Run the tool and Python regression tests:

```bash
scripts/smoke_test_environment.sh
python -m unittest discover -s test -p 'test_*.py' -v
```

Run the deterministic simulated minimap2/WhatsHap integration test:

```bash
test/simulated/run.sh
```

See `test/simulated/README.md` for the truth variants and test design.

## Running IGenotyper

The input PacBio BAM must have a SAMtools `.bai` index.

```bash
samtools index reads.bam
IG phase --sample SAMPLE --threads 8 reads.bam output
IG assembly --threads 8 output
IG detect output
```

## Usage

### `IG phase`

```text
IG phase [--rhesus] [--sample SAMPLE] [--threads THREADS] [--mem MEM]
         [--cluster] [--queue QUEUE] [--walltime WALLTIME] [--tmp TMP]
         [--data-dir DATA_DIR] [--input_vcf VCF] BAM OUTDIR
```

- `--sample SAMPLE`: sample name; default `sample`
- `--threads THREADS`: worker threads; default `1`
- `--mem MEM`: cluster memory request; default `20`
- `--cluster`: submit work through the optional Watson-IG cluster package
- `--queue QUEUE`: cluster queue; default `premium`
- `--walltime WALLTIME`: cluster walltime; default `24`
- `--tmp TMP`: temporary directory; default `<OUTDIR>/tmp`
- `--data-dir DATA_DIR`: directory containing `reference.fasta`
- `--input_vcf VCF`: use an existing phased VCF
- `--rhesus`: use the bundled rhesus resources

### `IG assembly`

```text
IG assembly [--rhesus] [--threads THREADS] [--mem MEM] [--cluster]
            [--queue QUEUE] [--walltime WALLTIME] [--data-dir DATA_DIR]
            OUTDIR
```

Assembly contigs are aligned exclusively with minimap2's `asm20` preset.

### `IG detect`

```text
IG detect [--rhesus] [--hom HOM] [--data-dir DATA_DIR] OUTDIR
```

### `IG alleles`

```text
IG alleles [--database DB] [--num_reads NUM_READS] [--data-dir DATA_DIR]
           OUTDIR
```

Run `IG <command> --help` for the descriptions of every option.

## Explanation of steps
### Phase
In the `phase` step, CCS reads are aligned with minimap2 and phased. Read group
annotations 1 and 2 correspond to haplotypes 1 and 2; read group 0 contains
unassigned reads. In IGV, these can be viewed by grouping alignments by read
group.

### Assemble
In the `assembly` step, Canu assembles each haplotype block. A script and output
directory are created for every region/haplotype block. With `--cluster`, these
independent jobs can be submitted through the optional cluster integration.

### Detect
In the `detect` step, SNVs, indels, structural variants, genes, and alleles are
genotyped. SNVs are written as VCF, indels and structural variants as BED, and
gene/allele calls as tab-delimited output.

## Output directories
alignments  alleles  assembly  logs  plots  preprocessed  report.html  tmp  variants
| Directories            | Description                                          |
|------------------------|------------------------------------------------------|
| `<output>/alignments`  | Alignments of CCS, subreads and contigs (phased and unphased) |
| `<output>/assembly`    | Assembly of IGH locus                                |
| `<output>/variants`    | SNVs, indels and SVs                                 |
| `<output>/alleles`     | Alleles in sample                                    |
| `<output>/logs`        | Log files with input parameters                   |
| `<output>/tmp`         | Temporary files. Could be deleted.                   |

## Output files
1. `alignments/`
    1. `ccs_to_ref*`: CCS reads aligned to reference
    2. `contigs_to_ref*`: All assembled contigs aligned to the reference
    3. `igh_contigs_to_ref*`: IGH assembled contigs aligned to igh reference
2. `assembly/`
    1. `contigs.fasta`: All assembled contigs
    2. `igh_contigs.fasta`: IGH assembled contigs
3. `alleles/`
    1. `assembly_alleles.bed`: Alleles extracted from the assembly for each gene
    2. `assembly_genes.fasta`: Fasta sequence from the assembly for each gene/allele
    3. `ccs_alleles.bed`: Alleles extracted from the CCS reads for each gene
    4. `ccs_genes.fasta`: Fasta sequence from the CCS reads for each gene/allele
4. `logs/`
    1. `gene_cov.txt`: Haplotype coverage for each gene
5. `variants/`
    1. `snvs_phased_from_ccs.vcf`: Phased SNVs detected from the CCS reads
    2. `snvs_assembly.vcf`: SNVS detected from the assembly
    3. `indel_assembly.bed`: Indels detected from the assembly
    4. `sv_assembly.bed`: SVs detected from the assembly
    5. `phased_blocks.txt`: Phased haplotype blocks


# Todo
1. Add IGL and IGK alleles to data/alleles.fasta

### Recovering interrupted or older runs

See [pipeline fixes and recovery](docs/pipeline-recovery.md) for the WhatsHap 2.8
read-name collision fix, completion records, legacy-output migration, assembly
platform behavior, and synthetic regression tests.

Assembly requires at least 20x mean depth across IG target bases. Lower-coverage
samples exit cleanly with `assembly/assembly_status.json` recording
`insufficient_coverage` and `retryable: false`; phasing outputs are preserved.
Use `IG assembly --coverage-bed IG_TARGETS.bed OUTDIR` for custom IG coordinates.
