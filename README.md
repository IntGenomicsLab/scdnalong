# scdnalong

> [!IMPORTANT]
> **scdnalong is archived and no longer maintained.** It has been superseded by
> [**LRSinglecell**](https://github.com/IntGenomicsLab/lrsinglecell), which processes
> single-cell long-read DNA *and* cDNA in one pipeline. Use LRSinglecell for new work.
> This repository stays available so that published results produced with scdnalong
> can be reproduced. The instructions below are for that purpose.

scdnalong is a Nextflow pipeline that processes single-cell long-read **DNA** data
from Oxford Nanopore sequencing of 10x Genomics barcoded libraries (the ATAC arm of
10x Multiome). Per sample, it identifies cell barcodes in the reads, assigns each
read to a cell, aligns to a reference genome, and marks PCR duplicates per cell. It
produces a cell-barcode-tagged, duplicate-marked BAM file and QC reports.

**Relation to other pipelines.** scdnalong was built from the
[nf-core](https://nf-co.re) template and borrows its QC subworkflows from
[nf-core/scnanoseq](https://github.com/nf-core/scnanoseq), but it is **not** an
nf-core pipeline and is **not** part of nf-core/scnanoseq (the `nf-core/` prefix that
appears in log messages comes from the template). nf-core/scnanoseq handles
single-cell long-read *RNA*; scdnalong handled the *DNA* arm. LRSinglecell now
combines both.

## Version used in publications

| What | Value |
| --- | --- |
| Release | [`0.0.1`](https://github.com/IntGenomicsLab/scdnalong/releases/tag/0.0.1) |
| Commit | `d600a80bdd54a8cb9f11b65f3eb7b7da10d8461e` |
| Nextflow | 24.10.4 (see [Nextflow version](#nextflow-version)) |

Always run with `-r 0.0.1`. The `dev` branch has changed since then (newer flexiplex
version, an extra filtering option, Slurm job arrays) and will not give identical
output.

## Pipeline steps

1. Concatenate all FASTQ files of a sample (`cat`).
2. Read QC: FastQC, NanoPlot, ToulligQC, `seqkit stats`, NanoComp.
3. **Barcode discovery**: [flexiplex](https://github.com/DavidsonGroup/flexiplex)
   (v1.01) counts candidate 16 bp cell barcodes, searched between the 10x ATAC
   flanking sequences (left `ACCGAGATCTACAC`, right
   `CGCGTCTGTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG`, `-f 8` edit distance on flanks).
4. **Barcode filtering**: `flexiplex-filter` keeps barcodes that are on the 10x
   whitelist and above the knee (inflection point) of the barcode rank plot.
5. **Read assignment**: reads are split into `--split_fastq_n` chunks
   (`seqkit split2 --by-part`), demultiplexed by flexiplex against the filtered
   barcode list in parallel, and merged again. The barcode is moved into the read name.
6. **Alignment**: minimap2 `-ax map-ont` (secondary alignments kept by default).
7. **Tagging**: [flexiformatter](https://github.com/ljwharbers/flexiformatter) (1.0.2)
   moves the barcode from the read name into the `CB` BAM tag.
8. **Duplicate marking**: Picard MarkDuplicates `--BARCODE_TAG CB`, so reads are only
   called duplicates within the same cell. Duplicates are flagged, not removed.
9. `samtools sort`/`index`/`stats`/`flagstat`/`idxstats`, NanoComp on the BAMs, and
   MultiQC.

All tools run from pinned BioContainers images; exact versions are in each
module's `container` line and are reported in
`<outdir>/pipeline_info/nf_core_pipeline_software_mqc_versions.yml` after a run.

## Requirements

- [Nextflow](https://www.nextflow.io/) **24.10.x** (not 25.10 or later, see below).
- Docker, Singularity or Apptainer.
- A reference genome FASTA (plain or `.gz`).
- The 10x Genomics ARC ATAC barcode whitelist `737K-arc-v1.txt` (shipped with
  Cell Ranger ARC as `cellranger-arc-*/lib/python/atac/barcodes/737K-arc-v1.txt.gz`;
  also bundled in LRSinglecell as
  [`assets/whitelist/cellranger_arc_atac.737K-arc-v1.txt.gz`](https://github.com/IntGenomicsLab/lrsinglecell/tree/main/assets/whitelist)).
  Decompress it before use.

### Nextflow version

`0.0.1` uses config syntax that the strict config parser of Nextflow 25.10+ rejects
(`Variable declarations cannot be mixed with config statements`). Pin Nextflow
24.10.4 for the run; this downloads it on first use and does not change your
installed version:

```bash
export NXF_VER=24.10.4
```

## Input: samplesheet

A CSV file with a header and **two** columns. One row per FASTQ file; rows with the
same `sample` are concatenated before processing.

```csv
sample,fastq_1
D0,/path/to/D0_ATAC_run1.fastq.gz
D0,/path/to/D0_ATAC_run2.fastq.gz
Q1,/path/to/Q1_ATAC.fastq.gz
```

| Column | Description |
| --- | --- |
| `sample` | Sample name, no spaces. Used as the output folder name. |
| `fastq_1` | Basecalled Nanopore reads, gzipped (`.fastq.gz` or `.fq.gz`). Local path or URL. |

There is no `fastq_2` column: Nanopore reads are single-end.

## Running

```bash
export NXF_VER=24.10.4

nextflow run IntGenomicsLab/scdnalong -r 0.0.1 \
    -profile singularity \
    --input samplesheet.csv \
    --outdir results \
    --fasta GRCh38.fasta \
    --whitelist 737K-arc-v1.txt \
    --barcode_format 10x_atac \
    --split_fastq_n 24
```

Use `-profile docker` or `-profile apptainer` as appropriate, plus your
institutional profile (or a `-c` config) for a cluster scheduler.

### Parameters

| Parameter | Required | Default | Description |
| --- | --- | --- | --- |
| `--input` | yes | | Samplesheet (see above). |
| `--outdir` | yes | `./results` | Output directory. |
| `--fasta` | yes | | Reference genome FASTA, plain or gzipped. |
| `--whitelist` | yes | | Uncompressed 10x barcode whitelist (`737K-arc-v1.txt`). |
| `--barcode_format` | yes | | Barcode layout. Only `10x_atac` is implemented. **If omitted, flexiplex silently uses its built-in 10x 3′ RNA layout**, which is wrong for these libraries. |
| `--split_fastq_n` | no | `1` | Number of chunks each sample's FASTQ is split into for parallel demultiplexing (not a number of reads). Does not change results. |
| `--save_secondary_alignment` | no | `true` | Keep secondary alignments from minimap2. |
| `--skip_dedup` | no | `false` | Skip Picard MarkDuplicates. |
| `--skip_qc` | no | `false` | Skip all QC steps. Individual QC tools: `--skip_fastqc`, `--skip_nanoplot`, `--skip_toulligqc`, `--skip_fastq_nanocomp`, `--skip_bam_nanocomp`, `--skip_multiqc`. |

`--whitelist`, `--barcode_format` and `--split_fastq_n` are missing from the
parameter schema, so Nextflow prints *"The following invalid input values have
been detected"* for them. This warning is harmless; the values are used.

### Small datasets

`flexiplex-filter` searches for the knee between barcode rank 50 and the rank at the
95th percentile of the barcode read counts. With very few barcodes (e.g. a heavily
subsampled test file) that window is empty and the step fails with
`ValueError: The lower bound must be greater than the upper bound`. `0.0.1` has no
option to switch the knee search off; this does not happen with real datasets.

## Output

```
results/
├── <sample>/
│   ├── bam/dedup/<sample>.dedup.sorted.bam(.bai)   # aligned, CB-tagged, duplicates flagged
│   ├── flexiplex/
│   │   ├── flexiplex_barcodes_counts.txt          # reads per candidate barcode
│   │   ├── <sample>_known_barcodes.txt            # barcodes kept (whitelist + knee)
│   │   └── <sample>.merged.fastq.gz               # demultiplexed reads
│   └── qc/                                        # fastqc, nanoplot, toulligqc, seqkit, samtools
├── batch_qcs/
│   ├── multiqc/{raw_qc,final_qc}/multiqc_report.html
│   └── nanocomp/{fastq,bam}/
└── pipeline_info/                                 # execution report, trace, software versions
```

In the BAM, `CB` holds the cell barcode. Duplicates carry SAM flag `0x400`; filter
them with `samtools view -F 0x400` if needed.

## Test run on public data

The pipeline was re-checked on 2026-11-09 with public data from
[IntGenomicsLab/test-datasets](https://github.com/IntGenomicsLab/test-datasets/tree/main/lrsinglecell)
(`sub_RPMI8402_DNA.fastq.gz`, 10x Multiome ATAC, and `chr21.fa`). On the `dev`
branch this is available as `-profile test`:

```bash
NXF_VER=24.10.4 nextflow run IntGenomicsLab/scdnalong -r dev -profile test,singularity --outdir results_test
```

The test file contains only ~140 barcodes, so the test profile turns off the knee
search (`--flexiplex_filter_args --no-inflection`, an option that exists on `dev`
only; see [Small datasets](#small-datasets)) and uses `assets/whitelist_test.txt`, the
101 barcodes of the 10x ARC ATAC whitelist that occur in the test reads. It also
disables the Slurm job array of `SEQKIT_SPLIT2` so it runs on a laptop or a single
node. A successful run ends with 410 aligned reads, all carrying a `CB` tag, across
101 cells.

## Credits

scdnalong was written by Luuk Harbers. QC subworkflows were adapted from
nf-core/scnanoseq (Austyn Trull, Lara Ianov). The pipeline skeleton comes from the
nf-core template.

## Citations

If you use results produced with scdnalong, please cite the tools listed in
[`CITATIONS.md`](CITATIONS.md), in particular flexiplex, minimap2, Picard and samtools,
and the nf-core framework:

> Ewels PA, Peltzer A, Fillinger S, et al. The nf-core framework for
> community-curated bioinformatics pipelines. *Nat Biotechnol.* 2020;38:276–278.
> doi: [10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x)
