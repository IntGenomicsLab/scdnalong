# scdnalong: Usage

> [!IMPORTANT]
> scdnalong is archived and superseded by
> [LRSinglecell](https://github.com/IntGenomicsLab/lrsinglecell). This page documents
> how to reproduce results made with scdnalong `0.0.1`.

The complete instructions (requirements, samplesheet, parameters, outputs and known
issues) are in the [README](../README.md). This page repeats the essentials.

## Quick start

```bash
export NXF_VER=24.10.4   # 0.0.1 does not parse on Nextflow 25.10+

nextflow run IntGenomicsLab/scdnalong -r 0.0.1 \
    -profile singularity \
    --input samplesheet.csv \
    --outdir results \
    --fasta GRCh38.fasta \
    --whitelist 737K-arc-v1.txt \
    --barcode_format 10x_atac \
    --split_fastq_n 24
```

## Samplesheet

Two columns, comma-separated, with a header. One row per gzipped Nanopore FASTQ file;
rows sharing a `sample` name are concatenated.

```csv
sample,fastq_1
D0,/path/to/D0_ATAC_run1.fastq.gz
D0,/path/to/D0_ATAC_run2.fastq.gz
Q1,/path/to/Q1_ATAC.fastq.gz
```

| Column | Description |
| --- | --- |
| `sample` | Sample name without spaces; becomes the output folder name. |
| `fastq_1` | Gzipped FASTQ (`.fastq.gz`/`.fq.gz`), local path or URL. |

An example is in [`assets/samplesheet.csv`](../assets/samplesheet.csv).

## Required parameters

| Parameter | Description |
| --- | --- |
| `--input` | The samplesheet. |
| `--outdir` | Output directory. |
| `--fasta` | Reference genome FASTA (plain or gzipped). |
| `--whitelist` | Uncompressed 10x ARC ATAC whitelist `737K-arc-v1.txt`. |
| `--barcode_format` | Must be `10x_atac` (the only layout implemented). |

See the README for optional parameters.

## Profiles

- `docker`, `singularity`, `apptainer`, `podman`, `conda`: software packaging.
- `test`: public test data (on the `dev` branch only; see README).

Combine with an institutional profile or `-c` config to run on a cluster.
