# scdnalong: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## dev - archived

scdnalong is superseded by [LRSinglecell](https://github.com/IntGenomicsLab/lrsinglecell)
and archived.

### `Changed`

- Replaced the nf-core template README and usage docs with real documentation:
  samplesheet format, parameters, outputs, the Nextflow version needed to run
  `0.0.1`, known issues, and the relation to nf-core/scnanoseq and LRSinglecell.
- `assets/samplesheet.csv` now shows the actual two-column format (`sample,fastq_1`).
- `-profile test` now uses public data from IntGenomicsLab/test-datasets instead of
  paths on the VSC cluster, and runs with the local executor.

### Changes on `dev` since `0.0.1` (not in the published version)

- flexiplex updated from 1.01 to 1.02.5.
- `--flexiplex_filter_args` to pass options to `flexiplex-filter`.
- `SEQKIT_SPLIT2` runs as a job array (`array 100`); requires a grid executor.
- FastQC/NanoComp module updates, MultiQC and minimap2 index fixes, MarkDuplicates
  temp dir.

## 0.0.1 - 2025-07-18

First release (commit `d600a80`). Version used for published results.
