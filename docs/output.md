# scdnalong: Output

> [!IMPORTANT]
> scdnalong is archived and superseded by
> [LRSinglecell](https://github.com/IntGenomicsLab/lrsinglecell). This page describes the
> output of scdnalong `0.0.1`.

All paths are relative to `--outdir`.

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

## Per-sample output

### Aligned reads: `<sample>/bam/dedup/`

- `<sample>.dedup.sorted.bam` and `.bai`: reads aligned with minimap2 (`-ax map-ont`), sorted and indexed.
- The `CB` tag holds the cell barcode of each read.
- PCR/optical duplicates are marked per cell (Picard MarkDuplicates `--BARCODE_TAG CB`). They carry SAM flag `0x400` and are **not removed**. Drop them with `samtools view -F 0x400` if needed.

The pipeline does not write a 10x-style `fragments.tsv.gz`. The BAM is the main result; split it per cell on the `CB` tag (e.g. `samtools split -d CB`) or convert it with your own tools.

### Barcodes: `<sample>/flexiplex/`

- `flexiplex_barcodes_counts.txt`: number of reads for each candidate 16 bp barcode found between the 10x ATAC flanks.
- `<sample>_known_barcodes.txt`: barcodes kept by `flexiplex-filter`, i.e. on the `--whitelist` and above the knee of the barcode-rank plot. This is the list of called cells.
- `<sample>.merged.fastq.gz`: reads that were assigned to a kept barcode, with the barcode in the read name.

### QC: `<sample>/qc/`

FastQC, NanoPlot, ToulligQC and `seqkit stats` on the raw reads, plus `samtools stats`, `flagstat` and `idxstats` on the final BAM.

## Batch output

- `batch_qcs/multiqc/raw_qc/multiqc_report.html` and `batch_qcs/multiqc/final_qc/multiqc_report.html`: MultiQC reports that collect the read-level and alignment-level QC of all samples.
- `batch_qcs/nanocomp/{fastq,bam}/`: NanoComp comparisons of read length and quality across samples.

## Pipeline information: `pipeline_info/`

- `execution_report_*.html`, `execution_timeline_*.html`, `execution_trace_*.txt`, `pipeline_dag_*.html`: Nextflow reports on runtime and resource use.
- `nf_core_pipeline_software_mqc_versions.yml`: the exact version of every tool used in the run.
- `params_*.json`: the parameters the pipeline was run with.
