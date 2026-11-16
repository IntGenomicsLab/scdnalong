# scdnalong: Citations

If you use results produced with scdnalong, please cite the scdnalong release you used
(0.0.1: [10.5281/zenodo.16099232](https://doi.org/10.5281/zenodo.16099232)) and the tools below.

## [nf-core](https://pubmed.ncbi.nlm.nih.gov/32055031/)

scdnalong is built on the nf-core template and is not part of nf-core.

> Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PubMed PMID: 32055031.

## [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)

> Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

## Pipeline tools

- [flexiplex](https://github.com/DavidsonGroup/flexiplex) (barcode discovery, filtering and demultiplexing)

  > Cheng O, Ling MH, Wang C, Wu S, Ritchie ME, Göke J, Amin N, Davidson NM. Flexiplex: a versatile demultiplexer and search tool for omics data. Bioinformatics. 2024;40(3):btae102. doi: 10.1093/bioinformatics/btae102.

- [flexiformatter](https://github.com/ljwharbers/flexiformatter) (moves the barcode into the `CB` tag)

  > Harbers L. flexiformatter. https://github.com/ljwharbers/flexiformatter

- [minimap2](https://github.com/lh3/minimap2) (alignment)

  > Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018;34(18):3094-3100. doi: 10.1093/bioinformatics/bty191.

- [Picard](https://broadinstitute.github.io/picard/) (MarkDuplicates)

  > Broad Institute. Picard Toolkit. https://broadinstitute.github.io/picard/

- [SAMtools](https://www.htslib.org/) (sort, index, stats)

  > Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. Twelve years of SAMtools and BCFtools. GigaScience. 2021;10(2):giab008. doi: 10.1093/gigascience/giab008.

- [SeqKit](https://bioinf.shenwei.me/seqkit/) (read statistics and splitting)

  > Shen W, Le S, Li Y, Hu F. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PLoS One. 2016;11(10):e0163962. doi: 10.1371/journal.pone.0163962.

- [NanoPlot](https://github.com/wdecoster/NanoPlot) and [NanoComp](https://github.com/wdecoster/nanocomp) (read QC)

  > De Coster W, Rademakers R. NanoPack2: population-scale evaluation of long-read sequencing data. Bioinformatics. 2023;39(5):btad311. doi: 10.1093/bioinformatics/btad311.

- [ToulligQC](https://github.com/GenomiqueENS/toulligQC) (read QC)

  > GenomiqueENS. ToulligQC: post-sequencing QC for Oxford Nanopore data. https://github.com/GenomiqueENS/toulligQC

- [pigz](https://zlib.net/pigz/) (compression)

  > Adler M. pigz: a parallel implementation of gzip. https://zlib.net/pigz/

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

  > Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data [Online].

- [MultiQC](https://pubmed.ncbi.nlm.nih.gov/27312411/)

  > Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. doi: 10.1093/bioinformatics/btw354. Epub 2016 Jun 16. PubMed PMID: 27312411; PubMed Central PMCID: PMC5039924.

## Software packaging/containerisation tools

- [Anaconda](https://anaconda.com)

  > Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.

- [Bioconda](https://pubmed.ncbi.nlm.nih.gov/29967506/)

  > Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018 Jul;15(7):475-476. doi: 10.1038/s41592-018-0046-7. PubMed PMID: 29967506.

- [BioContainers](https://pubmed.ncbi.nlm.nih.gov/28379341/)

  > da Veiga Leprevost F, Grüning B, Aflitos SA, Röst HL, Uszkoreit J, Barsnes H, Vaudel M, Moreno P, Gatto L, Weber J, Bai M, Jimenez RC, Sachsenberg T, Pfeuffer J, Alvarez RV, Griss J, Nesvizhskii AI, Perez-Riverol Y. BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics. 2017 Aug 15;33(16):2580-2582. doi: 10.1093/bioinformatics/btx192. PubMed PMID: 28379341; PubMed Central PMCID: PMC5870671.

- [Docker](https://dl.acm.org/doi/10.5555/2600239.2600241)

  > Merkel, D. (2014). Docker: lightweight linux containers for consistent development and deployment. Linux Journal, 2014(239), 2. doi: 10.5555/2600239.2600241.

- [Singularity](https://pubmed.ncbi.nlm.nih.gov/28494014/)

  > Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.
