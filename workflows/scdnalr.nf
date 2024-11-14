/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                            } from '../modules/nf-core/fastqc/main'
include { MULTIQC                           } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                  } from 'plugin/nf-schema'
include { paramsSummaryMultiqc              } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText            } from '../subworkflows/local/utils_nfcore_scdnalr_pipeline'

include { CAT_FASTQ as CAT_FASTQ_SAMPLE     } from '../modules/nf-core/cat/fastq/main'
include { SEQKIT_STATS as SEQKIT_STATS_PRE  } from '../modules/nf-core/seqkit/stats/main'
include { SEQKIT_STATS as SEQKIT_STATS_POST } from '../modules/nf-core/seqkit/stats/main'
include { NANOCOMP as NANOCOMP_FASTQ        } from '../modules/nf-core/nanocomp/main'
include { NANOCOMP as NANOCOMP_BAM          } from '../modules/nf-core/nanocomp/main'
include { MINIMAP2_ALIGN                    } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_INDEX                    } from '../modules/nf-core/minimap2/index/main'
include { FLEXIFORMATTER                    } from '../modules/local/flexiformatter/main'
include { PICARD_MARKDUPLICATES             } from '../modules/nf-core/picard/markduplicates/main' 

/*
 * Import subworkflows
*/
include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_PRE_FLEXIPLEX  } from '../subworkflows/local/toulligqc_nanoplot_fastqc'
include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_POST_FLEXIPLEX } from '../subworkflows/local/toulligqc_nanoplot_fastqc'
include { RUN_FLEXIPLEX                                             } from '../subworkflows/local/run_flexiplex'
include { PREPARE_REFERENCE_FILES                                   } from '../subworkflows/local/prepare_reference_files'
include { BAM_SORT_STATS_SAMTOOLS                                   } from '../subworkflows/nf-core/bam_sort_stats_samtools/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SCDNALR {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    
    //
    // MODULE: Combine fastqs from the same sample
    //
    CAT_FASTQ_SAMPLE ( ch_samplesheet )
        .reads
        .set { ch_cat_fastq }

    ch_versions = ch_versions.mix (CAT_FASTQ_SAMPLE.out.versions.first().ifEmpty(null))
    
    //
    // SUBWORKFLOW: Fastq QC with Nanoplot and FastQC - Pre Flexiplex
    // Credits for this subworkflow go to nf-core/scnanoseq developers
    ch_fastqc_multiqc_pre_flexiplex = Channel.empty()
    ch_seqkit_stats_pre = Channel.empty()
    if (!params.skip_qc){
        FASTQC_NANOPLOT_PRE_FLEXIPLEX (
            ch_cat_fastq,
            params.skip_nanoplot,
            params.skip_toulligqc,
            params.skip_fastqc
        )

        ch_fastqc_multiqc_pre_flexiplex = FASTQC_NANOPLOT_PRE_FLEXIPLEX.out.fastqc_multiqc.ifEmpty([])
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_PRE_FLEXIPLEX.out.nanoplot_version.first().ifEmpty(null))
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_PRE_FLEXIPLEX.out.toulligqc_version.first().ifEmpty(null))
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_PRE_FLEXIPLEX.out.fastqc_version.first().ifEmpty(null))
        
        SEQKIT_STATS_PRE (
            ch_cat_fastq
        )
        ch_seqkit_stats_pre = SEQKIT_STATS_PRE.out.stats
        ch_versions = ch_versions.mix(SEQKIT_STATS_PRE.out.versions.first().ifEmpty(null)) 
    }
    
    //
    // MODULE: NanoComp for FastQ files
    // Credits for this module go to nf-core/scnanoseq developers
    ch_cat_fastq
        .collect{it[1]}
        .map {
        [ [ 'id': 'nanocomp_fastq.' ] , it ]
        }
    ch_nanocomp_fastq_html = Channel.empty()
    ch_nanocomp_fastq_txt = Channel.empty()
    if (!params.skip_qc && !params.skip_fastq_nanocomp) {

        NANOCOMP_FASTQ (
            ch_cat_fastq
                .collect{it[1]}
                .map{
                    [ [ 'id': 'nanocomp_fastq.' ] , it ]
                }
        )

        ch_nanocomp_fastq_html = NANOCOMP_FASTQ.out.report_html
        ch_nanocomp_fastq_txt = NANOCOMP_FASTQ.out.stats_txt

        ch_versions = ch_versions.mix( NANOCOMP_FASTQ.out.versions )

    }
    

    
    //
    //
    // SUBWORKFLOW: RUN_FLEXIPLEX
    //
    RUN_FLEXIPLEX (
        ch_cat_fastq,
        params.whitelist
    )
    
    RUN_FLEXIPLEX.out.flexiplex_fastq
        .set { ch_flexiplex_fastq }
    
    ch_versions = ch_versions.mix(RUN_FLEXIPLEX.out.versions)
    
    //
    // SUBWORKFLOW: Fastq QC with Nanoplot and FastQC - post flexiplex
    // Credits for this subworkflow go to nf-core/scnanoseq developers
    ch_fastqc_multiqc_post_flexiplex = Channel.empty()
    ch_seqkit_stats_post = Channel.empty()
    if (!params.skip_qc){
        FASTQC_NANOPLOT_POST_FLEXIPLEX (
            ch_flexiplex_fastq,
            params.skip_nanoplot,
            params.skip_toulligqc,
            params.skip_fastqc
        )

        ch_fastqc_multiqc_post_flexiplex = FASTQC_NANOPLOT_POST_FLEXIPLEX.out.fastqc_multiqc.ifEmpty([])
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_FLEXIPLEX.out.nanoplot_version.first().ifEmpty(null))
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_FLEXIPLEX.out.toulligqc_version.first().ifEmpty(null))
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_FLEXIPLEX.out.fastqc_version.first().ifEmpty(null))

        SEQKIT_STATS_POST (
            ch_flexiplex_fastq
        )
        ch_seqkit_stats_post = SEQKIT_STATS_POST.out.stats
        ch_versions = ch_versions.mix(SEQKIT_STATS_POST.out.versions.first().ifEmpty(null))
    }
    
    //
    // SUBWORKFLOW: PREPARE_REFERENCE_FILES
    //
    
    PREPARE_REFERENCE_FILES (
        params.fasta
    )
    
    ch_fasta = PREPARE_REFERENCE_FILES.out.prepped_fasta
    ch_fai = PREPARE_REFERENCE_FILES.out.prepped_fai
    
    ch_versions = ch_versions.mix(PREPARE_REFERENCE_FILES.out.versions)
    
    //
    // MODULE: Run MINIMAP2_INDEX
    //
    
    // Create minimap2 index channel
    
    
    if (!params.skip_save_minimap2_index) {
        
        MINIMAP2_INDEX ( ch_fasta )
        ch_minimap_index = MINIMAP2_INDEX.out.index
        
        ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
    }

    //
    // MODULE: Run MINIMAP2_ALIGN TODO: add module specific options (no splicing)
    //
    MINIMAP2_ALIGN (
        ch_flexiplex_fastq,
        ch_minimap_index,
        true,
        'bai',
        "",
        ""
    )

    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
    MINIMAP2_ALIGN.out.bam 
        | set { ch_minimap_bam }
    
    //
    // MODULE: Run FLEXI_FORMATTER
    //
    FLEXIFORMATTER (
        ch_minimap_bam
    )
    ch_versions = ch_versions.mix(FLEXIFORMATTER.out.versions)
    FLEXIFORMATTER.out.bam
        | set { ch_tagged_bam }
    
    
    //
    // MODULE: MarkDuplicates
    //
    PICARD_MARKDUPLICATES ( 
        ch_tagged_bam,
        ch_fasta,
        ch_fai
    )
    
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)
    
    //
    // SUBWORKFLOW: BAM_SORT_STATS_SAMTOOLS
    // 
    BAM_SORT_STATS_SAMTOOLS (
        ch_tagged_bam,
        ch_fasta 
    )
    
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)
    
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )
    
    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )
    
    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
