//
// Creates gtfs to that add introns as features
//

include { PIGZ_UNCOMPRESS                           } from '../../modules/nf-core/pigz/uncompress/main'
include { PIGZ_COMPRESS                             } from '../../modules/nf-core/pigz/compress/main'
include { SEQKIT_SPLIT2                             } from '../../modules/nf-core/seqkit/split2/main'
include { FLEXIPLEX_DISCOVERY                       } from '../../modules/local/flexiplex/discovery/main'
include { FLEXIPLEX_FILTER                          } from '../../modules/local/flexiplex/filter/main'
include { FLEXIPLEX_ASSIGN                          } from '../../modules/local/flexiplex/assign/main'
include { CAT_FASTQ as CAT_FASTQ_SPLIT              } from '../../modules/nf-core/cat/fastq/main'


workflow RUN_FLEXIPLEX {
    take:
        reads
        whitelist

    main:
        ch_versions = Channel.empty()

		//
        // Check if reads are zipped
        //
        gzipped = reads.map { meta, fastq_list -> 
            def all_gzipped = fastq_list.every { it.endsWith('.gz') }
            def none_gzipped = fastq_list.every { !it.endsWith('.gz') }

            if (!all_gzipped && !none_gzipped) {
                throw new IllegalArgumentException("Error: Mixed gzipped and non-gzipped files in ${fastq_list}")
            }
            
            return all_gzipped
        }
        
        // Uncompress if gzipped
        if (gzipped){
            PIGZ_UNCOMPRESS( reads )

            ch_reads = PIGZ_UNCOMPRESS.out.file
            ch_versions = ch_versions.mix(PIGZ_UNCOMPRESS.out.versions)
        } else {
            ch_reads = reads
        }

        //
        // MODULE: Run flexiplex
        //
        FLEXIPLEX_DISCOVERY (
            ch_reads
    	)
        
        ch_versions = ch_versions.mix(FLEXIPLEX_DISCOVERY.out.versions)
        
        // 
        // MODULE: Filter flexiplex
        //
        
        FLEXIPLEX_FILTER (
            FLEXIPLEX_DISCOVERY.out.barcode_counts,
            whitelist
        )
        
        ch_versions = ch_versions.mix(FLEXIPLEX_FILTER.out.versions)
        
        //
        // MODULE: Run SEQKIT_SPLIT2
        //
        SEQKIT_SPLIT2 (
            ch_reads
        )
        
        // Transpose channel and add part to metadata
        SEQKIT_SPLIT2.out.reads
            | transpose
            | map { meta, reads ->
                part = (reads =~ /.*part_(\d+)\.fastq(?:\.gz)?$/)[0][1]
                newmap = [part: part]
                [meta + newmap, reads] }
            | set { ch_split_fastq }
                
        ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)
        
        // Merge the reads and barcodes channels
        ch_split_fastq 
            | combine(FLEXIPLEX_FILTER.out.barcodes)
            | map { meta, reads, meta2, barcodes -> { 
                meta.id == meta2.id ? [meta, reads, barcodes] : null }}
            | set { ch_split_fastq_barcode }
    
        
        // 
        // MODULE: Assign flexiplex
        //
        FLEXIPLEX_ASSIGN (
            ch_split_fastq_barcode,
        )
        
        ch_versions = ch_versions.mix(FLEXIPLEX_ASSIGN.out.versions)
        // Group by ID for CATFASTQ
        FLEXIPLEX_ASSIGN.out.reads
            | map { meta, reads ->
                [meta.subMap('id', 'single_end'), meta.part, reads] }
            | groupTuple
            | map { meta, part, reads -> [meta + [partcount: part.size()], reads] }
            | set { ch_grouped_flexiplex_fastq }
        
        //
        // MODULE: cat fastq
        //
        
        CAT_FASTQ_SPLIT (
            ch_grouped_flexiplex_fastq
        )
        
        ch_versions = ch_versions.mix(CAT_FASTQ_SPLIT.out.versions)
        
        //
        // MODULE: Compress Fastq
        //
        PIGZ_COMPRESS ( CAT_FASTQ_SPLIT.out.reads )
        
        
    emit:
        flexiplex_fastq = PIGZ_COMPRESS.out.archive
        versions = ch_versions
}