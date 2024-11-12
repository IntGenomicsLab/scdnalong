//
// Creates gtfs to that add introns as features
//

include { PIGZ_UNCOMPRESS          } from '../../modules/nf-core/pigz/uncompress/main'
include { PIGZ_COMPRESS            } from '../../modules/nf-core/pigz/compress/main'
include { FLEXIPLEX_DISCOVERY      } from '../../modules/local/flexiplex/discovery/main'
include { FLEXIPLEX_FILTER         } from '../../modules/local/flexiplex/filter/main'
include { FLEXIPLEX_ASSIGN         } from '../../modules/local/flexiplex/assign/main'

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
        // MODULE: Assign flexiplex
        //
        
        FLEXIPLEX_ASSIGN (
            ch_reads,
            FLEXIPLEX_FILTER.out.barcodes
        )
        
        ch_versions = ch_versions.mix(FLEXIPLEX_ASSIGN.out.versions)
        
        //
        // MODULE: Compress Fastq
        //
        PIGZ_COMPRESS ( FLEXIPLEX_ASSIGN.out.reads )
        
        
    emit:
        flexiplex_fastq = PIGZ_COMPRESS.out.archive
        versions = ch_versions
}