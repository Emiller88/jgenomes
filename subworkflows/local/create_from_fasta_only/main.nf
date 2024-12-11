include { BOWTIE2_BUILD                  } from '../../../modules/nf-core/bowtie2/build'
include { BOWTIE_BUILD as BOWTIE1_BUILD  } from '../../../modules/nf-core/bowtie/build'
include { BWAMEM2_INDEX                  } from '../../../modules/nf-core/bwamem2/index'
include { BWA_INDEX as BWAMEM1_INDEX     } from '../../../modules/nf-core/bwa/index'
include { DRAGMAP_HASHTABLE              } from '../../../modules/nf-core/dragmap/hashtable'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/gatk4/createsequencedictionary'
include { GAWK as BUILD_INTERVALS        } from '../../../modules/nf-core/gawk'
include { MSISENSORPRO_SCAN              } from '../../../modules/nf-core/msisensorpro/scan'
include { SAMTOOLS_FAIDX                 } from '../../../modules/nf-core/samtools/faidx'

workflow CREATE_FROM_FASTA_ONLY {
    take:
    fasta                        // channel: [meta, fasta]
    fasta_fai                    // channel: [meta, fasta_fai]
    run_bowtie1                  // boolean: true/false
    run_bowtie2                  // boolean: true/false
    run_bwamem1                  // boolean: true/false
    run_bwamem2                  // boolean: true/false
    run_createsequencedictionary // boolean: true/false
    run_dragmap                  // boolean: true/false
    run_faidx                    // boolean: true/false
    run_intervals                // boolean: true/false
    run_msisensorpro             // boolean: true/false
    run_sizes                    // boolean: true/false

    main:
    bowtie1_index = Channel.empty()
    bowtie2_index = Channel.empty()
    bwamem1_index = Channel.empty()
    bwamem2_index = Channel.empty()
    dragmap_hashmap = Channel.empty()
    fasta_dict = Channel.empty()
    fasta_fai = Channel.empty()
    fasta_sizes = Channel.empty()
    intervals_bed = Channel.empty()
    msisensorpro_list = Channel.empty()

    versions = Channel.empty()

    if (run_bowtie1) {
        fasta_bowtie1 = fasta.map { meta, fasta_ ->
            return meta.run_bowtie1 ? [meta, fasta_] : null
        }

        BOWTIE1_BUILD(fasta_bowtie1)

        bowtie1_index = BOWTIE1_BUILD.out.index
        versions = versions.mix(BOWTIE1_BUILD.out.versions)
    }

    if (run_bowtie2) {
        fasta_bowtie2 = fasta.map { meta, fasta_ ->
            return meta.run_bowtie2 ? [meta, fasta_] : null
        }

        BOWTIE2_BUILD(fasta_bowtie2)

        bowtie2_index = BOWTIE2_BUILD.out.index
        versions = versions.mix(BOWTIE2_BUILD.out.versions)
    }

    if (run_bwamem1) {
        fasta_bwamem1 = fasta.map { meta, fasta_ ->
            return meta.run_bwamem1 ? [meta, fasta_] : null
        }

        BWAMEM1_INDEX(fasta_bwamem1)

        bwamem1_index = BWAMEM1_INDEX.out.index
        versions = versions.mix(BWAMEM1_INDEX.out.versions)
    }

    if (run_bwamem2) {
        fasta_bwamem2 = fasta.map { meta, fasta_ ->
            return meta.run_bwamem2 ? [meta, fasta_] : null
        }

        BWAMEM2_INDEX(fasta_bwamem2)

        bwamem2_index = BWAMEM2_INDEX.out.index
        versions = versions.mix(BWAMEM2_INDEX.out.versions)
    }

    if (run_dragmap) {
        fasta_dragmap = fasta.map { meta, fasta_ ->
            return meta.run_dragmap ? [meta, fasta_] : null
        }

        DRAGMAP_HASHTABLE(fasta_dragmap)

        dragmap_hashmap = DRAGMAP_HASHTABLE.out.hashmap
        versions = versions.mix(DRAGMAP_HASHTABLE.out.versions)
    }

    if (run_createsequencedictionary) {
        fasta_gat4kdict = fasta.map { meta, fasta_ ->
            return meta.run_gat4kdict ? [meta, fasta_] : null
        }

        GATK4_CREATESEQUENCEDICTIONARY(fasta_gat4kdict)

        fasta_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict
        versions = versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    }

    if (run_faidx || run_intervals || run_sizes) {
        fasta_samtools = fasta.map { meta, fasta_ ->
            return meta.run_faidx ? [meta, fasta_] : null
        }

        SAMTOOLS_FAIDX(
            fasta_samtools,
            [[id: 'no_fai'], []],
            run_sizes
        )

        fasta_fai = fasta_fai.mix(SAMTOOLS_FAIDX.out.fai)
        fasta_sizes = SAMTOOLS_FAIDX.out.sizes
        versions = versions.mix(SAMTOOLS_FAIDX.out.versions)

        if (run_intervals) {
            fasta_fai_intervals = fasta_fai.map { meta, fasta_fai_ ->
                return meta.run_intervals ? [meta, fasta_fai_] : null
            }

            BUILD_INTERVALS(fasta_fai_intervals, [])
            intervals_bed = BUILD_INTERVALS.out.output
            versions = versions.mix(BUILD_INTERVALS.out.versions)
        }
    }

    if (run_msisensorpro) {
        fasta_msisensorpro = fasta.map { meta, fasta_ ->
            return meta.run_msisensorpro ? [meta, fasta_] : null
        }

        MSISENSORPRO_SCAN(fasta_msisensorpro)

        msisensorpro_list = MSISENSORPRO_SCAN.out.list
        versions = versions.mix(MSISENSORPRO_SCAN.out.versions)
    }

    emit:
    bowtie1_index     // channel: [meta, BowtieIndex/]
    bowtie2_index     // channel: [meta, Bowtie2Index/]
    bwamem1_index     // channel: [meta, BWAmemIndex/]
    bwamem2_index     // channel: [meta, BWAmem2memIndex/]
    dragmap_hashmap   // channel: [meta, DragmapHashtable/]
    fasta_dict        // channel: [meta, *.fa(sta).dict]
    fasta_fai         // channel: [meta, *.fa(sta).fai]
    fasta_sizes       // channel: [meta, *.fa(sta).sizes]
    intervals_bed     // channel: [meta, *.bed]
    msisensorpro_list // channel: [meta, *.list]
    versions          // channel: [versions.yml]
}
