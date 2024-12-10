include { BOWTIE_BUILD as BOWTIE1_BUILD } from '../../../modules/nf-core/bowtie/build'
include { BOWTIE2_BUILD                 } from '../../../modules/nf-core/bowtie2/build'
include { BWAMEM2_INDEX                 } from '../../../modules/nf-core/bwamem2/index'
include { BWA_INDEX as BWAMEM1_INDEX    } from '../../../modules/nf-core/bwa/index'
include { DRAGMAP_HASHTABLE             } from '../../../modules/nf-core/dragmap/hashtable'

workflow CREATE_ALIGN_INDEX {
    take:
    fasta       // channel: [meta, fasta]
    run_bowtie1 // boolean: true/false
    run_bowtie2 // boolean: true/false
    run_bwamem1 // boolean: true/false
    run_bwamem2 // boolean: true/false
    run_dragmap // boolean: true/false

    main:
    bowtie1_index = Channel.empty()
    bowtie2_index = Channel.empty()
    bwamem1_index = Channel.empty()
    bwamem2_index = Channel.empty()
    dragmap_hashmap = Channel.empty()

    versions = Channel.empty()

    if (run_bowtie1) {
        fasta_bowtie1 = fasta.map { meta, map_fasta ->
            return meta.run_bowtie1 ? [meta, map_fasta] : null
        }

        BOWTIE1_BUILD(fasta_bowtie1)

        bowtie1_index = BOWTIE1_BUILD.out.index
        versions = versions.mix(BOWTIE1_BUILD.out.versions)
    }

    if (run_bowtie2) {
        fasta_bowtie2 = fasta.map { meta, map_fasta ->
            return meta.run_bowtie2 ? [meta, map_fasta] : null
        }

        BOWTIE2_BUILD(fasta_bowtie2)

        bowtie2_index = BOWTIE2_BUILD.out.index
        versions = versions.mix(BOWTIE2_BUILD.out.versions)
    }

    if (run_bwamem1) {
        fasta_bwamem1 = fasta.map { meta, map_fasta ->
            return meta.run_bwamem1 ? [meta, map_fasta] : null
        }

        BWAMEM1_INDEX(fasta_bwamem1)

        bwamem1_index = BWAMEM1_INDEX.out.index
        versions = versions.mix(BWAMEM1_INDEX.out.versions)
    }

    if (run_bwamem2) {
        fasta_bwamem2 = fasta.map { meta, map_fasta ->
            return meta.run_bwamem2 ? [meta, map_fasta] : null
        }

        BWAMEM2_INDEX(fasta_bwamem2)

        bwamem2_index = BWAMEM2_INDEX.out.index
        versions = versions.mix(BWAMEM2_INDEX.out.versions)
    }

    if (run_dragmap) {
        fasta_dragmap = fasta.map { meta, map_fasta ->
            return meta.run_dragmap ? [meta, map_fasta] : null
        }

        DRAGMAP_HASHTABLE(fasta_dragmap)

        dragmap_hashmap = DRAGMAP_HASHTABLE.out.hashmap
        versions = versions.mix(DRAGMAP_HASHTABLE.out.versions)
    }

    emit:
    bowtie1_index   // channel: [meta, BowtieIndex/]
    bowtie2_index   // channel: [meta, Bowtie2Index/]
    bwamem1_index   // channel: [meta, BWAmemIndex/]
    bwamem2_index   // channel: [meta, BWAmem2memIndex/]
    dragmap_hashmap // channel: [meta, DragmapHashtable/]
    versions        // channel: [versions.yml]
}
