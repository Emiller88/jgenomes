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
        BOWTIE1_BUILD(fasta)

        bowtie1_index = BOWTIE1_BUILD.out.index
        versions = versions.mix(BOWTIE1_BUILD.out.versions)
    }

    if (run_bowtie2) {
        BOWTIE2_BUILD(fasta)

        bowtie2_index = BOWTIE2_BUILD.out.index
        versions = versions.mix(BOWTIE2_BUILD.out.versions)
    }

    if (run_bwamem1) {
        BWAMEM1_INDEX(fasta)

        bwamem1_index = BWAMEM1_INDEX.out.index
        versions = versions.mix(BWAMEM1_INDEX.out.versions)
    }

    if (run_bwamem2) {
        BWAMEM2_INDEX(fasta)

        bwamem2_index = BWAMEM2_INDEX.out.index
        versions = versions.mix(BWAMEM2_INDEX.out.versions)
    }

    if (run_dragmap) {
        DRAGMAP_HASHTABLE(fasta)

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
