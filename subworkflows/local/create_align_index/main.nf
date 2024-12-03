include { BOWTIE_BUILD as BOWTIE1_BUILD } from '../../../modules/nf-core/bowtie/build'
include { BOWTIE2_BUILD                 } from '../../../modules/nf-core/bowtie2/build'
include { BWAMEM2_INDEX                 } from '../../../modules/nf-core/bwamem2/index'
include { BWA_INDEX as BWAMEM1_INDEX    } from '../../../modules/nf-core/bwa/index'
include { DRAGMAP_HASHTABLE             } from '../../../modules/nf-core/dragmap/hashtable'

workflow CREATE_ALIGN_INDEX {
    take:
    fasta
    run_bowtie1
    run_bowtie2
    run_bwamem1
    run_bwamem2
    run_dragmap

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
    bowtie1_index
    bowtie2_index
    bwamem1_index
    bwamem2_index
    dragmap_hashmap
    versions
}
