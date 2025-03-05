include { BWAMEM2_INDEX                  } from '../../../modules/nf-core/bwamem2/index'
include { BWA_INDEX as BWAMEM1_INDEX     } from '../../../modules/nf-core/bwa/index'
include { DRAGMAP_HASHTABLE              } from '../../../modules/nf-core/dragmap/hashtable'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/gatk4/createsequencedictionary'
include { GAWK as BUILD_INTERVALS        } from '../../../modules/nf-core/gawk'
include { MSISENSORPRO_SCAN              } from '../../../modules/nf-core/msisensorpro/scan'
include { SAMTOOLS_FAIDX                 } from '../../../modules/nf-core/samtools/faidx'
include { TABIX_BGZIPTABIX               } from '../../../modules/nf-core/tabix/bgziptabix'
include { TABIX_TABIX                    } from '../../../modules/nf-core/tabix/tabix'

workflow PREPARE_GENOME_DNASEQ {
    take:
    fasta                        // channel: [meta, fasta]
    fasta_fai                    // channel: [meta, fasta_fai]
    vcf                          // channel: [meta, vcf]
    run_bwamem1                  // boolean: true/false
    run_bwamem2                  // boolean: true/false
    run_createsequencedictionary // boolean: true/false
    run_dragmap                  // boolean: true/false
    run_faidx                    // boolean: true/false
    run_intervals                // boolean: true/false
    run_msisensorpro             // boolean: true/false
    run_tabix                    // boolean: true/false

    main:
    bwamem1_index = Channel.empty()
    bwamem2_index = Channel.empty()
    dragmap_hashmap = Channel.empty()
    fasta_dict = Channel.empty()
    fasta_sizes = Channel.empty()
    intervals_bed = Channel.empty()
    msisensorpro_list = Channel.empty()
    vcf_gz = Channel.empty()
    vcf_tbi = Channel.empty()

    versions = Channel.empty()

    if (run_bwamem1) {
        // Do not run BWAMEM1_INDEX if the condition is false
        fasta_bwamem1 = fasta.map { meta, fasta_ -> meta.run_bwamem1 ? [meta, fasta_] : null }

        BWAMEM1_INDEX(fasta_bwamem1)

        bwamem1_index = BWAMEM1_INDEX.out.index
        versions = versions.mix(BWAMEM1_INDEX.out.versions)
    }

    if (run_bwamem2) {
        // Do not run BWAMEM2_INDEX if the condition is false
        fasta_bwamem2 = fasta.map { meta, fasta_ -> meta.run_bwamem2 ? [meta, fasta_] : null }

        BWAMEM2_INDEX(fasta_bwamem2)

        bwamem2_index = BWAMEM2_INDEX.out.index
        versions = versions.mix(BWAMEM2_INDEX.out.versions)
    }

    if (run_dragmap) {
        // Do not run DRAGMAP_HASHTABLE if the condition is false
        fasta_dragmap = fasta.map { meta, fasta_ -> meta.run_dragmap ? [meta, fasta_] : null }

        DRAGMAP_HASHTABLE(fasta_dragmap)

        dragmap_hashmap = DRAGMAP_HASHTABLE.out.hashmap
        versions = versions.mix(DRAGMAP_HASHTABLE.out.versions)
    }

    if (run_createsequencedictionary) {
        // Do not run GATK4_CREATESEQUENCEDICTIONARY if the condition is false
        fasta_gat4kdict = fasta.map { meta, fasta_ -> meta.run_createsequencedictionary ? [meta, fasta_] : null }

        GATK4_CREATESEQUENCEDICTIONARY(fasta_gat4kdict)

        fasta_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict
        versions = versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    }

    if (run_faidx || run_intervals) {

        if (run_faidx) {

            // Do not run SAMTOOLS_FAIDX if the condition is false
            fasta_samtools = fasta.map { meta, fasta_ -> meta.run_faidx ? [meta, fasta_] : null }

            // No need to generate sizes for DNAseq
            generate_sizes = false

            SAMTOOLS_FAIDX(
                fasta_samtools,
                [[id: 'no_fai'], []],
                generate_sizes,
            )

            fasta_fai = fasta_fai.mix(SAMTOOLS_FAIDX.out.fai)
            versions = versions.mix(SAMTOOLS_FAIDX.out.versions)
        }

        if (run_intervals) {
            // Do not run BUILD_INTERVALS if the condition is false
            fasta_fai_intervals = fasta_fai.map { meta, fasta_fai_ -> meta.run_intervals ? [meta, fasta_fai_] : null }

            BUILD_INTERVALS(fasta_fai_intervals, [])
            intervals_bed = BUILD_INTERVALS.out.output
            versions = versions.mix(BUILD_INTERVALS.out.versions)
        }
    }

    if (run_msisensorpro) {
        // Do not run MSISENSORPRO_SCAN if the condition is false
        fasta_msisensorpro = fasta.map { meta, fasta_ -> meta.run_msisensorpro ? [meta, fasta_] : null }

        MSISENSORPRO_SCAN(fasta_msisensorpro)

        msisensorpro_list = MSISENSORPRO_SCAN.out.list
        versions = versions.mix(MSISENSORPRO_SCAN.out.versions)
    }

    if (run_tabix) {
        // Do not run TABIX_TABIX if the condition is false
        vcf_tabix = vcf.map { meta, vcf_ -> meta.run_tabix ? [meta, vcf_] : null }

        // Do not run TABIX_BGZIPTABIX if the condition is false
        vcf_bgziptabix = vcf.map { meta, vcf_ -> meta.compress_vcf ? [meta, vcf_] : null }

        TABIX_BGZIPTABIX(vcf_bgziptabix)
        TABIX_TABIX(vcf_tabix)

        vcf_gz = TABIX_BGZIPTABIX.out.gz_tbi.map { meta, vcf_gz_, _vcf_tbi -> [meta, vcf_gz_] }
        vcf_tbi = TABIX_TABIX.out.tbi.mix(TABIX_BGZIPTABIX.out.gz_tbi.map { meta, _vcf_gz, vcf_tbi_ -> [meta, vcf_tbi_] })

        versions = versions.mix(TABIX_BGZIPTABIX.out.versions)
        versions = versions.mix(TABIX_TABIX.out.versions)
    }

    emit:
    bwamem1_index     // channel: [meta, BWAmemIndex/]
    bwamem2_index     // channel: [meta, BWAmem2memIndex/]
    dragmap_hashmap   // channel: [meta, DragmapHashtable/]
    fasta_dict        // channel: [meta, *.fa(sta).dict]
    fasta_fai         // channel: [meta, *.fa(sta).fai]
    intervals_bed     // channel: [meta, *.bed]
    msisensorpro_list // channel: [meta, *.list]
    vcf_gz            // channel: [meta, *.vcf.gz]
    vcf_tbi           // channel: [meta, *.vcf.gz.tbi]
    versions          // channel: [versions.yml]
}
