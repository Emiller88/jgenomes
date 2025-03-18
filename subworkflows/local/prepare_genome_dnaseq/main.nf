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
    intervals_bed = Channel.empty()
    msisensorpro_list = Channel.empty()
    vcf_gz = Channel.empty()
    vcf_tbi = Channel.empty()

    versions = Channel.empty()

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

    if (run_createsequencedictionary) {
        GATK4_CREATESEQUENCEDICTIONARY(fasta)

        fasta_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict
        versions = versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    }

    if (run_faidx || run_intervals) {

        if (run_faidx) {
            // Do not generate sizes for DNAseq
            generate_sizes = false

            SAMTOOLS_FAIDX(fasta.map { meta, fasta_ -> [meta, fasta_, []] }, generate_sizes)

            fasta_fai = fasta_fai.mix(SAMTOOLS_FAIDX.out.fai)
            versions = versions.mix(SAMTOOLS_FAIDX.out.versions)
        }

        if (run_intervals) {
            BUILD_INTERVALS(fasta_fai, [], false)
            intervals_bed = BUILD_INTERVALS.out.output
            versions = versions.mix(BUILD_INTERVALS.out.versions)
        }
    }

    if (run_msisensorpro) {
        MSISENSORPRO_SCAN(fasta)

        msisensorpro_list = MSISENSORPRO_SCAN.out.list
        versions = versions.mix(MSISENSORPRO_SCAN.out.versions)
    }

    if (run_tabix) {
        vcf_to_index = vcf.branch { _meta, vcf_ ->
            vcf_gz: vcf_.toString().endsWith('.gz')
            vcf: true
        }

        TABIX_BGZIPTABIX(vcf_to_index.vcf)
        TABIX_TABIX(vcf_to_index.vcf_gz)

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
