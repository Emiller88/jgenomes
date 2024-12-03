include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/gatk4/createsequencedictionary'
include { GAWK as BUILD_INTERVALS        } from '../../../modules/nf-core/gawk'
include { MSISENSORPRO_SCAN              } from '../../../modules/nf-core/msisensorpro/scan'
include { SAMTOOLS_FAIDX                 } from '../../../modules/nf-core/samtools/faidx'

workflow INDEX_FASTA {
    take:
    fasta
    input_fasta_fai
    run_createsequencedictionary
    run_faidx
    run_intervals
    run_msisensorpro
    run_sizes

    main:
    intervals_bed = Channel.empty()
    fasta_fai = Channel.empty()
    fasta_dict = Channel.empty()
    fasta_sizes = Channel.empty()
    msisensorpro_list = Channel.empty()

    versions = Channel.empty()

    if (run_createsequencedictionary) {
        GATK4_CREATESEQUENCEDICTIONARY(fasta)

        fasta_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict
        versions = versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    }

    if (run_faidx || run_intervals || run_sizes) {
        fasta_samtools = fasta.map { meta, map_fasta ->
            return meta.run_faidx ? [meta, map_fasta] : null
        }

        SAMTOOLS_FAIDX(
            fasta_samtools,
            [[id: 'no_fai'], []],
            run_sizes
        )

        fasta_fai = input_fasta_fai.mix(SAMTOOLS_FAIDX.out.fai)
        fasta_sizes = SAMTOOLS_FAIDX.out.sizes
        versions = versions.mix(SAMTOOLS_FAIDX.out.versions)

        if (run_intervals) {
            fasta_fai_intervals = fasta_fai.map { meta, map_fasta_fai ->
                return meta.intervals ? [meta, map_fasta_fai] : null
            }

            BUILD_INTERVALS(fasta_fai_intervals, [])
            intervals_bed = BUILD_INTERVALS.out.output
            versions = versions.mix(BUILD_INTERVALS.out.versions)
        }
    }

    if (run_msisensorpro) {
        MSISENSORPRO_SCAN(fasta)

        msisensorpro_list = MSISENSORPRO_SCAN.out.list
        versions = versions.mix(MSISENSORPRO_SCAN.out.versions)
    }

    emit:
    intervals_bed
    fasta_fai
    fasta_dict
    fasta_sizes
    msisensorpro_list
    versions
}
