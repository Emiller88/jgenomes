include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/gatk4/createsequencedictionary'
include { GAWK as BUILD_INTERVALS        } from '../../../modules/nf-core/gawk'
include { MSISENSORPRO_SCAN              } from '../../../modules/nf-core/msisensorpro/scan'
include { SAMTOOLS_FAIDX                 } from '../../../modules/nf-core/samtools/faidx'

workflow INDEX_FASTA {
    take:
    fasta                        // channel: [meta, fasta]
    input_fasta_fai              // channel: [meta, fasta_fai]
    run_createsequencedictionary // boolean: true/false
    run_faidx                    // boolean: true/false
    run_intervals                // boolean: true/false
    run_msisensorpro             // boolean: true/false
    run_sizes                    // boolean: true/false

    main:
    intervals_bed = Channel.empty()
    fasta_fai = Channel.empty()
    fasta_dict = Channel.empty()
    fasta_sizes = Channel.empty()
    msisensorpro_list = Channel.empty()

    versions = Channel.empty()

    if (run_createsequencedictionary) {
        fasta_gat4kdict = fasta.map { meta, map_fasta ->
            return meta.run_gat4kdict ? [meta, map_fasta] : null
        }

        GATK4_CREATESEQUENCEDICTIONARY(fasta_gat4kdict)

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
                return meta.run_intervals ? [meta, map_fasta_fai] : null
            }

            BUILD_INTERVALS(fasta_fai_intervals, [])
            intervals_bed = BUILD_INTERVALS.out.output
            versions = versions.mix(BUILD_INTERVALS.out.versions)
        }
    }

    if (run_msisensorpro) {
        fasta_msisensorpro = fasta.map { meta, map_fasta ->
            return meta.run_msisensorpro ? [meta, map_fasta] : null
        }

        MSISENSORPRO_SCAN(fasta_msisensorpro)

        msisensorpro_list = MSISENSORPRO_SCAN.out.list
        versions = versions.mix(MSISENSORPRO_SCAN.out.versions)
    }

    emit:
    fasta_dict        // channel: [meta, *.fa(sta).dict]
    fasta_fai         // channel: [meta, *.fa(sta).fai]
    fasta_sizes       // channel: [meta, *.fa(sta).sizes]
    intervals_bed     // channel: [meta, *.bed]
    msisensorpro_list // channel: [meta, *.list]
    versions          // channel: [versions.yml]
}
