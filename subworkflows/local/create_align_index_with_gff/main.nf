include { GFFREAD                                               } from '../../../modules/nf-core/gffread'
include { HISAT2_BUILD                                          } from '../../../modules/nf-core/hisat2/build'
include { HISAT2_EXTRACTSPLICESITES                             } from '../../../modules/nf-core/hisat2/extractsplicesites'
include { KALLISTO_INDEX                                        } from '../../../modules/nf-core/kallisto/index'
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA       } from '../../../modules/nf-core/rsem/preparereference'
include { RSEM_PREPAREREFERENCE as RSEM_PREPAREREFERENCE_GENOME } from '../../../modules/nf-core/rsem/preparereference'
include { SALMON_INDEX                                          } from '../../../modules/nf-core/salmon/index'
include { STAR_GENOMEGENERATE                                   } from '../../../modules/nf-core/star/genomegenerate'

workflow CREATE_ALIGN_INDEX_WITH_GFF {
    take:
    fasta                          // channel: [meta, fasta]
    input_gff                      // channel: [meta, gff]
    input_gtf                      // channel: [meta, gtf]
    input_splice_sites             // channel: [meta, splice_sites]
    input_transcript_fasta         // channel: [meta, transcript_fasta]
    run_hisat2                     // boolean: true/false
    run_hisat2_extractsplicesites  // boolean: true/false
    run_kallisto                   // boolean: true/false
    run_rsem                       // boolean: true/false
    run_rsem_make_transcript_fasta // boolean: true/false
    run_salmon                     // boolean: true/false
    run_star                       // boolean: true/false

    main:
    gtf = Channel.empty()
    hisat2_index = Channel.empty()
    kallisto_index = Channel.empty()
    rsem_index = Channel.empty()
    salmon_index = Channel.empty()
    splice_sites = Channel.empty()
    star_index = Channel.empty()
    transcript_fasta = Channel.empty()

    versions = Channel.empty()

    if (run_hisat2 || run_kallisto || run_rsem || run_rsem_make_transcript_fasta || run_salmon || run_star) {

        GFFREAD(
            input_gff,
            []
        )

        versions = versions.mix(GFFREAD.out.versions)

        gtf = input_gtf
            .mix(GFFREAD.out.gtf)
            .groupTuple()
            .map { meta, file ->
                return file[1] ? [meta, file[1]] : [meta, file]
            }

        if (run_hisat2 || run_hisat2_extractsplicesites) {
            gtf_hisat2 = gtf.map { meta, map_gtf ->
                return meta.run_hisat2 ? [meta, map_gtf] : null
            }

            HISAT2_EXTRACTSPLICESITES(gtf_hisat2)

            versions = versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)
            splice_sites = input_splice_sites.mix(HISAT2_EXTRACTSPLICESITES.out.txt)

            if (run_hisat2) {
                HISAT2_BUILD(
                    fasta,
                    gtf,
                    splice_sites
                )

                hisat2_index = HISAT2_BUILD.out.index

                versions = versions.mix(HISAT2_BUILD.out.versions)
            }
        }

        if (run_kallisto || run_rsem_make_transcript_fasta || run_salmon) {
            fasta_make_transcripts_fasta = fasta.map { meta, map_fasta ->
                return meta.run_rsem_make_transcript_fasta ? [meta, map_fasta] : null
            }

            MAKE_TRANSCRIPTS_FASTA(
                fasta_make_transcripts_fasta,
                gtf
            )
            versions = versions.mix(MAKE_TRANSCRIPTS_FASTA.out.versions)

            transcript_fasta = input_transcript_fasta.mix(MAKE_TRANSCRIPTS_FASTA.out.transcript_fasta)

            if (run_kallisto) {
                KALLISTO_INDEX(transcript_fasta)

                kallisto_index = KALLISTO_INDEX.out.index
                versions = versions.mix(KALLISTO_INDEX.out.versions)
            }

            if (run_salmon) {
                SALMON_INDEX(
                    fasta,
                    transcript_fasta
                )

                salmon_index = SALMON_INDEX.out.index
                versions = versions.mix(SALMON_INDEX.out.versions)
            }
        }

        if (run_rsem) {
            RSEM_PREPAREREFERENCE_GENOME(fasta, gtf)

            rsem_index = RSEM_PREPAREREFERENCE_GENOME.out.index
            versions = versions.mix(RSEM_PREPAREREFERENCE_GENOME.out.versions)
        }

        if (run_star) {
            STAR_GENOMEGENERATE(fasta, gtf)

            star_index = STAR_GENOMEGENERATE.out.index
            versions = versions.mix(STAR_GENOMEGENERATE.out.versions)
        }
    }

    emit:
    gtf              // channel: [meta, gtf]
    hisat2_index     // channel: [meta, Hisat2Index/]
    kallisto_index   // channel: [meta, KallistoIndex]
    rsem_index       // channel: [meta, RSEMIndex/]
    salmon_index     // channel: [meta, SalmonIndex/]
    splice_sites     // channel: [meta, *.splice_sites.txt]
    star_index       // channel: [meta, STARIndex/]
    transcript_fasta // channel: [meta, *.transcripts.fasta]
    versions         // channel: [versions.yml]
}
