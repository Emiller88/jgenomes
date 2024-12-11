include { GFFREAD                                               } from '../../../modules/nf-core/gffread'
include { HISAT2_BUILD                                          } from '../../../modules/nf-core/hisat2/build'
include { HISAT2_EXTRACTSPLICESITES                             } from '../../../modules/nf-core/hisat2/extractsplicesites'
include { KALLISTO_INDEX                                        } from '../../../modules/nf-core/kallisto/index'
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA       } from '../../../modules/nf-core/rsem/preparereference'
include { RSEM_PREPAREREFERENCE as RSEM_PREPAREREFERENCE_GENOME } from '../../../modules/nf-core/rsem/preparereference'
include { SALMON_INDEX                                          } from '../../../modules/nf-core/salmon/index'
include { STAR_GENOMEGENERATE                                   } from '../../../modules/nf-core/star/genomegenerate'

workflow CREATE_FROM_FASTA_AND_ANNOTATION {
    take:
    fasta                          // channel: [meta, fasta]
    gff                            // channel: [meta, gff]
    gtf                            // channel: [meta, gtf]
    splice_sites                   // channel: [meta, splice_sites]
    transcript_fasta               // channel: [meta, transcript_fasta]
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
        gff_gffread = gff.map { meta, _gff ->
            return meta.run_gffread ? [meta, _gff] : null
        }

        GFFREAD(
            gff_gffread,
            []
        )

        versions = versions.mix(GFFREAD.out.versions)

        gtf = gtf
            .mix(GFFREAD.out.gtf)
            .groupTuple()
            .map { meta, _gtf ->
                return _gtf[1] ? [meta, _gtf[1]] : [meta, _gtf]
            }

        if (run_hisat2 || run_hisat2_extractsplicesites) {
            gtf_hisat2 = gtf.map { meta, _gtf ->
                return meta.run_hisat2 ? [meta, _gtf] : null
            }

            HISAT2_EXTRACTSPLICESITES(gtf_hisat2)

            versions = versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)
            splice_sites = splice_sites.mix(HISAT2_EXTRACTSPLICESITES.out.txt)

            if (run_hisat2) {
                fasta_hisat2 = fasta.map { meta, _fasta ->
                    return meta.run_hisat2 ? [meta, _fasta] : null
                }

                HISAT2_BUILD(
                    fasta_hisat2,
                    gtf,
                    splice_sites
                )

                hisat2_index = HISAT2_BUILD.out.index

                versions = versions.mix(HISAT2_BUILD.out.versions)
            }
        }

        if (run_kallisto || run_rsem_make_transcript_fasta || run_salmon) {
            fasta_make_transcripts_fasta = fasta.map { meta, _fasta ->
                return meta.run_rsem_make_transcript_fasta ? [meta, _fasta] : null
            }

            MAKE_TRANSCRIPTS_FASTA(
                fasta_make_transcripts_fasta,
                gtf
            )
            versions = versions.mix(MAKE_TRANSCRIPTS_FASTA.out.versions)

            transcript_fasta = transcript_fasta.mix(MAKE_TRANSCRIPTS_FASTA.out.transcript_fasta)

            if (run_kallisto) {
                transcript_fasta_kallisto = transcript_fasta.map { meta, _transcript_fasta ->
                    return meta.run_kallisto ? [meta, _transcript_fasta] : null
                }

                KALLISTO_INDEX(transcript_fasta_kallisto)

                kallisto_index = KALLISTO_INDEX.out.index
                versions = versions.mix(KALLISTO_INDEX.out.versions)
            }

            if (run_salmon) {
                fasta_salmon = fasta.map { meta, _fasta ->
                    return meta.run_salmon ? [meta, _fasta] : null
                }

                transcript_fasta_salmon = transcript_fasta.map { meta, _transcript_fasta ->
                    return meta.run_salmon ? [meta, _transcript_fasta] : null
                }

                SALMON_INDEX(
                    fasta_salmon,
                    transcript_fasta_salmon
                )

                salmon_index = SALMON_INDEX.out.index
                versions = versions.mix(SALMON_INDEX.out.versions)
            }
        }

        if (run_rsem) {
            fasta_rsem = fasta.map { meta, _fasta ->
                return meta.run_rsem ? [meta, _fasta] : null
            }

            RSEM_PREPAREREFERENCE_GENOME(fasta_rsem, gtf)

            rsem_index = RSEM_PREPAREREFERENCE_GENOME.out.index
            versions = versions.mix(RSEM_PREPAREREFERENCE_GENOME.out.versions)
        }

        if (run_star) {
            fasta_star = fasta.map { meta, _fasta ->
                return meta.run_star ? [meta, _fasta] : null
            }


            STAR_GENOMEGENERATE(fasta_star, gtf)

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
