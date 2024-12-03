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
    fasta
    gff
    gtf
    splice_sites
    transcript_fasta
    run_hisat2
    run_kallisto
    run_rsem
    run_rsem_make_transcript_fasta
    run_salmon
    run_star

    main:
    gff_gtf = Channel.empty()
    hisat2_index = Channel.empty()
    hisat2_splice_sites = Channel.empty()
    kallisto_index = Channel.empty()
    rsem_index = Channel.empty()
    rsem_transcript_fasta = Channel.empty()
    salmon_index = Channel.empty()
    star_index = Channel.empty()

    versions = Channel.empty()

    if (run_hisat2 || run_kallisto || run_rsem || run_rsem_make_transcript_fasta || run_salmon || run_star) {

        gff_gffread = gff
            .join(gtf)
            .map { meta, input_gff, input_gtf ->
                return input_gtf ? null : [meta, input_gff]
            }

        GFFREAD(
            gff_gffread.map { meta, file ->
                return file ? [meta, file] : null
            },
            []
        )

        gff_gtf = gtf
            .mix(GFFREAD.out.gtf)
            .groupTuple()
            .map { meta, file ->
                return file[1] ? [meta, file[1]] : [meta, file]
            }
    }

    if (run_hisat2) {
        // TODO: be smarter about input assets
        //   Here we either return an empty channel if we have a splice_sites so that HISAT2_EXTRACTSPLICESITES is not triggered
        //   Or we return the provided gtf so that HISAT2_EXTRACTSPLICESITES is run
        ch_gtf_hisat2 = gff_gtf
            .join(splice_sites)
            .groupTuple()
            .map { meta, input_gtf, input_splice_sites ->
                return input_splice_sites[0][0] ? null : [meta, input_gtf]
            }

        HISAT2_EXTRACTSPLICESITES(ch_gtf_hisat2)
        versions = versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)

        // TODO: be smarter about input assets
        //   Here we either mix+GT an empty channel (either no output or no input splice_sites) with the splice_sites return splice_sites
        //   And we filter out the empty value
        hisat2_splice_sites = splice_sites
            .mix(HISAT2_EXTRACTSPLICESITES.out.txt)
            .groupTuple()
            .map { meta, txt ->
                return txt[1] ? [meta, txt[1]] : [meta, txt]
            }

        if (run_hisat2) {
            HISAT2_BUILD(
                fasta,
                gff_gtf,
                hisat2_splice_sites
            )

            hisat2_index = HISAT2_BUILD.out.index
            versions = versions.mix(HISAT2_BUILD.out.versions)
        }
    }

    if (run_kallisto || run_rsem_make_transcript_fasta || run_salmon) {
        // TODO: be smarter about input assets
        //   Here we either return an empty channel if we have a transcript_fasta so that MAKE_TRANSCRIPTS_FASTA is not triggered
        //   Or we return the provided gtf so that MAKE_TRANSCRIPTS_FASTA is run
        gtf_rsem = gff_gtf
            .join(transcript_fasta)
            .groupTuple()
            .map { meta, input_gtf, input_transcript_fasta ->
                return input_transcript_fasta[0][0] ? null : [meta, input_gtf]
            }

        MAKE_TRANSCRIPTS_FASTA(
            fasta,
            gtf_rsem
        )
        versions = versions.mix(MAKE_TRANSCRIPTS_FASTA.out.versions)

        // TODO: be smarter about input assets
        //   Here we either mix+GT an empty channel (either no output or no input transcript_fasta) with the transcript_fasta return transcript_fasta
        //   And we filter out the empty value
        rsem_transcript_fasta = transcript_fasta
            .mix(MAKE_TRANSCRIPTS_FASTA.out.transcript_fasta)
            .groupTuple()
            .map { meta, txt ->
                return txt[1] ? [meta, txt[1]] : [meta, txt]
            }

        if (run_kallisto) {
            KALLISTO_INDEX(rsem_transcript_fasta)

            kallisto_index = KALLISTO_INDEX.out.index
            versions = versions.mix(KALLISTO_INDEX.out.versions)
        }

        if (run_salmon) {
            SALMON_INDEX(
                fasta,
                rsem_transcript_fasta
            )

            salmon_index = SALMON_INDEX.out.index
            versions = versions.mix(SALMON_INDEX.out.versions)
        }
    }

    if (run_rsem) {
        RSEM_PREPAREREFERENCE_GENOME(fasta, gff_gtf)

        rsem_index = RSEM_PREPAREREFERENCE_GENOME.out.index
        versions = versions.mix(RSEM_PREPAREREFERENCE_GENOME.out.versions)
    }

    if (run_star) {
        STAR_GENOMEGENERATE(fasta, gff_gtf)

        star_index = STAR_GENOMEGENERATE.out.index
        versions = versions.mix(STAR_GENOMEGENERATE.out.versions)
    }

    emit:
    gff_gtf
    hisat2_index
    hisat2_splice_sites
    kallisto_index
    rsem_index
    rsem_transcript_fasta
    salmon_index
    star_index
    versions
}
