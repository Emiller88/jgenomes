include { BOWTIE2_BUILD                                         } from '../../../modules/nf-core/bowtie2/build'
include { BOWTIE_BUILD as BOWTIE1_BUILD                         } from '../../../modules/nf-core/bowtie/build'
include { GFFREAD                                               } from '../../../modules/nf-core/gffread'
include { HISAT2_BUILD                                          } from '../../../modules/nf-core/hisat2/build'
include { HISAT2_EXTRACTSPLICESITES                             } from '../../../modules/nf-core/hisat2/extractsplicesites'
include { KALLISTO_INDEX                                        } from '../../../modules/nf-core/kallisto/index'
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA       } from '../../../modules/nf-core/rsem/preparereference'
include { RSEM_PREPAREREFERENCE as RSEM_PREPAREREFERENCE_GENOME } from '../../../modules/nf-core/rsem/preparereference'
include { SALMON_INDEX                                          } from '../../../modules/nf-core/salmon/index'
include { SAMTOOLS_FAIDX                                        } from '../../../modules/nf-core/samtools/faidx'
include { STAR_GENOMEGENERATE                                   } from '../../../modules/nf-core/star/genomegenerate'

workflow PREPARE_GENOME_RNASEQ {
    take:
    fasta                          // channel: [meta, fasta]
    fasta_fai                      // channel: [meta, fasta_fai]
    gff                            // channel: [meta, gff]
    gtf                            // channel: [meta, gtf]
    splice_sites                   // channel: [meta, splice_sites]
    transcript_fasta               // channel: [meta, transcript_fasta]
    run_bowtie1                    // boolean: true/false
    run_bowtie2                    // boolean: true/false
    run_faidx                      // boolean: true/false
    run_hisat2                     // boolean: true/false
    run_hisat2_extractsplicesites  // boolean: true/false
    run_kallisto                   // boolean: true/false
    run_rsem                       // boolean: true/false
    run_rsem_make_transcript_fasta // boolean: true/false
    run_salmon                     // boolean: true/false
    run_star                       // boolean: true/false
    run_sizes                      // boolean: true/false

    main:
    bowtie1_index = Channel.empty()
    bowtie2_index = Channel.empty()
    hisat2_index = Channel.empty()
    kallisto_index = Channel.empty()
    rsem_index = Channel.empty()
    salmon_index = Channel.empty()
    star_index = Channel.empty()
    fasta_sizes = Channel.empty()

    versions = Channel.empty()

    if (run_bowtie1) {
        // Do not run BOWTIE1_BUILD if the condition is false
        fasta_bowtie1 = fasta.map { meta, fasta_ -> meta.run_bowtie1 ? [meta, fasta_] : null }

        BOWTIE1_BUILD(fasta_bowtie1)

        bowtie1_index = BOWTIE1_BUILD.out.index
        versions = versions.mix(BOWTIE1_BUILD.out.versions)
    }

    if (run_bowtie2) {
        // Do not run BOWTIE2_BUILD if the condition is false
        fasta_bowtie2 = fasta.map { meta, fasta_ -> meta.run_bowtie2 ? [meta, fasta_] : null }

        BOWTIE2_BUILD(fasta_bowtie2)

        bowtie2_index = BOWTIE2_BUILD.out.index
        versions = versions.mix(BOWTIE2_BUILD.out.versions)
    }

    if (run_faidx || run_sizes) {
        // Do not run SAMTOOLS_FAIDX if the condition is false
        fasta_samtools = fasta.map { meta, fasta_ -> meta.run_faidx ? [meta, fasta_] : null }

        SAMTOOLS_FAIDX(
            fasta_samtools,
            [[id: 'no_fai'], []],
            run_sizes,
        )

        fasta_fai = fasta_fai.mix(SAMTOOLS_FAIDX.out.fai)
        fasta_sizes = SAMTOOLS_FAIDX.out.sizes
        versions = versions.mix(SAMTOOLS_FAIDX.out.versions)
    }

    if (run_hisat2 || run_kallisto || run_rsem || run_rsem_make_transcript_fasta || run_salmon || run_star) {
        // Do not run GFFREAD if the condition is false
        gff_gffread = gff.map { meta, gff_ -> meta.run_gffread ? [meta, gff_] : null }

        GFFREAD(
            gff_gffread,
            [],
        )

        versions = versions.mix(GFFREAD.out.versions)

        gtf = gtf
            .mix(GFFREAD.out.gtf)
            .groupTuple()
            .map { meta, gtf_ -> gtf_[1] ? [meta, gtf_[1]] : [meta, gtf_] }

        if (run_hisat2 || run_hisat2_extractsplicesites) {
            // Do not run HISAT2_EXTRACTSPLICESITES if the condition is false
            gtf_hisat2 = gtf.map { meta, gtf_ -> meta.run_hisat2 ? [meta, gtf_] : null }

            HISAT2_EXTRACTSPLICESITES(gtf_hisat2)

            versions = versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)
            splice_sites = splice_sites.mix(HISAT2_EXTRACTSPLICESITES.out.txt)

            if (run_hisat2) {
                // Do not run HISAT2_BUILD if the condition is false
                fasta_hisat2 = fasta.map { meta, fasta_ -> meta.run_hisat2 ? [meta, fasta_] : null }

                HISAT2_BUILD(
                    fasta_hisat2,
                    gtf,
                    splice_sites,
                )

                hisat2_index = HISAT2_BUILD.out.index

                versions = versions.mix(HISAT2_BUILD.out.versions)
            }
        }

        if (run_kallisto || run_rsem_make_transcript_fasta || run_salmon) {
            // Do not run MAKE_TRANSCRIPTS_FASTA if the condition is false
            fasta_make_transcripts_fasta = fasta.map { meta, fasta_ -> meta.run_rsem_make_transcript_fasta ? [meta, fasta_] : null }

            MAKE_TRANSCRIPTS_FASTA(
                fasta_make_transcripts_fasta,
                gtf,
            )
            versions = versions.mix(MAKE_TRANSCRIPTS_FASTA.out.versions)

            transcript_fasta = transcript_fasta.mix(MAKE_TRANSCRIPTS_FASTA.out.transcript_fasta)

            if (run_kallisto) {
                // Do not run KALLISTO_INDEX if the condition is false
                transcript_fasta_kallisto = transcript_fasta.map { meta, transcript_fasta_ -> meta.run_kallisto ? [meta, transcript_fasta_] : null }

                KALLISTO_INDEX(transcript_fasta_kallisto)

                kallisto_index = KALLISTO_INDEX.out.index
                versions = versions.mix(KALLISTO_INDEX.out.versions)
            }

            if (run_salmon) {
                // Do not run SALMON_INDEX if the condition is false
                fasta_salmon = fasta.map { meta, fasta_ -> meta.run_salmon ? [meta, fasta_] : null }

                // Do not run SALMON_INDEX if the condition is false
                transcript_fasta_salmon = transcript_fasta.map { meta, transcript_fasta_ -> meta.run_salmon ? [meta, transcript_fasta_] : null }

                SALMON_INDEX(
                    fasta_salmon,
                    transcript_fasta_salmon,
                )

                salmon_index = SALMON_INDEX.out.index
                versions = versions.mix(SALMON_INDEX.out.versions)
            }
        }

        if (run_rsem) {
            // Do not run RSEM_PREPAREREFERENCE_GENOME if the condition is false
            fasta_rsem = fasta.map { meta, fasta_ -> meta.run_rsem ? [meta, fasta_] : null }

            RSEM_PREPAREREFERENCE_GENOME(fasta_rsem, gtf)

            rsem_index = RSEM_PREPAREREFERENCE_GENOME.out.index
            versions = versions.mix(RSEM_PREPAREREFERENCE_GENOME.out.versions)
        }

        if (run_star) {
            // Do not run STAR_GENOMEGENERATE if the condition is false
            fasta_star = fasta.map { meta, fasta_ -> meta.run_star ? [meta, fasta_] : null }

            STAR_GENOMEGENERATE(fasta_star, gtf)

            star_index = STAR_GENOMEGENERATE.out.index
            versions = versions.mix(STAR_GENOMEGENERATE.out.versions)
        }
    }

    emit:
    bowtie1_index    // channel: [meta, BowtieIndex/]
    bowtie2_index    // channel: [meta, Bowtie2Index/]
    fasta_fai        // channel: [meta, *.fa(sta).fai]
    fasta_sizes      // channel: [meta, *.fa(sta).sizes]
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
