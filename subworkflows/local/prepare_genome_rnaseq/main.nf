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
    run_sizes                      // boolean: true/false
    run_star                       // boolean: true/false

    main:
    def join_by_meta_id = { channel1, channel2, channel3 = null ->
        if (channel3) {
            channel1
                .map { meta, content1_ -> [meta.id, content1_, meta] }
                .join(channel2.map { meta, content2_ -> [meta.id, content2_, meta] })
                .join(channel3.map { meta, content3_ -> [meta.id, content3_, meta] })
                .map { _id, content1_, meta1, content2_, meta2, content3_, meta3 -> [meta3 + meta2 + meta1, content1_, content2_, content3_] }
        }
        else {
            channel1
                .map { meta, content1_ -> [meta.id, content1_, meta] }
                .join(channel2.map { meta, content2_ -> [meta.id, content2_, meta] })
                .map { _id, content1_, meta1, content2_, meta2 -> [meta2 + meta1, content1_, content2_] }
        }
    }

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
        BOWTIE1_BUILD(fasta)

        bowtie1_index = BOWTIE1_BUILD.out.index
        versions = versions.mix(BOWTIE1_BUILD.out.versions)
    }

    if (run_bowtie2) {
        BOWTIE2_BUILD(fasta)

        bowtie2_index = BOWTIE2_BUILD.out.index
        versions = versions.mix(BOWTIE2_BUILD.out.versions)
    }

    if (run_faidx || run_sizes) {
        SAMTOOLS_FAIDX(fasta.map { meta, fasta_ -> [meta, fasta_, []] }, run_sizes)

        fasta_fai = fasta_fai.mix(SAMTOOLS_FAIDX.out.fai)
        fasta_sizes = SAMTOOLS_FAIDX.out.sizes
        versions = versions.mix(SAMTOOLS_FAIDX.out.versions)
    }

    if (run_hisat2 || run_kallisto || run_rsem || run_rsem_make_transcript_fasta || run_salmon || run_star) {
        GFFREAD(gff.map { meta, gff_ -> [meta, [], gff_] })

        versions = versions.mix(GFFREAD.out.versions)

        gtf = gtf
            .mix(GFFREAD.out.gtf)
            .groupTuple()
            .map { meta, gtf_ -> gtf_[1] ? [meta, gtf_[1]] : [meta, gtf_[0]] }

        if (run_hisat2 || run_hisat2_extractsplicesites) {
            HISAT2_EXTRACTSPLICESITES(gtf)

            versions = versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)
            splice_sites = splice_sites.mix(HISAT2_EXTRACTSPLICESITES.out.txt)

            if (run_hisat2) {
                HISAT2_BUILD(join_by_meta_id(fasta, gtf, splice_sites))

                hisat2_index = HISAT2_BUILD.out.index

                versions = versions.mix(HISAT2_BUILD.out.versions)
            }
        }

        if (run_kallisto || run_rsem_make_transcript_fasta || run_salmon) {
            MAKE_TRANSCRIPTS_FASTA(join_by_meta_id(fasta, gtf))

            versions = versions.mix(MAKE_TRANSCRIPTS_FASTA.out.versions)

            transcript_fasta = transcript_fasta.mix(MAKE_TRANSCRIPTS_FASTA.out.transcript_fasta)

            if (run_kallisto) {
                KALLISTO_INDEX(transcript_fasta)

                kallisto_index = KALLISTO_INDEX.out.index
                versions = versions.mix(KALLISTO_INDEX.out.versions)
            }

            if (run_salmon) {
                SALMON_INDEX(join_by_meta_id(fasta, transcript_fasta))

                salmon_index = SALMON_INDEX.out.index
                versions = versions.mix(SALMON_INDEX.out.versions)
            }
        }

        if (run_rsem) {
            RSEM_PREPAREREFERENCE_GENOME(join_by_meta_id(fasta, gtf))

            rsem_index = RSEM_PREPAREREFERENCE_GENOME.out.index
            versions = versions.mix(RSEM_PREPAREREFERENCE_GENOME.out.versions)
        }

        if (run_star) {
            STAR_GENOMEGENERATE(join_by_meta_id(fasta, gtf))

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
