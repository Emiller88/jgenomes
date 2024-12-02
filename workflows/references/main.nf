include { BOWTIE2_BUILD                                         } from '../../modules/nf-core/bowtie2/build'
include { BOWTIE_BUILD as BOWTIE1_BUILD                         } from '../../modules/nf-core/bowtie/build'
include { BWAMEM2_INDEX                                         } from '../../modules/nf-core/bwamem2/index'
include { BWA_INDEX as BWAMEM1_INDEX                            } from '../../modules/nf-core/bwa/index'
include { DRAGMAP_HASHTABLE                                     } from '../../modules/nf-core/dragmap/hashtable'
include { GATK4_CREATESEQUENCEDICTIONARY                        } from '../../modules/nf-core/gatk4/createsequencedictionary'
include { GAWK as BUILD_INTERVALS                               } from '../../modules/nf-core/gawk'
include { GFFREAD                                               } from '../../modules/nf-core/gffread'
include { HISAT2_BUILD                                          } from '../../modules/nf-core/hisat2/build'
include { HISAT2_EXTRACTSPLICESITES                             } from '../../modules/nf-core/hisat2/extractsplicesites'
include { KALLISTO_INDEX                                        } from '../../modules/nf-core/kallisto/index'
include { MSISENSORPRO_SCAN                                     } from '../../modules/nf-core/msisensorpro/scan'
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA       } from '../../modules/nf-core/rsem/preparereference'
include { RSEM_PREPAREREFERENCE as RSEM_PREPAREREFERENCE_GENOME } from '../../modules/nf-core/rsem/preparereference'
include { SALMON_INDEX                                          } from '../../modules/nf-core/salmon/index'
include { SAMTOOLS_FAIDX                                        } from '../../modules/nf-core/samtools/faidx'
include { STAR_GENOMEGENERATE                                   } from '../../modules/nf-core/star/genomegenerate'
include { TABIX_TABIX                                           } from '../../modules/nf-core/tabix/tabix'

// include { BBMAP_BBSPLIT                                         } from '../../modules/nf-core/bbmap/bbsplit'
// include { CUSTOM_CATADDITIONALFASTA                             } from '../../modules/nf-core/custom/catadditionalfasta'
// include { SORTMERNA as SORTMERNA_INDEX                          } from '../../modules/nf-core/sortmerna'

workflow REFERENCES {
    take:
    reference // fasta, gff, gtf, splice_sites, transcript_fasta
    tools     // bowtie|bowtie2|bwamem1|bwamem2|createsequencedictionary|dragmap|faidx|gffread|intervals|hisat2|hisat2_extractsplicesites|kallisto|msisensorpro|rsem|rsem_make_transcripts_fasta|salmon|star|tabix

    main:
    ch_bowtie1 = Channel.empty()
    ch_bowtie2 = Channel.empty()
    ch_fasta_fai = Channel.empty()
    ch_gff_gtf = Channel.empty()
    ch_hisat2 = Channel.empty()
    ch_hisat2_splice_sites = Channel.empty()
    ch_intervals_bed = Channel.empty()
    ch_kallisto = Channel.empty()
    ch_msisensorpro = Channel.empty()
    ch_rsem = Channel.empty()
    ch_rsem_transcript_fasta = Channel.empty()
    ch_salmon = Channel.empty()
    ch_sizes = Channel.empty()
    ch_star = Channel.empty()
    ch_vcf_tbi = Channel.empty()
    versions = Channel.empty()

    input = reference.multiMap { meta, intervals_bed, fasta, fasta_dict, fasta_fai, fasta_sizes, gff, gtf, splice_sites, transcript_fasta, vcf, readme, bed12, mito_name, macs_gsize ->
        fasta: [meta, fasta]
        fasta_dict: [meta, fasta_dict]
        fasta_fai: [meta, fasta_fai]
        fasta_sizes: [meta, fasta_sizes]
        gff: [meta, gff]
        gtf: [meta, gtf]
        splice_sites: [meta, splice_sites]
        transcript_fasta: [meta, transcript_fasta]
        readme: [meta, readme]
        bed12: [meta, bed12]
        mito_name: [meta, mito_name]
        macs_gsize: [meta, macs_gsize]
        intervals_bed: [meta, intervals_bed]
        bwamem1_fasta: tools.contains('bwamem1') && fasta ? [meta, file(fasta)] : [[:], []]
        bwamem2_fasta: tools.contains('bwamem2') && fasta ? [meta, file(fasta)] : [[:], []]
        fasta_msisensorpro: tools.contains('msisensorpro') && fasta ? [meta, file(fasta)] : [[:], []]
        createsequencedictionary_fasta: tools.contains('createsequencedictionary') && fasta ? [meta, file(fasta)] : [[:], []]
        dragmap_fasta: tools.contains('dragmap') && fasta ? [meta, file(fasta)] : [[:], []]
        fasta_samtools: ((tools.contains('faidx') || tools.contains('sizes')) && !(fasta_fai || fasta_sizes) && fasta) || (tools.contains('intervals') && !(fasta_fai || intervals_bed) && fasta) ? [meta, file(fasta)] : [[:], []]
        gff_gffread: !gtf && gff && (tools.contains('gffread') || tools.contains('hisat2') || tools.contains('kallisto') || tools.contains('rsem') || tools.contains('salmon') || tools.contains('star')) ? [meta, file(gff)] : [[:], []]
        vcf: tools.contains('tabix') && vcf ? [meta, file(vcf)] : [[:], []]
    }
    // I should be able to output null instead of `[[:], []] and have that registered as an empty channel and not trigger downstream processes
    // but not working currently

    if (tools && tools.split(',').contains('bowtie1')) {
        BOWTIE1_BUILD(
            input.fasta.map { meta, file ->
                return file ? [meta, file] : null
            }
        )

        ch_bowtie1 = BOWTIE1_BUILD.out.index
        versions = versions.mix(BOWTIE1_BUILD.out.versions)
    }

    if (tools && tools.split(',').contains('bowtie2')) {
        BOWTIE2_BUILD(
            input.fasta.map { meta, file ->
                return file ? [meta, file] : null
            }
        )

        ch_bowtie2 = BOWTIE2_BUILD.out.index
        versions = versions.mix(BOWTIE2_BUILD.out.versions)
    }

    // the whole map -> null should be removed once I managed to make it work properly
    BWAMEM1_INDEX(
        input.bwamem1_fasta.map { meta, file ->
            return file ? [meta, file] : null
        }
    )
    ch_bwamem1 = BWAMEM1_INDEX.out.index

    BWAMEM2_INDEX(
        input.bwamem2_fasta.map { meta, file ->
            return file ? [meta, file] : null
        }
    )
    ch_bwamem2 = BWAMEM2_INDEX.out.index

    DRAGMAP_HASHTABLE(
        input.dragmap_fasta.map { meta, file ->
            return file ? [meta, file] : null
        }
    )
    ch_dragmap = DRAGMAP_HASHTABLE.out.hashmap

    GATK4_CREATESEQUENCEDICTIONARY(
        input.createsequencedictionary_fasta.map { meta, file ->
            return file ? [meta, file] : null
        }
    )
    ch_fasta_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict

    SAMTOOLS_FAIDX(
        input.fasta_samtools.map { meta, file ->
            return file ? [meta, file] : null
        },
        [[id: 'no_fai'], []],
        tools.contains('sizes')
    )

    // TODO: be smarter about input assets
    //   Here we either mix+GT an empty channel (either no output or no input faidx) with the faidx return faidx
    //   And we filter out the empty value
    ch_fasta_fai = input.fasta_fai
        .mix(SAMTOOLS_FAIDX.out.fai)
        .groupTuple()
        .map { meta, file ->
            return file[1] ? [meta, file[1]] : [meta, file]
        }

    // TODO: be smarter about input assets
    //   Here we either mix+GT an empty channel (either no output or no input sizes) with the sizes return sizes
    //   And we filter out the empty value
    ch_sizes = input.fasta_sizes
        .mix(SAMTOOLS_FAIDX.out.sizes)
        .groupTuple()
        .map { meta, file ->
            return file[1] ? [meta, file[1]] : [meta, file]
        }

    ch_fasta_fai_intervals_bed = input.intervals_bed
        .mix(ch_fasta_fai)
        .groupTuple()
        .map { meta, file ->
            return !tools.contains('intervals') ? null : file[0][0] ? null : file.flatten() ? [meta, file.flatten()] : [meta, file]
        }

    BUILD_INTERVALS(ch_fasta_fai_intervals_bed, [])
    ch_intervals_bed = BUILD_INTERVALS.out.output

    TABIX_TABIX(
        input.vcf.map { meta, file ->
            return file ? [meta, file] : null
        }.transpose()
    )
    ch_vcf_tbi = TABIX_TABIX.out.tbi

    GFFREAD(
        input.gff_gffread.map { meta, file ->
            return file ? [meta, file] : null
        },
        []
    )

    ch_gff_gtf = input.gtf
        .mix(GFFREAD.out.gtf)
        .groupTuple()
        .map { meta, file ->
            return file[1] ? [meta, file[1]] : [meta, file]
        }

    if (tools.contains('hisat2')) {
        // TODO: be smarter about input assets
        //   Here we either return an empty channel if we have a splice_sites so that HISAT2_EXTRACTSPLICESITES is not triggered
        //   Or we return the provided gtf so that HISAT2_EXTRACTSPLICESITES is run
        ch_gtf_hisat2 = ch_gff_gtf
            .join(input.splice_sites)
            .groupTuple()
            .map { meta, gtf, splice_sites ->
                return splice_sites[0][0] ? null : [meta, gtf]
            }

        HISAT2_EXTRACTSPLICESITES(ch_gtf_hisat2)
        versions = versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)

        // TODO: be smarter about input assets
        //   Here we either mix+GT an empty channel (either no output or no input splice_sites) with the splice_sites return splice_sites
        //   And we filter out the empty value
        ch_hisat2_splice_sites = input.splice_sites
            .mix(HISAT2_EXTRACTSPLICESITES.out.txt)
            .groupTuple()
            .map { meta, txt ->
                return txt[1] ? [meta, txt[1]] : [meta, txt]
            }

        if (tools && tools.split(',').contains('hisat2')) {
            HISAT2_BUILD(
                input.fasta.map { meta, file ->
                    return file ? [meta, file] : null
                },
                ch_gff_gtf,
                ch_hisat2_splice_sites
            )

            ch_hisat2 = HISAT2_BUILD.out.index
            versions = versions.mix(HISAT2_BUILD.out.versions)
        }
    }

    if (tools.contains('kallisto') || tools.contains('rsem_make_transcript_fasta') || tools.contains('salmon')) {
        // TODO: be smarter about input assets
        //   Here we either return an empty channel if we have a transcript_fasta so that MAKE_TRANSCRIPTS_FASTA is not triggered
        //   Or we return the provided gtf so that MAKE_TRANSCRIPTS_FASTA is run
        ch_gtf_rsem = ch_gff_gtf
            .join(input.transcript_fasta)
            .groupTuple()
            .map { meta, gtf, transcript_fasta ->
                return transcript_fasta[0][0] ? null : [meta, gtf]
            }

        MAKE_TRANSCRIPTS_FASTA(
            input.fasta.map { meta, file ->
                return file ? [meta, file] : null
            },
            ch_gtf_rsem
        )
        versions = versions.mix(MAKE_TRANSCRIPTS_FASTA.out.versions)

        // TODO: be smarter about input assets
        //   Here we either mix+GT an empty channel (either no output or no input transcript_fasta) with the transcript_fasta return transcript_fasta
        //   And we filter out the empty value
        ch_rsem_transcript_fasta = input.transcript_fasta
            .mix(MAKE_TRANSCRIPTS_FASTA.out.transcript_fasta)
            .groupTuple()
            .map { meta, txt ->
                return txt[1] ? [meta, txt[1]] : [meta, txt]
            }

        if (tools.contains('kallisto')) {
            KALLISTO_INDEX(ch_rsem_transcript_fasta)

            ch_kallisto = KALLISTO_INDEX.out.index
            versions = versions.mix(KALLISTO_INDEX.out.versions)
        }

        if (tools.contains('salmon')) {
            SALMON_INDEX(
                input.fasta.map { meta, file ->
                    return file ? [meta, file] : null
                },
                ch_rsem_transcript_fasta
            )

            ch_salmon = SALMON_INDEX.out.index
            versions = versions.mix(SALMON_INDEX.out.versions)
        }
    }

    if (tools.contains('msisensorpro')) {
        MSISENSORPRO_SCAN(
            input.fasta_msisensorpro.map { meta, file ->
                return file ? [meta, file] : null
            }
        )

        ch_msisensorpro = MSISENSORPRO_SCAN.out.list
        versions = versions.mix(MSISENSORPRO_SCAN.out.versions)
    }

    if (tools && tools.split(',').contains('rsem')) {
        RSEM_PREPAREREFERENCE_GENOME(
            input.fasta.map { meta, file ->
                return file ? [meta, file] : null
            },
            ch_gff_gtf
        )

        ch_rsem = RSEM_PREPAREREFERENCE_GENOME.out.index
        versions = versions.mix(RSEM_PREPAREREFERENCE_GENOME.out.versions)
    }

    if (tools.contains('star')) {
        STAR_GENOMEGENERATE(
            input.fasta.map { meta, file ->
                return file ? [meta, file] : null
            },
            ch_gff_gtf
        )

        ch_star = STAR_GENOMEGENERATE.out.index
        versions = versions.mix(STAR_GENOMEGENERATE.out.versions)
    }

    // versions
    versions = versions.mix(BUILD_INTERVALS.out.versions)
    versions = versions.mix(BWAMEM1_INDEX.out.versions)
    versions = versions.mix(BWAMEM2_INDEX.out.versions)
    versions = versions.mix(DRAGMAP_HASHTABLE.out.versions)
    versions = versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    versions = versions.mix(GFFREAD.out.versions)
    versions = versions.mix(SAMTOOLS_FAIDX.out.versions)
    versions = versions.mix(TABIX_TABIX.out.versions)

    // input fasta
    ch_fasta = input.fasta

    emit:
    bowtie1               = ch_bowtie1
    bowtie2               = ch_bowtie2
    bwamem1               = ch_bwamem1
    bwamem2               = ch_bwamem2
    dragmap               = ch_dragmap
    fasta                 = ch_fasta
    fasta_dict            = ch_fasta_dict
    fasta_fai             = ch_fasta_fai
    gff_gtf               = ch_gff_gtf
    hisat2                = ch_hisat2
    hisat2_splice_sites   = ch_hisat2_splice_sites
    intervals_bed         = ch_intervals_bed
    kallisto              = ch_kallisto
    msisensorpro          = ch_msisensorpro
    rsem                  = ch_rsem
    rsem_transcript_fasta = ch_rsem_transcript_fasta
    salmon                = ch_salmon
    sizes                 = ch_sizes
    star                  = ch_star
    vcf_tbi               = ch_vcf_tbi
    versions              = versions

    publish:
    ch_bowtie1 >> 'bowtie1'
    ch_bowtie2 >> 'bowtie2'
    ch_bwamem1 >> 'bwamem1'
    ch_bwamem2 >> 'bwamem2'
    ch_dragmap >> 'dragmap'
    ch_fasta >> 'fasta'
    ch_fasta_dict >> 'fasta_dict'
    ch_fasta_fai >> 'fasta_fai'
    ch_gff_gtf >> 'gffread'
    ch_hisat2 >> 'hisat2'
    ch_hisat2_splice_sites >> 'hisat2'
    ch_intervals_bed >> 'intervals'
    ch_kallisto >> 'kallisto'
    ch_msisensorpro >> 'msisensorpro'
    ch_rsem >> 'rsem'
    ch_rsem_transcript_fasta >> 'make'
    ch_salmon >> 'salmon'
    ch_sizes >> 'fasta_sizes'
    ch_star >> 'star'
    ch_vcf_tbi >> 'tabix_vcf_tbi'
}
