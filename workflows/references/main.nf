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

// include { BBMAP_BBSPLIT                                         } from '../../modules/nf-core/bbmap/bbsplit'
// include { CUSTOM_CATADDITIONALFASTA                             } from '../../modules/nf-core/custom/catadditionalfasta'
// include { SORTMERNA as SORTMERNA_INDEX                          } from '../../modules/nf-core/sortmerna'

workflow REFERENCES {
    take:
    reference // fasta, gff, gtf, splice_sites, transcript_fasta
    tools     // bowtie|bowtie2|bwamem1|bwamem2|createsequencedictionary|dragmap|faidx|gffread|intervals|hisat2|hisat2_extractsplicesites|kallisto|msisensorpro|rsem|rsem_make_transcripts_fasta|salmon|star

    main:
    bed_intervals = Channel.empty()
    bowtie1 = Channel.empty()
    bowtie2 = Channel.empty()
    faidx = Channel.empty()
    gffread = Channel.empty()
    hisat2 = Channel.empty()
    hisat2_splice_sites = Channel.empty()
    kallisto = Channel.empty()
    msisensorpro = Channel.empty()
    rsem = Channel.empty()
    rsem_transcript_fasta = Channel.empty()
    salmon = Channel.empty()
    sizes = Channel.empty()
    star = Channel.empty()
    versions = Channel.empty()

    input = reference.multiMap { meta, intervals, fasta, fasta_dict, fasta_fai, fasta_sizes, gff, gtf, splice_sites, transcript_fasta, readme, bed12, mito_name, macs_gsize ->
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
        intervals: [meta, intervals]
        bwamem1_fasta: tools.contains('bwamem1') ? [meta, fasta] : [[:], []]
        bwamem2_fasta: tools.contains('bwamem2') ? [meta, fasta] : [[:], []]
        createsequencedictionary_fasta: tools.contains('createsequencedictionary') ? [meta, fasta] : [[:], []]
        dragmap_fasta: tools.contains('dragmap') ? [meta, fasta] : [[:], []]
        fasta_samtools: ((tools.contains('faidx') || tools.contains('sizes')) && !(fasta_fai || fasta_sizes)) || (tools.contains('intervals') && !(fasta_fai || intervals)) ? [meta, fasta] : [[:], []]
    }
    // I should be able to output null instead of `[[:], []] and have that registered as an empty channel and not trigger downstream processes
    // but not working currently

    if (tools && tools.split(',').contains('bowtie1')) {
        BOWTIE1_BUILD(input.fasta)

        bowtie1 = BOWTIE1_BUILD.out.index
        versions = versions.mix(BOWTIE1_BUILD.out.versions)
    }

    if (tools && tools.split(',').contains('bowtie2')) {
        BOWTIE2_BUILD(input.fasta)

        bowtie2 = BOWTIE2_BUILD.out.index
        versions = versions.mix(BOWTIE2_BUILD.out.versions)
    }

    // the whole map -> null should be removed once I managed to make it work properly
    BWAMEM1_INDEX(
        input.bwamem1_fasta.map { meta, file ->
            return file ? [meta, file] : null
        }
    )

    BWAMEM2_INDEX(
        input.bwamem2_fasta.map { meta, file ->
            return file ? [meta, file] : null
        }
    )

    DRAGMAP_HASHTABLE(
        input.dragmap_fasta.map { meta, file ->
            return file ? [meta, file] : null
        }
    )

    GATK4_CREATESEQUENCEDICTIONARY(
        input.createsequencedictionary_fasta.map { meta, file ->
            return file ? [meta, file] : null
        }
    )

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
    faidx = input.fasta_fai
        .mix(SAMTOOLS_FAIDX.out.fai)
        .groupTuple()
        .map { meta, file ->
            return file[1] ? [meta, file[1]] : [meta, file]
        }

    // TODO: be smarter about input assets
    //   Here we either mix+GT an empty channel (either no output or no input sizes) with the sizes return sizes
    //   And we filter out the empty value
    sizes = input.fasta_sizes
        .mix(SAMTOOLS_FAIDX.out.sizes)
        .groupTuple()
        .map { meta, file ->
            return file[1] ? [meta, file[1]] : [meta, file]
        }

    faidx_intervals = faidx
        .mix(input.intervals)
        .groupTuple()
        .map { meta, file ->
            return file[1] ? null : [meta, file]
        }

    BUILD_INTERVALS(faidx_intervals, [])

    if (tools.contains('gffread')) {
        GFFREAD(input.gff, [])

        gffread = GFFREAD.out.gtf.map { it[1] }
        versions = versions.mix(GFFREAD.out.versions)
    }

    if (tools.contains('hisat2')) {
        // TODO: be smarter about input assets
        //   Here we either return an empty channel if we have a splice_sites so that HISAT2_EXTRACTSPLICESITES is not triggered
        //   Or we return the provided gtf so that HISAT2_EXTRACTSPLICESITES is run
        gtf_hisat2 = input.gtf
            .join(input.splice_sites)
            .groupTuple()
            .map { meta, gtf, splice_sites ->
                return splice_sites[0][0] ? null : [meta, gtf]
            }

        HISAT2_EXTRACTSPLICESITES(gtf_hisat2)
        versions = versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)

        // TODO: be smarter about input assets
        //   Here we either mix+GT an empty channel (either no output or no input splice_sites) with the splice_sites return splice_sites
        //   And we filter out the empty value
        hisat2_extractsplicesites = input.splice_sites
            .mix(HISAT2_EXTRACTSPLICESITES.out.txt)
            .groupTuple()
            .map { meta, txt ->
                return txt[1] ? [meta, txt[1]] : [meta, txt]
            }

        if (tools && tools.split(',').contains('hisat2')) {
            HISAT2_BUILD(input.fasta, input.gtf, hisat2_extractsplicesites)

            hisat2 = HISAT2_BUILD.out.index
            versions = versions.mix(HISAT2_BUILD.out.versions)
        }
    }

    if (tools.contains('kallisto') || tools.contains('rsem_make_transcript_fasta') || tools.contains('salmon')) {
        // TODO: be smarter about input assets
        //   Here we either return an empty channel if we have a transcript_fasta so that MAKE_TRANSCRIPTS_FASTA is not triggered
        //   Or we return the provided gtf so that MAKE_TRANSCRIPTS_FASTA is run
        gtf_rsem = input.gtf
            .join(input.transcript_fasta)
            .groupTuple()
            .map { meta, gtf, transcript_fasta ->
                return transcript_fasta[0][0] ? null : [meta, gtf]
            }

        MAKE_TRANSCRIPTS_FASTA(input.fasta, gtf_rsem)
        versions = versions.mix(MAKE_TRANSCRIPTS_FASTA.out.versions)

        // TODO: be smarter about input assets
        //   Here we either mix+GT an empty channel (either no output or no input transcript_fasta) with the transcript_fasta return transcript_fasta
        //   And we filter out the empty value
        rsem_transcript_fasta = input.transcript_fasta
            .mix(MAKE_TRANSCRIPTS_FASTA.out.transcript_fasta)
            .groupTuple()
            .map { meta, txt ->
                return txt[1] ? [meta, txt[1]] : [meta, txt]
            }

        if (tools.contains('kallisto')) {
            KALLISTO_INDEX(rsem_transcript_fasta)

            kallisto = KALLISTO_INDEX.out.index
            versions = versions.mix(KALLISTO_INDEX.out.versions)
        }

        if (tools.contains('salmon')) {
            SALMON_INDEX(input.fasta, rsem_transcript_fasta)

            salmon = SALMON_INDEX.out.index
            versions = versions.mix(SALMON_INDEX.out.versions)
        }
    }

    if (tools.contains('msisensorpro')) {
        MSISENSORPRO_SCAN(input.fasta)

        msisensorpro = MSISENSORPRO_SCAN.out.list
        versions = versions.mix(MSISENSORPRO_SCAN.out.versions)
    }

    if (tools && tools.split(',').contains('rsem')) {
        RSEM_PREPAREREFERENCE_GENOME(input.fasta, input.gtf)

        rsem = RSEM_PREPAREREFERENCE_GENOME.out.index
        versions = versions.mix(RSEM_PREPAREREFERENCE_GENOME.out.versions)
    }

    if (tools.contains('star')) {
        STAR_GENOMEGENERATE(input.fasta, input.gtf)

        star = STAR_GENOMEGENERATE.out.index
        versions = versions.mix(STAR_GENOMEGENERATE.out.versions)
    }

    versions = versions.mix(BUILD_INTERVALS.out.versions)
    versions = versions.mix(BWAMEM1_INDEX.out.versions)
    versions = versions.mix(BWAMEM2_INDEX.out.versions)
    versions = versions.mix(DRAGMAP_HASHTABLE.out.versions)
    versions = versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    versions = versions.mix(SAMTOOLS_FAIDX.out.versions)

    emit:
    bowtie1               = bowtie1
    bowtie2               = bowtie2
    bwamem1               = BWAMEM1_INDEX.out.index
    bwamem2               = BWAMEM2_INDEX.out.index
    dict                  = GATK4_CREATESEQUENCEDICTIONARY.out.dict
    dragmap               = DRAGMAP_HASHTABLE.out.hashmap
    faidx                 = faidx
    fasta                 = input.fasta
    gffread               = gffread
    hisat2                = hisat2
    hisat2_splice_sites   = hisat2_splice_sites
    intervals             = bed_intervals
    kallisto              = kallisto
    msisensorpro          = msisensorpro
    rsem                  = rsem
    rsem_transcript_fasta = rsem_transcript_fasta
    salmon                = salmon
    sizes                 = sizes
    star                  = star
    versions              = versions
}
