include { BBMAP_BBSPLIT                                         } from '../../modules/nf-core/bbmap/bbsplit'
include { BOWTIE2_BUILD                                         } from '../../modules/nf-core/bowtie2/build'
include { BOWTIE_BUILD as BOWTIE1_BUILD                         } from '../../modules/nf-core/bowtie/build'
include { BWAMEM2_INDEX                                         } from '../../modules/nf-core/bwamem2/index'
include { BWA_INDEX as BWAMEM1_INDEX                            } from '../../modules/nf-core/bwa/index'
include { CUSTOM_CATADDITIONALFASTA                             } from '../../modules/nf-core/custom/catadditionalfasta'
include { DRAGMAP_HASHTABLE                                     } from '../../modules/nf-core/dragmap/hashtable'
include { GATK4_CREATESEQUENCEDICTIONARY                        } from '../../modules/nf-core/gatk4/createsequencedictionary'
include { GFFREAD                                               } from '../../modules/nf-core/gffread'
include { HISAT2_BUILD                                          } from '../../modules/nf-core/hisat2/build'
include { HISAT2_EXTRACTSPLICESITES                             } from '../../modules/nf-core/hisat2/extractsplicesites'
include { KALLISTO_INDEX                                        } from '../../modules/nf-core/kallisto/index'
include { MSISENSORPRO_SCAN                                     } from '../../modules/nf-core/msisensorpro/scan'
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA       } from '../../modules/nf-core/rsem/preparereference'
include { RSEM_PREPAREREFERENCE as RSEM_PREPAREREFERENCE_GENOME } from '../../modules/nf-core/rsem/preparereference'
include { SALMON_INDEX                                          } from '../../modules/nf-core/salmon/index'
include { SAMTOOLS_FAIDX                                        } from '../../modules/nf-core/samtools/faidx'
include { SORTMERNA as SORTMERNA_INDEX                          } from '../../modules/nf-core/sortmerna'
include { STAR_GENOMEGENERATE                                   } from '../../modules/nf-core/star/genomegenerate'

workflow REFERENCES {
    take:
    reference // fasta, gff, gtf, splice_sites, transcript_fasta
    tools     // bowtie|bowtie2|bwamem1|bwamem2|createsequencedictionary|dragmap|faidx|gffread|hisat2|hisat2_extractsplicesites|kallisto|msisensorpro|rsem|rsem_make_transcripts_fasta|salmon|star

    main:
    bowtie1 = Channel.empty()
    bowtie2 = Channel.empty()
    bwamem1 = Channel.empty()
    bwamem2 = Channel.empty()
    dict = Channel.empty()
    dragmap = Channel.empty()
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

    input = reference.multiMap { meta, fasta, fasta_dict, fasta_fai, fasta_sizes, gff, gtf, splice_sites, transcript_fasta, readme, bed12, mito_name, macs_gsize ->
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
    }

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

    if (tools && tools.split(',').contains('bwamem1')) {
        BWAMEM1_INDEX(input.fasta)

        bwamem1 = BWAMEM1_INDEX.out.index
        versions = versions.mix(BWAMEM1_INDEX.out.versions)
    }

    if (tools && tools.split(',').contains('bwamem2')) {
        BWAMEM2_INDEX(input.fasta)

        bwamem2 = BWAMEM2_INDEX.out.index
        versions = versions.mix(BWAMEM2_INDEX.out.versions)
    }

    if (tools && tools.split(',').contains('dragmap')) {
        DRAGMAP_HASHTABLE(input.fasta)

        dragmap = DRAGMAP_HASHTABLE.out.hashmap
        versions = versions.mix(DRAGMAP_HASHTABLE.out.versions)
    }

    if (tools && tools.split(',').contains('createsequencedictionary')) {
        GATK4_CREATESEQUENCEDICTIONARY(input.fasta)

        dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict
        versions = versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    }

    if (tools && (tools.split(',').contains('faidx') || tools.split(',').contains('sizes'))) {
        // TODO: be smarter about input assets
        //   Here we either return an empty channel if we have a fai and sizes so that SAMTOOLS_FAIDX is not triggered
        //   Or we return the provided fasta so that SAMTOOLS_FAIDX is run
        fasta_samtools = input.fasta
            .join(input.fasta_fai)
            .join(input.fasta_sizes)
            .groupTuple()
            .map { meta, fasta, fai, sizes ->
                return fai[0][0] && sizes[0][0] ? null : [meta, fasta]
            }

        generate_sizes = tools.split(',').contains('sizes')
        SAMTOOLS_FAIDX(fasta_samtools, [[id: 'no_fai'], []], generate_sizes)


        // TODO: be smarter about input assets
        //   Here we either mix+GT an empty channel (either no output or no input faidx) with the faidx return faidx
        //   And we filter out the empty value
        faidx = input.fasta_fai
            .mix(SAMTOOLS_FAIDX.out.fai)
            .groupTuple()
            .map { meta, fai ->
                return fai[1] ? [meta, fai[1]] : [meta, fai]
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

        versions = versions.mix(SAMTOOLS_FAIDX.out.versions)
    }

    if (tools && tools.split(',').contains('gffread')) {
        GFFREAD(input.gff, [])

        gffread = GFFREAD.out.gtf.map { it[1] }
        versions = versions.mix(GFFREAD.out.versions)
    }

    if (tools && (tools.split(',').contains('hisat2') || tools.split(',').contains('hisat2_extractsplicesites'))) {
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

        if (tools.split(',').contains('hisat2')) {
            HISAT2_BUILD(input.fasta, input.gtf, hisat2_extractsplicesites)

            hisat2 = HISAT2_BUILD.out.index
            versions = versions.mix(HISAT2_BUILD.out.versions)
        }
    }

    if (tools && tools.split(',').contains('kallisto') || tools.split(',').contains('rsem_make_transcript_fasta') || tools.split(',').contains('salmon')) {
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

        if (tools.split(',').contains('kallisto')) {
            KALLISTO_INDEX(rsem_transcript_fasta)

            kallisto = KALLISTO_INDEX.out.index
            versions = versions.mix(KALLISTO_INDEX.out.versions)
        }

        if (tools.split(',').contains('salmon')) {
            SALMON_INDEX(input.fasta, rsem_transcript_fasta)

            salmon = SALMON_INDEX.out.index
            versions = versions.mix(SALMON_INDEX.out.versions)
        }
    }

    if (tools && tools.split(',').contains('msisensorpro')) {
        MSISENSORPRO_SCAN(input.fasta)

        msisensorpro = MSISENSORPRO_SCAN.out.list
        versions = versions.mix(MSISENSORPRO_SCAN.out.versions)
    }

    if (tools && tools.split(',').contains('rsem')) {
        RSEM_PREPAREREFERENCE_GENOME(input.fasta, input.gtf)

        rsem = RSEM_PREPAREREFERENCE_GENOME.out.index
        versions = versions.mix(RSEM_PREPAREREFERENCE_GENOME.out.versions)
    }

    if (tools && tools.split(',').contains('star')) {
        STAR_GENOMEGENERATE(input.fasta, input.gtf)

        star = STAR_GENOMEGENERATE.out.index
        versions = versions.mix(STAR_GENOMEGENERATE.out.versions)
    }

    emit:
    bowtie1
    bowtie2
    bwamem1
    bwamem2
    dict
    dragmap
    faidx
    gffread
    hisat2
    hisat2_splice_sites
    kallisto
    msisensorpro
    rsem
    rsem_transcript_fasta
    salmon
    sizes
    star
    versions
}
