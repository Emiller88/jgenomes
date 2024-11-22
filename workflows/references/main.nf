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
    bowtie1             = Channel.empty()
    bowtie2             = Channel.empty()
    bwamem1             = Channel.empty()
    bwamem2             = Channel.empty()
    dict                = Channel.empty()
    dragmap             = Channel.empty()
    faidx               = Channel.empty()
    gffread             = Channel.empty()
    hisat2              = Channel.empty()
    hisat2_splice_sites = Channel.empty()
    kallisto            = Channel.empty()
    msisensorpro        = Channel.empty()
    rsem                = Channel.empty()
    salmon              = Channel.empty()
    sizes               = Channel.empty()
    star                = Channel.empty()
    transcript_fasta    = Channel.empty()
    versions            = Channel.empty()

    input = reference.multiMap{ meta, fasta, gff, gtf, splice_sites, transcript_fasta, readme, bed, mito, sizes ->
        fasta:            [meta, fasta]
        gff:              [meta, gff]
        gtf:              [meta, gtf]
        splice_sites:     [meta, splice_sites]
        transcript_fasta: [meta, transcript_fasta]
        readme:           [meta, readme]
    }

    if (tools && tools.split(',').contains('bowtie1')) {
        BOWTIE1_BUILD(input.fasta)

        bowtie1 = BOWTIE1_BUILD.out.index
        versions = versions.mix(BOWTIE1_BUILD.out.versions.first())
    }

    if (tools && tools.split(',').contains('bowtie2')) {
        BOWTIE2_BUILD(input.fasta)

        bowtie2 = BOWTIE2_BUILD.out.index
        versions = versions.mix(BOWTIE2_BUILD.out.versions.first())
    }

    if (tools && tools.split(',').contains('bwamem1')) {
        BWAMEM1_INDEX(input.fasta)

        bwamem1 = BWAMEM1_INDEX.out.index
        versions = versions.mix(BWAMEM1_INDEX.out.versions.first())
    }

    if (tools && tools.split(',').contains('bwamem2')) {
        BWAMEM2_INDEX(input.fasta)

        bwamem2 = BWAMEM2_INDEX.out.index
        versions = versions.mix(BWAMEM2_INDEX.out.versions.first())
    }

    if (tools && tools.split(',').contains('dragmap')) {
        DRAGMAP_HASHTABLE(input.fasta)

        dragmap = DRAGMAP_HASHTABLE.out.hashmap
        versions = versions.mix(DRAGMAP_HASHTABLE.out.versions.first())
    }

    if (tools && tools.split(',').contains('createsequencedictionary')) {
        GATK4_CREATESEQUENCEDICTIONARY(input.fasta)

        dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict
        versions = versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions.first())
    }

    if (tools && tools.split(',').contains('faidx')) {
        generate_sizes = tools.split(',').contains('sizes')
        SAMTOOLS_FAIDX(input.fasta, [ [ id:'no_fai' ], [] ], generate_sizes)

        faidx = SAMTOOLS_FAIDX.out.fai.first()
        sizes = SAMTOOLS_FAIDX.out.sizes.first()
        versions = versions.mix(SAMTOOLS_FAIDX.out.versions.first())
    }

    if (tools && tools.split(',').contains('gffread')) {
        GFFREAD(input.gff, [])

        gtf      = GFFREAD.out.gtf.map{ it[1] }
        versions = versions.mix(GFFREAD.out.versions)
    }

    if (tools && (tools.split(',').contains('hisat2') || tools.split(',').contains('hisat2_extractsplicesites')) ) {
        // TODO: be smarter about input assets
        //   Here we either return an empty channel if we have a splice_sites so that HISAT2_EXTRACTSPLICESITES is not triggered
        //   Or we return the provided gtf so that HISAT2_EXTRACTSPLICESITES is run
        gtf = input.gtf.join(input.splice_sites).groupTuple()
            .map{ meta, gtf, splice_sites -> return splice_sites[0][0] ? null : [meta, gtf] }

        HISAT2_EXTRACTSPLICESITES(gtf)
        versions = versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions.first())

        // TODO: be smarter about input assets
        //   Here we either mix+GT an empty channel (either no output or no input splice_sites) with the splice_sites return splice_sites
        //   And we filter out the empty value
        hisat2_extractsplicesites = input.splice_sites.mix(HISAT2_EXTRACTSPLICESITES.out.txt).groupTuple()
            .map{ meta, txt -> return txt[1] ? [meta, txt[1]] : [meta, txt] }

        if (tools.split(',').contains('hisat2')) {
            HISAT2_BUILD(input.fasta, input.gtf, hisat2_extractsplicesites)

            hisat2 = HISAT2_BUILD.out.index
            versions = versions.mix(HISAT2_BUILD.out.versions.first())
        }
    }

    if (tools && (tools.split(',').contains('kallisto')) || tools.split(',').contains('rsem_make_transcript_fasta') || tools.split(',').contains('salmon')) {
        // TODO: be smarter about input assets
        //   Here we either return an empty channel if we have a splice_sites so that MAKE_TRANSCRIPTS_FASTA is not triggered
        //   Or we return the provided gtf so that MAKE_TRANSCRIPTS_FASTA is run
        gtf = input.gtf.join(input.transcript_fasta).groupTuple()
            .map{ meta, gtf, transcript_fasta -> return transcript_fasta[0][0] ? null : [meta, gtf] }

        MAKE_TRANSCRIPTS_FASTA(input.fasta, gtf)
        versions = versions.mix(MAKE_TRANSCRIPTS_FASTA.out.versions.first())

        // TODO: be smarter about input assets
        //   Here we either mix+GT an empty channel (either no output or no input splice_sites) with the splice_sites return splice_sites
        //   And we filter out the empty value
        transcript_fasta = input.transcript_fasta.mix(MAKE_TRANSCRIPTS_FASTA.out.transcript_fasta).groupTuple()
            .map{ meta, txt -> return txt[1] ? [meta, txt[1]] : [meta, txt] }

        if (tools.split(',').contains('kallisto')) {
            KALLISTO_INDEX(transcript_fasta)

            kallisto = KALLISTO_INDEX.out.index
            versions = versions.mix(KALLISTO_INDEX.out.versions.first())
        }

        if (tools.split(',').contains('salmon')) {
            SALMON_INDEX(input.fasta, transcript_fasta)

            salmon = SALMON_INDEX.out.index
            versions = versions.mix(SALMON_INDEX.out.versions.first())
        }

    }

    if (tools && tools.split(',').contains('msisensorpro')) {
        MSISENSORPRO_SCAN(input.fasta)

        msisensorpro = MSISENSORPRO_SCAN.out.list
        versions = versions.mix(MSISENSORPRO_SCAN.out.versions.first())
    }

    if (tools && tools.split(',').contains('rsem')) {
        RSEM_PREPAREREFERENCE_GENOME(input.fasta, input.gtf)

        rsem = RSEM_PREPAREREFERENCE_GENOME.out.index
        versions = versions.mix(RSEM_PREPAREREFERENCE_GENOME.out.versions.first())
    }

    if (tools && tools.split(',').contains('star')) {
        STAR_GENOMEGENERATE(input.fasta, input.gtf)

        star = STAR_GENOMEGENERATE.out.index
        versions = versions.mix(STAR_GENOMEGENERATE.out.versions.first())
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
    salmon
    sizes
    star
    transcript_fasta
    versions
}
