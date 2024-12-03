include { CREATE_ALIGN_INDEX          } from '../../subworkflows/local/create_align_index'
include { CREATE_ALIGN_INDEX_WITH_GFF } from '../../subworkflows/local/create_align_index_with_gff'
include { INDEX_FASTA                 } from '../../subworkflows/local/index_fasta'
include { INDEX_VCF                   } from '../../subworkflows/local/index_vcf'
include { SAMPLESHEET_TO_CHANNEL      } from '../../subworkflows/local/samplesheet_to_channel'

workflow REFERENCES {
    take:
    reference // fasta, gff, gtf, splice_sites, transcript_fasta, vcf
    tools     // bowtie|bowtie2|bwamem1|bwamem2|createsequencedictionary|dragmap|faidx|gffread|intervals|hisat2|hisat2_extractsplicesites|kallisto|msisensorpro|rsem|rsem_make_transcripts_fasta|salmon|star|tabix

    main:
    versions = Channel.empty()

    SAMPLESHEET_TO_CHANNEL(reference)

    ch_intervals_bed = SAMPLESHEET_TO_CHANNEL.out.intervals_bed
    ch_fasta = SAMPLESHEET_TO_CHANNEL.out.fasta
    ch_fasta_dict = SAMPLESHEET_TO_CHANNEL.out.fasta_dict
    ch_fasta_fai = SAMPLESHEET_TO_CHANNEL.out.fasta_fai
    ch_fasta_sizes = SAMPLESHEET_TO_CHANNEL.out.fasta_sizes
    ch_gff = SAMPLESHEET_TO_CHANNEL.out.gff
    ch_gtf = SAMPLESHEET_TO_CHANNEL.out.gtf
    ch_splice_sites = SAMPLESHEET_TO_CHANNEL.out.splice_sites
    ch_transcript_fasta = SAMPLESHEET_TO_CHANNEL.out.transcript_fasta
    ch_vcf = SAMPLESHEET_TO_CHANNEL.out.vcf

    CREATE_ALIGN_INDEX(
        ch_fasta,
        tools.split(',').contains('bowtie1'),
        tools.split(',').contains('bowtie2'),
        tools.split(',').contains('bwamem1'),
        tools.split(',').contains('bwamem2'),
        tools.split(',').contains('dragmap')
    )

    CREATE_ALIGN_INDEX_WITH_GFF(
        ch_fasta,
        ch_gff,
        ch_gtf,
        ch_splice_sites,
        ch_transcript_fasta,
        tools.split(',').contains('hisat2'),
        tools.split(',').contains('hisat2_extractsplicesites'),
        tools.split(',').contains('kallisto'),
        tools.split(',').contains('rsem'),
        tools.split(',').contains('rsem_make_transcript_fasta'),
        tools.split(',').contains('salmon'),
        tools.split(',').contains('star')
    )

    INDEX_FASTA(
        ch_fasta,
        ch_fasta_fai,
        tools.split(',').contains('createsequencedictionary'),
        tools.split(',').contains('faidx'),
        tools.split(',').contains('intervals'),
        tools.split(',').contains('msisensorpro'),
        tools.split(',').contains('sizes')
    )

    INDEX_VCF(
        ch_vcf,
        tools.split(',').contains('tabix')
    )

    ch_bowtie1 = CREATE_ALIGN_INDEX.out.bowtie1_index
    ch_bowtie2 = CREATE_ALIGN_INDEX.out.bowtie2_index
    ch_bwamem1 = CREATE_ALIGN_INDEX.out.bwamem1_index
    ch_bwamem2 = CREATE_ALIGN_INDEX.out.bwamem2_index
    ch_dragmap = CREATE_ALIGN_INDEX.out.dragmap_hashmap

    ch_gff_gtf = CREATE_ALIGN_INDEX_WITH_GFF.out.gff_gtf
    ch_hisat2 = CREATE_ALIGN_INDEX_WITH_GFF.out.hisat2_index
    ch_splice_sites = ch_splice_sites.mix(CREATE_ALIGN_INDEX_WITH_GFF.out.hisat2_splice_sites)
    ch_kallisto = CREATE_ALIGN_INDEX_WITH_GFF.out.kallisto_index
    ch_rsem = CREATE_ALIGN_INDEX_WITH_GFF.out.rsem_index
    ch_transcript_fasta = ch_transcript_fasta.mix(CREATE_ALIGN_INDEX_WITH_GFF.out.rsem_transcript_fasta)
    ch_salmon = CREATE_ALIGN_INDEX_WITH_GFF.out.salmon_index
    ch_star = CREATE_ALIGN_INDEX_WITH_GFF.out.star_index
    ch_star = CREATE_ALIGN_INDEX_WITH_GFF.out.star_index

    ch_fasta_dict = ch_fasta_dict.mix(INDEX_FASTA.out.fasta_dict)
    ch_fasta_fai = ch_fasta_fai.mix(INDEX_FASTA.out.fasta_fai)
    ch_intervals_bed = ch_intervals_bed.mix(INDEX_FASTA.out.intervals_bed)
    ch_fasta_sizes = ch_fasta_sizes.mix(INDEX_FASTA.out.fasta_sizes)
    ch_msisensorpro = INDEX_FASTA.out.msisensorpro_list

    ch_vcf_tbi = INDEX_VCF.out.vcf_tbi

    versions = versions.mix(CREATE_ALIGN_INDEX.out.versions)
    versions = versions.mix(CREATE_ALIGN_INDEX_WITH_GFF.out.versions)
    versions = versions.mix(INDEX_FASTA.out.versions)
    versions = versions.mix(INDEX_VCF.out.versions)

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
    hisat2_splice_sites   = ch_splice_sites
    intervals_bed         = ch_intervals_bed
    kallisto              = ch_kallisto
    msisensorpro          = ch_msisensorpro
    rsem                  = ch_rsem
    rsem_transcript_fasta = ch_transcript_fasta
    salmon                = ch_salmon
    sizes                 = ch_fasta_sizes
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
    ch_splice_sites >> 'hisat2'
    ch_intervals_bed >> 'intervals'
    ch_kallisto >> 'kallisto'
    ch_msisensorpro >> 'msisensorpro'
    ch_rsem >> 'rsem'
    ch_transcript_fasta >> 'make'
    ch_salmon >> 'salmon'
    ch_fasta_sizes >> 'fasta_sizes'
    ch_star >> 'star'
    ch_vcf_tbi >> 'vcf_tbi'
}
