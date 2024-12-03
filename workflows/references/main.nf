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

    intervals_bed = SAMPLESHEET_TO_CHANNEL.out.intervals_bed
    fasta = SAMPLESHEET_TO_CHANNEL.out.fasta
    fasta_dict = SAMPLESHEET_TO_CHANNEL.out.fasta_dict
    fasta_fai = SAMPLESHEET_TO_CHANNEL.out.fasta_fai
    fasta_sizes = SAMPLESHEET_TO_CHANNEL.out.fasta_sizes
    gff = SAMPLESHEET_TO_CHANNEL.out.gff
    gtf = SAMPLESHEET_TO_CHANNEL.out.gtf
    splice_sites = SAMPLESHEET_TO_CHANNEL.out.splice_sites
    transcript_fasta = SAMPLESHEET_TO_CHANNEL.out.transcript_fasta
    vcf = SAMPLESHEET_TO_CHANNEL.out.vcf

    CREATE_ALIGN_INDEX(
        fasta,
        tools.split(',').contains('bowtie1'),
        tools.split(',').contains('bowtie2'),
        tools.split(',').contains('bwamem1'),
        tools.split(',').contains('bwamem2'),
        tools.split(',').contains('dragmap')
    )

    CREATE_ALIGN_INDEX_WITH_GFF(
        fasta,
        gff,
        gtf,
        splice_sites,
        transcript_fasta,
        tools.split(',').contains('hisat2'),
        tools.split(',').contains('hisat2_extractsplicesites'),
        tools.split(',').contains('kallisto'),
        tools.split(',').contains('rsem'),
        tools.split(',').contains('rsem_make_transcript_fasta'),
        tools.split(',').contains('salmon'),
        tools.split(',').contains('star')
    )

    INDEX_FASTA(
        fasta,
        fasta_fai,
        tools.split(',').contains('createsequencedictionary'),
        tools.split(',').contains('faidx'),
        tools.split(',').contains('intervals'),
        tools.split(',').contains('msisensorpro'),
        tools.split(',').contains('sizes')
    )

    INDEX_VCF(
        vcf,
        tools.split(',').contains('tabix')
    )

    bowtie1_index = CREATE_ALIGN_INDEX.out.bowtie1_index
    bowtie2_index = CREATE_ALIGN_INDEX.out.bowtie2_index
    bwamem1_index = CREATE_ALIGN_INDEX.out.bwamem1_index
    bwamem2_index = CREATE_ALIGN_INDEX.out.bwamem2_index
    dragmap_hashmap = CREATE_ALIGN_INDEX.out.dragmap_hashmap

    gtf = gtf.mix(CREATE_ALIGN_INDEX_WITH_GFF.out.gtf)
    hisat2_index = CREATE_ALIGN_INDEX_WITH_GFF.out.hisat2_index
    splice_sites = splice_sites.mix(CREATE_ALIGN_INDEX_WITH_GFF.out.splice_sites)
    kallisto_index = CREATE_ALIGN_INDEX_WITH_GFF.out.kallisto_index
    rsem_index = CREATE_ALIGN_INDEX_WITH_GFF.out.rsem_index
    transcript_fasta = transcript_fasta.mix(CREATE_ALIGN_INDEX_WITH_GFF.out.transcript_fasta)
    salmon_index = CREATE_ALIGN_INDEX_WITH_GFF.out.salmon_index
    star_index = CREATE_ALIGN_INDEX_WITH_GFF.out.star_index

    fasta_dict = fasta_dict.mix(INDEX_FASTA.out.fasta_dict)
    fasta_fai = fasta_fai.mix(INDEX_FASTA.out.fasta_fai)
    intervals_bed = intervals_bed.mix(INDEX_FASTA.out.intervals_bed)
    fasta_sizes = fasta_sizes.mix(INDEX_FASTA.out.fasta_sizes)
    msisensorpro_list = INDEX_FASTA.out.msisensorpro_list

    vcf_tbi = INDEX_VCF.out.vcf_tbi

    versions = versions.mix(CREATE_ALIGN_INDEX.out.versions)
    versions = versions.mix(CREATE_ALIGN_INDEX_WITH_GFF.out.versions)
    versions = versions.mix(INDEX_FASTA.out.versions)
    versions = versions.mix(INDEX_VCF.out.versions)

    emit:
    bowtie1_index     // channel: [meta, BowtieIndex/]
    bowtie2_index     // channel: [meta, Bowtie2Index/]
    bwamem1_index     // channel: [meta, BWAmemIndex/]
    bwamem2_index     // channel: [meta, BWAmem2memIndex/]
    dragmap_hashmap   // channel: [meta, DragmapHashtable/]
    fasta             // channel: [meta, *.fa(sta)]
    fasta_dict        // channel: [meta, *.fa(sta).dict]
    fasta_fai         // channel: [meta, *.fa(sta).fai]
    fasta_sizes       // channel: [meta, *.fa(sta).sizes]
    gff               // channel: [meta, gtf]
    hisat2_index      // channel: [meta, Hisat2Index/]
    intervals_bed     // channel: [meta, *.bed]
    kallisto_index    // channel: [meta, KallistoIndex]
    msisensorpro_list // channel: [meta, *.list]
    rsem_index        // channel: [meta, RSEMIndex/]
    salmon_index      // channel: [meta, SalmonIndex/]
    splice_sites      // channel: [meta, *.splice_sites.txt]
    star_index        // channel: [meta, STARIndex/]
    transcript_fasta  // channel: [meta, *.transcripts.fasta]
    vcf_tbi           // channel: [meta, *.vcf.tbi]
    versions          // channel: [versions.yml]

    publish:
    bowtie1_index >> 'bowtie1_index'
    bowtie2_index >> 'bowtie2_index'
    bwamem1_index >> 'bwamem1_index'
    bwamem2_index >> 'bwamem2_index'
    dragmap_hashmap >> 'dragmap_hashmap'
    fasta >> 'fasta'
    fasta_dict >> 'fasta_dict'
    fasta_fai >> 'fasta_fai'
    fasta_sizes >> 'fasta_sizes'
    gff >> 'gff'
    hisat2_index >> 'hisat2_index'
    intervals_bed >> 'intervals_bed'
    kallisto_index >> 'kallisto_index'
    msisensorpro_list >> 'msisensorpro_list'
    rsem_index >> 'rsem_index'
    salmon_index >> 'salmon_index'
    splice_sites >> 'splice_sites'
    star_index >> 'star_index'
    transcript_fasta >> 'transcript_fasta'
    vcf_tbi >> 'vcf_tbi'
}
