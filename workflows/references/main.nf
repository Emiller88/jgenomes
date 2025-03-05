include { PREPARE_GENOME_DNASEQ } from '../../subworkflows/local/prepare_genome_dnaseq'
include { PREPARE_GENOME_RNASEQ } from '../../subworkflows/local/prepare_genome_rnaseq'

workflow REFERENCES {
    take:
    ascat_alleles
    ascat_loci
    ascat_loci_gc
    ascat_loci_rt
    chr_dir
    fasta
    fasta_dict
    fasta_fai
    fasta_sizes
    gff
    gtf
    intervals_bed
    splice_sites
    transcript_fasta
    vcf
    tools            // List: Can contain any combination of tools of the list of available tools, or just no_tools

    main:

    versions = Channel.empty()

    // Create references for dnaseq based pipelines such as nf-core/sarek
    PREPARE_GENOME_DNASEQ(
        fasta,
        fasta_fai,
        vcf,
        tools.split(',').contains('bwamem1'),
        tools.split(',').contains('bwamem2'),
        tools.split(',').contains('createsequencedictionary'),
        tools.split(',').contains('dragmap'),
        tools.split(',').contains('faidx'),
        tools.split(',').contains('intervals'),
        tools.split(',').contains('msisensorpro'),
        tools.split(',').contains('tabix'),
    )

    fasta_fai = fasta_fai.mix(PREPARE_GENOME_DNASEQ.out.fasta_fai)

    // Create references for rnaseq based pipelines such as nf-core/riboseq, nf-core/rnaseq, nf-core/rnavar
    PREPARE_GENOME_RNASEQ(
        fasta,
        fasta_fai,
        gff,
        gtf,
        splice_sites,
        transcript_fasta,
        tools.split(',').contains('bowtie1'),
        tools.split(',').contains('bowtie2'),
        tools.split(',').contains('faidx'),
        tools.split(',').contains('hisat2'),
        tools.split(',').contains('hisat2_extractsplicesites'),
        tools.split(',').contains('kallisto'),
        tools.split(',').contains('rsem'),
        tools.split(',').contains('rsem_make_transcript_fasta'),
        tools.split(',').contains('salmon'),
        tools.split(',').contains('sizes'),
        tools.split(',').contains('star'),
    )

    // This works with a mixture of input and computed references
    fasta_dict = fasta_dict.mix(PREPARE_GENOME_DNASEQ.out.fasta_dict)
    fasta_sizes = fasta_sizes.mix(PREPARE_GENOME_RNASEQ.out.fasta_sizes)
    gtf = gtf.mix(PREPARE_GENOME_RNASEQ.out.gtf)
    intervals_bed = intervals_bed.mix(PREPARE_GENOME_DNASEQ.out.intervals_bed)
    splice_sites = splice_sites.mix(PREPARE_GENOME_RNASEQ.out.splice_sites)
    transcript_fasta = transcript_fasta.mix(PREPARE_GENOME_RNASEQ.out.transcript_fasta)

    // TODO: This does not work YET with a mixture of input and computed references
    bowtie1_index = PREPARE_GENOME_RNASEQ.out.bowtie1_index
    bowtie2_index = PREPARE_GENOME_RNASEQ.out.bowtie2_index
    bwamem1_index = PREPARE_GENOME_DNASEQ.out.bwamem1_index
    bwamem2_index = PREPARE_GENOME_DNASEQ.out.bwamem2_index
    dragmap_hashmap = PREPARE_GENOME_DNASEQ.out.dragmap_hashmap
    hisat2_index = PREPARE_GENOME_RNASEQ.out.hisat2_index
    kallisto_index = PREPARE_GENOME_RNASEQ.out.kallisto_index
    msisensorpro_list = PREPARE_GENOME_DNASEQ.out.msisensorpro_list
    rsem_index = PREPARE_GENOME_RNASEQ.out.rsem_index
    salmon_index = PREPARE_GENOME_RNASEQ.out.salmon_index
    star_index = PREPARE_GENOME_RNASEQ.out.star_index
    vcf_tbi = PREPARE_GENOME_DNASEQ.out.vcf_tbi

    // TODO: Refactor this with topics
    versions = versions.mix(PREPARE_GENOME_DNASEQ.out.versions)
    versions = versions.mix(PREPARE_GENOME_RNASEQ.out.versions)

    reference = Channel
        .empty()
        .mix(
            bowtie1_index.map { meta, reference_ -> [meta + [file: 'bowtie1_index'], reference_] },
            bowtie2_index.map { meta, reference_ -> [meta + [file: 'bowtie2_index'], reference_] },
            bwamem1_index.map { meta, reference_ -> [meta + [file: 'bwamem1_index'], reference_] },
            bwamem2_index.map { meta, reference_ -> [meta + [file: 'bwamem2_index'], reference_] },
            dragmap_hashmap.map { meta, reference_ -> [meta + [file: 'dragmap_hashmap'], reference_] },
            fasta.map { meta, reference_ -> [meta + [file: 'fasta'], reference_] },
            fasta_dict.map { meta, reference_ -> [meta + [file: 'fasta_dict'], reference_] },
            fasta_fai.map { meta, reference_ -> [meta + [file: 'fasta_fai'], reference_] },
            fasta_sizes.map { meta, reference_ -> [meta + [file: 'fasta_sizes'], reference_] },
            gff.map { meta, reference_ -> [meta + [file: 'gff'], reference_] },
            gtf.map { meta, reference_ -> [meta + [file: 'gtf'], reference_] },
            hisat2_index.map { meta, reference_ -> [meta + [file: 'hisat2_index'], reference_] },
            intervals_bed.map { meta, reference_ -> [meta + [file: 'intervals_bed'], reference_] },
            kallisto_index.map { meta, reference_ -> [meta + [file: 'kallisto_index'], reference_] },
            msisensorpro_list.map { meta, reference_ -> [meta + [file: 'msisensorpro_list'], reference_] },
            rsem_index.map { meta, reference_ -> [meta + [file: 'rsem_index'], reference_] },
            salmon_index.map { meta, reference_ -> [meta + [file: 'salmon_index'], reference_] },
            splice_sites.map { meta, reference_ -> [meta + [file: 'splice_sites'], reference_] },
            star_index.map { meta, reference_ -> [meta + [file: 'star_index'], reference_] },
            transcript_fasta.map { meta, reference_ -> [meta + [file: 'transcript_fasta'], reference_] },
            vcf.map { meta, reference_ -> [meta + [file: "${meta.type}_vcf"], reference_] },
            vcf_tbi.map { meta, reference_ -> [meta + [file: "${meta.type}_vcf_tbi"], reference_] },
        )

    emit:
    reference // channel: [meta, *]
    versions  // channel: [versions.yml]

    publish:
    reference >> 'reference'
}
