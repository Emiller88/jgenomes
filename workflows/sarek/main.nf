include { BWA_INDEX as BWAMEM1_INDEX } from '../../modules/nf-core/bwa/index/main'
include { BWAMEM2_INDEX } from '../../modules/nf-core/bwamem2/index/main'
include { DRAGMAP_HASHTABLE } from '../../modules/nf-core/dragmap/hashtable/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/gatk4/createsequencedictionary/main'
include { MSISENSORPRO_SCAN } from '../../modules/nf-core/msisensorpro/scan/main'
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main'


workflow SAREK {
    take:
    reference // fasta, gtf

    main:
    versions = Channel.empty()
    reference
        .multiMap { meta, fasta, gtf, gff, bed, readme, mito, size ->
            fasta: tuple(meta, fasta)
            gtf:   tuple(meta, gtf)
            bed:   tuple(meta, bed)
        }
        .set { input }

    BWAMEM1_INDEX(input.fasta)
    BWAMEM2_INDEX(input.fasta)
    DRAGMAP_HASHTABLE(input.fasta)

    GATK4_CREATESEQUENCEDICTIONARY(input.fasta)
    MSISENSORPRO_SCAN(input.fasta)
    // FIXME
    // SAMTOOLS_FAIDX(input.fasta, [ [ id: input.meta.id ], [] ] )

    versions = versions.mix(BWAMEM1_INDEX.out.versions)
    versions = versions.mix(BWAMEM2_INDEX.out.versions)
    versions = versions.mix(DRAGMAP_HASHTABLE.out.versions)
    versions = versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    versions = versions.mix(MSISENSORPRO_SCAN.out.versions)

    emit:
    bwa = BWAMEM1_INDEX.out.index.map{ meta, index -> [index] }.collect()
    bwamem2 = BWAMEM2_INDEX.out.index.map{ meta, index -> [index] }.collect()
    versions
}
