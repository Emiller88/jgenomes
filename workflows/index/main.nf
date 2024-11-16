include { BOWTIE_BUILD  } from "./../../modules/nf-core/bowtie/build/main"
include { BOWTIE2_BUILD } from "./../../modules/nf-core/bowtie2/build/main"
include { HOMER_UNIQMAP } from '../../subworkflows/local/homer_uniqmap/main.nf'

workflow INDEX {
    take:
    reference // fasta, gtf

    main:
    versions = Channel.empty()
    reference
        .multiMap { meta, fasta, gtf, bed, readme, mito, size ->
            fasta: tuple(meta, fasta)
            gtf: tuple(meta, gtf)
            bed: tuple(meta, bed)
        }
        .set { input }

    // FIXME Can't access container
    // BOWTIE_BUILD ( input.fasta )
    // version = versions.mix(BOWTIE_BUILD.out.versions)
    BOWTIE2_BUILD(input.fasta)
    versions = versions.mix(BOWTIE2_BUILD.out.versions)

    HOMER_UNIQMAP(input.fasta)

    emit:
    bowtie2_index = BOWTIE2_BUILD.out.index
    versions
}
