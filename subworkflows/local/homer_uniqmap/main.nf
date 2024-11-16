include { HOMER_GETMAPPABLEREGIONS } from '../../../modules/local/homer/getmappableregions/main'
include { HOMER_CREATEUNIQMAP      } from '../../../modules/local/homer/createuniqmap/main'

workflow HOMER_UNIQMAP {
    take:
    fasta //    file: genome.fa

    main:
    ch_versions = Channel.empty()

    // Split FASTA by chromosome
    split_fastas = fasta
        .splitFasta(by: 1, file: true)
        .toSortedList()

    // Generate mappable regions
    HOMER_GETMAPPABLEREGIONS(
        split_fastas,
        1000000000,
        50
    )
    ch_versions = ch_versions.mix(HOMER_GETMAPPABLEREGIONS.out.versions)

    // Create uniqmap directory
    HOMER_CREATEUNIQMAP(
        HOMER_GETMAPPABLEREGIONS.out.txt
    )
    ch_versions = ch_versions.mix(HOMER_CREATEUNIQMAP.out.versions)

    emit:
    uniqmap = HOMER_CREATEUNIQMAP.out.uniqmap_dir
}
