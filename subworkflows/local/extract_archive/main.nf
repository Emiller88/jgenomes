//
// extract reference genome files
//

include { GUNZIP } from '../../../modules/nf-core/gunzip'
include { UNTAR  } from '../../../modules/nf-core/untar'
include { UNZIP  } from '../../../modules/nf-core/unzip'

workflow EXTRACT_ARCHIVE {
    take:
    extract_gz
    extract_tar
    extract_zip

    main:
    versions = Channel.empty()

    // extract reference
    GUNZIP(extract_gz)
    UNTAR(extract_tar)
    UNZIP(extract_zip)

    extracted = Channel
        .empty()
        .mix(
            GUNZIP.out.gunzip,
            UNTAR.out.untar,
            UNZIP.out.unzipped_archive,
        )

    versions = versions.mix(GUNZIP.out.versions)
    versions = versions.mix(UNTAR.out.versions)
    versions = versions.mix(UNZIP.out.versions)

    emit:
    extracted // channel: [ meta, extracted_archive ]
    versions  // channel: [ versions.yml ]
}
