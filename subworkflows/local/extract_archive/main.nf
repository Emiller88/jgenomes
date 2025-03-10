//
// extract archive from any format
//

include { GUNZIP } from '../../../modules/nf-core/gunzip'
include { UNTAR  } from '../../../modules/nf-core/untar'
include { UNZIP  } from '../../../modules/nf-core/unzip'

workflow EXTRACT_ARCHIVE {
    take:
    archive

    main:
    versions = Channel.empty()

    archive_to_extract = archive.branch { _meta, archive_ ->
        tar: archive_.toString().endsWith('.tar.gz')
        gz: archive_.toString().endsWith('.gz')
        zip: archive_.toString().endsWith('.zip')
        non_assigned: true
    }

    // This is a confidence check
    not_extracted = archive_to_extract.non_assigned
    not_extracted.view { log.warn("Archive not in the expected format: " + it) }

    // extract archive
    GUNZIP(archive_to_extract.gz)
    UNTAR(archive_to_extract.tar)
    UNZIP(archive_to_extract.zip)

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
    extracted     // channel: [ meta, extracted_archive ]
    not_extracted // channel: [ meta, not_recognized_archive ]
    versions      // channel: [ versions.yml ]
}
