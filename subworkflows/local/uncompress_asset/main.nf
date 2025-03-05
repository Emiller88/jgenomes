//
// Uncompress reference genome files
//

include { GUNZIP } from '../../../modules/nf-core/gunzip'
include { UNTAR  } from '../../../modules/nf-core/untar'
include { UNZIP  } from '../../../modules/nf-core/unzip'

workflow UNCOMPRESS_ASSET {
    take:
    extract_gz
    extract_tar
    extract_zip

    main:
    versions = Channel.empty()

    // Do not run GUNZIP, UNTAR, UNZIP if the condition is false
    // extract_gz = extract_gz.map { meta, extract_gz_ -> meta.decompress_gz ? [meta, extract_gz_] : null }
    // extract_tar = extract_tar.map { meta, extract_tar_ -> meta.decompress_tar ? [meta, extract_tar_] : null }
    // extract_zip = extract_zip.map { meta, extract_zip_ -> meta.decompress_zip ? [meta, extract_zip_] : null }

    // Uncompress the assets
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
