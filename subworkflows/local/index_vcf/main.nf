include { TABIX_BGZIPTABIX } from '../../../modules/nf-core/tabix/bgziptabix'
include { TABIX_TABIX      } from '../../../modules/nf-core/tabix/tabix'

workflow INDEX_VCF {
    take:
    vcf       // channel: [meta, vcf]
    run_tabix // boolean: true/false

    main:
    vcf_gz = Channel.empty()
    vcf_tbi = Channel.empty()
    versions = Channel.empty()

    if (run_tabix) {
        // Do not run TABIX_TABIX if the condition is false
        vcf_tabix = vcf.map { meta, vcf_ -> meta.run_tabix ? [meta, vcf_] : null }

        // Do not run TABIX_BGZIPTABIX if the condition is false
        vcf_bgziptabix = vcf.map { meta, vcf_ -> meta.compress_vcf ? [meta, vcf_] : null }

        TABIX_BGZIPTABIX(vcf_bgziptabix)
        TABIX_TABIX(vcf_tabix)

        vcf_gz = TABIX_BGZIPTABIX.out.gz_tbi.map { meta, vcf_gz_, _vcf_tbi -> [meta, vcf_gz_] }
        vcf_tbi = TABIX_TABIX.out.tbi.mix(TABIX_BGZIPTABIX.out.gz_tbi.map { meta, _vcf_gz, vcf_tbi_ -> [meta, vcf_tbi_] })

        versions = versions.mix(TABIX_BGZIPTABIX.out.versions)
        versions = versions.mix(TABIX_TABIX.out.versions)
    }

    emit:
    vcf_gz   // channel: [meta, *.vcf.gz]
    vcf_tbi  // channel: [meta, *.vcf.gz.tbi]
    versions // channel: [versions.yml]
}
