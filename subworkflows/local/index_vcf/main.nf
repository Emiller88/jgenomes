include { TABIX_TABIX } from '../../../modules/nf-core/tabix/tabix'

workflow INDEX_VCF {
    take:
    vcf       // channel: [meta, vcf]
    run_tabix // boolean: true/false

    main:
    vcf_tbi = Channel.empty()
    versions = Channel.empty()


    if (run_tabix) {
        vcf_tabix = vcf.map { meta, vcf_ ->
            return meta.run_tabix ? [meta, vcf_] : null
        }

        TABIX_TABIX(vcf_tabix)

        vcf_tbi = TABIX_TABIX.out.tbi
        versions = TABIX_TABIX.out.versions
    }

    emit:
    vcf_tbi  // channel: [meta, *.vcf.tbi]
    versions // channel: [versions.yml]
}
