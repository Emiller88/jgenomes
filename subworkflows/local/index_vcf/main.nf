include { TABIX_TABIX } from '../../../modules/nf-core/tabix/tabix'

workflow INDEX_VCF {
    take:
    vcf       // channel: [meta, vcf]
    run_tabix // boolean: true/false

    main:
    vcf_tbi = Channel.empty()
    versions = Channel.empty()


    if (run_tabix) {
        TABIX_TABIX(vcf)

        vcf_tbi = TABIX_TABIX.out.tbi
        versions = TABIX_TABIX.out.versions
    }

    emit:
    vcf_tbi
    versions
}
