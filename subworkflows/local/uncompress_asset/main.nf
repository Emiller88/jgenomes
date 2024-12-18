//
// Uncompress reference genome files
//

include { GUNZIP as GUNZIP_FASTA } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF   } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GFF   } from '../../../modules/nf-core/gunzip'

workflow UNCOMPRESS_ASSET {
    take:
    fasta // channel: [meta, fasta]
    gff   // channel: [meta, gff]
    gtf   // channel: [meta, gtf]

    main:
    versions = Channel.empty()

    fasta = fasta.map { meta, fasta_ -> meta.decompress_fasta ? [meta, fasta_] : null }
    gff = gff.map { meta, gff_ -> meta.decompress_gff ? [meta, gff_] : null }
    gtf = gtf.map { meta, gtf_ -> meta.decompress_gtf ? [meta, gtf_] : null }

    GUNZIP_FASTA(fasta)
    GUNZIP_GFF(gff)
    GUNZIP_GTF(gtf)

    fasta = GUNZIP_FASTA.out.gunzip
    gff = GUNZIP_GFF.out.gunzip
    gtf = GUNZIP_GTF.out.gunzip

    versions = versions.mix(GUNZIP_FASTA.out.versions)
    versions = versions.mix(GUNZIP_GFF.out.versions)
    versions = versions.mix(GUNZIP_GTF.out.versions)

    emit:
    fasta    // channel: [ meta, genome.fasta]
    gff      // channel: [ meta, genome.gff]
    gtf      // channel: [ meta, genome.gtf]
    versions // channel: [ versions.yml ]
}
