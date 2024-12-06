//
// Uncompress reference genome files
//

include { GUNZIP as GUNZIP_FASTA } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF   } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GFF   } from '../../../modules/nf-core/gunzip'

workflow UNCOMPRESS_REFERENCES {
    take:
    fasta_input // channel: [meta, fasta]
    gff_input   // channel: [meta, gff]
    gtf_input   // channel: [meta, gtf]

    main:
    versions = Channel.empty()

    fasta_input = fasta_input.map { meta, fasta_map ->
        meta.decompress_fasta ? [meta, fasta_map] : null
    }
    gff_input = gff_input.map { meta, gff_map ->
        meta.decompress_gff ? [meta, gff_map] : null
    }
    gtf_input = gtf_input.map { meta, gtf_map ->
        meta.decompress_gtf ? [meta, gtf_map] : null
    }

    GUNZIP_FASTA(fasta_input)
    GUNZIP_GFF(gff_input)
    GUNZIP_GTF(gtf_input)

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
