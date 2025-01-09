//
// Uncompress reference genome files
//

include { GUNZIP as GUNZIP_FASTA } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GFF   } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF   } from '../../../modules/nf-core/gunzip'
include { UNTAR as UNTAR_CHR_DIR } from '../../../modules/nf-core/untar'
include { UNZIP as UNZIP_ALLELES } from '../../../modules/nf-core/unzip'
include { UNZIP as UNZIP_GC      } from '../../../modules/nf-core/unzip'
include { UNZIP as UNZIP_LOCI    } from '../../../modules/nf-core/unzip'
include { UNZIP as UNZIP_RT      } from '../../../modules/nf-core/unzip'

workflow UNCOMPRESS_ASSET {
    take:
    ascat_alleles // channel: [meta, ascat_alleles]
    ascat_loci    // channel: [meta, ascat_loci]
    ascat_loci_gc // channel: [meta, ascat_loci_gc]
    ascat_loci_rt // channel: [meta, ascat_loci_rt]
    chr_dir       // channel: [meta, chr_dir]
    fasta         // channel: [meta, fasta]
    gff           // channel: [meta, gff]
    gtf           // channel: [meta, gtf]

    main:
    versions = Channel.empty()

    // Do not run GUNZIP_*, UNTAR_*, UNZIP_* if the condition is false
    ascat_alleles = ascat_alleles.map { meta, ascat_alleles_ -> meta.decompress_ascat_alleles ? [meta, ascat_alleles_] : null }
    ascat_loci = ascat_loci.map { meta, ascat_loci_ -> meta.decompress_ascat_loci ? [meta, ascat_loci_] : null }
    ascat_loci_gc = ascat_loci_gc.map { meta, ascat_loci_gc_ -> meta.decompress_ascat_loci_gc ? [meta, ascat_loci_gc_] : null }
    ascat_loci_rt = ascat_loci_rt.map { meta, ascat_loci_rt_ -> meta.decompress_ascat_loci_rt ? [meta, ascat_loci_rt_] : null }
    chr_dir = chr_dir.map { meta, chr_dir_ -> meta.decompress_chr_dir ? [meta, chr_dir_] : null }
    fasta = fasta.map { meta, fasta_ -> meta.decompress_fasta ? [meta, fasta_] : null }
    gff = gff.map { meta, gff_ -> meta.decompress_gff ? [meta, gff_] : null }
    gtf = gtf.map { meta, gtf_ -> meta.decompress_gtf ? [meta, gtf_] : null }

    // Uncompress the assets
    GUNZIP_FASTA(fasta)
    GUNZIP_GFF(gff)
    GUNZIP_GTF(gtf)
    UNTAR_CHR_DIR(chr_dir)
    UNZIP_ALLELES(ascat_alleles)
    UNZIP_GC(ascat_loci_gc)
    UNZIP_LOCI(ascat_loci)
    UNZIP_RT(ascat_loci_rt)

    ascat_alleles = UNZIP_ALLELES.out.unzipped_archive
    ascat_loci = UNZIP_LOCI.out.unzipped_archive
    ascat_loci_gc = UNZIP_GC.out.unzipped_archive
    ascat_loci_rt = UNZIP_RT.out.unzipped_archive
    chr_dir = UNTAR_CHR_DIR.out.untar
    fasta = GUNZIP_FASTA.out.gunzip
    gff = GUNZIP_GFF.out.gunzip
    gtf = GUNZIP_GTF.out.gunzip

    versions = versions.mix(GUNZIP_FASTA.out.versions)
    versions = versions.mix(GUNZIP_GFF.out.versions)
    versions = versions.mix(GUNZIP_GTF.out.versions)
    versions = versions.mix(UNTAR_CHR_DIR.out.versions)
    versions = versions.mix(UNZIP_ALLELES.out.versions)
    versions = versions.mix(UNZIP_GC.out.versions)
    versions = versions.mix(UNZIP_LOCI.out.versions)
    versions = versions.mix(UNZIP_RT.out.versions)

    emit:
    ascat_alleles // channel: [ meta, ascat_alleles ]
    ascat_loci    // channel: [ meta, ascat_loci ]
    ascat_loci_gc // channel: [ meta, ascat_loci_gc ]
    ascat_loci_rt // channel: [ meta, ascat_loci_rt ]
    chr_dir       // channel: [ meta, chr_dir ]
    fasta         // channel: [ meta, genome.fasta]
    gff           // channel: [ meta, genome.gff]
    gtf           // channel: [ meta, genome.gtf]
    versions      // channel: [ versions.yml ]
}
