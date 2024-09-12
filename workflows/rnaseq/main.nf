include { GUNZIP as GUNZIP_GTF } from '../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GFF } from '../../modules/nf-core/gunzip'

include { GFFREAD } from '../../modules/nf-core/gffread'

include { STAR_GENOMEGENERATE } from "../../modules/nf-core/star/genomegenerate/main"
include { HISAT2_EXTRACTSPLICESITES } from "../../modules/nf-core/hisat2/extractsplicesites"
include { HISAT2_BUILD } from "../../modules/nf-core/hisat2/build"
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA } from "../../modules/nf-core/rsem/preparereference"
include { SALMON_INDEX } from "../../modules/nf-core/salmon/index"
include { KALLISTO_INDEX } from "../../modules/nf-core/kallisto/index"
include { RSEM_PREPAREREFERENCE as RSEM_PREPAREREFERENCE_GENOME } from "../../modules/nf-core/rsem/preparereference"

workflow RNASEQ {
    take:
    reference // fasta, gtf

    main:
    reference
        .multiMap { meta, fasta, gtf, gff, bed, readme, mito, size ->
            fasta: tuple(meta, fasta)
            gtf:   tuple(meta, gtf)
            gff:   tuple(meta, gff)
            bed:   tuple(meta, bed)
        }
        .set { input }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (input.gtf || input.gff) {
        if (input.gtf) {
            // TODO Maybe gunzip? Might not need to since there's no downstream analysis
            ch_gtf = input.gtf
        } else if (input.gff) {
            // TODO Add fasta cause why not
            ch_gtf = GFFREAD ( input.gff, [] ).input.gtf.map { it[1] }
        }
    }

    STAR_GENOMEGENERATE ( input.fasta, ch_gtf )

    ch_splicesites = HISAT2_EXTRACTSPLICESITES ( ch_gtf ).txt.map { it[1] }
    HISAT2_BUILD ( input.fasta, ch_gtf, ch_splicesites.map { [ [:], it ] } )

    ch_transcript_fasta = MAKE_TRANSCRIPTS_FASTA ( input.fasta, ch_gtf ).transcript_fasta

    SALMON_INDEX ( input.fasta, ch_transcript_fasta )

    KALLISTO_INDEX ( ch_transcript_fasta.map{[ [:], it]} )

    RSEM_PREPAREREFERENCE_GENOME ( input.fasta, ch_gtf )

    emit:
    star_index = STAR_GENOMEGENERATE.out.index
    hisat2_index = HISAT2_BUILD.out.index
    transcript_fasta = ch_transcript_fasta
    salmon_index = SALMON_INDEX.out.index
    kallisto_index = KALLISTO_INDEX.out.index
    rsem_index = RSEM_PREPAREREFERENCE_GENOME.out.index
}
