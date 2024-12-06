workflow SAMPLESHEET_TO_CHANNEL {
    take:
    reference // channel: [meta, intervals_bed, fasta, fasta_dict, fasta_fai, fasta_sizes, gff, gtf, splice_sites, transcript_fasta, vcf, readme, bed12, mito_name, macs_gsize]

    main:

    intervals_bed = reference.map { meta, input_intervals_bed, input_fasta, input_fasta_dict, input_fasta_fai, input_fasta_sizes, input_gff, input_gtf, input_splice_sites, input_transcript_fasta, input_vcf, input_readme, input_bed12, input_mito_name, input_macs_gsize ->
        return input_intervals_bed ? [meta, input_intervals_bed] : null
    }

    fasta = reference.map { meta, input_intervals_bed, input_fasta, input_fasta_dict, input_fasta_fai, input_fasta_sizes, input_gff, input_gtf, input_splice_sites, input_transcript_fasta, input_vcf, input_readme, input_bed12, input_mito_name, input_macs_gsize ->
        return input_fasta ? [meta + [decompress_fasta: input_fasta.endsWith('.gz') ?: false] + [run_faidx: input_fasta_fai && input_fasta_sizes ? false : true] + [run_intervals: input_intervals_bed ? false : true] + [run_rsem_make_transcript_fasta: input_transcript_fasta ? false : true], input_fasta] : null
    }

    fasta_dict = reference.map { meta, input_intervals_bed, input_fasta, input_fasta_dict, input_fasta_fai, input_fasta_sizes, input_gff, input_gtf, input_splice_sites, input_transcript_fasta, input_vcf, input_readme, input_bed12, input_mito_name, input_macs_gsize ->
        return input_fasta_dict ? [meta, input_fasta_dict] : null
    }

    fasta_fai = reference.map { meta, input_intervals_bed, input_fasta, input_fasta_dict, input_fasta_fai, input_fasta_sizes, input_gff, input_gtf, input_splice_sites, input_transcript_fasta, input_vcf, input_readme, input_bed12, input_mito_name, input_macs_gsize ->
        return input_fasta_fai ? [meta + [run_intervals: input_intervals_bed ? false : true], input_fasta_fai] : null
    }

    fasta_sizes = reference.map { meta, input_intervals_bed, input_fasta, input_fasta_dict, input_fasta_fai, input_fasta_sizes, input_gff, input_gtf, input_splice_sites, input_transcript_fasta, input_vcf, input_readme, input_bed12, input_mito_name, input_macs_gsize ->
        return input_fasta_sizes ? [meta, input_fasta_sizes] : null
    }

    gff = reference.map { meta, input_intervals_bed, input_fasta, input_fasta_dict, input_fasta_fai, input_fasta_sizes, input_gff, input_gtf, input_splice_sites, input_transcript_fasta, input_vcf, input_readme, input_bed12, input_mito_name, input_macs_gsize ->
        return input_gff && !input_gtf ? [meta + [decompress_gff: input_gff.endsWith('.gz') ?: false] + [run_hisat2: input_splice_sites ? false : true], input_gff] : null
    }

    gtf = reference.map { meta, input_intervals_bed, input_fasta, input_fasta_dict, input_fasta_fai, input_fasta_sizes, input_gff, input_gtf, input_splice_sites, input_transcript_fasta, input_vcf, input_readme, input_bed12, input_mito_name, input_macs_gsize ->
        return input_gtf ? [meta + [decompress_gtf: input_gtf.endsWith('.gz') ?: false] + [run_hisat2: input_splice_sites ? false : true], input_gtf] : null
    }

    splice_sites = reference.map { meta, input_intervals_bed, input_fasta, input_fasta_dict, input_fasta_fai, input_fasta_sizes, input_gff, input_gtf, input_splice_sites, input_transcript_fasta, input_vcf, input_readme, input_bed12, input_mito_name, input_macs_gsize ->
        return input_splice_sites ? [meta, input_splice_sites] : null
    }

    transcript_fasta = reference.map { meta, input_intervals_bed, input_fasta, input_fasta_dict, input_fasta_fai, input_fasta_sizes, input_gff, input_gtf, input_splice_sites, input_transcript_fasta, input_vcf, input_readme, input_bed12, input_mito_name, input_macs_gsize ->
        return input_transcript_fasta ? [meta, input_transcript_fasta] : null
    }

    vcf = reference
        .map { meta, input_intervals_bed, input_fasta, input_fasta_dict, input_fasta_fai, input_fasta_sizes, input_gff, input_gtf, input_splice_sites, input_transcript_fasta, input_vcf, input_readme, input_bed12, input_mito_name, input_macs_gsize ->
            return input_vcf ? [meta, file(input_vcf)] : null
        }
        .transpose()

    emit:
    intervals_bed
    fasta
    fasta_dict
    fasta_fai
    fasta_sizes
    gff
    gtf
    splice_sites
    transcript_fasta
    vcf
}
