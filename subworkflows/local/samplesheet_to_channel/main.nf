workflow SAMPLESHEET_TO_CHANNEL {
    take:
    reference // channel: [meta, bed12, bowtie1_index, bowtie2_index, bwamem1_index, bwamem2_index, dragmap_hashtable, fasta, fasta_dict, fasta_fai, fasta_sizes, gff, gtf, hisat2_index, intervals_bed, kallisto_index, macs_gsize, mito_name, msisensorpro_list, readme, rsem_index, salmon_index, splice_sites, star_index, transcript_fasta, vcf]

    main:

    intervals_bed = reference.map { meta, input_bed12, input_bowtie1_index, input_bowtie2_index, input_bwamem1_index, input_bwamem2_index, input_dragmap_hashtable, input_fasta, input_fasta_dict, input_fasta_fai, input_fasta_sizes, input_gff, input_gtf, input_hisat2_index, input_intervals_bed, input_kallisto_index, input_macs_gsize, input_mito_name, input_msisensorpro_list, input_readme, input_rsem_index, input_salmon_index, input_splice_sites, input_star_index, input_transcript_fasta, input_vcf ->
        return input_intervals_bed ? [meta, input_intervals_bed] : null
    }

    fasta = reference.map { meta, input_bed12, input_bowtie1_index, input_bowtie2_index, input_bwamem1_index, input_bwamem2_index, input_dragmap_hashtable, input_fasta, input_fasta_dict, input_fasta_fai, input_fasta_sizes, input_gff, input_gtf, input_hisat2_index, input_intervals_bed, input_kallisto_index, input_macs_gsize, input_mito_name, input_msisensorpro_list, input_readme, input_rsem_index, input_salmon_index, input_splice_sites, input_star_index, input_transcript_fasta, input_vcf ->
        return input_fasta ? [meta + [decompress_fasta: input_fasta.endsWith('.gz') ?: false] + [run_bowtie1: input_bowtie1_index ? false : true] + [run_bowtie2: input_bowtie2_index ? false : true] + [run_bwamem1: input_bwamem1_index ? false : true] + [run_bwamem2: input_bwamem2_index ? false : true] + [run_dragmap: input_dragmap_hashtable ? false : true] + [run_faidx: input_fasta_fai && input_fasta_sizes ? false : true] + [run_gatkdict: input_fasta_dict ? false : true] + [run_hisat2: input_hisat2_index ? false : true] + [run_intervals: input_intervals_bed ? false : true] + [run_kallisto: input_kallisto_index ? false : true] + [run_msisenpro: input_msisensorpro_list ? false : true] + [run_rsem: input_rsem_index ? false : true] + [run_rsem_make_transcript_fasta: input_transcript_fasta ? false : true] + [run_salmon: input_salmon_index ? false : true] + [run_star: input_star_index ? false : true], input_fasta] : null
    }

    fasta_dict = reference.map { meta, input_bed12, input_bowtie1_index, input_bowtie2_index, input_bwamem1_index, input_bwamem2_index, input_dragmap_hashtable, input_fasta, input_fasta_dict, input_fasta_fai, input_fasta_sizes, input_gff, input_gtf, input_hisat2_index, input_intervals_bed, input_kallisto_index, input_macs_gsize, input_mito_name, input_msisensorpro_list, input_readme, input_rsem_index, input_salmon_index, input_splice_sites, input_star_index, input_transcript_fasta, input_vcf ->
        return input_fasta_dict ? [meta, input_fasta_dict] : null
    }

    fasta_fai = reference.map { meta, input_bed12, input_bowtie1_index, input_bowtie2_index, input_bwamem1_index, input_bwamem2_index, input_dragmap_hashtable, input_fasta, input_fasta_dict, input_fasta_fai, input_fasta_sizes, input_gff, input_gtf, input_hisat2_index, input_intervals_bed, input_kallisto_index, input_macs_gsize, input_mito_name, input_msisensorpro_list, input_readme, input_rsem_index, input_salmon_index, input_splice_sites, input_star_index, input_transcript_fasta, input_vcf ->
        return input_fasta_fai ? [meta + [run_intervals: input_intervals_bed ? false : true], input_fasta_fai] : null
    }

    fasta_sizes = reference.map { meta, input_bed12, input_bowtie1_index, input_bowtie2_index, input_bwamem1_index, input_bwamem2_index, input_dragmap_hashtable, input_fasta, input_fasta_dict, input_fasta_fai, input_fasta_sizes, input_gff, input_gtf, input_hisat2_index, input_intervals_bed, input_kallisto_index, input_macs_gsize, input_mito_name, input_msisensorpro_list, input_readme, input_rsem_index, input_salmon_index, input_splice_sites, input_star_index, input_transcript_fasta, input_vcf ->
        return input_fasta_sizes ? [meta, input_fasta_sizes] : null
    }

    gff = reference.map { meta, input_bed12, input_bowtie1_index, input_bowtie2_index, input_bwamem1_index, input_bwamem2_index, input_dragmap_hashtable, input_fasta, input_fasta_dict, input_fasta_fai, input_fasta_sizes, input_gff, input_gtf, input_hisat2_index, input_intervals_bed, input_kallisto_index, input_macs_gsize, input_mito_name, input_msisensorpro_list, input_readme, input_rsem_index, input_salmon_index, input_splice_sites, input_star_index, input_transcript_fasta, input_vcf ->
        return input_gff && !input_gtf ? [meta + [decompress_gff: input_gff.endsWith('.gz') ?: false] + [run_gffread: input_fasta ?: false] + [run_hisat2: input_splice_sites ? false : true], input_gff] : null
    }

    gtf = reference.map { meta, input_bed12, input_bowtie1_index, input_bowtie2_index, input_bwamem1_index, input_bwamem2_index, input_dragmap_hashtable, input_fasta, input_fasta_dict, input_fasta_fai, input_fasta_sizes, input_gff, input_gtf, input_hisat2_index, input_intervals_bed, input_kallisto_index, input_macs_gsize, input_mito_name, input_msisensorpro_list, input_readme, input_rsem_index, input_salmon_index, input_splice_sites, input_star_index, input_transcript_fasta, input_vcf ->
        return input_gtf ? [meta + [decompress_gtf: input_gtf.endsWith('.gz') ?: false] + [run_hisat2: input_splice_sites ? false : true], input_gtf] : null
    }

    splice_sites = reference.map { meta, input_bed12, input_bowtie1_index, input_bowtie2_index, input_bwamem1_index, input_bwamem2_index, input_dragmap_hashtable, input_fasta, input_fasta_dict, input_fasta_fai, input_fasta_sizes, input_gff, input_gtf, input_hisat2_index, input_intervals_bed, input_kallisto_index, input_macs_gsize, input_mito_name, input_msisensorpro_list, input_readme, input_rsem_index, input_salmon_index, input_splice_sites, input_star_index, input_transcript_fasta, input_vcf ->
        return input_splice_sites ? [meta, input_splice_sites] : null
    }

    transcript_fasta = reference.map { meta, input_bed12, input_bowtie1_index, input_bowtie2_index, input_bwamem1_index, input_bwamem2_index, input_dragmap_hashtable, input_fasta, input_fasta_dict, input_fasta_fai, input_fasta_sizes, input_gff, input_gtf, input_hisat2_index, input_intervals_bed, input_kallisto_index, input_macs_gsize, input_mito_name, input_msisensorpro_list, input_readme, input_rsem_index, input_salmon_index, input_splice_sites, input_star_index, input_transcript_fasta, input_vcf ->
        return input_transcript_fasta ? [meta + [decompress_fasta: input_fasta.endsWith('.gz') ?: false] + [run_hisat2: input_hisat2_index ? false : true] + [run_kallisto: input_kallisto_index ? false : true] + [run_rsem: input_rsem_index ? false : true] + [run_salmon: input_salmon_index ? false : true] + [run_star: input_star_index ? false : true], input_transcript_fasta] : null
    }

    vcf = reference
        .map { meta, input_bed12, input_bowtie1_index, input_bowtie2_index, input_bwamem1_index, input_bwamem2_index, input_dragmap_hashtable, input_fasta, input_fasta_dict, input_fasta_fai, input_fasta_sizes, input_gff, input_gtf, input_hisat2_index, input_intervals_bed, input_kallisto_index, input_macs_gsize, input_mito_name, input_msisensorpro_list, input_readme, input_rsem_index, input_salmon_index, input_splice_sites, input_star_index, input_transcript_fasta, input_vcf ->
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
