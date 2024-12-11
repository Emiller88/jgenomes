workflow ASSET_TO_CHANNEL {
    take:
    asset // channel: [meta, fasta]

    main:

    // All the files and meta data are contained in the meta map (except for fasta)
    // They are extracted out of the meta map in their own channel in this subworkflow
    // When adding a new field in the assets/schema_input.json, also add it in the meta map
    // And in this scrip, add a new map operation and a new output corresponding to this input

    intervals_bed = asset.map { meta, _fasta ->
        return meta.intervals_bed ? [meta, meta.intervals_bed] : null
    }

    fasta = asset.map { meta, fasta ->
        return fasta ? [meta + [decompress_fasta: fasta.endsWith('.gz') ?: false] + [run_bowtie1: meta.bowtie1_index ? false : true] + [run_bowtie2: meta.bowtie2_index ? false : true] + [run_bwamem1: meta.bwamem1_index ? false : true] + [run_bwamem2: meta.bwamem2_index ? false : true] + [run_dragmap: meta.dragmap_hashtable ? false : true] + [run_faidx: meta.fasta_fai && meta.fasta_sizes ? false : true] + [run_gatkdict: meta.fasta_dict ? false : true] + [run_hisat2: meta.hisat2_index ? false : true] + [run_intervals: meta.intervals_bed ? false : true] + [run_kallisto: meta.kallisto_index ? false : true] + [run_msisenpro: meta.msisensorpro_list ? false : true] + [run_rsem: meta.rsem_index ? false : true] + [run_rsem_make_transcript_fasta: meta.transcript_fasta ? false : true] + [run_salmon: meta.salmon_index ? false : true] + [run_star: meta.star_index ? false : true], fasta] : null
    }

    fasta_dict = asset.map { meta, _fasta ->
        return meta.fasta_dict ? [meta, meta.fasta_dict] : null
    }

    fasta_fai = asset.map { meta, _fasta ->
        return meta.fasta_fai ? [meta + [run_intervals: meta.intervals_bed ? false : true], meta.fasta_fai] : null
    }

    fasta_sizes = asset.map { meta, _fasta ->
        return meta.fasta_sizes ? [meta, meta.fasta_sizes] : null
    }

    gff = asset.map { meta, fasta_gff ->
        return meta.gff && !meta.gtf ? [meta + [decompress_gff: meta.gff.endsWith('.gz') ?: false] + [run_gffread: fasta_gff ?: false] + [run_hisat2: meta.splice_sites ? false : true], meta.gff] : null
    }

    gtf = asset.map { meta, _fasta ->
        return meta.gtf ? [meta + [decompress_gtf: meta.gtf.endsWith('.gz') ?: false] + [run_hisat2: meta.splice_sites ? false : true], meta.gtf] : null
    }

    splice_sites = asset.map { meta, _fasta ->
        return meta.splice_sites ? [meta, meta.splice_sites] : null
    }

    transcript_fasta = asset.map { meta, _fasta ->
        return meta.transcript_fasta ? [meta + [run_hisat2: meta.hisat2_index ? false : true] + [run_kallisto: meta.kallisto_index ? false : true] + [run_rsem: meta.rsem_index ? false : true] + [run_salmon: meta.salmon_index ? false : true] + [run_star: meta.star_index ? false : true], meta.transcript_fasta] : null
    }

    // Using transpose here because we want to catch vcf with globs in the path because of nf-core/Sarek
    vcf = asset
        .map { meta, _fasta ->
            return meta.vcf && !meta.vcf_tbi ? [meta, file(meta.vcf)] : null
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
