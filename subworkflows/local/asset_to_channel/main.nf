workflow ASSET_TO_CHANNEL {
    take:
    asset // channel: [meta, fasta]

    main:
    // All the files and meta data are contained in the meta map (except for fasta)
    // They are extracted out of the meta map in their own channel in this subworkflow
    // When adding a new field in the assets/schema_input.json, also add it in the meta map
    // And in this scrip, add a new branch and a new output corresponding to this input
    // And in the emit, add the new output to the channel

    // Only keep the actual meta data in the meta map
    // Add a field here if it is a relevant meta data
    def reduce = { meta -> meta.subMap(['genome', 'id', 'source', 'source_vcf', 'source_version', 'species']) }

    intervals_bed_branch = asset.branch { meta, _fasta ->
        file: meta.intervals_bed
        return [reduce(meta), meta.intervals_bed]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }
    intervals_bed = intervals_bed_branch.file

    fasta_branch = asset.branch { meta, fasta_ ->
        file: fasta_
        // If ends with .gz, decompress it
        def meta_extra = [decompress_fasta: fasta_.endsWith('.gz') ?: false]
        // If any of the asset exists, then adding run_tools to false and skip the asset creation from the fasta file
        meta_extra += [run_bowtie1: meta.bowtie1_index ? false : true]
        meta_extra += [run_bowtie2: meta.bowtie2_index ? false : true]
        meta_extra += [run_bwamem1: meta.bwamem1_index ? false : true]
        meta_extra += [run_bwamem2: meta.bwamem2_index ? false : true]
        meta_extra += [run_dragmap: meta.dragmap_hashtable ? false : true]
        meta_extra += [run_faidx: meta.fasta_fai && meta.fasta_sizes ? false : true]
        meta_extra += [run_gatkdict: meta.fasta_dict ? false : true]
        meta_extra += [run_hisat2: meta.hisat2_index ? false : true]
        meta_extra += [run_intervals: meta.intervals_bed ? false : true]
        meta_extra += [run_kallisto: meta.kallisto_index ? false : true]
        meta_extra += [run_msisenpro: meta.msisensorpro_list ? false : true]
        meta_extra += [run_rsem: meta.rsem_index ? false : true]
        meta_extra += [run_rsem_make_transcript_fasta: meta.transcript_fasta ? false : true]
        meta_extra += [run_salmon: meta.salmon_index ? false : true]
        meta_extra += [run_star: meta.star_index ? false : true]
        return [reduce(meta) + meta_extra, fasta_]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }
    fasta = fasta_branch.file


    fasta_dict_branch = asset.branch { meta, _fasta ->
        file: meta.fasta_dict
        return [reduce(meta), meta.fasta_dict]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }
    fasta_dict = fasta_dict_branch.file


    // If we have intervals_bed, then we don't need to run faidx
    fasta_fai_branch = asset.branch { meta, _fasta ->
        file: meta.fasta_fai
        // If we have intervals_bed, then we don't need to run faidx
        def meta_extra = [run_intervals: meta.intervals_bed ? false : true]
        return [reduce(meta) + meta_extra, meta.fasta_fai]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }
    fasta_fai = fasta_fai_branch.file


    fasta_sizes_branch = asset.branch { meta, _fasta ->
        file: meta.fasta_sizes
        return [reduce(meta), meta.fasta_sizes]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }
    fasta_sizes = fasta_sizes_branch.file


    gff_branch = asset.branch { meta, fasta_ ->
        file: meta.gff
        // If ends with .gz, decompress it
        def meta_extra = [decompress_gff: meta.gff.endsWith('.gz') ?: false]
        // If any of the asset exists, then adding run_tools to false and skip the asset creation from the annotation derived file
        // (gff, gtf or transcript_fasta)
        meta_extra += [run_gffread: fasta_ && !meta.gtf ?: false]
        meta_extra += [run_hisat2: meta.splice_sites ? false : true]
        return [reduce(meta) + meta_extra, meta.gff]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }
    gff = gff_branch.file


    gtf_branch = asset.branch { meta, _fasta ->
        file: meta.gtf
        // If ends with .gz, decompress it
        def meta_extra = [decompress_gtf: meta.gtf.endsWith('.gz') ?: false]
        // If any of the asset exists, then adding run_tools to false and skip the asset creation from the annotation derived file
        // (gff, gtf or transcript_fasta)
        meta_extra += [run_hisat2: meta.splice_sites ? false : true]
        return [reduce(meta) + meta_extra, meta.gtf]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }
    gtf = gtf_branch.file


    splice_sites_branch = asset.branch { meta, _fasta ->
        file: meta.splice_sites
        return [reduce(meta), meta.splice_sites]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }
    splice_sites = splice_sites_branch.file

    transcript_fasta_branch = asset.branch { meta, _fasta ->
        file: meta.transcript_fasta
        // If any of the asset exists, then adding run_tools to false and skip the asset creation from the annotation derived file
        // (gff, gtf or transcript_fasta)
        def meta_extra = [run_hisat2: meta.hisat2_index ? false : true]
        meta_extra += [run_kallisto: meta.kallisto_index ? false : true]
        meta_extra += [run_rsem: meta.rsem_index ? false : true]
        meta_extra += [run_salmon: meta.salmon_index ? false : true]
        meta_extra += [run_star: meta.star_index ? false : true]
        return [reduce(meta) + meta_extra, meta.transcript_fasta]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }
    transcript_fasta = transcript_fasta_branch.file


    vcf_branch = asset.branch { meta, _fasta ->
        file: meta.vcf
        // If we already have the vcf_tbi, then we don't need to index the vcf
        def meta_extra = [run_tabix: meta.vcf_tbi || meta.vcf.endsWith('.vcf') ? false : true]
        meta_extra += [compress_vcf: meta.vcf.endsWith('.vcf') ?: false]
        // return a file, because we can catch globs this way, but it create issues with publishing
        return [reduce(meta) + meta_extra, file(meta.vcf)]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }
    // Using transpose here because we want to catch vcf with globs in the path because of nf-core/Sarek
    vcf = vcf_branch.file.transpose()

    emit:
    intervals_bed    // channel: [meta, *.bed]
    fasta            // channel: [meta, *.f(ast|n)?a]
    fasta_dict       // channel: [meta, *.f(ast|n)?a.dict]
    fasta_fai        // channel: [meta, *.f(ast|n)?a.fai]
    fasta_sizes      // channel: [meta, *.f(ast|n)?a.sizes]
    gff              // channel: [meta, gff]
    gtf              // channel: [meta, gtf]
    splice_sites     // channel: [meta, *.splice_sites.txt]
    transcript_fasta // channel: [meta, *.transcripts.fasta]
    vcf              // channel: [meta, *.vcf.gz]
}
