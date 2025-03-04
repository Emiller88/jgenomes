workflow ASSET_TO_CHANNEL {
    take:
    asset // channel: [meta, fasta]
    tools // List: Can contain any combination of tools of the list of available tools, or just no_tools

    main:
    // All the files and meta data are contained in the meta map (except for fasta)
    // They are extracted out of the meta map in their own channel in this subworkflow
    // When adding a new field in the assets/schema_input.json, also add it in the meta map
    // And in this scrip, add a new branch and a new output corresponding to this input
    // And in the emit, add the new output to the channel

    // Only keep the actual meta data in the meta map
    // Add a field here if it is a relevant meta data
    def reduce = { meta -> meta.subMap(['genome', 'id', 'source', 'source_version', 'species']) }

    ascat_alleles_branch = asset.branch { meta, _readme ->
        file: meta.ascat_alleles
        return [reduce(meta), meta.ascat_alleles]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }
    ascat_alleles = ascat_alleles_branch.file

    ascat_loci_branch = asset.branch { meta, _readme ->
        file: meta.ascat_loci
        return [reduce(meta), meta.ascat_loci]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }
    ascat_loci = ascat_loci_branch.file

    ascat_loci_gc_branch = asset.branch { meta, _readme ->
        file: meta.ascat_loci_gc
        return [reduce(meta), meta.ascat_loci_gc]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }
    ascat_loci_gc = ascat_loci_gc_branch.file

    ascat_loci_rt_branch = asset.branch { meta, _readme ->
        file: meta.ascat_loci_rt
        return [reduce(meta), meta.ascat_loci_rt]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }
    ascat_loci_rt = ascat_loci_rt_branch.file

    chr_dir_branch = asset.branch { meta, _readme ->
        file: meta.chr_dir
        return [reduce(meta), meta.chr_dir]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }
    chr_dir = chr_dir_branch.file

    intervals_bed_branch = asset.branch { meta, _readme ->
        file: meta.intervals_bed
        return [reduce(meta), meta.intervals_bed]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }
    intervals_bed = intervals_bed_branch.file

    fasta_branch = asset.branch { meta, _readme ->
        file: meta.fasta
        // If any of the asset exists, then adding run_tools to false and skip the asset creation from the fasta file
        def meta_extra = [run_bowtie1: meta.bowtie1_index ? false : true]
        meta_extra += [run_bowtie2: meta.bowtie2_index ? false : true]
        meta_extra += [run_bwamem1: meta.bwamem1_index ? false : true]
        meta_extra += [run_bwamem2: meta.bwamem2_index ? false : true]
        meta_extra += [run_createsequencedictionary: meta.fasta_dict ? false : true]
        meta_extra += [run_dragmap: meta.dragmap_hashtable ? false : true]
        meta_extra += [run_faidx: meta.fasta_fai && (meta.fasta_sizes || !tools.contains('sizes')) ? false : true]
        meta_extra += [run_hisat2: meta.hisat2_index ? false : true]
        meta_extra += [run_intervals: meta.intervals_bed ? false : true]
        meta_extra += [run_kallisto: meta.kallisto_index ? false : true]
        meta_extra += [run_msisenpro: meta.msisensorpro_list ? false : true]
        meta_extra += [run_rsem: meta.rsem_index ? false : true]
        meta_extra += [run_rsem_make_transcript_fasta: meta.transcript_fasta ? false : true]
        meta_extra += [run_salmon: meta.salmon_index ? false : true]
        meta_extra += [run_star: meta.star_index ? false : true]
        return [reduce(meta) + meta_extra, meta.fasta]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }
    fasta = fasta_branch.file

    fasta_dict_branch = asset.branch { meta, _readme ->
        file: meta.fasta_dict
        return [reduce(meta), meta.fasta_dict]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }
    fasta_dict = fasta_dict_branch.file

    // If we have intervals_bed, then we don't need to run faidx
    fasta_fai_branch = asset.branch { meta, _readme ->
        file: meta.fasta_fai
        // If we have intervals_bed, then we don't need to run faidx
        def meta_extra = [run_intervals: meta.intervals_bed ? false : true]
        return [reduce(meta) + meta_extra, meta.fasta_fai]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }
    fasta_fai = fasta_fai_branch.file

    fasta_sizes_branch = asset.branch { meta, _readme ->
        file: meta.fasta_sizes
        return [reduce(meta), meta.fasta_sizes]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }
    fasta_sizes = fasta_sizes_branch.file

    gff_branch = asset.branch { meta, _readme ->
        file: meta.gff
        // If any of the asset exists, then adding run_tools to false and skip the asset creation from the annotation derived file
        // (gff, gtf or transcript_fasta)
        def meta_extra = [run_gffread: meta.fasta && !meta.gtf ?: false]
        meta_extra += [run_hisat2: meta.splice_sites ? false : true]
        return [reduce(meta) + meta_extra, meta.gff]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }
    gff = gff_branch.file

    gtf_branch = asset.branch { meta, _readme ->
        file: meta.gtf
        // If any of the asset exists, then adding run_tools to false and skip the asset creation from the annotation derived file
        // (gff, gtf or transcript_fasta)
        def meta_extra = [run_hisat2: meta.splice_sites ? false : true]
        return [reduce(meta) + meta_extra, meta.gtf]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }
    gtf = gtf_branch.file

    splice_sites_branch = asset.branch { meta, _readme ->
        file: meta.splice_sites
        return [reduce(meta), meta.splice_sites]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }
    splice_sites = splice_sites_branch.file

    transcript_fasta_branch = asset.branch { meta, _readme ->
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

    // HANDLING OF VCF

    dbsnp_branch = asset.branch { meta, _readme ->
        file: meta.vcf_dbsnp_vcf

        // If we already have the vcf_tbi, then we don't need to index the vcf
        def meta_extra = [run_tabix: meta.vcf_dbsnp_vcf_tbi || meta.vcf_dbsnp_vcf.endsWith('.vcf') ? false : true]
        meta_extra += [compress_vcf: meta.vcf_dbsnp_vcf.endsWith('.vcf') ?: false]
        meta_extra += [type: 'dbsnp', source_vcf: meta.vcf_dbsnp_vcf_source]
        return [reduce(meta) + meta_extra, meta.vcf_dbsnp_vcf.contains('{') ? file(meta.vcf_dbsnp_vcf) : meta.vcf_dbsnp_vcf]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }

    germline_resource_branch = asset.branch { meta, _readme ->
        file: meta.vcf_germline_resource_vcf
        // If we already have the vcf_tbi, then we don't need to index the vcf
        def meta_extra = [run_tabix: meta.vcf_germline_resource_vcf_tbi || meta.vcf_germline_resource_vcf.endsWith('.vcf') ? false : true]
        meta_extra += [compress_vcf: meta.vcf_germline_resource_vcf.endsWith('.vcf') ?: false]
        meta_extra += [type: 'germline_resource', source_vcf: meta.vcf_germline_resource_vcf_source]
        return [reduce(meta) + meta_extra, meta.vcf_germline_resource_vcf]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }

    known_indels_branch = asset.branch { meta, _readme ->
        file: meta.vcf_known_indels_vcf
        // If we already have the vcf_tbi, then we don't need to index the vcf
        def meta_extra = [run_tabix: meta.vcf_known_indels_vcf_tbi || meta.vcf_known_indels_vcf.endsWith('.vcf') ? false : true]
        meta_extra += [compress_vcf: meta.vcf_known_indels_vcf.endsWith('.vcf') ?: false]
        meta_extra += [type: 'known_indels', source_vcf: meta.vcf_known_indels_vcf_source]
        return [reduce(meta) + meta_extra, meta.vcf_known_indels_vcf]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }

    known_snps_branch = asset.branch { meta, _readme ->
        file: meta.vcf_known_snps_vcf
        // If we already have the vcf_tbi, then we don't need to index the vcf
        def meta_extra = [run_tabix: meta.vcf_known_snps_vcf_tbi || meta.vcf_known_snps_vcf.endsWith('.vcf') ? false : true]
        meta_extra += [compress_vcf: meta.vcf_known_snps_vcf.endsWith('.vcf') ?: false]
        meta_extra += [type: 'known_snps', source_vcf: meta.vcf_known_snps_vcf_source]
        return [reduce(meta) + meta_extra, meta.vcf_known_snps_vcf]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }

    pon_branch = asset.branch { meta, _readme ->
        file: meta.vcf_pon_vcf
        // If we already have the vcf_tbi, then we don't need to index the vcf
        def meta_extra = [run_tabix: meta.vcf_pon_vcf_tbi || meta.vcf_pon_vcf.endsWith('.vcf') ? false : true]
        meta_extra += [compress_vcf: meta.vcf_pon_vcf.endsWith('.vcf') ?: false]
        meta_extra += [type: 'pon', source_vcf: meta.vcf_pon_vcf_source]
        return [reduce(meta) + meta_extra, meta.vcf_pon_vcf]
        other: true
        // If the asset doesn't exist, then we return nothing
        return null
    }

    vcf = Channel.empty().mix(dbsnp_branch.file, germline_resource_branch.file, known_indels_branch.file, known_snps_branch.file, pon_branch.file).transpose()

    emit:
    ascat_alleles    // channel: [meta, *.ascat_alleles.txt]
    ascat_loci       // channel: [meta, *.ascat_loci.txt]
    ascat_loci_gc    // channel: [meta, *.ascat_loci_gc.txt]
    ascat_loci_rt    // channel: [meta, *.ascat_loci_rt.txt]
    chr_dir          // channel: [meta, *.chr_dir]
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
