workflow ASSET_TO_CHANNEL {
    take:
    asset // channel: [meta, fasta]

    main:

    // All the files and meta data are contained in the meta map (except for fasta)
    // They are extracted out of the meta map in their own channel in this subworkflow
    // When adding a new field in the assets/schema_input.json, also add it in the meta map
    // And in this scrip, add a new map operation and a new output corresponding to this input

    // If any of the asset does not exist, then we return null
    // That way, the channel will be empty and does not trigger anything

    def reduce = { meta -> meta.subMap(['genome', 'id', 'source', 'source_vcf', 'source_version', 'species']) }

    intervals_bed_branch = asset.branch { meta, _fasta ->
        file: meta.intervals_bed
        return [reduce(meta), file(meta.intervals_bed)]
        other: true
        return null
    }
    intervals_bed = intervals_bed_branch.file


    // If ends with .gz, decompress it
    // If any of the asset exists, then adding run_tools to false and skip the asset creation from the fasta file
    fasta_branch = asset.branch { meta, fasta_ ->
        file: fasta_
        return [reduce(meta) + [decompress_fasta: fasta_.endsWith('.gz') ?: false] + [run_bowtie1: meta.bowtie1_index ? false : true] + [run_bowtie2: meta.bowtie2_index ? false : true] + [run_bwamem1: meta.bwamem1_index ? false : true] + [run_bwamem2: meta.bwamem2_index ? false : true] + [run_dragmap: meta.dragmap_hashtable ? false : true] + [run_faidx: meta.fasta_fai && meta.fasta_sizes ? false : true] + [run_gatkdict: meta.fasta_dict ? false : true] + [run_hisat2: meta.hisat2_index ? false : true] + [run_intervals: meta.intervals_bed ? false : true] + [run_kallisto: meta.kallisto_index ? false : true] + [run_msisenpro: meta.msisensorpro_list ? false : true] + [run_rsem: meta.rsem_index ? false : true] + [run_rsem_make_transcript_fasta: meta.transcript_fasta ? false : true] + [run_salmon: meta.salmon_index ? false : true] + [run_star: meta.star_index ? false : true], file(fasta_)]
        other: true
        return null
    }
    fasta = fasta_branch.file


    fasta_dict_branch = asset.branch { meta, _fasta ->
        file: meta.fasta_dict
        return [reduce(meta), file(meta.fasta_dict)]
        other: true
        return null
    }
    fasta_dict = fasta_dict_branch.file


    // If we have intervals_bed, then we don't need to run faidx
    fasta_fai_branch = asset.branch { meta, _fasta ->
        file: meta.fasta_fai
        return [reduce(meta) + [run_intervals: meta.intervals_bed ? false : true], file(meta.fasta_fai)]
        other: true
        return null
    }
    fasta_fai = fasta_fai_branch.file


    fasta_sizes_branch = asset.branch { meta, _fasta ->
        file: meta.fasta_sizes
        return [reduce(meta), file(meta.fasta_sizes)]
        other: true
        return null
    }
    fasta_sizes = fasta_sizes_branch.file


    // If ends with .gz, decompress it
    // If any of the asset exists, then adding run_tools to false and skip the asset creation from the annotation derived file (gff, gtf or transcript_fasta)
    gff_branch = asset.branch { meta, fasta_ ->
        file: meta.gff
        return [reduce(meta) + [decompress_gff: meta.gff.endsWith('.gz') ?: false] + [run_gffread: fasta_ && !meta.gtf ?: false] + [run_hisat2: meta.splice_sites ? false : true], file(meta.gff)]
        other: true
        return null
    }
    gff = gff_branch.file


    // If ends with .gz, decompress it
    // If any of the asset exists, then adding run_tools to false and skip the asset creation from the annotation derived file (gff, gtf or transcript_fasta)
    gtf_branch = asset.branch { meta, _fasta ->
        file: meta.gtf
        return [reduce(meta) + [decompress_gtf: meta.gtf.endsWith('.gz') ?: false] + [run_hisat2: meta.splice_sites ? false : true], file(meta.gtf)]
        other: true
        return null
    }
    gtf = gtf_branch.file


    splice_sites_branch = asset.branch { meta, _fasta ->
        file: meta.splice_sites
        return [reduce(meta), file(meta.splice_sites)]
        other: true
        return null
    }
    splice_sites = splice_sites_branch.file


    // If any of the asset exists, then adding run_tools to false and skip the asset creation from the annotation derived file (gff, gtf or transcript_fasta)
    transcript_fasta_branch = asset.branch { meta, _fasta ->
        file: meta.transcript_fasta
        return [reduce(meta) + [run_hisat2: meta.hisat2_index ? false : true] + [run_kallisto: meta.kallisto_index ? false : true] + [run_rsem: meta.rsem_index ? false : true] + [run_salmon: meta.salmon_index ? false : true] + [run_star: meta.star_index ? false : true], file(meta.transcript_fasta)]
        other: true
        return null
    }
    transcript_fasta = transcript_fasta_branch.file


    // Using transpose here because we want to catch vcf with globs in the path because of nf-core/Sarek
    // return a file, because we can catch globs this way, but it create issues with publishing
    // If we already have the vcf_tbi, then we don't need to index the vcf
    vcf_branch = asset.branch { meta, _fasta ->
        file: meta.vcf
        return [reduce(meta) + [run_tabix: meta.vcf_tbi ? false : true], file(meta.vcf)]
        other: true
        return null
    }
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
