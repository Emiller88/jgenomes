include { ASSET_TO_CHANNEL                 } from '../../subworkflows/local/asset_to_channel'
include { CREATE_FROM_FASTA_AND_ANNOTATION } from '../../subworkflows/local/create_from_fasta_and_annotation'
include { CREATE_FROM_FASTA_ONLY           } from '../../subworkflows/local/create_from_fasta_only'
include { INDEX_VCF                        } from '../../subworkflows/local/index_vcf'
include { UNCOMPRESS_ASSET                 } from '../../subworkflows/local/uncompress_asset'

workflow REFERENCES {
    take:
    asset // Channel: [meta, fasta]
    tools // List: Can contain any combination of tools of the list of available tools, or just no_tools

    main:
    versions = Channel.empty()

    // Create channels from the input file(s)
    // Channels are empty when no assets are corresponding
    ASSET_TO_CHANNEL(asset, tools)

    intervals_bed = ASSET_TO_CHANNEL.out.intervals_bed
    fasta_dict = ASSET_TO_CHANNEL.out.fasta_dict
    fasta_fai = ASSET_TO_CHANNEL.out.fasta_fai
    fasta_sizes = ASSET_TO_CHANNEL.out.fasta_sizes
    splice_sites = ASSET_TO_CHANNEL.out.splice_sites
    transcript_fasta = ASSET_TO_CHANNEL.out.transcript_fasta
    vcf = ASSET_TO_CHANNEL.out.vcf

    // Assess if assets needs to be uncompress or not
    // (We do not uncompress VCFs)

    ascat_alleles_input = ASSET_TO_CHANNEL.out.ascat_alleles
        .map { meta, ascat_alleles_ -> [meta + [file: 'ascat_alleles'], ascat_alleles_] }
        .branch { _meta, archive_ ->
            extract_tar: archive_.endsWith('.tar.gz')
            extract_gz: archive_.endsWith('.gz')
            extract_zip: archive_.endsWith('.zip')
            other: true
        }

    ascat_loci_input = ASSET_TO_CHANNEL.out.ascat_loci
        .map { meta, ascat_loci_ -> [meta + [file: 'ascat_loci'], ascat_loci_] }
        .branch { _meta, archive_ ->
            extract_tar: archive_.endsWith('.tar.gz')
            extract_gz: archive_.endsWith('.gz')
            extract_zip: archive_.endsWith('.zip')
            other: true
        }

    ascat_loci_gc_input = ASSET_TO_CHANNEL.out.ascat_loci_gc
        .map { meta, ascat_loci_gc_ -> [meta + [file: 'ascat_loci_gc'], ascat_loci_gc_] }
        .branch { _meta, archive_ ->
            extract_tar: archive_.endsWith('.tar.gz')
            extract_gz: archive_.endsWith('.gz')
            extract_zip: archive_.endsWith('.zip')
            other: true
        }

    ascat_loci_rt_input = ASSET_TO_CHANNEL.out.ascat_loci_rt
        .map { meta, ascat_loci_rt_ -> [meta + [file: 'ascat_loci_rt'], ascat_loci_rt_] }
        .branch { _meta, archive_ ->
            extract_tar: archive_.endsWith('.tar.gz')
            extract_gz: archive_.endsWith('.gz')
            extract_zip: archive_.endsWith('.zip')
            other: true
        }

    chr_dir_input = ASSET_TO_CHANNEL.out.chr_dir
        .map { meta, chr_dir_ -> [meta + [file: 'chr_dir'], chr_dir_] }
        .branch { _meta, archive_ ->
            extract_tar: archive_.endsWith('.tar.gz')
            extract_gz: archive_.endsWith('.gz')
            extract_zip: archive_.endsWith('.zip')
            other: true
        }

    fasta_input = ASSET_TO_CHANNEL.out.fasta
        .map { meta, fasta_ -> [meta + [file: 'fasta'], fasta_] }
        .branch { _meta, archive_ ->
            extract_tar: archive_.endsWith('.tar.gz')
            extract_gz: archive_.endsWith('.gz')
            extract_zip: archive_.endsWith('.zip')
            other: true
        }

    gff_input = ASSET_TO_CHANNEL.out.gff
        .map { meta, gff_ -> [meta + [file: 'gff'], gff_] }
        .branch { _meta, archive_ ->
            extract_tar: archive_.endsWith('.tar.gz')
            extract_gz: archive_.endsWith('.gz')
            extract_zip: archive_.endsWith('.zip')
            other: true
        }

    gtf_input = ASSET_TO_CHANNEL.out.gtf
        .map { meta, gtf_ -> [meta + [file: 'gtf'], gtf_] }
        .branch { _meta, archive_ ->
            extract_tar: archive_.endsWith('.tar.gz')
            extract_gz: archive_.endsWith('.gz')
            extract_zip: archive_.endsWith('.zip')
            other: true
        }

    files_to_extract_gz = Channel
        .empty()
        .mix(
            ascat_alleles_input.extract_gz,
            ascat_loci_input.extract_gz,
            ascat_loci_gc_input.extract_gz,
            ascat_loci_rt_input.extract_gz,
            chr_dir_input.extract_gz,
            fasta_input.extract_gz,
            gff_input.extract_gz,
            gtf_input.extract_gz,
        )

    files_to_extract_tar = Channel
        .empty()
        .mix(
            ascat_alleles_input.extract_tar,
            ascat_loci_input.extract_tar,
            ascat_loci_gc_input.extract_tar,
            ascat_loci_rt_input.extract_tar,
            chr_dir_input.extract_tar,
            fasta_input.extract_tar,
            gff_input.extract_tar,
            gtf_input.extract_tar,
        )

    files_to_extract_zip = Channel
        .empty()
        .mix(
            ascat_alleles_input.extract_zip,
            ascat_loci_input.extract_zip,
            ascat_loci_gc_input.extract_zip,
            ascat_loci_rt_input.extract_zip,
            chr_dir_input.extract_zip,
            fasta_input.extract_zip,
            gff_input.extract_zip,
            gtf_input.extract_zip,
        )

    // Uncompress any assets that need to be
    UNCOMPRESS_ASSET(
        files_to_extract_gz,
        files_to_extract_tar,
        files_to_extract_zip,
    )

    extracted_asset = UNCOMPRESS_ASSET.out.extracted.branch { meta_, _extracted_asset ->
        ascat_alleles: meta_.file == 'ascat_alleles'
        ascat_loci: meta_.file == 'ascat_loci'
        ascat_loci_gc: meta_.file == 'ascat_loci_gc'
        ascat_loci_rt: meta_.file == 'ascat_loci_rt'
        chr_dir: meta_.file == 'chr_dir'
        fasta: meta_.file == 'fasta'
        gff: meta_.file == 'gff'
        gtf: meta_.file == 'gtf'
        other: true
    }

    extracted_asset.other.view { "Non assigned extracted asset: " + it }


    // This covers a mixture of compressed and uncompressed assets
    ascat_alleles = ascat_alleles_input.other.mix(extracted_asset.ascat_alleles)
    ascat_loci = ascat_loci_input.other.mix(extracted_asset.ascat_loci)
    ascat_loci_gc = ascat_loci_gc_input.other.mix(extracted_asset.ascat_loci_gc)
    ascat_loci_rt = ascat_loci_rt_input.other.mix(extracted_asset.ascat_loci_rt)
    chr_dir = chr_dir_input.other.mix(extracted_asset.chr_dir)
    fasta = fasta_input.other.mix(extracted_asset.fasta)
    gff = gff_input.other.mix(extracted_asset.gff)
    gtf = gtf_input.other.mix(extracted_asset.gtf)

    // Create reference assets from fasta only
    CREATE_FROM_FASTA_ONLY(
        fasta,
        fasta_fai,
        tools.split(',').contains('bowtie1'),
        tools.split(',').contains('bowtie2'),
        tools.split(',').contains('bwamem1'),
        tools.split(',').contains('bwamem2'),
        tools.split(',').contains('createsequencedictionary'),
        tools.split(',').contains('dragmap'),
        tools.split(',').contains('faidx'),
        tools.split(',').contains('intervals'),
        tools.split(',').contains('msisensorpro'),
        tools.split(',').contains('sizes'),
    )

    // Create reference assets from fasta and gene annotation (so gff, gtf and transcript_fasta)
    CREATE_FROM_FASTA_AND_ANNOTATION(
        fasta,
        gff,
        gtf,
        splice_sites,
        transcript_fasta,
        tools.split(',').contains('hisat2'),
        tools.split(',').contains('hisat2_extractsplicesites'),
        tools.split(',').contains('kallisto'),
        tools.split(',').contains('rsem'),
        tools.split(',').contains('rsem_make_transcript_fasta'),
        tools.split(',').contains('salmon'),
        tools.split(',').contains('star'),
    )

    // Index VCF
    INDEX_VCF(
        vcf,
        tools.split(',').contains('tabix'),
    )

    // This works with a mixture of input and computed assets
    fasta_dict = fasta_dict.mix(CREATE_FROM_FASTA_ONLY.out.fasta_dict)
    fasta_fai = fasta_fai.mix(CREATE_FROM_FASTA_ONLY.out.fasta_fai)
    fasta_sizes = fasta_sizes.mix(CREATE_FROM_FASTA_ONLY.out.fasta_sizes)
    gtf = gtf.mix(CREATE_FROM_FASTA_AND_ANNOTATION.out.gtf)
    intervals_bed = intervals_bed.mix(CREATE_FROM_FASTA_ONLY.out.intervals_bed)
    splice_sites = splice_sites.mix(CREATE_FROM_FASTA_AND_ANNOTATION.out.splice_sites)
    transcript_fasta = transcript_fasta.mix(CREATE_FROM_FASTA_AND_ANNOTATION.out.transcript_fasta)

    // TODO: This does not work YET with a mixture of input and computed assets
    bowtie1_index = CREATE_FROM_FASTA_ONLY.out.bowtie1_index
    bowtie2_index = CREATE_FROM_FASTA_ONLY.out.bowtie2_index
    bwamem1_index = CREATE_FROM_FASTA_ONLY.out.bwamem1_index
    bwamem2_index = CREATE_FROM_FASTA_ONLY.out.bwamem2_index
    dragmap_hashmap = CREATE_FROM_FASTA_ONLY.out.dragmap_hashmap
    hisat2_index = CREATE_FROM_FASTA_AND_ANNOTATION.out.hisat2_index
    kallisto_index = CREATE_FROM_FASTA_AND_ANNOTATION.out.kallisto_index
    msisensorpro_list = CREATE_FROM_FASTA_ONLY.out.msisensorpro_list
    rsem_index = CREATE_FROM_FASTA_AND_ANNOTATION.out.rsem_index
    salmon_index = CREATE_FROM_FASTA_AND_ANNOTATION.out.salmon_index
    star_index = CREATE_FROM_FASTA_AND_ANNOTATION.out.star_index
    vcf_tbi = INDEX_VCF.out.vcf_tbi

    // TODO: Refactor this with topics
    versions = versions.mix(CREATE_FROM_FASTA_ONLY.out.versions)
    versions = versions.mix(CREATE_FROM_FASTA_AND_ANNOTATION.out.versions)
    versions = versions.mix(INDEX_VCF.out.versions)
    versions = versions.mix(UNCOMPRESS_ASSET.out.versions)

    reference = Channel
        .empty()
        .mix(
            bowtie1_index.map { meta, reference_ -> [meta + [file: 'bowtie1_index'], reference_] },
            bowtie2_index.map { meta, reference_ -> [meta + [file: 'bowtie2_index'], reference_] },
            bwamem1_index.map { meta, reference_ -> [meta + [file: 'bwamem1_index'], reference_] },
            bwamem2_index.map { meta, reference_ -> [meta + [file: 'bwamem2_index'], reference_] },
            dragmap_hashmap.map { meta, reference_ -> [meta + [file: 'dragmap_hashmap'], reference_] },
            fasta.map { meta, reference_ -> [meta + [file: 'fasta'], reference_] },
            fasta_dict.map { meta, reference_ -> [meta + [file: 'fasta_dict'], reference_] },
            fasta_fai.map { meta, reference_ -> [meta + [file: 'fasta_fai'], reference_] },
            fasta_sizes.map { meta, reference_ -> [meta + [file: 'fasta_sizes'], reference_] },
            gff.map { meta, reference_ -> [meta + [file: 'gff'], reference_] },
            gtf.map { meta, reference_ -> [meta + [file: 'gtf'], reference_] },
            hisat2_index.map { meta, reference_ -> [meta + [file: 'hisat2_index'], reference_] },
            intervals_bed.map { meta, reference_ -> [meta + [file: 'intervals_bed'], reference_] },
            kallisto_index.map { meta, reference_ -> [meta + [file: 'kallisto_index'], reference_] },
            msisensorpro_list.map { meta, reference_ -> [meta + [file: 'msisensorpro_list'], reference_] },
            rsem_index.map { meta, reference_ -> [meta + [file: 'rsem_index'], reference_] },
            salmon_index.map { meta, reference_ -> [meta + [file: 'salmon_index'], reference_] },
            splice_sites.map { meta, reference_ -> [meta + [file: 'splice_sites'], reference_] },
            star_index.map { meta, reference_ -> [meta + [file: 'star_index'], reference_] },
            transcript_fasta.map { meta, reference_ -> [meta + [file: 'transcript_fasta'], reference_] },
            vcf.map { meta, reference_ -> [meta + [file: "${meta.type}_vcf"], reference_] },
            vcf_tbi.map { meta, reference_ -> [meta + [file: "${meta.type}_vcf_tbi"], reference_] },
        )

    emit:
    reference // channel: [meta, *]
    versions  // channel: [versions.yml]

    publish:
    reference >> 'reference'
}
