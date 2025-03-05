include { CREATE_FROM_FASTA_AND_ANNOTATION } from '../../subworkflows/local/create_from_fasta_and_annotation'
include { CREATE_FROM_FASTA_ONLY           } from '../../subworkflows/local/create_from_fasta_only'
include { EXTRACT_REFERENCE                } from '../../subworkflows/local/extract_reference'
include { INDEX_VCF                        } from '../../subworkflows/local/index_vcf'
include { YAML_TO_CHANNEL                  } from '../../subworkflows/local/yaml_to_channel'

workflow REFERENCES {
    take:
    asset // Channel: [meta, fasta]
    tools // List: Can contain any combination of tools of the list of available tools, or just no_tools

    main:
    def need_extract = { channel, type ->
        channel
            .map { meta, reference_ -> [meta + [reference: type], reference_] }
            .branch { _meta, reference_ ->
                tar: reference_.endsWith('.tar.gz')
                gz: reference_.endsWith('.gz')
                zip: reference_.endsWith('.zip')
                other: true
            }
    }

    versions = Channel.empty()

    // Create channels from the input file(s)
    // Channels are empty when no reference corresponds
    YAML_TO_CHANNEL(asset, tools)

    intervals_bed = YAML_TO_CHANNEL.out.intervals_bed
    fasta_dict = YAML_TO_CHANNEL.out.fasta_dict
    fasta_fai = YAML_TO_CHANNEL.out.fasta_fai
    fasta_sizes = YAML_TO_CHANNEL.out.fasta_sizes
    splice_sites = YAML_TO_CHANNEL.out.splice_sites
    transcript_fasta = YAML_TO_CHANNEL.out.transcript_fasta
    vcf = YAML_TO_CHANNEL.out.vcf

    // Check if reference needs to be extracted
    // Copy the reference type to the meta
    // (We do not uncompress VCFs)
    ascat_alleles_input = need_extract(YAML_TO_CHANNEL.out.ascat_alleles, 'ascat_alleles')
    ascat_loci_input = need_extract(YAML_TO_CHANNEL.out.ascat_loci, 'ascat_loci')
    ascat_loci_gc_input = need_extract(YAML_TO_CHANNEL.out.ascat_loci_gc, 'ascat_loci_gc')
    ascat_loci_rt_input = need_extract(YAML_TO_CHANNEL.out.ascat_loci_rt, 'ascat_loci_rt')
    chr_dir_input = need_extract(YAML_TO_CHANNEL.out.chr_dir, 'chr_dir')
    fasta_input = need_extract(YAML_TO_CHANNEL.out.fasta, 'fasta')
    gff_input = need_extract(YAML_TO_CHANNEL.out.gff, 'gff')
    gtf_input = need_extract(YAML_TO_CHANNEL.out.gtf, 'gtf')

    extract_gz = Channel
        .empty()
        .mix(
            ascat_alleles_input.gz,
            ascat_loci_input.gz,
            ascat_loci_gc_input.gz,
            ascat_loci_rt_input.gz,
            chr_dir_input.gz,
            fasta_input.gz,
            gff_input.gz,
            gtf_input.gz,
        )

    extract_tar = Channel
        .empty()
        .mix(
            ascat_alleles_input.tar,
            ascat_loci_input.tar,
            ascat_loci_gc_input.tar,
            ascat_loci_rt_input.tar,
            chr_dir_input.tar,
            fasta_input.tar,
            gff_input.tar,
            gtf_input.tar,
        )

    extract_zip = Channel
        .empty()
        .mix(
            ascat_alleles_input.zip,
            ascat_loci_input.zip,
            ascat_loci_gc_input.zip,
            ascat_loci_rt_input.zip,
            chr_dir_input.zip,
            fasta_input.zip,
            gff_input.zip,
            gtf_input.zip,
        )

    // Uncompress any assets that need to be
    // From any archive format
    EXTRACT_REFERENCE(
        extract_gz,
        extract_tar,
        extract_zip,
    )

    extracted_asset = EXTRACT_REFERENCE.out.extracted.branch { meta_, _extracted_asset ->
        ascat_alleles: meta_.reference == 'ascat_alleles'
        ascat_loci: meta_.reference == 'ascat_loci'
        ascat_loci_gc: meta_.reference == 'ascat_loci_gc'
        ascat_loci_rt: meta_.reference == 'ascat_loci_rt'
        chr_dir: meta_.reference == 'chr_dir'
        fasta: meta_.reference == 'fasta'
        gff: meta_.reference == 'gff'
        gtf: meta_.reference == 'gtf'
        other: true
    }

    // This is a sanity check
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

    // This works with a mixture of input and computed references
    fasta_dict = fasta_dict.mix(CREATE_FROM_FASTA_ONLY.out.fasta_dict)
    fasta_fai = fasta_fai.mix(CREATE_FROM_FASTA_ONLY.out.fasta_fai)
    fasta_sizes = fasta_sizes.mix(CREATE_FROM_FASTA_ONLY.out.fasta_sizes)
    gtf = gtf.mix(CREATE_FROM_FASTA_AND_ANNOTATION.out.gtf)
    intervals_bed = intervals_bed.mix(CREATE_FROM_FASTA_ONLY.out.intervals_bed)
    splice_sites = splice_sites.mix(CREATE_FROM_FASTA_AND_ANNOTATION.out.splice_sites)
    transcript_fasta = transcript_fasta.mix(CREATE_FROM_FASTA_AND_ANNOTATION.out.transcript_fasta)

    // TODO: This does not work YET with a mixture of input and computed references
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
    versions = versions.mix(CREATE_FROM_FASTA_AND_ANNOTATION.out.versions)
    versions = versions.mix(CREATE_FROM_FASTA_ONLY.out.versions)
    versions = versions.mix(EXTRACT_REFERENCE.out.versions)
    versions = versions.mix(INDEX_VCF.out.versions)

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
