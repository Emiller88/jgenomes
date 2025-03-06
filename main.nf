#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/references
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/references
    Website: https://nf-co.re/references
    Slack  : https://nfcore.slack.com/channels/references
----------------------------------------------------------------------------------------
*/

nextflow.preview.output = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap        } from 'plugin/nf-schema'
include { methodsDescriptionText  } from './subworkflows/local/utils_nfcore_references_pipeline'
include { paramsSummaryMultiqc    } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML  } from './subworkflows/nf-core/utils_nfcore_pipeline'

include { MULTIQC                 } from './modules/nf-core/multiqc'

include { EXTRACT_ARCHIVE         } from './subworkflows/local/extract_archive'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_references_pipeline'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_references_pipeline'
include { YAML_TO_CHANNEL         } from './subworkflows/local/yaml_to_channel'

include { REFERENCES              } from "./workflows/references"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    main:
    // SUBWORKFLOW: Run initialisation tasks
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        args,
        params.outdir,
        params.input,
    )

    // WORKFLOW: Run main workflow
    if (!params.tools) {
        log.warn("No tools specified")
    }

    NFCORE_REFERENCES(PIPELINE_INITIALISATION.out.references, params.tools ?: "no_tools")

    ch_multiqc_files = Channel.empty()

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(softwareVersionsToYAML(NFCORE_REFERENCES.out.versions).collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_references_software_mqc_versions.yml', sort: true, newLine: true))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    multiqc_data = MULTIQC.out.data
    multiqc_plots = MULTIQC.out.plots
    multiqc_report = MULTIQC.out.report

    // SUBWORKFLOW: Run completion tasks
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        MULTIQC.out.report.toList(),
    )

    publish:
    multiqc_data >> 'multiqc'
    multiqc_plots >> 'multiqc'
    multiqc_report >> 'multiqc'
}

output {
    'multiqc' {
        path "multiqc"
    }
    'reference' {
        path { meta, _file ->
            { file ->
                if (meta.file == "bowtie1_index") {
                    "${meta.species}/${meta.source}/${meta.genome}/Sequence/BowtieIndex/version1.3.1"
                }
                else if (meta.file == "bowtie2_index") {
                    "${meta.species}/${meta.source}/${meta.genome}/Sequence/Bowtie2Index/version2.5.2"
                }
                else if (meta.file == "bwamem1_index") {
                    "${meta.species}/${meta.source}/${meta.genome}/Sequence/BWAIndex/version0.7.18"
                }
                else if (meta.file == "bwamem2_index") {
                    "${meta.species}/${meta.source}/${meta.genome}/Sequence/BWAmem2Index/version2.2.1"
                }
                else if (meta.file == "dragmap_hashmap") {
                    "${meta.species}/${meta.source}/${meta.genome}/Sequence/dragmap/version1.2.1"
                }
                else if (meta.file == "fasta") {
                    "${meta.species}/${meta.source}/${meta.genome}/Sequence/WholeGenomeFasta/${file}"
                }
                else if (meta.file == "fasta_dict") {
                    "${meta.species}/${meta.source}/${meta.genome}/Sequence/WholeGenomeFasta/${file}"
                }
                else if (meta.file == "fasta_fai") {
                    "${meta.species}/${meta.source}/${meta.genome}/Sequence/WholeGenomeFasta/${file}"
                }
                else if (meta.file == "fasta_sizes") {
                    "${meta.species}/${meta.source}/${meta.genome}/Sequence/WholeGenomeFasta/${file}"
                }
                else if (meta.file == "gff") {
                    "${meta.species}/${meta.source}/${meta.genome}/Annotation/Genes/${file}"
                }
                else if (meta.file == "gtf") {
                    "${meta.species}/${meta.source}/${meta.genome}/Annotation/Genes/${file}"
                }
                else if (meta.file == "hisat2_index") {
                    meta.source_version == "unknown"
                        ? "${meta.species}/${meta.source}/${meta.genome}/Sequence/Hisat2Index/version2.2.1"
                        : "${meta.species}/${meta.source}/${meta.genome}/Sequence/Hisat2Index/${meta.source_version}/version2.2.1"
                }
                else if (meta.file == "intervals_bed") {
                    "${meta.species}/${meta.source}/${meta.genome}/Annotation/intervals/${file}"
                }
                else if (meta.file == "kallisto_index") {
                    meta.source_version == "unknown"
                        ? "${meta.species}/${meta.source}/${meta.genome}/Sequence/KallistoIndex/version0.51.1/${file}"
                        : "${meta.species}/${meta.source}/${meta.genome}/Sequence/KallistoIndex/${meta.source_version}/version0.51.1/${file}"
                }
                else if (meta.file == "msisensorpro_list") {
                    "${meta.species}/${meta.source}/${meta.genome}/Annotation/msisensorpro/${file}"
                }
                else if (meta.file == "rsem_index") {
                    meta.source_version == "unknown"
                        ? "${meta.species}/${meta.source}/${meta.genome}/Sequence/RSEMIndex/version1.3.1/"
                        : "${meta.species}/${meta.source}/${meta.genome}/Sequence/RSEMIndex/${meta.source_version}/version1.3.1/"
                }
                else if (meta.file == "salmon_index") {
                    meta.source_version == "unknown"
                        ? "${meta.species}/${meta.source}/${meta.genome}/Sequence/SalmonIndex/version1.10.3/"
                        : "${meta.species}/${meta.source}/${meta.genome}/Sequence/SalmonIndex/${meta.source_version}/version1.10.3/"
                }
                else if (meta.file == "splice_sites") {
                    "${meta.species}/${meta.source}/${meta.genome}/Sequence/SpliceSites/${file}"
                }
                else if (meta.file == "star_index") {
                    meta.source_version == "unknown"
                        ? "${meta.species}/${meta.source}/${meta.genome}/Sequence/STARIndex/version2.7.11b/"
                        : "${meta.species}/${meta.source}/${meta.genome}/Sequence/STARIndex/${meta.source_version}/version2.7.11b/"
                }
                else if (meta.file == "transcript_fasta") {
                    "${meta.species}/${meta.source}/${meta.genome}/Sequence/TranscriptFasta/${file}"
                }
                else if (meta.file == "${meta.type}_vcf") {
                    "${meta.species}/${meta.source}/${meta.genome}/Annotation/${meta.source_vcf}/${file}"
                }
                else if (meta.file == "${meta.type}_vcf_tbi") {
                    "${meta.species}/${meta.source}/${meta.genome}/Annotation/${meta.source_vcf}/${file}"
                }
                else {
                    null
                }
            }
        }

        index {
            path "index.json"
            mapper { meta, reference -> ["${meta.file}:${reference}"] }
        }
    }
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Build references depending on type of asset and the tools specified
//
workflow NFCORE_REFERENCES {
    take:
    references
    tools      // list of tools to use to build references

    main:
    // Helper closure to check if a reference needs to be extracted
    // Add the reference type to the meta
    // Depending on the extension, return the appropriate channel
    def need_extract = { channel, type ->
        channel
            .map { meta, reference_ -> [meta + [reference: type], reference_] }
            .branch { _meta, reference_ ->
                to_extract: reference_.toString().endsWith('.gz') || reference_.toString().endsWith('.zip')
                not_extracted: true
            }
    }

    YAML_TO_CHANNEL(references, params.tools ?: "no_tools")

    // References that need to be extracted
    // (VCFs are not extracted)
    ascat_alleles_input = need_extract(YAML_TO_CHANNEL.out.ascat_alleles, 'ascat_alleles')
    ascat_loci_input = need_extract(YAML_TO_CHANNEL.out.ascat_loci, 'ascat_loci')
    ascat_loci_gc_input = need_extract(YAML_TO_CHANNEL.out.ascat_loci_gc, 'ascat_loci_gc')
    ascat_loci_rt_input = need_extract(YAML_TO_CHANNEL.out.ascat_loci_rt, 'ascat_loci_rt')
    chr_dir_input = need_extract(YAML_TO_CHANNEL.out.chr_dir, 'chr_dir')
    fasta_input = need_extract(YAML_TO_CHANNEL.out.fasta, 'fasta')
    gff_input = need_extract(YAML_TO_CHANNEL.out.gff, 'gff')
    gtf_input = need_extract(YAML_TO_CHANNEL.out.gtf, 'gtf')

    // gather all archived references
    archive_to_extract = Channel
        .empty()
        .mix(
            ascat_alleles_input.to_extract,
            ascat_loci_input.to_extract,
            ascat_loci_gc_input.to_extract,
            ascat_loci_rt_input.to_extract,
            chr_dir_input.to_extract,
            fasta_input.to_extract,
            gff_input.to_extract,
            gtf_input.to_extract,
        )

    // Extract references from any archive format
    EXTRACT_ARCHIVE(
        archive_to_extract
    )

    // return to the appropriate channels
    extracted_asset = EXTRACT_ARCHIVE.out.extracted.branch { meta_, _extracted_asset ->
        ascat_alleles: meta_.reference == 'ascat_alleles'
        ascat_loci: meta_.reference == 'ascat_loci'
        ascat_loci_gc: meta_.reference == 'ascat_loci_gc'
        ascat_loci_rt: meta_.reference == 'ascat_loci_rt'
        chr_dir: meta_.reference == 'chr_dir'
        fasta: meta_.reference == 'fasta'
        gff: meta_.reference == 'gff'
        gtf: meta_.reference == 'gtf'
        non_assigned: true
    }

    // This is a confidence check
    extracted_asset.non_assigned.view { log.warn("Non assigned extracted asset: " + it) }

    // WORKFLOW: Run pipeline
    // Mix the references that were extracted with the references that did not need to be extracted
    // Some references are not extracted because they are usually not stored in an archived format
    // TODO: check if more references need to be extracted
    REFERENCES(
        ascat_alleles_input.not_extracted.mix(extracted_asset.ascat_alleles),
        ascat_loci_input.not_extracted.mix(extracted_asset.ascat_loci),
        ascat_loci_gc_input.not_extracted.mix(extracted_asset.ascat_loci_gc),
        ascat_loci_rt_input.not_extracted.mix(extracted_asset.ascat_loci_rt),
        chr_dir_input.not_extracted.mix(extracted_asset.chr_dir),
        fasta_input.not_extracted.mix(extracted_asset.fasta),
        YAML_TO_CHANNEL.out.fasta_dict,
        YAML_TO_CHANNEL.out.fasta_fai,
        YAML_TO_CHANNEL.out.fasta_sizes,
        gff_input.not_extracted.mix(extracted_asset.gff),
        gtf_input.not_extracted.mix(extracted_asset.gtf),
        YAML_TO_CHANNEL.out.intervals_bed,
        YAML_TO_CHANNEL.out.splice_sites,
        YAML_TO_CHANNEL.out.transcript_fasta,
        YAML_TO_CHANNEL.out.vcf,
        tools,
    )

    emit:
    reference = REFERENCES.out.reference
    versions  = REFERENCES.out.versions.mix(EXTRACT_ARCHIVE.out.versions)
}
