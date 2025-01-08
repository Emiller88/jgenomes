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

include { MULTIQC                 } from './modules/nf-core/multiqc/main'

include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_references_pipeline'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_references_pipeline'

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

    NFCORE_REFERENCES(PIPELINE_INITIALISATION.out.asset, params.tools ?: "no_tools")

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
        ch_multiqc_logo.getVal(),
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
                else if (meta.file == "vcf") {
                    "${meta.species}/${meta.source}/${meta.genome}/Annotation/${meta.source_vcf}/${file}"
                }
                else if (meta.file == "vcf_tbi") {
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
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_REFERENCES {
    take:
    input // channel: asset reference yml file read in from --input
    tools // list of tools to use to build references

    main:
    // WORKFLOW: Run pipeline
    REFERENCES(input, tools)

    emit:
    reference = REFERENCES.out.reference
    versions  = REFERENCES.out.versions
}
