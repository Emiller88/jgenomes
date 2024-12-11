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

include { REFERENCES              } from "./workflows/references/main"
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
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
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
        ch_multiqc_logo.toList(),
        [],
        []
    )

    multiqc_data = MULTIQC.out.data
    multiqc_plots = MULTIQC.out.plots
    multiqc_report = MULTIQC.out.report

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        MULTIQC.out.report.toList()
    )

    publish:
    multiqc_data >> 'multiqc'
    multiqc_plots >> 'multiqc'
    multiqc_report >> 'multiqc'
}

output {
    'bowtie1_index' {
        path { meta, _bowtie1_index -> { _file -> "${meta.species}/${meta.source}/${meta.genome}/Sequence/BowtieIndex/version1.3.1" } }
    }
    'bowtie2_index' {
        path { meta, _bowtie2_index -> { _file -> "${meta.species}/${meta.source}/${meta.genome}/Sequence/Bowtie2Index/version2.5.2" } }
    }
    'bwamem1_index' {
        path { meta, _bwamem1_index -> { _file -> "${meta.species}/${meta.source}/${meta.genome}/Sequence/BWAIndex/version0.7.18" } }
    }
    'bwamem2_index' {
        path { meta, _bwamem2_index -> { _file -> "${meta.species}/${meta.source}/${meta.genome}/Sequence/BWAmem2Index/version2.2.1" } }
    }
    'dragmap_hashmap' {
        path { meta, _index -> { _file -> "${meta.species}/${meta.source}/${meta.genome}/Sequence/dragmap/version1.2.1" } }
    }
    'fasta' {
        path { meta, _fasta -> { file -> "${meta.species}/${meta.source}/${meta.genome}/Sequence/WholeGenomeFasta/${file}" } }
    }
    'fasta_dict' {
        path { meta, _fasta_dict -> { file -> "${meta.species}/${meta.source}/${meta.genome}/Sequence/WholeGenomeFasta/${file}" } }
    }
    'fasta_fai' {
        path { meta, _fasta_fai -> { file -> "${meta.species}/${meta.source}/${meta.genome}/Sequence/WholeGenomeFasta/${file}" } }
    }
    'fasta_sizes' {
        path { meta, _fasta_sizes -> { file -> "${meta.species}/${meta.source}/${meta.genome}/Sequence/WholeGenomeFasta/${file}" } }
    }
    'gff' {
        path { meta, _gff -> { file -> "${meta.species}/${meta.source}/${meta.genome}/Annotation/Genes/${file}" } }
    }
    'gtf' {
        path { meta, _gtf -> { file -> "${meta.species}/${meta.source}/${meta.genome}/Annotation/Genes/${file}" } }
    }
    'hisat2_index' {
        path { meta, _hisat2_index ->
            { _file ->
                meta.source_version == "unknown"
                    ? "${meta.species}/${meta.source}/${meta.genome}/Sequence/Hisat2Index/version2.2.1"
                    : "${meta.species}/${meta.source}/${meta.genome}/Sequence/Hisat2Index/${meta.source_version}/version2.2.1"
            }
        }
    }
    'intervals_bed' {
        path { meta, _intervals_bed -> { file -> "${meta.species}/${meta.source}/${meta.genome}/Annotation/intervals/${file}" } }
    }
    'kallisto_index' {
        path { meta, _kallisto_index ->
            { _file ->
                meta.source_version == "unknown"
                    ? "${meta.species}/${meta.source}/${meta.genome}/Sequence/KallistoIndex/version0.51.1/${file}"
                    : "${meta.species}/${meta.source}/${meta.genome}/Sequence/KallistoIndex/${meta.source_version}/version0.51.1/${file}"
            }
        }
    }
    'msisensorpro_list' {
        path { meta, _index -> { file -> "${meta.species}/${meta.source}/${meta.genome}/Annotation/msisensorpro/${file}" } }
    }
    'multiqc_data' {
        path { _multiqc_data -> { _file -> "multiqc/multiqc_data" } }
    }
    'multiqc_plots' {
        path { _multiqc_plots -> { _file -> "multiqc/multiqc_plots" } }
    }
    'multiqc_report' {
        path { _multiqc_report -> { _file -> "multiqc/multiqc_report" } }
    }
    'rsem_index' {
        path { meta, _rsem_index ->
            { _file ->
                meta.source_version == "unknown"
                    ? "${meta.species}/${meta.source}/${meta.genome}/Sequence/RSEMIndex/version1.3.1"
                    : "${meta.species}/${meta.source}/${meta.genome}/Sequence/RSEMIndex/${meta.source_version}/version1.3.1"
            }
        }
    }
    'salmon_index' {
        path { meta, _salmon_index ->
            { _file ->
                meta.source_version == "unknown"
                    ? "${meta.species}/${meta.source}/${meta.genome}/Sequence/SalmonIndex/version1.10.3"
                    : "${meta.species}/${meta.source}/${meta.genome}/Sequence/SalmonIndex/${meta.source_version}/version1.10.3"
            }
        }
    }
    'splice_sites' {
        path { meta, _splice_sites -> { file -> "${meta.species}/${meta.source}/${meta.genome}/Sequence/SpliceSites/${file}" } }
    }
    'star_index' {
        path { meta, _star_index ->
            { _file ->
                meta.source_version == "unknown"
                    ? "${meta.species}/${meta.source}/${meta.genome}/Sequence/STARIndex/version2.7.11b"
                    : "${meta.species}/${meta.source}/${meta.genome}/Sequence/STARIndex/${meta.source_version}/version2.7.11b"
            }
        }
    }
    'transcript_fasta' {
        path { meta, _transcript_fasta -> { file -> "${meta.species}/${meta.source}/${meta.genome}/Sequence/TranscriptFasta/${file}" } }
    }
    'vcf_tbi' {
        path { meta, _vcf_tbi -> { file -> "${meta.species}/${meta.source}/${meta.genome}/Annotation/${meta.source_vcf}/${file}" } }
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
    bowtie1_index     = REFERENCES.out.bowtie1_index
    bowtie2_index     = REFERENCES.out.bowtie2_index
    bwamem1_index     = REFERENCES.out.bwamem1_index
    bwamem2_index     = REFERENCES.out.bwamem2_index
    dragmap_hashmap   = REFERENCES.out.dragmap_hashmap
    fasta             = REFERENCES.out.fasta
    fasta_dict        = REFERENCES.out.fasta_dict
    fasta_fai         = REFERENCES.out.fasta_fai
    fasta_sizes       = REFERENCES.out.fasta_sizes
    gtf               = REFERENCES.out.gtf
    hisat2_index      = REFERENCES.out.hisat2_index
    splice_sites      = REFERENCES.out.splice_sites
    intervals_bed     = REFERENCES.out.intervals_bed
    kallisto_index    = REFERENCES.out.kallisto_index
    msisensorpro_list = REFERENCES.out.msisensorpro_list
    rsem_index        = REFERENCES.out.rsem_index
    transcript_fasta  = REFERENCES.out.transcript_fasta
    salmon_index      = REFERENCES.out.salmon_index
    star_index        = REFERENCES.out.star_index
    vcf_tbi           = REFERENCES.out.vcf_tbi
    versions          = REFERENCES.out.versions
}
