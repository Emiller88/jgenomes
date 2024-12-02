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
        log.error("No tools specified")
        error("EXIT: No tools specified")
    }

    NFCORE_REFERENCES(PIPELINE_INITIALISATION.out.samplesheet, params.tools)

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
    'bowtie1' {
        path 'bowtie1'
    }
    'bowtie2' {
        path 'bowtie2'
    }
    'bwamem1' {
        path { meta, index -> { file -> "${meta.species}/${meta.source}/${meta.id}/Sequence/BWAIndex/" } }
    }
    'bwamem2' {
        path { meta, index -> { file -> "${meta.species}/${meta.source}/${meta.id}/Sequence/BWAmem2Index/" } }
    }
    'dragmap' {
        path { meta, index -> { file -> "${meta.species}/${meta.source}/${meta.id}/Sequence/dragmap/" } }
    }
    'fasta' {
        path { meta, fasta -> { file -> "${meta.species}/${meta.source}/${meta.id}/Sequence/WholeGenomeFasta/${file}" } }
    }
    'fasta_dict' {
        path { meta, dict -> { file -> "${meta.species}/${meta.source}/${meta.id}/Sequence/WholeGenomeFasta/${file}" } }
    }
    'fasta_fai' {
        path { meta, fai -> { file -> "${meta.species}/${meta.source}/${meta.id}/Sequence/WholeGenomeFasta/${file}" } }
    }
    'fasta_sizes' {
        path 'fasta_sizes'
    }
    'gffread' {
        path 'gffread'
    }
    'hisat2' {
        path 'hisat2'
    }
    'intervals' {
        path { meta, intervals -> { file -> "${meta.species}/${meta.source}/${meta.id}/Annotation/intervals/${file}" } }
    }
    'kallisto' {
        path 'kallisto'
    }
    'make' {
        path 'make'
    }
    'msisensorpro' {
        path { meta, index -> { file -> "${meta.species}/${meta.source}/${meta.id}/Annotation/msisensorpro/${file}" } }
    }
    'multiqc_data' {
        path 'multiqc'
    }
    'multiqc_plots' {
        path 'multiqc'
    }
    'multiqc_report' {
        path 'multiqc'
    }
    'rsem' {
        path 'rsem'
    }
    'salmon' {
        path 'salmon'
    }
    'star' {
        path 'star'
    }
    'tabix_vcf_tbi' {
        path { meta, tbi -> { file -> "${meta.species}/${meta.source}/${meta.id}/Annotation/${meta.source_vcf}/${file}" } }
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
    bowtie1               = REFERENCES.out.bowtie1
    bowtie2               = REFERENCES.out.bowtie2
    bwamem1               = REFERENCES.out.bwamem1
    bwamem2               = REFERENCES.out.bwamem2
    dragmap               = REFERENCES.out.dragmap
    fasta                 = REFERENCES.out.fasta
    fasta_dict            = REFERENCES.out.fasta_dict
    fasta_fai             = REFERENCES.out.fasta_fai
    gffread               = REFERENCES.out.gff_gtf
    hisat2                = REFERENCES.out.hisat2
    hisat2_splice_sites   = REFERENCES.out.hisat2_splice_sites
    intervals             = REFERENCES.out.intervals_bed
    kallisto              = REFERENCES.out.kallisto
    msisensorpro          = REFERENCES.out.msisensorpro
    rsem                  = REFERENCES.out.rsem
    rsem_transcript_fasta = REFERENCES.out.rsem_transcript_fasta
    salmon                = REFERENCES.out.salmon
    sizes                 = REFERENCES.out.sizes
    star                  = REFERENCES.out.star
    vcf_tbi               = REFERENCES.out.vcf_tbi
    versions              = REFERENCES.out.versions
}
