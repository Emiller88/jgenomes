<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-references_logo_dark.png">
    <img alt="nf-core/references" src="docs/images/nf-core-references_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/references/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/references/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/references/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/references/actions/workflows/linting.yml)
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/references/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.14576225-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.14576225)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/references)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23references-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/references)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)
[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/references** is a bioinformatics pipeline that build references, for multiple use cases.

It is primarily designed to build references for common organisms and store it on [AWS iGenomes](https://github.com/ewels/AWS-iGenomes/).

From a fasta file, it will be able to build the following references:

- Bowtie1 index
- Bowtie2 index
- BWA-MEM index
- BWA-MEM2 index
- DRAGMAP hashtable
- Fasta dictionary (with GATK4)
- Fasta fai (with SAMtools)
- Fasta sizes (with SAMtools)
- Fasta intervals bed (with GATK4)
- MSIsensor-pro list

With an additional annotation file describing the genes (either GFF3 or GTF), it will be able to build the following references:

- GTF (from GFF3 with GFFREAD)
- HISAT2 index
- Kallisto index
- RSEM index
- Salmon index
- Splice sites (with HISAT2)
- STAR index
- Transcript fasta (with RSEM)

With a vcf file, it will compress it, if it was not already compressed, and tabix index it.

## Assets

Assets are stored in [references-assets](https://github.com/nf-core/references-assets).

## Running

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow.Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

`asset.yml`:

```yml
- genome: GRCh38_chr21
  fasta: "https://raw.githubusercontent.com/nf-core/test-datasets/references/references/GRCh38_chr21/GRCh38_chr21.fa"
  gtf: "https://raw.githubusercontent.com/nf-core/test-datasets/references/references/GRCh38_chr21/GRCh38_chr21.gtf"
  source_version: "CUSTOM"
  readme: "https://raw.githubusercontent.com/nf-core/test-datasets/references/references/GRCh38_chr21/README.md"
  source: "nf-core/references"
  source_vcf: "GATK_BUNDLE"
  species: "Homo_sapiens"
  vcf: "https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/vcf/dbsnp_146.hg38.vcf.gz"
```

Each line represents a source for building a reference, a reference already built or a metadata.

Now, you can run the pipeline using:

```bash
nextflow run nf-core/references \
   -profile <docker/singularity/.../institute> \
   --input asset.yml \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/references/usage) and the [parameter documentation](https://nf-co.re/references/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/references/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/references/output).

## Credits

nf-core/references was originally written by [Maxime U Garcia](https://github.com/maxulysse) | [Edmund Miller](https://github.com/edmundmiller) | [Phil Ewels](https://github.com/ewels).

We thank the following people for their extensive assistance in the development of this pipeline:

- [Adam Talbot](https://github.com/adamrtalbot)
- [Friederike Hanssen](https://github.com/FriederikeHanssen)
- [Harshil Patel](https://github.com/drpatelh)
- [Jonathan Manning](https://github.com/pinin4fjords)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#references` channel](https://nfcore.slack.com/channels/references) (you can join with [this invite](https://nf-co.re/join/slack)).

### How to hack on it

0. Have docker, and Nextflow installed
1. `nextflow run main.nf`

### Some thoughts on reference building

- We could use the glob and if you just drop a fasta in s3 bucket it'll get picked up and new resources built
  - Could take this a step further and make it a little config file that has the fasta, gtf, genome_size etc.
- How do we avoid rebuilding? Ideally we should build once on a new minor release of an aligner/reference. IMO kinda low priority because the main cost is going to be egress, not compute.
- How much effort is too much effort?
  - Should it be as easy as adding a file on s3?
    - No that shouldn't be a requirement, should be able to link to a reference externally (A "source of truth" ie an FTP link), and the workflow will build the references
    - So like mulled biocontainers, just make a PR to the samplesheet and boom new reference in the s3 bucket if it's approved?

### Roadmap

PoC for v1.0:

- Replace aws-igenomes
  - bwa, bowtie2, star, bismark need to be built
  - fasta, gtf, bed12, mito_name, macs_gsize, blacklist, copied over

Other nice things to have:

- Building our test-datasets
- Down-sampling for a unified genomics test dataset creation, (Thinking about viralitegration/rnaseq/wgs) and spiking in test cases of interest (Specific variants for example)

## Citations

If you use nf-core/references for your analysis, please cite it using the following doi: [10.5281/zenodo.14576225](https://doi.org/10.5281/zenodo.14576225)

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
