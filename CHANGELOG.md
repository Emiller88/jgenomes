# nf-core/references: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## dev

### Added

- [72](https://github.com/nf-core/references/pull/72) - Compress vcf files when they are not already compressed
- [74](https://github.com/nf-core/references/pull/74) - Add DOI
- [76](https://github.com/nf-core/references/pull/76) - Add untar and unzip

### Changed

- [73](https://github.com/nf-core/references/pull/73) - Update modules
- [74](https://github.com/nf-core/references/pull/74) - Back to dev
- [75](https://github.com/nf-core/references/pull/75) - Template update for nf-core/tools v3.1.1
- [78](https://github.com/nf-core/references/pull/78) - Rename input to asset
- [79](https://github.com/nf-core/references/pull/79) - Improve VCF handling and refactor schema for VCFs
- [80](https://github.com/nf-core/references/pull/80) - Improve VCF publishing in index
- [86](https://github.com/nf-core/references/pull/86) - Template update for nf-core/tools v3.2.0
- [87](https://github.com/nf-core/references/pull/87) - Automatic nf-test sharding

### Fixed

- [77](https://github.com/nf-core/references/pull/77) - Fix createsequencedictionary usage
- [81](https://github.com/nf-core/references/pull/81) - Fix fasta_fai unnecessary regeneration

### Dependencies

| modules                           | old version | new version |
| --------------------------------- | ----------- | ----------- |
| gunzip                            | 1.1         | 1.13        |
| samtools (in star/genomegenerate) | 1.20        | 1.21        |
| untar                             |             | 1.34        |
| unzip                             |             | 16.02       |

> NB: Dependency has been updated if both old and new version information is present.
> NB: Dependency has been added if just the new version information is present.
> NB: Dependency has been removed if new version information isn't present.

### Subworkflows

### Deprecated

## [0.1](https://github.com/nf-core/references/releases/tag/0.1) - Tar Tarasque

Initial pre-release of nf-core/references, created with the [nf-core](https://nf-co.re/) template.
Tar is a dark grey color ( #383838), and the Tarasque is a legendary dragon from the South of France.

### Added

- [5](https://github.com/nf-core/references/pull/5) - Brainstorm about input
- [14](https://github.com/nf-core/references/pull/14) - nf-core/rnaseq references
- [17](https://github.com/nf-core/references/pull/17) - nf-core/sarek references
- [22](https://github.com/nf-core/references/pull/22) - Test pipeline output with nf-test for sarek
- [23](https://github.com/nf-core/references/pull/23) - Tools selection
- [24](https://github.com/nf-core/references/pull/24) - Test pipeline output with nf-test for rnaseq
- [27](https://github.com/nf-core/references/pull/27) - Add faidx, gffread, sizes generation
- [28](https://github.com/nf-core/references/pull/28) - Add hisat2 generation
- [29](https://github.com/nf-core/references/pull/29) - Add hisat2_extract_sites generation
- [31](https://github.com/nf-core/references/pull/31) - Add rsem, rsem_make_transcript_fasta generation
- [32](https://github.com/nf-core/references/pull/32) - Add kallisto, salmon generation
- [37](https://github.com/nf-core/references/pull/37) - Add tabix tbi generation
- [42](https://github.com/nf-core/references/pull/42) - Add multiple references build tests
- [43](https://github.com/nf-core/references/pull/43) - Add support for multiple known_indels_vcf via s3 globs
- [47](https://github.com/nf-core/references/pull/47) - Add fasta assets for files in igenomes
- [51](https://github.com/nf-core/references/pull/51) - Add fasta_fai assets for files in igenomes
- [52](https://github.com/nf-core/references/pull/52) - Add abundantsequences_fasta, bismark_index, bowtie1_index, bowtie2_index, bwamem1_index, bwamem2_index, chrom_info, chromosomes_fasta, dragmap_hashtable, fasta_dict, genes_bed, genes_refflat, genes_refgene, genome_size_xml, gtf, hairpin_fasta, mature_fasta, readme, source, source_vcf, species, star_index and vcf assets for files in igenomes
- [56](https://github.com/nf-core/references/pull/56) - Add fields for bowtie1_index, bowtie2_index, bwamem1_index, bwamem2_index, dragmap_hashtable, hisat2_index, kallisto_index, msisensorpro_list, rsem_index, salmon_index, star_index, vcf_tbi in assets
- [56](https://github.com/nf-core/references/pull/56) - Add new params: kallisto_make_unique to use the --make-unique option for kallisto
- [56](https://github.com/nf-core/references/pull/56) - New file assets/genomes/Caenorhabditis_elegans/NCBI/WBcel235_updated.yml, build from assets/genomes/Caenorhabditis_elegans/NCBI/WBcel235.yml
- [62](https://github.com/nf-core/references/pull/62) - Added comments to the code
- [66](https://github.com/nf-core/references/pull/66) - Output index

### Changed

- [21](https://github.com/nf-core/references/pull/21) - Template update for nf-core/tools v3.0.2
- [23](https://github.com/nf-core/references/pull/23) - Merge all scripts into one
- [31](https://github.com/nf-core/references/pull/31) - Test refactor: hisat2 and rsem have their own tests
- [33](https://github.com/nf-core/references/pull/33) - default.yml asset file is now in the test-dataset repo
- [35](https://github.com/nf-core/references/pull/35) - Use new output system
- [35](https://github.com/nf-core/references/pull/35) - Samtools + intervals are separate from Sarek tests now
- [35](https://github.com/nf-core/references/pull/35) - Rename bed_intervals to intervals_bed
- [35](https://github.com/nf-core/references/pull/35) - Rename rnaseq tests to tools tests
- [41](https://github.com/nf-core/references/pull/41) - Better sarek tests
- [41](https://github.com/nf-core/references/pull/41) - Better publishing for sarek related files
- [43](https://github.com/nf-core/references/pull/43) - Fasta is no longer a required asset
- [48](https://github.com/nf-core/references/pull/48) - Simplify VCF tabix index generation and related assets
- [48](https://github.com/nf-core/references/pull/48) - Code refactoring (new subworfklows for each type of operations)
- [49](https://github.com/nf-core/references/pull/49) - Better publishing for all files
- [53](https://github.com/nf-core/references/pull/53) - Better publishing for all aligner indexes
- [56](https://github.com/nf-core/references/pull/56) - reference_version -> source_version
- [59](https://github.com/nf-core/references/pull/59) - Simplify input structure
- [60](https://github.com/nf-core/references/pull/60) - Just 2 shards in CI instead of 5
- [61](https://github.com/nf-core/references/pull/61) - Move assets in nf-core/references-assets
- [62](https://github.com/nf-core/references/pull/62) - samplesheet -> asset
- [62](https://github.com/nf-core/references/pull/62) - Refactor and simplify codebase
- [63](https://github.com/nf-core/references/pull/63) - Unpack gff even when gtf is present
- [64](https://github.com/nf-core/references/pull/64) - Improve documentation
- [68](https://github.com/nf-core/references/pull/68) - Minor refactoring
- [69](https://github.com/nf-core/references/pull/69) - Better comments
- [70](https://github.com/nf-core/references/pull/70) - Prepare release 0.1

### Fixed

- [19](https://github.com/nf-core/references/pull/19) - Use nf-core TEMPLATE
- [23](https://github.com/nf-core/references/pull/23) - No generation of bowtie2 index for sarek
- [30](https://github.com/nf-core/references/pull/30) - Deal with existing splice_sites
- [33](https://github.com/nf-core/references/pull/33) - Deal with existing faidx, sizes
- [36](https://github.com/nf-core/references/pull/36) - Deal intervals generation
- [39](https://github.com/nf-core/references/pull/39) - Fix gtf generation and dependencies
- [50](https://github.com/nf-core/references/pull/50) - Minimal JAVA is 17
- [51](https://github.com/nf-core/references/pull/51) - Fix missing fasta assets for GATK build
- [56](https://github.com/nf-core/references/pull/56) - Add new logic for skip creation of existing assets
- [62](https://github.com/nf-core/references/pull/62) - Remove failure when no tools are specified, because one might just want to unpack assets

### Dependencies

### Subworkflows

| old name                    | new name                         |
| --------------------------- | -------------------------------- |
| SAMPLESHEET_TO_CHANNEL      | ASSET_TO_CHANNEL                 |
| CREATE_ALIGN_INDEX_WITH_GFF | CREATE_FROM_FASTA_AND_ANNOTATION |
| CREATE_ALIGN_INDEX          | CREATE_FROM_FASTA_ONLY           |
| INDEX_FASTA                 | CREATE_FROM_FASTA_ONLY           |
| UNCOMPRESS_REFERENCES       | UNCOMPRESS_ASSET                 |

> NB: Subworkflow has been updated if both old and new name is present.
> NB: Subworkflow has been added if just the new name is present.
> NB: Subworkflow has been removed if new name isn't present.
> NB: Subworkflows have been merged if several subworkflows have the same new name.

### Deprecated
