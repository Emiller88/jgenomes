# nf-core/references: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of nf-core/references, created with the [nf-core](https://nf-co.re/) template.

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

### Fixed

- [19](https://github.com/nf-core/references/pull/19) - Use nf-core TEMPLATE
- [23](https://github.com/nf-core/references/pull/23) - No generation of bowtie2 index for sarek
- [30](https://github.com/nf-core/references/pull/30) - Deal with existing splice_sites
- [33](https://github.com/nf-core/references/pull/33) - Deal with existing faidx, sizes
- [36](https://github.com/nf-core/references/pull/36) - Deal intervals generation
- [39](https://github.com/nf-core/references/pull/39) - Fix gtf generation and dependencies
- [50](https://github.com/nf-core/references/pull/50) - Minimal JAVA is 17

### Dependencies

### Deprecated
