# Brainstorming

## Generate

- md5 checksums (validate downloads if possible)

## Track within the pipeline

- software_versions
- copy of command.sh (or just save Nextflow report?)
- Asset input paths
- Show skipped reference types if already existed
- Allow appending to the readme (treat like changelog), in case new asset types added

## Strategy

When adding a new asset, build for the latest reference versions only. Do all genomes.
Optionally backfill old releases on demand if specifically triggered.
