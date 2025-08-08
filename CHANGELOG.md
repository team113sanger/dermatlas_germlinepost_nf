# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [0.3.3] - 2025-08-08
### Changed 
- Quality of life improvements. Correction of file names and locations to better mirror the manual pipeline.
### Added
- Added template assets for fetching by Dermanager.

## [0.3.2] - 2025-07-15
### Changed 
- Altering paths and defaults for new Dermatlas resources dir 

## [0.3.1] - 2024-05-08
### Added 
- Fixed publishing paths for post-processing steps to be consistent with old manual pipeline

## [0.3.0] - 2024-04-01
### Added 
- Refactor to use the new Germline post-processing steps (unified MAF generation with somatic variant pipeline, updated oncoplots).
- Make publication of intermediate files optional with `publish_intermediates` parameter.

### Fixed
- Patching an issue with re-entrancy that can occur when failing after VCF chrom spitting 

## [0.2.7] - 2024-02-12
### Fixed
- Patching an issue with post-process-only set to true where "input name collision" was declared. May related to double declaration of `vcf_ch` in `post_process_only.nf`.
### Added
- Improved testing coverage and end-to-end pipeline test

## [0.2.6] - 2024-01-08
### Changed
- Incorporate GERMLINE 0.5.1 container as default and upstream change in Cosmic filtering should now be reflected here

## [0.2.5] - 2024-01-07
### Changed
- Updated default order of arguments for VEP. 

## [0.2.4] - 2024-11-18
### Changed
- Removed `--no_stats` flag from vep runs in order to prevent an observed VEP bug caused
by mutliallelic variants (see [here](https://github.com/Ensembl/ensembl-vep/issues/1013) and [here](https://github.com/Ensembl/ensembl-vep/issues/818))

## [0.2.3] - 2024-11-08
### Changed
- Updated the output path for the oncoplotting and summary steps (summ_tabs)

## [0.2.2] - 2024-11-06
### Added
- Some updates to docs and config files for better running on Farm for dermatlas

## [0.2.1] - 2024-10-30
### Added
- Add `--protein` flag to Dermatlas and FUR Vep annotation process.

## [0.2.0] - 2024-10-02
### Added
- Workflow steps for generating variant calls from Dermatlas bam files. Triggerable via post-process only = False. New parameters to support these additional steps.

### Changed 
- The way that custom annotation files are specified to vep so that the process can be made generic. Now use a single custom files param (with ; seperation) and a custom args flag rather than specifying in the vep config, which was opinionated. 

## [0.1.1] - 2024-09-27
### Added
- Configuration edits required for Farm22. Modifying resource allocation, cpus and publish dirs 
### Fixed
- Fix an issue where sample map path would cause failures 

## [0.1.0] - 2024-09-25
- Initial release of the dermatlas germline post-processing pipeline for user-testing. Tested on farm22 up until tsv conversion
