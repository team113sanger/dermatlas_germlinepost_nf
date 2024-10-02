# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2024-10-02
### Added
- Workflow steps for generating variant calls from dermatlas bams. Triggerable via post-process only = False

### Changed 
- The way that custom annotation files are specified to vep so that the process can be made generic. Now use a single custom files param (with ; seperation) and a custom args flag rather than specifying in the vep config, which was opinionated. 

## [0.1.1] - 2024-09-27
### Added
- Configuration edits required for Farm22. Modifying resource allocation, cpus and publish dirs 
### Fixed
- Fix an issue where sample map path would cause failures 

## [0.1.0] - 2024-09-25
- Initial release of the dermatlas germline post-processing pipeline for user-testing. Tested on farm22 up until tsv conversion
