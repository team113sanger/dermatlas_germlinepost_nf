# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Keywords

As of version 0.4.0 the following *keywords* are used at the start of each
changelog entry to indicate the impact of the change:

- **REPRODUCIBILITY** - a change to the pipeline's scientific processing that
  may cause the same input data to produce different scientific outputs or
  results, including changes to algorithms, tolerances, randomisation,
  scientific functionality, or output formats.
- **ROBUSTNESS** - a fix or improvement to the pipeline's scientific
  functionality that improves correctness, reliability, or the range of inputs
  that can be processed, without intentionally changing the scientific results
  of an equivalent successful analysis.
- **INTEGRATION** - a change to how the pipeline integrates with other systems
  or infrastructure, without changing its scientific processing or results.

## [Unreleased]

## [0.4.0] - 2026-09-29
### Added
- **INTEGRATION** - run reporting. `workflow.onComplete` calls `Utils.reportRun`
  (`lib/Utils.groovy`, shared byte-for-byte with the other Dermatlas pipelines), which
  records the run in the Dermatlas website's analysis log via `dermatlas-http cohort
  analysis-log` (>= 0.6.1) and posts a Slack message. Each is opt-in/opt-out through
  `DERMATLAS_WEBSITE_LOGGING` / `DERMATLAS_SLACK_NOTIFICATIONS`, reads its values
  (`COHORT_SLUG`, `SAMPLE_LIST_VERSION_FILE`, `SELF_DESCRIBING_API`,
  `SLACK_WEBHOOK_URL`) from the environment, never fires on a stub run, and never
  changes the pipeline's exit status. New params: `analysis_pipeline_slug`
  (`germline_pipe`), `is_stub`, `trace_file`.
- **INTEGRATION** - `run_germline.sh` is rebuilt from the Dermatlas launcher template
  (`dermatlas_rnafusions_nf` 0.4.15): it sources `source_me.sh` (`SOURCE_ME`, or
  `"none"` plus the MANUAL ENVIRONMENT OVERRIDES block for git-clone runs), checks every
  variable the config needs before `nextflow run` starts, reports a failed launch to
  stderr and (opted in) Slack, holds an exclusive `flock` on
  `${PROJECT_DIR}/germline_pipe/.lock` for the life of the run, and writes
  `.completed_successfully` / `.completed_with_error` at exit. See "Reclaiming disk
  space" in the README.
- **INTEGRATION** - after a successful run the launcher writes
  `stats/resource-stats-<RUN_ID>.txt` (wall time, work-dir bytes and inodes), reports
  it to the website with `dermatlas-http cohort analysis-workdir-stats` (module-loaded
  via `DERMATLAS_HTTP_MODULE`, default `dermatlas-http`; requires dermatlas-web-client
  >= 0.6.2; best effort), and deletes the run's work directory unless
  `DERMATLAS_CLEANUP_WORK_DIR=false`. A failed or killed run always keeps it.
- **INTEGRATION** - one launcher-owned `RUN_ID` names the run's nextflow logs
  (`logs/nextflow-run-<RUN_ID>.log`, `logs/nextflow-pull-<RUN_ID>.log`), execution trace
  (new) and execution report, both under `${PROJECT_DIR}/germline_pipe/traces/`.
  `nextflow run` uses a per-revision clone (`clones/<REVISION>`) and a pinned
  singularity cache.
- **INTEGRATION** - `.update-version.sh` sets the version in every file that records it
  (`assets/run_germline.sh`, `docs/source/conf.py`, `nextflow.config`); see "Cutting a
  release" in the README.

### Changed
- **INTEGRATION** - **Breaking:** `germline_variants.config` takes `tsv_file` from
  `DNA_GERMLINE_NORMAL_MANIFEST`, the normal manifest dermanager now generates (replacing
  `germline_normal_select.R`), and `outdir` from `${ANALYSIS_DIR}/germline`, and
  `run_germline.sh` requires both. A `source_me.sh` that does not export
  `DNA_GERMLINE_NORMAL_MANIFEST` fails the launch with the variable named.
- **INTEGRATION** - **Breaking:** the launcher reads its config from
  `commands/germline_pipe/germline_variants.config` (was `commands/germline_variants.config`)
  and runs in `${PROJECT_DIR}/germline_pipe` (was `germline_pipeline`), matching
  dermanager's slug. A run started under the old directory cannot be `-resume`d from the
  new one.
- **INTEGRATION** - **Breaking:** a run killed by `bkill` or an LSF limit now exits
  `128+n` and records `.completed_with_error`; a second submission while a run holds the
  pipeline directory exits 75 without touching it.
- **INTEGRATION** - the post-run `chmod -R ug+rw` of the results now runs inside the
  launcher's exit handler on success only, and its failure is a `NOTE:`, not a failed job.
- **INTEGRATION** - `.github/workflows/publish-assets.yml` is replaced with the reference
  copy: the rolling `main-latest` / `develop-latest` tags are created once and never
  moved (only the attached bundle is replaced), which kept breaking `git hf release
  finish`. The repository is GitHub-primary; README, docs and `manifest.homePage` no
  longer point at GitLab.

### Fixed
- **REPRODUCIBILITY** - the `farm22` profile's `reference_genome` now matches the asset
  config (`references/germline/genome.fa`). Managed runs already used this genome; a
  run relying on the profile alone previously used
  `resources/ascat/GRCh38_full_analysis_set_plus_decoy_hla.fa`.
- **INTEGRATION** - the execution report is written beside the trace instead of the
  launch directory, and the dead top-level `tracedir` setting is removed.
- **INTEGRATION** - `docs/source/conf.py` carried release 0.3.3; the user docs no longer
  embed a stale copy of the launcher or config.

## [0.3.6] - 2026-08-27
### Added
- `.github/workflows/publish-assets.yml` publishes `assets/` to GitHub Releases as
  `projectify_asset_bundle.tar.gz` (and a `.sha256` of it) on every push to `main` and
  `develop` - as the rolling `main-latest` and `develop-latest` pre-releases - and on
  every `X.Y.Z` tag. `dermanager projectify` fetches assets from those release URLs
  instead of the GitHub API, which needs no token and is not rate limited. See
  "Asset release bundles" in the README.

## [0.3.5] - 2026-03-27
### Fixed
- Fixed VEP annotation parameters in `vep_annotation.nf` for `homo_sapiens` to ensure match with the ones used for somatic calling pipeline. 

## [0.3.4] - 2025-10-06
### Changed 
- Moved to documentation via sphinx and bundled with repo rather than confleunce

### Added
- Asset updates for multi-pipline running setup as default
- Documentation updates to reflect multi-pipeline


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
