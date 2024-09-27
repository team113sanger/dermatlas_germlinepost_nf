# dermatlas_germlinepost_nf

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.04.5-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

dermatlas_germlinepost_nf is a bioinformatics pipeline written in [Nextflow](http://www.nextflow.io) for generating and/or performing post-processing of germline variants generated with GATK on cohorts of tumors within the Dermatlas project. 

## Pipeline summary

In brief, the pipeline takes a set samples that have been pre-processed by the Dermatlas ingestion pipeline and then:
- Generates a GenomicsDB datastore for joint calling germline variants
- Creates index files required by GATK for processing your genome of interest
- Generates per-chromosome variant call files for the cohort
- Merges those per-chrom files into a single cohort VCF and indexes it
- Selects, marks and filters SNPs and Indels as per GATK 
- Annotates the final variant sets with VEP
- Reformats and then summarises the data to produce germline oncoplots and tables.

## Inputs 

Inputs will depend on whether you are runnning in post-processing mode or end-to-end. 
Only post-processing currently supported. Inputs can be split into those which are cohort dependent and independent.

### Cohort-dependent variables
`study_id`: Prefix number for the cohort
`geno_vcf`: a path to a set of .vcf files in a project directory. Note: the pipeline assumes that corresponding index files have been pre-generated and are co-located with vcf and you should use a ** glob match to recursively collect all bamfiles in the directory
`sample_map`: path to a tab delimited file containing Sample IDs and the vcf files that they correspond to 
`outdir`: path to the where you would like the pipeline to output results

### Cohort-independent variables
Reference files that are reused across pipeline executions have been placed within the pipeline's default `nextflow.config` file to simplify configuration. These can be ommited from setup. Behind the scences though, the following reference files are required for a run: 

- `reference_genome`: path to a reference genome file
- `baitset`: path to a `.bed` file describing the analysed genomic regions
- `vep_cache`: path to the release directory that contains a vep cache 
- `vep_config`: path to a file containing vep options. See `assets/vep_config.ini`
- `gnomad_file`: path to a gnomad annotation file to use in VEP 
- `dbsnp_file`: path to a DBSNP annotation file to use in VEP 
- `clinvar_file`: path to a DBSNP annotation file to use in VEP 
- `cosmic_file`: path to a COSMIC annotation file to use in VEP  
- `nih_germline_resource`: path to file containing the information of the set of genes used by the NHS for [germline cancer predisposition diagnosis - prepared by mdc1@sanger.ac.uk](https://gitlab.internal.sanger.ac.uk/DERMATLAS/resources/national_genomic_test_germline_cancer_genes/-/tree/0.1.0?ref_type=tags)
- `cancer_gene_census_resoruce`: Cancer gene Census list of genes form COSMIC v97 
- `flag_genes`: path to a list of [FLAG](https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-014-0064-y#Sec11) genes, frequently mutated in normal exomes.

Default reference file values supplied within the `nextflow.config` file can be overided by adding them to the params `.json` file. An example complete params file `tests/test_data/test_params.json` is supplied within this repo for demonstation.

## Usage 

The recommended way to launch this pipeline is using a wrapper script (e.g. `bsub < my_wrapper.sh`) that submits nextflow as a job and records the version (**e.g.** `-r 0.1.0`)  and the `.json` parameter file supplied for a run.

An example wrapper script:
```
#!/bin/bash
#BSUB -q normal
#BSUB -G team113
#BSUB -R "select[mem>8000] rusage[mem=8000] span[hosts=1]"
#BSUB -M 8000
#BSUB -oo nf_out.o
#BSUB -eo nf_out.e

PARAMS_FILE="/lustre/scratch125/casm/team113da/users/jb63/nf_germline_testing/params.json"

# Load module dependencies
module load nextflow-23.10.0
module load /software/modules/ISG/singularity/3.11.4
module load /software/team113/modules/modulefiles/tw/0.6.2

# Create a nextflow job that will spawn other jobs

nextflow run 'https://gitlab.internal.sanger.ac.uk/DERMATLAS/analysis-methods/dermatlas_germlinepost_nf' \
-r 0.1.0 \
-params-file $PARAMS_FILE \
-c nextflow.config \
-profile farm22 
```


When running the pipeline for the first time on the farm you will need to provide credentials to pull singularity containers from the team113 sanger gitlab. These should be provided as environment variables:
`SINGULARITY_DOCKER_USERNAME`=userid@sanger.ac.uk
`SINGULARITY_DOCKER_PASSWORD`=YOUR_GITLAB_LOGIN_PASSWORD

You can fix these variables to load by default by adding the following lines to your `~/.bashrc` file
```
export SINGULARITY_DOCKER_USERNAME=userid@sanger.ac.uk
export SINGULARITY_DOCKER_PASSWORD=YOUR_GITLAB_LOGIN_PASSWORD
```

The pipeline can configured to run on either Sanger OpenStack secure-lustre instances or farm22 by changing the profile speicified:
`-profile secure_lustre` or `-profile farm22`. 

## Pipeline visualisation 
Created using nextflow's in-built visualitation features.

```mermaid
flowchart TB
    
```

## Testing

This pipeline has been developed with the [nf-test](http://nf-test.com) testing framework. Unit tests and small test data are provided within the pipeline `test` subdirectory. A snapshot has been taken of the outputs of most steps in the pipeline to help detect regressions when editing. You can run all tests on openstack with:

```
nf-test test 
```
and individual tests with:
```
nf-test test tests/modules/ascat_exomes.nf.test
```

For faster testing of the flow of data through the pipeline **without running any of the tools involved**, stubs have been provided to mock the results of each succesful step.
```
nextflow run main.nf \
-params-file params.json \
-c tests/nextflow.config \
--stub-run
```


