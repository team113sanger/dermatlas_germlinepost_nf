# dermatlas_germlinepost_nf

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.04.5-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

dermatlas_germlinepost_nf is a bioinformatics pipeline written in [Nextflow](http://www.nextflow.io) for generating and/or performing post-processing of germline variants generated with GATK on cohorts of tumors within the Dermatlas project. 

## Pipeline summary

In brief, the pipeline takes a set samples that have been pre-processed by the Dermatlas ingestion pipeline and then:
- Optional: Prepares samples for calling with GATK haplotype caller
- Generates a GenomicsDB datastore for joint calling germline variants
- Creates index files required by GATK for processing your genome of interest
- Generates per-chromosome variant call files for the cohort
- Merges those per-chrom files into a single cohort VCF and indexes it
- Selects, marks and filters SNPs and Indels as per GATK 
- Annotates the final variant sets with VEP
- Reformats and then summarises the data to produce germline oncoplots and tables.

## Inputs 

Inputs will depend on whether you are runnning in post-processing mode or end-to-end. Inputs can also be split into those which are cohort dependent and independent.

### Cohort-dependent variables
- `study_id`: prefix string to be applied to cohort-level summary files
- `outdir`: path to the where you would like the pipeline to output results
- `post_process_only`: logical determining whether to run post processing (VCF-> oncoplot) or the end-to-end (BAM -> oncoplot) germline analysis. 

**If true, the following inputs are required:**
- `geno_vcf`: a path to a set of .vcf files in a project directory. **Note: the pipeline assumes that corresponding index files have been pre-generated and are co-located with vcf and you should use a ** glob match to recursively collect all bamfiles in the directory**
- `sample_map`: path to a tab delimited file containing Sample IDs and the vcf files that they correspond to. Please see `tests/testdata/sample_map.tsv` for an example

**If false, the following inputs are required:**
- `tsv_file`: a manifest containing sample ids, associated bam files and their indexes. Please see `tests/testdata/manifest.tsv` for an example



### Cohort-independent variables
Reference files that are reused across pipeline executions have been placed within the pipeline's default `nextflow.config` file to simplify configuration. These can be ommited from setup. Behind the scences though, the following reference files are required for a run: 
- `chrom_list`: path to a text file containing ordered chrosome names See `assets/grch38_chromosome.txt`
- `reference_genome`: path to a reference genome file
- `baitset`: path to a `.bed` file describing the analysed genomic regions
- `vep_cache`: path to the release directory that contains a vep cache 
- `custom_files`: path to a set of annotation file to use in VEP (seperated by semi-colons)
- `custom_args`: path to a set of arguments to use with each custom file in VEP (seperated by semi-colons)
- `species`: VEP parameter, specifying the species being analysed (string)
- `db_name`: VEP parameter, specifying the ensembl data package version and corresponding db
- `assembly`: VEP parameter, specifying the reference genome build for the run.
- `summarise_results`: logical (whether to apply Dermatlas post processing into tables and figure)
**If true, the following inputs are required:**
- `nih_germline_resource`: path to file containing the information of the set of genes used by the NHS for [germline cancer predisposition diagnosis - prepared by mdc1@sanger.ac.uk](https://gitlab.internal.sanger.ac.uk/DERMATLAS/resources/national_genomic_test_germline_cancer_genes/-/tree/0.1.0?ref_type=tags)
- `cancer_gene_census_resoruce`: Cancer gene Census list of genes form COSMIC v97 
- `flag_genes`: path to a list of [FLAG](https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-014-0064-y#Sec11) genes, frequently mutated in normal exomes.
- `publish_intermediates`: logical (whether to publish large intermediate files (BAM and CRAMs)to the output directory)
- `alternative_transcripts`: path to a file containing Ensembl transcripts where we wish to modify the canonical transcript for accurate variant reporting.


Default values for reference files are supplied within the `nextflow.config` file and can be overided by adding them to the params `.json` file. An example complete params file `tests/test_data/test_params.json` is supplied within this repository for demonstation.

## Usage 

The recommended way to launch this pipeline on Sangers HPC is using a wrapper script (e.g. `bsub < my_wrapper.sh`) that submits nextflow as a job and records the version (**e.g.** `-r 0.1.1`)  and the `.json` parameter file supplied for a run.

An example wrapper script:
```
#!/bin/bash
#BSUB -q oversubscribed
#BSUB -G team113-grp
#BSUB -R "select[mem>8000] rusage[mem=8000] span[hosts=1]"
#BSUB -M 8000
#BSUB -oo logs/germline_variant_calling_%J.o
#BSUB -eo logs/germline_variant_calling_%J.e

PARAMS_FILE="/lustre/scratch125/casm/team113da/users/jb63/nf_germline_testing/params.json"

# Load module dependencies
module load nextflow-23.10.0
module load /software/modules/ISG/singularity/3.11.4

# Create a nextflow job that will spawn other jobs

nextflow run 'https://gitlab.internal.sanger.ac.uk/DERMATLAS/analysis-methods/dermatlas_germlinepost_nf' \
-r 0.3.0 \
-params-file $PARAMS_FILE \
-profile farm22 
```

When running the pipeline for the first time on the farm you will need to provide credentials to pull singularity containers from the team113 sanger gitlab. You should be able to do this by running
```
module load singularity/3.11.4 
singularity remote login --username $(whoami) docker://gitlab-registry.internal.sanger.ac.uk
```

The pipeline can configured to run on either Sanger OpenStack secure-lustre instances or farm22 by changing the profile speicified:
`-profile secure_lustre` or `-profile farm22`. 

## Pipeline visualisation 
Created using nextflow's in-built visualitation features.

```mermaid
flowchart TB
    subgraph " "
    v0["Channel.of"]
    v3["Channel.of"]
    v6["Channel.fromPath"]
    v12["ref_genome"]
    v14["Channel.fromPath"]
    v24["baitset"]
    v42["baitset"]
    v44["baitset"]
    v46["vep_cache"]
    v47["ref_genome"]
    v48["species"]
    v49["assembly_name"]
    v50["db_version"]
    v54["baitset"]
    v56["baitset"]
    v58["vep_cache"]
    v59["ref_genome"]
    v60["species"]
    v61["assembly_name"]
    v62["db_version"]
    v65["NIH_GERMLINE_TSV"]
    v66["CANCER_GENE_CENSUS"]
    v67["FLAG_GENES"]
    v71["NIH_GERMLINE_TSV"]
    v72["CANCER_GENE_CENSUS"]
    end
    v13([CREATE_DICT])
    subgraph GERMLINE
    subgraph NF_DEEPVARIANT
    v18([sort_cram])
    v19([markDuplicates])
    v21([coord_sort_cram])
    v22([bam_to_cram])
    v25([gatk_haplotypecaller])
    v15(( ))
    v26(( ))
    v27(( ))
    v29(( ))
    end
    end
    subgraph " "
    v20[" "]
    v23[" "]
    v28[" "]
    v69[" "]
    v70[" "]
    v74[" "]
    v75[" "]
    end
    v32([GENERATE_GENOMICS_DB])
    v33([GATK_GVCF_PER_CHROM])
    v38([MERGE_COHORT_VCF])
    v39([INDEX_COHORT_VCF])
    subgraph PROCESS_SNPS
    v41([SELECT_VARIANTS])
    v43([MARK_VARIANTS])
    v45([FILTER_VARIANTS])
    v51([ANNOTATE_VARIANTS])
    v52([CONVERT_TO_TSV])
    v1(( ))
    v4(( ))
    end
    subgraph PROCESS_INDELS
    v53([SELECT_VARIANTS])
    v55([MARK_VARIANTS])
    v57([FILTER_VARIANTS])
    v63([ANNOTATE_VARIANTS])
    v64([CONVERT_TO_TSV])
    end
    v68([COMBINED_SUMMARY])
    v73([CONVERT_TO_MAF])
    v7(( ))
    v34(( ))
    v40(( ))
    v0 --> v1
    v3 --> v4
    v6 --> v7
    v12 --> v13
    v13 --> v19
    v13 --> v22
    v13 --> v25
    v13 --> v33
    v14 --> v15
    v15 --> v18
    v18 --> v19
    v19 --> v21
    v19 --> v20
    v21 --> v22
    v21 --> v25
    v22 --> v23
    v24 --> v25
    v25 --> v26
    v25 --> v27
    v25 --> v29
    v27 --> v28
    v7 --> v32
    v26 --> v32
    v29 --> v32
    v32 --> v33
    v7 --> v33
    v33 --> v34
    v34 --> v38
    v38 --> v39
    v39 --> v40
    v40 --> v41
    v41 --> v43
    v42 --> v43
    v43 --> v45
    v44 --> v45
    v45 --> v51
    v46 --> v51
    v47 --> v51
    v48 --> v51
    v49 --> v51
    v50 --> v51
    v1 --> v51
    v4 --> v51
    v51 --> v52
    v52 --> v68
    v40 --> v53
    v53 --> v55
    v54 --> v55
    v55 --> v57
    v56 --> v57
    v57 --> v63
    v58 --> v63
    v59 --> v63
    v60 --> v63
    v61 --> v63
    v62 --> v63
    v1 --> v63
    v4 --> v63
    v63 --> v64
    v64 --> v68
    v65 --> v68
    v66 --> v68
    v67 --> v68
    v68 --> v73
    v68 --> v70
    v68 --> v69
    v71 --> v73
    v72 --> v73
    v73 --> v75
    v73 --> v74
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


