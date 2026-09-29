# Nextflow: Germline variant calling pipeline

Germline variant calling and post-processing for DERMATLAS can be run mostly with a single nextflow pipeline in a largely "set-and-forget" manner to reproduce the manual steps detailed in [DERMATLAS - Germline calling with GATK - for WES - using Nextflow Tower](https://confluence.sanger.ac.uk/x/BJOeB). This document contains an SOP for configuring and running the pipeline. For a more detailed explanation of the pipeline, the inputs, steps and requirements for running can be found within the pipeline project [README](https://github.com/team113sanger/dermatlas_germlinepost_nf/blob/develop/README.md)

## Workflow Overview

1. Generate the normal input table
2. Generating the cohort config file
3. Running the pipeline
   - i) From BAMS
   - ii) From VCFs
4. Make a release folder
5. Cleanup the intermediate BAM files created by the pipeline

## Workflow Steps

### 1. Generate Input Table

The pipeline runs on a `.tsv` file detailing the normal samples to run, with a header and the following columns:

| sample   | object                                                                                                                                   | object_index                                                                                                                                 |
|:---------|:-----------------------------------------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------|
| PD42171b | /lustre/scratch124/casm/team113/projects/5534_Landscape_sebaceous_tumours_GRCh38_Remap_germline/BAMS/PD42171b.sample.dupmarked.bam | /lustre/scratch124/casm/team113/projects/5534_Landscape_sebaceous_tumours_GRCh38_Remap_germline/BAMS/PD42171b.sample.dupmarked.bam.bai |

Provided that you have set up your project with [dermanager](https://confluence.sanger.ac.uk/x/fwBTCQ), this table is generated for you from the cohort's retained, usable matched tumour-normal pairs (one normal per patient, pointed at the staged BAM), and `source_me.sh` exports its path as `DNA_GERMLINE_NORMAL_MANIFEST`:

```
metadata/{study_id}_{canapps_id}-normal_one_per_patient_matched_germline.tsv
```

It replaces the file previously made by hand with `germline_normal_select.R`. Only matched pairs are used, so no in-silico normal can reach a germline call.

### 2. Generating the cohort config file

The nextflow pipeline's config file encodes all of the options and inputs we might want to pass to the pipeline. dermanager unpacks it, with the wrapper script that launches the pipeline, into `commands/germline_pipe/`:

- `commands/germline_pipe/germline_variants.config`
- `commands/germline_pipe/run_germline.sh`

Both are the files in the pipeline's [`assets/`](https://github.com/team113sanger/dermatlas_germlinepost_nf/tree/develop/assets) directory. The config reads its cohort-dependent inputs from the environment set by the project `source_me.sh`, so for most runs nothing needs editing:

- the study ID (used in labelling output files) - `${STUDY}`
- the normal samples `.tsv` file (Step 1) - `${DNA_GERMLINE_NORMAL_MANIFEST}`
- the output directory to publish results into - `${ANALYSIS_DIR}/germline`

The other parameters mostly select which steps are included in a pipeline run and point at reference files. For convenience of maintaining the pipeline in a way that in can be run on or off farm22, all the reference files used by the pipeline are duplicated in `/lustre/scratch127/casm/teams/team113//secure-lustre/projects/dermatlas/resources/dermatlas`. These are direct copies of the resources directory you might find in other dermatlas PUs

### 3. Running the pipeline

:::{important}
**Different entry points**

This pipeline has two entry points: one starting from the raw sample BAMs for a cohort and one starting from the VCFs that have been produced by GATK haplotype caller. This is partly to help speed things up (so that you can avoid the computationally intensive steps at the start of the pipeline if you have already run without the need for a cached run of the pipeline.

As you might be aware, nextflow has a helpful cache-ing feature which keeps a record of which steps have been run and skips them. This ordinarily works very smoothly but there is a race condition which prevents caching working properly for this pipeline (calling of variants by chromosome) and breaks things when you try to rerun. 

This bug has been resolved by preventing cacheing in later versions of the pipeline (0.3.3+) but if you have a run that complains in this way,  please see the ii) Running from VCFs section
:::

```
WARN: [DERMATLAS_GERMLINE:GATK_GVCF_PER_CHROM (25)] Unable to resume cached task -- See log file for details
WARN: [DERMATLAS_GERMLINE:GATK_GVCF_PER_CHROM (24)] Unable to resume cached task -- See log file for details
```

#### **i) From bams**

Provided that you have set up your project with dermanager, launch the pipeline from the project directory with:

```bash
cd $PROJECT_DIR
bsub -e logs/germline.e -o logs/germline.o < commands/germline_pipe/run_germline.sh
```

The bsub magic at the start of the wrapper script will send a nextflow "master job", that looks after all other jobs to the oversubscribed queue (where it can live in peace running for a long period without fear of termination). Nextflow will shortly start submitting jobs on your behalf to the relevant queues.

The wrapper sources `source_me.sh`, checks the environment it needs, and runs the pipeline in `${PROJECT_DIR}/germline_pipe/`. Its `REVISION` selects which version of the pipeline to run; see the project [GitHub releases](https://github.com/team113sanger/dermatlas_germlinepost_nf/releases) for the latest versions. After a successful run it deletes the run's `work/` directory and makes `analysis/germline` group read-writable. Website logging, Slack notifications and the work-directory cleanup are each controlled by a toggle.

If you haven't initialised your project with dermanager, see "Without the website" in the [README](https://github.com/team113sanger/dermatlas_germlinepost_nf/blob/develop/README.md) for the environment the wrapper needs and how to provide it.

#### ii) From VCFs

If for some reason you aren't able to relaunch the pipeline with a cache - then you might want to run only the later steps of the pipeline. This is fairly straightforward to do by altering your germline config file.

You need only make three edits.

- Change post\_process\_only to TRUE  in the germline config file
- Provide a sample map file (detailing the links between VCF files and sample PD ids)
- Provide a path to the new genotype vcfs.

Here is what the updated config file should look like:

```
params {
    study_id = "${STUDY}"
    tsv_file = "${DNA_GERMLINE_NORMAL_MANIFEST}"
    outdir = "${ANALYSIS_DIR}/germline"
    geno_vcf = "${ANALYSIS_DIR}/germline/gatk_haplotypecaller/**.vcf.gz"
    sample_map = "${ANALYSIS_DIR}/germline/sample_map.txt"
    chrom_list = "${baseDir}/assets/grch38_chromosome.txt"
    post_process_only = true
    summarise_results = true
    samples_to_process = -1
    run_mode = "sort_inputs"
    run_coord_sort_cram = true
    run_deepvariant = false
    run_haplotypecaller = true
    run_markDuplicates = true
    baitset = "/lustre/scratch127/casm/projects/dermatlas/resources/baitset/GRCh38_WES5_canonical_pad100.merged.bed"
    reference_genome = "/lustre/scratch127/casm/projects/dermatlas/references/germline/genome.fa"
    vep_cache = "/lustre/scratch127/casm/projects/dermatlas/references/vep/cache/103"
    custom_files = "/lustre/scratch127/casm/projects/dermatlas/references/vep/cosmic/v97/CosmicV97Coding_Noncoding.normal.counts.vcf.gz{,.tbi};/lustre/scratch127/casm/projects/dermatlas/references/vep/clinvar/20230121/clinvar_20230121.chr.canonical.vcf.gz{,.tbi};/lustre/scratch127/casm/projects/dermatlas/references/vep/dbsnp/155/dbSNP155.GRCh38.GCF_000001405.39.mod.vcf.gz{,.tbi};/lustre/scratch127/casm/projects/dermatlas/references/vep/gnomad/v3.1.2/gnomad.genomes.v3.1.2.short.vcf.gz{,.tbi}"
    custom_args = "CosmicV97Coding_Noncoding.normal.counts.vcf.gz,Cosmic,vcf,exact,0,CNT;clinvar_20230121.chr.canonical.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT;dbSNP155.GRCh38.GCF_000001405.39.mod.vcf.gz,dbSNP,vcf,exact,0;gnomad.genomes.v3.1.2.short.vcf.gz,gnomAD,vcf,exact,0,FLAG,AF"
    nih_germline_resource = "/lustre/scratch127/casm/projects/dermatlas/resources/germline/national_genomic_test_germline_cancer_genes/output/Cancer_national_genomic_test_directory_v7.2_June_2023_gene_smv_summary.tsv"
    cancer_gene_census_resource = "/lustre/scratch127/casm/projects/dermatlas/resources/COSMIC/cancer_gene_census.v97.genes.tsv"
    flag_genes = "/lustre/scratch127/casm/projects/dermatlas/resources/germline/FLAG_genes_maftools.tsv"
    species = "homo_sapiens"
    filter_col = "gnomAD_AF"
    db_version = "103"
    assembly = "GRCh38"
    samples_to_process = -1
    publish_intermediates = false
    alternative_transcripts = "/lustre/scratch127/casm/projects/dermatlas/resources/ensembl/dermatlas_noncanonical_transcripts_ens103.v2.tsv"
  
}
```

After you have made the edit, submit a new run like so:

```
bsub -e logs/germline_vcf.e -o logs/germline_vcf.o < commands/germline_pipe/run_germline.sh
```

### Troubleshooting problem nextflow runs:

 There are several reasons the gemline pipeline might fail including bugs in the pipeline; issues with LSF; or misconfiguration.  In most cases (especially when you suspect a farm/ LSF failure), simply re-submitting the pipeline with

```
bsub -e logs/germline_vcf.e -o logs/germline_vcf.o < commands/germline_pipe/run_germline.sh
```

will trigger the nextflow `-resume` directive and the pipeline will pick up where it left off. A failed run always keeps its work directory, so it can be resumed; a successful one has its work directory deleted unless `DERMATLAS_CLEANUP_WORK_DIR=false` was set.

It is often worth taking a glance at the pipeline logs (`<YOUR_PROJECT_DIR>/analysis/logs/germline_calling_%J.o`) to follow and see what's going on, especially if jobs have failed.

When jobs fail, nextflow will provide the path to the directory a failed job was run in. I'd recommend inspecting the files in here with `ls -la` and printing some of the log files for the job with

```
cat .command.err
cat .command.out
cat .command.sh

```

:::{important}
**Multiple runs**

The wrapper runs each project's pipeline in its own directory (`${PROJECT_DIR}/germline_pipe`), so different cohorts can run in parallel. A second submission for the same project while one is still running fails immediately with exit code 75 and names the run holding the directory: wait for it to finish, or kill it, before resubmitting.

:::

### 4. Make a variant release

To generate a variant release containing all the summary tables, filtered variants, oncoplots, MAF file for the cohort, and readme, we use the `make_germ_variant_release.sh` script that lives in the GERMLINE codebase

**Usage:**

```bash
 bash make_germ_variant_release.sh --help
Usage: make_germ_variant_release.sh PROJECTDIR STUDY RELEASE
Description: Script to generate a release directory at <PROJECTDIR>/analysis/germline/releasev<RELEASE> 
Arguments:
  PROJECTDIR        Project directory full path 
  STUDY             Sequencescape STUDY sequencing ID
  RELEASE           Release number for germline  

Example:
  make_germ_variant_release.sh /My/DERMATLAS/PROJECTDIR 6674 2 
```

**For running this script you will require a:**

- **STUDY**: Sequencing study ID
- **PROJECTDIR**: Project directory
- **FINALJOINTDIR:** Full path to the directory where  the Final VCF files are
- **VERSION**: Release version you are creating (Default 1)

**Example running command:**


```bash
# Navigate into the project directory and source the project environmental variables
source source_me.sh
VERSION=1

# Configure your environment 
source ${PROJECTDIR}/scripts/germline/source_me.sh
 
#VCF Directories
FINALJOINTDIR=${PROJECTDIR}/analysis/germline/Final_joint_call

cd ${FINALJOINTDIR} 

# Make a release 
bash ${PROJECTDIR}/scripts/germline/scripts/make_germ_variant_release.sh ${PROJECTDIR:?unset} ${STUDY:?unset} ${VERSION:?unset}

```

#### **Expected Outputs:**

The code above will generate a release directory which has the following outputs:

Files and directories expected from the **revised SOP**:

```
tree -L 2 releasev1
releasev1
├── README.md
└── sum_files
    ├── 6674_germline.keep.maf
    ├── results
    ├── results_noflags_high
    └── results_noflags_modhigh

```

Useful links: 

- GATK Workflow details for joint genotype calling for cohorts <https://gatk.broadinstitute.org/hc/en-us/articles/360035890411-Calling-variants-on-cohorts-of-samples-using-the-HaplotypeCaller-in-GVCF-mode>
- This is a quick overview from GATK documentation on how to apply the workflow in practice. For more details, see the [Best Practices workflows](https://gatk.zendesk.com/hc/en-us/articles/360035894751) documentation.
