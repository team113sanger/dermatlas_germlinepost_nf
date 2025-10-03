# Nextflow: Germline variant calling pipeline

Germline variant calling and post-processing for DERMATLAS can be run mostly with a single nextflow pipeline in a largely "set-and-forget" manner to reproduce the manual steps detailed in [DERMATLAS - Germline calling with GATK - for WES - using Nextflow Tower](https://confluence.sanger.ac.uk/x/BJOeB). This document contains an overview of how to configure and run the pipeline. For a more detailed explanation of the pipeline, the inputs, steps and requirements for running can be found within the pipeline project [README](https://gitlab.internal.sanger.ac.uk/DERMATLAS/analysis-methods/dermatlas_germlinepost_nf/-/blob/develop/README.md?ref_type=heads)

## Workflow Overview

1. Generate the normal input table
2. Generating the cohort config file
3. Running the pipeline
   - i) From BAMS
   - ii) From VCFs
4. Make a release folder
5. Cleanup the intermediate BAM files created by the pipeline

## Workflow Steps

:::{tab-set}

:::{tab-item} 1. Generate Input Table

Generating the input table for samples to run germline calling on a study works essentially in the same way as when manually running of the pipeline. It requires a `.tsv` file detailing the normal samples to run, with the following columns:

| sample   | object                                                                                                                                   | object index                                                                                                                                 |
|:---------|:-----------------------------------------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------|
| PD42171b | /lustre/scratch124/casm/team113/projects/5534\_Landscape\_sebaceous\_tumours\_GRCh38\_Remap\_germline/BAMS/PD42171b.sample.dupmarked.bam | /lustre/scratch124/casm/team113/projects/5534\_Landscape\_sebaceous\_tumours\_GRCh38\_Remap\_germline/BAMS/PD42171b.sample.dupmarked.bam.bai |

```

The easiest means to create this table is to run the `germline_normal_select.R` script from GERMLINE which will populate these fields from your matched tumour normal pairs. Provided that you have installed analysis methods with dermanager, this can be accomplished by:

Navigate into your project directory:
```bash
cd $PROJECT_DIR
```

Then get ready to run the germline_normal_select.R script like so:

```bash
# Setup project environmental variables

source source_me.sh

# Setup germline analysis environment variables 
source ${PROJECTDIR}/scripts/germline/source_me.sh

# Run the script to get the final table format for input to nf_deepvariant
Rscript ${PROJECTDIR}/scripts/germline/scripts/germline_normal_select.R \
--study_id ${STUDY} \
--bam_dir ${PROJECTDIR}/bams \
--sample_pairs ${PROJECTDIR}/metadata/${STUDY}_${PROJECT}-one_tumour_per_patient_matched.tsv \
--outdir ${PROJECTDIR}/metadata
```

:::{note} 

If Rejected samples are available

If DNA samples in the cohort need to be rejected from sample list creation then you can use the update the "rejected DNA samples" table  (`${PROJECTDIR}/metadata/rejected\_DNA\_samples.txt`) and ignore the new IDs by running the command with the `–remove\_list` parameter:



```bash
Rscript ${PROJECTDIR}/scripts/GERMLINE/scripts/germline_normal_select.R \
--study_id ${STUDY:?unset} \
--bam_dir ${PROJECTDIR:?unset}/bams \
--sample_pairs ${PROJECTDIR:?unset}/metadata/${STUDY:?unset}_${PROJECT:?unset}-one_tumour_per_patient_matched.tsv \
--remove_list ${PROJECTDIR:?unset}/metadata/rejected_DNA_samples.txt \
--outdir ${PROJECTDIR:?unset}/metadata 
```

This will generate a file called : \*\_**normal\_one\_per\_patient\_matched\_selected\_germl\_samples.tsv**

:::

:::


:::{tab-item} 2. Generating the cohort config file

The pipeline's config file encodes all of the options and inputs we might want to pass to the pipeline. For newer versions of Dermanager `commands/germline.config`

For most pipeline runs there are only **3** parameters you'll need to change to get things going:

- The study ID (used in labelling output files)
- The path to the normal samples `.tsv`  file (generated in Step 1)
- The output directory to publish results into

There is a large set of other parameters specified within this file but won't normally need changing. These other parameters mostly modify which steps are included in a pipeline run and paths to reference files. For convenience of maintaining the pipeline in a way that in can be run on or off farm22, all the reference files used by the pipeline are stored in

`/lustre/scratch124/casm/team113/secure-lustre/resources/dermatlas. T`hese are direct copies of the resources directory you might find in other dermatlas PUs

**germline.config**

```
params {
    study_id = "${STUDY}"
    tsv_file = "${PROJECT_DIR}/metadata/${STUDY}_normal_one_per_patient_matched_selected_germl_samples.tsv"
    outdir = "${PROJECT_DIR}/analysis/germline"
    chrom_list = "${baseDir}/assets/grch38_chromosome.txt"
    post_process_only = false
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

:::

:::{tab-item} 3. Running the pipeline

:::{important}
Different entry points

This pipeline has two entry points: one starting from the raw sample BAMs for a cohort and one starting from the VCFs that have been produced by GATK haplotype caller. This is partly to help speed things up (so that you can avoid the computationally intensive steps at the start of the pipeline if you have already run without the need for a cached run of the pipeline.

As you might be aware, nextflow has a helpful cache-ing feature which keeps a record of which steps have been run and skips them. This ordinarily works very smoothly but there is a race condition which prevents caching working properly for this pipeline (calling of variants by chromosome) and breaks things when you try to rerun. 

This bug has been resolved by preventing cacheing in later versions of the pipeline (0.3.3+) but if you have a run that complains in this way,  please see the ii) Running from VCFs section
:::

```
WARN: [DERMATLAS_GERMLINE:GATK_GVCF_PER_CHROM (25)] Unable to resume cached task -- See log file for details
WARN: [DERMATLAS_GERMLINE:GATK_GVCF_PER_CHROM (24)] Unable to resume cached task -- See log file for details
```

**i) From bams**

Once you have your inputs and are authenticated you can prepare to launch the pipeline by modifying and save this wrapper script. You will need to update the desitination of the params file + desired log file locations. 

In this script the "`-r"`  option specifies which version of the pipeline you'd like to run. Normally you should select the latest version.

**Example file:**

**run\_germline\_calling.sh**

```
#!/bin/bash
#BSUB -q oversubscribed
#BSUB -G team113-grp
#BSUB -R "select[mem>8000] rusage[mem=8000] span[hosts=1]"
#BSUB -M 8000
#BSUB -oo <CHANGE_ME_TO_PROJECT_DIR>/analysis/logs/gemline_calling_%J.o
#BSUB -eo <CHANGE_ME_TO_PROJECT_DIR>/analysis/logs/gemline_calling_%J.e
set -euo pipefail

source source_me.sh
export CONFIG="${PROJECT_DIR}/commands/germline_variants.config"
export REVISION="0.3.2"


# Load module dependencies
module load nextflow-23.10.0
module load /software/modules/ISG/singularity/3.11.4
 
nextflow pull 'https://github.com/team113sanger/dermatlas_germlinepost_nf' 

nextflow run 'https://github.com/team113sanger/dermatlas_germlinepost_nf' \
-resume \
-r "${REVISION}" \
-c "${CONFIG}" \
-profile farm22
```

If you called the script `run_germline.sh` then you'll now be able to submit:

```
bsub < run_germline.sh
```

The bsub magic at the start of the wrapper script will send a nextflow "master job", which looks after all other jobs to the oversubscribed queue (where it can live in peace running for a long period without fear of termination). Nextflow will shortly start submitting jobs on your behalf to the relevant queues

### ii) From VCFs

If for some reason you aren't able to relaunch the pipeline with a cache - then you might want to run only the later steps of the pipeline. This is fairly straightforward to do by altering your germline config file.

You need only make three edits.

- Change post\_process\_only to TRUE  in the germline config file
- Provide a sample map file (detailing the links between VCF files and sample PD ids)
- Provide a path to the new genotype vcfs.

Here is what the updated config file should look like:

```
params {
    study_id = "${STUDY}"
    tsv_file = "${PROJECT_DIR}/metadata/${STUDY}_normal_one_per_patient_matched_selected_germl_samples.tsv"
    outdir = "${PROJECT_DIR}/analysis/germline"
    geno_vcf = "${PROJECT_DIR}/analysis/germline/gatk_haplotypecaller/**.vcf.gz"
	sample_map = "${PROJECT_DIR}/analysis/germline/sample_map.txt"
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
bsub < run_germline.sh
```

:::

::::

### Troubleshooting problem nextflow runs:

 There are several reasons the gemline pipeline might fail including bugs in the pipeline; issues with LSF; or misconfiguration.  In most cases (especially when you suspect a farm/ LSF failure), simply re-submitting the pipeline with

```
bsub < run_germline.sh
```

will trigger the nextflow `-resume` directive and the pipeline will pick up where it left off.

It is often worth taking a glance at the pipeline logs (<YOUR\_PROJECT\_DIR>/analysis/logs/gemline\_calling\_%J.o) to follow and see what's going on, especially if things have failed/

When jobs fail, nextflow will provide the path to the directory a failed job was run in. I'd recommend inspecting the files in here with `ls -la` and printing some of the log files for the job with

```
cat .command.err
cat .command.out
cat .command.sh

```

:::{important}
Multiple runs

Nextflow is able to keep track of past runs by creating a .nextflow directory in the current location and stores intermediate files in a work. If you want to run the same pipeline but on different cohorts (e.g. hidradenomas and hidradenocarcionmas) in parallel, please ensure that you launch each instance of the pipeline in a seperate directory - otherwise nextflow can't keep track of what is going on an report errors about "nextflow lock files "

:::

## **Make a variant release**

---

To generate a variant release containing all the summary tables, filtered variants, oncoplots, MAF file for the cohort, and readme,we use the `make\_germ\_variant\_release.sh` script that lives in the GERMLINE analysis method

**Usage:**

```bash
 bash make_germ_variant_release.sh --help
Usage: make_germ_variant_release.sh PROJECTDIR STUDY RELEASE
Description: Script to generate a release directory at <PROJECTDIR>/analysis/germline/releasev<RELEASE> 
Arguments:
  PROJECTDIR        Project directory full path 
  STUDY             Sequencescape STUDY sequencing ID
  RELEASE           Relase number for germline  

Example:
  make_germ_variant_release.sh /My/DERMATLAS/PROJECTDIR 6674 2 
```

#### **For running this script you will require the:**

- **STUDY**: Sequencing study ID
- **PROJECTDIR**: Project directory
- **FINALJOINTDIR:** Full path to the directory where  the Final VCF files are
- VERSION: Release version you are creating (Default 1)

**Example running command:**

**[]()**It is required to run  **source\_me.sh** prior to the job execution 

**Convert VCF to MAFs and plot**

```bash
# Navigate into the project directory and source the project environmental variables
source source_me.sh
VERSION=1

# Configure your environment 
source ${PROJECTDIR}/scripts/germline/source_me.sh
 
#VCF Directories
FINALJOINTDIR=${PROJECTDIR}/analysis/germline/Final_joint_call

cd ${FINALJOINTDIR} 

#MAke a release 
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

Useful links : 

- GATK Workflow details for joint genotype calling for cohorts <https://gatk.broadinstitute.org/hc/en-us/articles/360035890411-Calling-variants-on-cohorts-of-samples-using-the-HaplotypeCaller-in-GVCF-mode>
- This is a quick overview from GATK documentation on how to apply the workflow in practice. For more details, see the [Best Practices workflows](https://gatk.zendesk.com/hc/en-us/articles/360035894751) documentation.
