#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process MERGE_COHORT_VCF {
    container "file:/software/team113/dermatlas/singularity_images/gatk__4.2.6.1.sif"
    input: 
    tuple val(meta), path(vcf_file_list)

    output: 
    tuple val(meta), path("*cohort_raw_genotypegvcf.vcf.gz"), emit: cohort_vcf

    script:
    """
    gatk --java-options "-Xmx8g -Xms8g" GatherVcfs \
    -I ${vcf_file_list}
    -O "${meta.study}_cohort_raw_genotypegvcf.vcf.gz"
    """
}

process INDEX_COHORT_VCF {
    container "file:/software/team113/dermatlas/singularity_images/gatk__4.2.6.1.sif"

    input: 
    tuple val(meta), path(cohort_vcf)
    
    output:
    tuple val(meta), path(cohort_vcf), path("*.tbi"), emit: indexed_cohort_vcf

    
    script:
    """
    gatk --java-options "-Xmx8g -Xms8g" IndexFeatureFile \
    -I ${cohort_vcf}
    """

}

process SELECT_SNP_VARIANTS {
    input: 
    tuple val(meta), path(raw_genotpye_vcf),  path(raw_genotpye_index) 

    output:
    tuple val(meta), path("*_cohort_raw_snps.vcf.gz"), emit: raw_variants
    
    script:
    """
    gatk --java-options "-Xmx8g -Xms8g" SelectVariants \
     -V ${raw_genotpye_vcf} \
    -select-type SNP \
    -O "${meta.study}_cohort_raw_snps.vcf.gz"
    """

}

process SELECT_INDEL_VARIANTS {
    input: 
    tuple val(meta), path(raw_genotpye_vcf),  path(raw_genotpye_index) 

    output:
    tuple val(meta), path("*_cohort_indel_raw.vcf.gz"), emit: raw_variants
    
    script:
    """
    gatk --java-options "-Xmx8g -Xms8g" SelectVariants \
     -V ${raw_genotpye_vcf} \
    -select-type INDEL \
    -O "${meta.study}_cohort_indel_raw.vcf.gz"
    """

}



process MARK_SNP_VARIANTS {
    input: 
    tuple val(meta), path(raw_genotpye_vcf),  path(raw_genotpye_index) 
    path(baitset)

    output:
    tuple val(meta), path("*_cohort_raw_snps.vcf.gz"), emit: marked_variants
    script:
    """
    gatk --java-options "-Xmx8g -Xms8g" VariantFiltration \
    -V ${raw_genotpye_vcf} \
    -L ${baitset} \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O "${meta.study}_cohort_snps.marked.vcf.gz"
    """

}


process MARK_INDEL_VARIANTS {
    input: 
    tuple val(meta), path(raw_genotpye_vcf),  path(raw_genotpye_index) 
    path(baitset)

    output:
    tuple val(meta), path("*_cohort_indel.marked.vcf.gz"), emit: marked_variants
    script:
    """
    gatk --java-options "-Xmx8g -Xms8g" VariantFiltration \
    -V ${raw_genotpye_vcf} \
    -L ${baitset} \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O "${meta.study}_cohort_indel.marked.vcf.gz"
    """

}

process FILTER_SNP_VARIANTS {
    module "bcftools"
    input: 
    tuple val(meta), path(marked_genotpye_vcf),  path(raw_genotpye_index) 
    path(baitset)

    output:
    tuple val(meta), path("*_cohort_snps.marked.target.pass.vcf.gz"), path("*_cohort_snps.marked.target.pass.vcf.gz.tbi"), emit: filtered_variants
    script:
    """
    bcftools view -f PASS ${marked_genotpye_vcf} \
    -Oz -o "${meta.study}_cohort_snps.marked.target.pass.vcf.gz" \
    -T ${baitset} 
    tabix -p vcf "${meta.study}_cohort_snps.marked.target.pass.vcf.gz"
    """

}


process FILTER_INDEL_VARIANTS {
    module "bcftools"
    input: 
    tuple val(meta), path(marked_genotpye_vcf),  path(raw_genotpye_index) 
    path(baitset)

    output:
    tuple val(meta), path("*_cohort_indels.marked.target.pass.vcf.gz"), path("*_cohort_snps.marked.target.pass.vcf.gz.tbi"), emit: filtered_variants
    script:
    """
    bcftools view -f PASS ${marked_genotpye_vcf} \
    -Oz -o "${meta.study}_cohort_indel.marked.target.pass.vcf.gz" \
    -T ${baitset} 
    tabix -p vcf "${meta.study}_cohort_indel.marked.target.pass.vcf.gz"
    """

}

process ANNOTATE_VARIANTS {
    module "ensembl_vep/103.1"
    input:
    tuple val(meta), path(vcf_file)
    path(gnomadfile)
    path(clinvarfile)
    path(dbsnpfile)
    path(cosmicfile)
    val(ens_ver)
    val(species)
    val(assembly)
    path(vep_cache)
    
    output: 
    tuple val(meta), path("*vep.vcf.gz"), emit: vep_annotation
    script: 
    def custom = "--custom $gnomadfile,gnomAD,vcf,exact,0,FLAG,AF --custom $clinvarfile,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT --custom $dbsnpfile,dbSNP,vcf,exact,0, --custom $cosmicfile,COSMIC,vcf,exact,0,CNT"
    def outfname = meta.id + "vep.vcf.gz"
    """
    vep -i ${vcf_file} \
    --cache_version $ens_ver \
    -t SO \
    --format vcf \
    -o $outfname \
    --cache \
    --dir $vep_cache \
    --buffer 20000 \
    --species $species \
    --offline \
    --symbol \
    --biotype --vcf --sift s \
    --no_stats \
    --assembly $assembly \
    --flag_pick_allele \
    --canonical \
    --hgvs \
    --shift_hgvs 1 \
    --fasta $ref \
    $custom \
    --compress_output bgzip \
    --mane --numbers --polyphen p \
    --domain --transcript_version \
    --show_ref_allele --fork 4
    """
}

process CONVERT_TO_TSV {
    input: 
    tuple val(meta), path(vep_vcf)

    output: 
    tuple val(meta), path("*.marked.target.pass.vep.tsv.gz")
    
    script:
    """
    zcat ${vep_vcf} | gatk_germline_full_vcf2table.v2.pl -> "${meta.study}_cohort_snps.marked.target.pass.vep.tsv"
    gzip "${meta.study}_cohort_snps.marked.target.pass.vep.tsv"
    """
}


process COMBINED_SUMMARY{
    publishDir "results", mode: 'copy'

    input: 
    tuple val(meta), path(VEP_SNP_TSVGZ), path(VEP_INDEL_TSVGZ)
    path(NIH_GERMLINE_TSV)
    path(CANCER_GENE_CENSUS)
    path(FLAG_GENES)


    output:
    path("*.tsv")
    path("*.xlsx")


    script:
    """
    get_XSLX_andfilt_COSMIC_and_clinvar_PATHVars.R \
    --study_id ${meta.study} \
    --input_snv_tsv ${VEP_SNP_TSVGZ} \
    --input_indel_tsv ${VEP_INDEL_TSVGZ} \
    --germ_pred_tsv ${NIH_GERMLINE_TSV} \
    --cgc_tsv ${CANCER_GENE_CENSUS} \
    --flags_tsv ${FLAG_GENES} \
    --outdir "."
    """
}

process CONVERT_TO_MAF {
    publishDir "results", mode: 'copy'

    input:
    tuple val(meta), path(SNP_INDEL_GNAMDF_TSVGZ), path(NIH_GERMLINE_TSV), path(CANCER_GENE_CENSUS)

    output: 
    path "*cohort_snps.indel.marked.target.pass.vep.gnmadF.maf.gz", emit: maf
    path("oncplots/*"), emit: oncoplots
    script: 
    """
    convert_germ_sumtsv_to_MAF_plot.R \
    --study_id ${meta.study} \
    --input_germ_tsv ${SNP_INDEL_GNAMDF_TSVGZ} \
    --germ_pred_tsv ${NIH_GERMLINE_TSV} \
    --cgc_tsv ${CANCER_GENE_CENSUS} \
    --outdir ${SUMTABDIR}
    """

}


workflow POST_PROCESS_VARIANTS {
    MERGE_COHORT_VCF()
    INDEX_COHORT_VCF(MERGE_COHORT_VCF.out.cohort_vcf)
    SELECT_SNP_VARIANTS(INDEX_COHORT_VCF.out.indexed_cohort_vcf)
    SELECT_INDEL_VARIANTS(INDEX_COHORT_VCF.out.indexed_cohort_vcf)

    MARK_SNP_VARIANTS(SELECT_SNP_VARIANTS.out.raw_variants)
    MARK_INDEL_VARIANTS(SELECT_INDEL_VARIANTS.out.raw_variants)

    FILTER_SNP_VARIANTS(MARK_SNP_VARIANTS.out.marked_variants)
    FILTER_INDEL_VARIANTS(MARK_INDEL_VARIANTS.out.marked_variants)
    
    ANNOTATE_VARIANTS(FILTER_SNP_VARIANTS.out.filtered_variants)
    //ANNOTATE_VARIANTS(MARK_INDEL_VARIANTS.out.filtered_variants)

    CONVERT_TO_TSV(ANNOTATE_VARIANTS.out.vep_annotation)
}