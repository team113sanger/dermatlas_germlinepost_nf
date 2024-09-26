

process CONVERT_TO_TSV {
    container "gitlab-registry.internal.sanger.ac.uk/dermatlas/analysis-methods/germline:0.5.0"
    publishDir "${params.outdir}/Final_joint_call", mode: "copy"
    input: 
    tuple val(meta), path(vep_vcf), path(vep_index)

    output: 
    tuple val(meta), path("*.marked.target.pass.vep.tsv.gz"), emit: tsv_file
    
    script:
    """
    zcat ${vep_vcf} | /opt/repo/scripts/gatk_germline_full_vcf2table.v2.pl -> "${meta.study_id}${meta.suffix}.marked.target.pass.vep.tsv"
    gzip "${meta.study_id}${meta.suffix}.marked.target.pass.vep.tsv"
    """
}


process COMBINED_SUMMARY{
    container "gitlab-registry.internal.sanger.ac.uk/dermatlas/analysis-methods/germline:0.5.0"
    publishDir "results/Final_joint_call/sumtabs", mode: 'copy'

    input: 
    tuple val(meta), path(VEP_SNP_TSVGZ)
    tuple val(meta), path(VEP_INDEL_TSVGZ)
    path(NIH_GERMLINE_TSV)
    path(CANCER_GENE_CENSUS)
    path(FLAG_GENES)


    output:
    tuple val(meta), path("*.marked.target.pass.vep.gnmadF.tsv.gz"), emit: outfile
    tuple val(meta), path("*.tsv.gz"), emit: summary
    tuple val(meta), path("*.xlsx"), emit: xlsheets


    script:
    """
    /opt/repo/scripts/get_XSLX_andfilt_COSMIC_and_clinvar_PATHVars.R \
    --study_id ${meta.study_id} \
    --input_snv_tsv ${VEP_SNP_TSVGZ} \
    --input_indel_tsv ${VEP_INDEL_TSVGZ} \
    --germ_pred_tsv ${NIH_GERMLINE_TSV} \
    --cgc_tsv ${CANCER_GENE_CENSUS} \
    --flags_tsv ${FLAG_GENES} \
    --outdir "."
    """
}

process CONVERT_TO_MAF {
    publishDir "results/Final_joint_call/sumtabs", mode: 'copy'
    container "gitlab-registry.internal.sanger.ac.uk/dermatlas/analysis-methods/germline:0.5.0"


    input:
    tuple val(meta), path(SNP_INDEL_GNAMDF_TSVGZ)
    path(NIH_GERMLINE_TSV)
    path(CANCER_GENE_CENSUS)

    output: 
    path ("*.gnmadF.maf.gz"), emit: maf
    path("oncplots/*"), emit: oncoplots
    script: 
    """
    /opt/repo/scripts/convert_germ_sumtsv_to_MAF_plot.R \
    --study_id ${meta.study_id} \
    --input_germ_tsv ${SNP_INDEL_GNAMDF_TSVGZ} \
    --germ_pred_tsv ${NIH_GERMLINE_TSV} \
    --cgc_tsv ${CANCER_GENE_CENSUS} \
    --outdir "."
    """

}
