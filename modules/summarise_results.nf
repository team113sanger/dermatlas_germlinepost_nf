

process CONVERT_TO_TSV {
    errorStrategy "ignore"
    container "gitlab-registry.internal.sanger.ac.uk/dermatlas/analysis-methods/germline:0.5.2"
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
    stub:
    """
    echo stub > ${meta.study_id}${meta.suffix}.marked.target.pass.vep.tsv.gz
    """
}
process FILTER_AND_ONCOPLOT {
    publishDir "${params.outdir}/Final_joint_call/sumtabs", mode: 'copy'
    container "gitlab-registry.internal.sanger.ac.uk/dermatlas/analysis-methods/germline/feature/upgrade_rdeps:79f4c529"

    input:
    tuple val(meta), path(maf)
    path(NIH_GERMLINE_TSV)
    path(CANCER_GENE_CENSUS)
    path(FLAG_GENES)
    path(samples_list)

    output: 
    path ("*.gnmadF.maf.gz"), emit: maf
    path("maftools_plots/*"), emit: oncoplots
    script: 
    """
    Rscript /opt/repo/scripts/filter_germline_mafs.R \
    --input_maf ${maf} \
    --germ_pred_tsv ${NIH_GERMLINE_TSV} \
    --flags_tsv ${FLAG_GENES} \
    --cgc_tsv ${CANCER_GENE_CENSUS} \
    --prefix ${meta.study_id}_germline \
    --samples ${samples_list} \
    --outdir "."
    """
    stub:
    """
    echo stub > demo.cohort_snps.indel.marked.target.pass.vep.gnmadF.maf.gz
    mkdir maftools_plots
    echo stub > maftools_plots/example.pdf
    """

}
