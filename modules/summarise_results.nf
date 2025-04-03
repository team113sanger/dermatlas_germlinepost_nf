process FILTER_AND_ONCOPLOT {
    publishDir "${params.outdir}/Final_joint_call/sum_files", mode: 'copy'
    container "gitlab-registry.internal.sanger.ac.uk/dermatlas/analysis-methods/germline:0.6.0"
    memory '16.GB'
    
    input:
    tuple val(meta), path(maf)
    path(NIH_GERMLINE_TSV)
    path(CANCER_GENE_CENSUS)
    path(FLAG_GENES)
    path(samples_list)

    output: 
    path("results*/*"), emit: oncoplots
    
    script: 
    """
    Rscript /opt/repo/scripts/filter_germline_mafs.R \
    --input_maf ${maf} \
    --outdir "." \
    --germ_pred_tsv ${NIH_GERMLINE_TSV} \
    --cgc_tsv ${CANCER_GENE_CENSUS} \
    --flags_tsv ${FLAG_GENES} \
    --prefix ${meta.study_id}_germline \
    --samples ${samples_list}
    """
    stub:
    """
    """

}
