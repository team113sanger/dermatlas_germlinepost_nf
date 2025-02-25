process VCF_TO_MAF {
    container "gitlab-registry.internal.sanger.ac.uk/dermatlas/analysis-methods/maf:0.6.2"
    input: 
        tuple val(meta), path(STUDY_COHORT_SNP_VCFS), path(snp_index)
        tuple val(indel_meta), path(STUDY_COHORT_INDEL_VCFS), path(indel_index)
        path(sample_list)
        val(assembly)
        val(filter_col)
        path(alt_transcripts)
    output: 
        tuple val(meta), path("*.maf"), emit: maf

    script:
    def keep_filter =  "--keep --af_col ${filter_col}"
    def pass_filter =  "--pass --canonical --exclude_noncoding"
    def args 
    if (meta.variant_set == 'keep'){
        args = keep_filter
    } else if (meta.variant_set == 'pass'){
        args = pass_filter
    } else {
        error "Unknown variant_set: ${meta.variant_set}"
    }
    """
    /opt/repo/reformat_vcf2maf.pl \
    --build ${assembly} \
    --tum_vaf 0.1 \
    --vaf_filter_type all \
    --transcripts ${alt_transcripts} \
    --indel_filter \
    --keep_multi \
    ${args} \
    --vcflist ${sample_list} > ${meta.study_id}_germline.${meta.variant_set}.maf
    """
}
