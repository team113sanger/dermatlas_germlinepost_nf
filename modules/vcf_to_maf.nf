process VCF_TO_KEEP_MAF{
    container "gitlab-registry.internal.sanger.ac.uk/dermatlas/analysis-methods/maf:0.6.2"


input:    
    tuple val(meta), path(STUDY_COHORT_VCFS)
    path(sample_list)
    val(assembly)
    val(filter_col)
    path(alt_transcripts)
output: 
    tuple val(meta), path("*_germline.keep.maf"), emit: maf

script:
"""
/opt/repo/reformat_vcf2maf.pl \
--build ${assembly} \
--tum_vaf 0.1 \
--vaf_filter_type all \
--transcripts ${alt_transcripts} \
--indel_filter \
--keep_multi \
--keep \
--af_col ${filter_col} \
--vcflist ${sample_list}
"""
}


process VCF_TO_PASS_MAF{
    container "gitlab-registry.internal.sanger.ac.uk/dermatlas/analysis-methods/maf:0.6.2"
    input:
        tuple val(meta), path(STUDY_COHORT_VCFS)
        path(sample_list)
        path(alt_transcripts)
        val(assembly)
        val(filter_col)
    output: 
        tuple val(meta), path("*_germline.pass.maf"), emit: maf

    script:
    """
    /opt/repo/reformat_vcf2maf.pl \
    --build ${assembly} \
    --tum_vaf 0.1 \
    --vaf_filter_type all \
    --transcripts ${alt_transcripts} \
    --indel_filter \
    --keep_multi \
    --pass \
    --canonical \
    --exclude_noncoding \
    --vcflist ${sample_list}
    """
}