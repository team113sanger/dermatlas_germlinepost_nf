process deepvariant {
    tag "$sample"
    publishDir "${params.outdir}/deepvariant/gvcf",
        mode: "${params.publish_dir_mode}",
        overwrite: true,
        pattern: "*dv.g.vcf.gz*"
    publishDir "${params.outdir}/deepvariant/vcf",
        mode: "${params.publish_dir_mode}",
        overwrite: true,
        pattern: "*dv.vcf.gz*"
     
    input:
    tuple val(sample), path(cram_file_sorted_dups_coord), path(cram_file_sorted_dups_coord_index)
    path(ref_genome)
    path(deepvariant_bed)
    output:
    tuple path("${cram_file_sorted_dups_coord}.dv.vcf.gz"), path("${cram_file_sorted_dups_coord}.dv.vcf.gz.tbi"), path("${cram_file_sorted_dups_coord}.dv.g.vcf.gz"), path("${cram_file_sorted_dups_coord}.dv.g.vcf.gz.tbi")

    script:
    """ 
    /opt/deepvariant/bin/run_deepvariant \\
    --ref=$ref_genome \\
    --regions=${params.bed_dir}/${params.bed_file_deepvariant} \\
    --reads=${cram_file_sorted_dups_coord} \\
    --output_vcf=${cram_file_sorted_dups_coord}.dv.vcf.gz \\
    --output_gvcf=${cram_file_sorted_dups_coord}.dv.g.vcf.gz \\
    --model_type=WES \\
    --customized_model=/opt/models/ukb_wes/model.ckpt-22236 \\
    --intermediate_results_dir ./tmp \\
    --num_shards=1
    """
}
