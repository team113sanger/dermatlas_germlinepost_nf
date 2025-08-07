process BAM_TO_CRAM {
    tag "$sample"
    if (params.publish_intermediates){
    publishDir path: "${params.outdir}/bam_to_cram/",
	mode: "${params.publish_dir_mode}",
	overwrite: "true"
    }

    when:
    params.run_coord_sort_cram
     
    input:
    tuple val(sample), path(cram_file_sorted_dups_coord), path(cram_file_sorted_dups_coord_index)
    tuple path(ref_genome), path(ref_genome_dict), path(reference_idx)

    output:
    tuple val(sample), path("${sample}.sorted.dups.coord.cram"), path("${sample}.sorted.dups.coord.cram.crai"), emit: processed_sample_cram_crai

    script:
    """ 
    samtools view ${cram_file_sorted_dups_coord} \\
    -o ${sample}.sorted.dups.coord.cram -O CRAM \\
    -T $ref_genome

    samtools index ${sample}.sorted.dups.coord.cram
    """
    stub:
    """
    echo stub > ${sample}.sorted.dups.coord.cram
    echo stub > ${sample}.sorted.dups.coord.cram.crai
    """
}

