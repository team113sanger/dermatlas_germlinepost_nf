process markDuplicates {
    tag "$sample"
    container "broadinstitute/gatk:4.2.6.1"
    label "gatk_steps"
    publishDir path: "${params.outdir}/markDuplicates/",
	mode: "${params.publish_dir_mode}",
	overwrite: "true"

    when:
    params.run_markDuplicates
     
    input:
    tuple val(sample), path(cram_file_sorted)
    tuple path(ref_genome), path(ref_genome_dict), path(reference_idx)


    output:
    tuple val(sample), path("${cram_file_sorted}.dups"), emit: markdup_sample_cram
    path("${cram_file_sorted}.dups")

    script:
    """ 
    gatk --java-options \"-Xms4g -Xmx4g  -XX:+UseSerialGC\" \\
      MarkDuplicates \\
      -I ${cram_file_sorted} \\
      -O ${cram_file_sorted}.dups \\
      --METRICS_FILE ${cram_file_sorted}.dups.metrics \\
      -R $ref_genome \\
      --VALIDATION_STRINGENCY SILENT \\
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \\
      --ASSUME_SORT_ORDER \"queryname\" \\
      --CLEAR_DT \"true\" \\
      --ADD_PG_TAG_TO_READS false \\
      --TMP_DIR ./tmp
    """
    stub: 
    """
    echo stub > ${cram_file_sorted}.dups
    """
}

