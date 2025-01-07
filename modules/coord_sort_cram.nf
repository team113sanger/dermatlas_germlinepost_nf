process coord_sort_cram {
    tag "$sample"
    publishDir path: "${params.outdir}/coord_sort_cram/",
	mode: "${params.publish_dir_mode}",
	overwrite: "true"

    when:
    params.run_coord_sort_cram
     
    input:
    tuple val(sample), path(cram_file_sorted_dups)

    output:
    tuple val(sample), path("${sample}.sorted.dups.coord.bam"), path("${sample}.sorted.dups.coord.bam.bai"), emit: markdup_sample_cram_crai

    script:
    """ 
    sambamba sort -p -m 7GB --tmpdir ./tmp ${cram_file_sorted_dups} -o ${sample}.sorted.dups.coord.bam
    sambamba index ${sample}.sorted.dups.coord.bam
    """
    stub: 
    """
    echo stub > "${sample}.sorted.dups.coord.bam"
    echo stub > "${sample}.sorted.dups.coord.bam.bai"
    """
}

