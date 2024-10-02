
process sort_cram {
    container "quay.io/biocontainers/sambamba:0.6.5--0"
    tag "$sample"
    publishDir path: "${params.outdir}/sort_cram/",
	mode: "${params.publish_dir_mode}",
	overwrite: "true"

    input:
    tuple val(sample), path(cram_file), path(cram_file_index)

    output:
    tuple val(sample), path("${cram_file}.sorted"), emit: sorted_sample_cram

    script:
    """ 
    sambamba sort -p -m 7GB -n \
    --tmpdir ./tmp /dev/stdin -o ${cram_file}.sorted \
    ${cram_file}
    """
    
    stub:
    """
    echo stub > ${cram_file}.sorted
    """
}

