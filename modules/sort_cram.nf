
process sort_cram {
    container "quay.io/biocontainers/sambamba:0.6.5--0"
    tag "$sample"
    if (params.publish_intermediates){
    publishDir path: "${params.outdir}/sort_cram/",
	mode: "${params.publish_dir_mode}",
	overwrite: "true"
    }

    input:
    tuple val(sample), path(cram_file), path(cram_file_index)

    output:
    tuple val(sample), path("${cram_file}.sorted"), emit: sorted_sample_cram

    script:
    """ 
    sambamba sort ${cram_file} -p -m 7GB -n \
    --tmpdir ./tmp /dev/stdin -o ${cram_file}.sorted \
   
    """
    
    stub:
    """
    echo stub > ${cram_file}.sorted
    """
}

