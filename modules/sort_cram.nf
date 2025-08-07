
process sort_cram {
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
    sambamba sort ${cram_file} -p -m ${task.memory.toGiga}GB -n \
    --tmpdir ./tmp /dev/stdin -o ${cram_file}.sorted \
   
    """
    
    stub:
    """
    echo stub > ${cram_file}.sorted
    """
}

