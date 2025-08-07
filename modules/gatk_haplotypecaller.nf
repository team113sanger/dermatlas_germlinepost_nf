process GATK_HAPLOTYPECALLER {
    tag "$meta"
    container "broadinstitute/gatk:4.2.6.1"
    label "gatk_steps"
    publishDir path: "${params.outdir}/gatk_haplotypecaller/",
	mode: "${params.publish_dir_mode}",
	overwrite: "true",
	pattern: "*gatk.g.vcf.gz*"
    
    when:
    params.run_haplotypecaller
     
    input:
    tuple val(meta), path(cram_file_sorted_dups_coord), path(cram_file_sorted_dups_coord_index)
    tuple path(reference_genome), path(ref_genome_dict), path(reference_idx)
    path(baitset)

    output:
    tuple val(meta), path("${cram_file_sorted_dups_coord}.gatk.g.vcf.gz"), path("${cram_file_sorted_dups_coord}.gatk.g.vcf.gz.tbi")

    script:
    """ 
    /gatk/gatk --java-options \"-Xms6g -Xmx6g -XX:+UseSerialGC\" \\
    HaplotypeCaller \\
    -I ${cram_file_sorted_dups_coord} \\
    -O ${cram_file_sorted_dups_coord}.gatk.g.vcf.gz \\
    -R $reference_genome \\
    -L $baitset \\
    -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \\
    -ERC GVCF \\
    -G StandardAnnotation \\
    -G StandardHCAnnotation \\
    -G AS_StandardAnnotation
    """
    stub:
    """
    echo stub > ${cram_file_sorted_dups_coord}.gatk.g.vcf.gz
    echo stub > ${cram_file_sorted_dups_coord}.gatk.g.vcf.gz.tbi
    """
    }

