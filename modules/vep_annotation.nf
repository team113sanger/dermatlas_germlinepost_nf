

process ANNOTATE_VARIANTS {
    publishDir "${params.outdir}/Final_joint_call", mode: "copy"
    container "ensemblorg/ensembl-vep:release_103.1"
    input:
    tuple val(meta), path(vcf_file), path(vcf_index)
    path(vep_cache)
    path(vep_config)
    path(gnomadfile)
    path(dbsnpfile)
    path(clinvarfile)
    path(cosmicfile)
    path(ref_genome)
    
    output: 
    tuple val(meta), path("*vep.vcf.gz"),path("*vep.vcf.gz.tbi"), emit: vep_annotation
    script: 
    def outfname = "${vcf_file}".replace(".vcf.gz", "") + ".vep.vcf.gz"
    """
    vep -i ${vcf_file} \
    --dir ${vep_cache} \
    --config ${vep_config} \
    --output_file ${outfname}
    tabix -p vcf ${outfname}
    """
    stub: 
    """
    echo stub > item.vep.vcf.gz
    echo stub > item.vep.vcf.gz.tbi
    """
    }