

process ANNOTATE_VARIANTS {
    publishDir "${params.outdir}/Final_joint_call", mode: "copy"

    input:
    tuple val(meta), path(vcf_file), path(vcf_index)
    path(vep_cache)
    path(vep_config)
    path(custom_files)
    val(custom_args)
    path(ref_genome)
    
    output: 
    tuple val(meta), path("*vep.vcf.gz"),path("*vep.vcf.gz.tbi"), emit: vep_annotation
    script: 
    def outfname = "${vcf_file}".replace(".vcf.gz", "") + ".vep.vcf.gz"
    """
    vep -i ${vcf_file} \
    --output_file ${outfname} \
    --dir_cache ${vep_cache} \
    --fasta ${ref_genome} \
    --cache \
    --offline \
    $custom_args \
    --cache \
    --db_version ${db_version} \
    -t SO \
    --format vcf \
    --buffer_size 20000 \ 
    --species {species_name} \
    --offline \
    --symbol \
    --biotype \
    --vcf \
    --sift s \
    --no_stats \
    --assembly {assembly_name} \
    --flag_pick_allele_gene \
    --canonical \
    --hgvs \
    --shift_hgvs 1 \
    --fasta {reference_fasta} \
    --compress_output bgzip \
    --mane \
    --numbers \
    --fork 4 \
    --domains
    tabix -p vcf ${outfname}
    """
    stub: 
    """
    echo stub > item.vep.vcf.gz
    echo stub > item.vep.vcf.gz.tbi
    """
    }


process ANNOTATE_FUR_VARIANTS {
    publishDir "${params.outdir}/Final_joint_call", mode: "copy"

    input:
    tuple val(meta), path(vcf_file), path(vcf_index)
    path(vep_cache)
    path(custom_files)
    val(custom_args)
    path(ref_genome)
    val(species)
    val(assembly_name)
    val(db_version)
    
    output: 
    tuple val(meta), path("*vep.vcf.gz"),path("*vep.vcf.gz.tbi"), emit: vep_annotation
    script: 
    def outfname = "${vcf_file}".replace(".vcf.gz", "") + ".vep.vcf.gz"
    """
    vep -i ${vcf_file} \
    --output_file ${outfname} \
    --cache \
    --dir ${vep_cache} \
    --fasta ${ref_genome} \
    --db_version ${db_version} \
    --species ${species} \
    --assembly ${assembly_name} \
    --offline \
    $custom_args \
    -t SO \
    --format vcf \
    --buffer_size 20000 \
    --offline \
    --symbol \
    --biotype \
    --vcf \
    --sift s \
    --no_stats \
    --flag_pick_allele_gene \
    --canonical \
    --hgvs \
    --shift_hgvs 1 \
    --compress_output bgzip \
    --mane \
    --numbers \
    --fork 4 \
    --domains
    tabix -p vcf ${outfname}
    """
    stub: 
    """
    echo stub > item.vep.vcf.gz
    echo stub > item.vep.vcf.gz.tbi
    """
    }
