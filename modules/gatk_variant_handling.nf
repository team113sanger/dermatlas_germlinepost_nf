
process GATK_GVCF_PER_CHROM {
    tag {CHR[0]}
    label "gatk_steps"
    cache false
    container "broadinstitute/gatk:4.2.6.1"
    publishDir "${params.outdir}/vcf/Raw_joint_call", mode: "copy"
    input: 
    tuple val(meta), path(GENDB)
    tuple path(ref_genome_file), path(ref_genome_dict), path(reference_idx)
    each CHR

    output: 
    tuple val(CHR), path("*vcf.gz"), emit: chrom_vcf

    script:
    def chrom = CHR[0]
    """
    gatk --java-options "-Xmx4g -Xms4g" GenotypeGVCFs \
       -R "${ref_genome_file}" \
       -V gendb://${GENDB} \
       -L ${chrom} \
       -O ${chrom}.vcf.gz \
       --genomicsdb-shared-posixfs-optimizations true 
    """

    stub:
    def chrom = CHR[0]
    """
    echo stub > ${chrom}.vcf.gz
    """

}

process MERGE_COHORT_VCF {
    container "broadinstitute/gatk:4.2.6.1"
    label "gatk_steps"
    publishDir "${params.outdir}/vcf/Raw_joint_call", mode: "copy"
    input: 
    tuple val(meta), path(vcf_file_list)

    output: 
    tuple val(meta), path("*cohort_raw_genotypegvcf.vcf.gz"), emit: cohort_vcf

    script:
    """
    gatk --java-options "-Xmx8g -Xms8g" GatherVcfs \
    ${vcf_file_list.collect{ f -> "-I ${f}" }.join(' ')} \
    -O "${meta.study_id}_cohort_raw_genotypegvcf.vcf.gz"
    """
    stub:
    """
    echo stub > TBC_cohort_raw_genotypegvcf.vcf.gz
    """
}

process INDEX_COHORT_VCF {
    publishDir "${params.outdir}/vcf/Raw_joint_call", mode: "copy"
    label "gatk_steps"
    container "broadinstitute/gatk:4.2.6.1"

    input: 
    tuple val(meta), path(cohort_vcf)
    
    output:
    tuple val(meta), path(cohort_vcf), path("*.tbi"), emit: indexed_cohort_vcf

    
    script:
    """
    gatk --java-options "-Xmx8g -Xms8g" IndexFeatureFile \
    -I ${cohort_vcf}
    """
    stub:
    """
    echo stub > TBC.vcf.gz
    echo stub > TBC.vcf.gz.tbi
    """

}

process SELECT_VARIANTS {
    container "broadinstitute/gatk:4.2.6.1"
    label "gatk_steps"
    publishDir "${params.outdir}/vcf/Final_joint_call", mode: "copy"
    input: 
    tuple val(meta), path(raw_genotpye_vcf), path(raw_genotpye_index) 

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: raw_variants
    
    script:
    """
    gatk --java-options "-Xmx8g -Xms8g" SelectVariants \
     -V ${raw_genotpye_vcf} \
    -select-type ${meta.variant_type} \
    -O "${meta.study_id}${meta.suffix}.vcf.gz"
    """

    stub: 
    """
    echo stub > TBC.vcf.gz
    echo stub > TBC_index.vcf.gz.tbi
    """

}



process MARK_VARIANTS {
    container "broadinstitute/gatk:4.2.6.1"
    label "gatk_steps"
    publishDir "${params.outdir}/vcf/Final_joint_call", mode: "copy"
    input: 
    tuple val(meta), path(raw_genotpye_vcf), path(raw_genotpye_index)
    path(baitset)

    output:
    tuple val(meta), path("*marked.vcf.gz"), path("*marked.vcf.gz.tbi"), emit: marked_variants
    
    script:
    def clean_suffix = meta.suffix.replaceAll('_raw_', '_').replaceAll('_raw', '')
    def snp_filter  =  '-filter "QD < 2.0" --filter-name "QD2" ' +
                       '-filter "QUAL < 30.0" --filter-name "QUAL30" ' +
                       '-filter "SOR > 3.0" --filter-name "SOR3" ' +
                       '-filter "FS > 60.0" --filter-name "FS60" ' +
                       '-filter "MQ < 40.0" --filter-name "MQ40" ' +
                       '-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" ' +
                       '-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"' 
    def indel_filter = '-filter "QD < 2.0" --filter-name "QD2" ' +
                       '-filter "QUAL < 30.0" --filter-name "QUAL30" ' +
                       '-filter "FS > 200.0" --filter-name "FS200" ' +
                       '-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"'
        // Assign the appropriate filter based on meta.variant
    def filter_to_apply
    if (meta.variant_type == 'SNP') {
        filter_to_apply = snp_filter
    } else if (meta.variant_type == 'INDEL') {
        filter_to_apply = indel_filter
    } else {
        error "Unknown variant type: ${meta.variant_type}"
    }
    """
    gatk --java-options "-Xmx8g -Xms8g" VariantFiltration \
    -V ${raw_genotpye_vcf} \
    -L ${baitset} \
     ${filter_to_apply} \
    -O "${meta.study_id}${clean_suffix}.marked.vcf.gz"
    """
    stub:
    """
    echo stub > TBC.marked.vcf.gz
    echo stub > TBC.marked.vcf.gz.tbi
    """

}


process FILTER_VARIANTS {
    publishDir "${params.outdir}/vcf/Final_joint_call", mode: "copy"
    
    input: 
    tuple val(meta), path(marked_genotpye_vcf), path(marked_genotpye_index)
    path(baitset)

    output:
    tuple val(meta), path("*.target.pass.vcf.gz"), path("*.target.pass.vcf.gz.tbi"), emit: filtered_variants
    script:
    def clean_suffix = meta.suffix.replaceAll('_raw_', '_').replaceAll('_raw', '')
    """
    bcftools view -f PASS ${marked_genotpye_vcf} \
    -Oz -o "${meta.study_id}${clean_suffix}.marked.target.pass.vcf.gz" \
    -T ${baitset} 
    tabix -p vcf "${meta.study_id}${clean_suffix}.marked.target.pass.vcf.gz"
    """
    stub:
    """
    echo stub > demo.marked.target.pass.vcf.gz
    echo stub > demo.marked.target.pass.vcf.gz.tbi
    """

}
