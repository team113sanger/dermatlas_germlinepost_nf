
process GENERATE_GENOMICS_DB {
    container "broadinstitute/gatk:4.2.6.1"
    publishDir "${params.outdir}"
    
    input: 
    path(SAMPLE_MAP)
    val(chroms)
    tuple val(meta), path(genotype_vcfs)

    output: 
    tuple val(meta), path("*data"), emit: genomicsdb

    script:
    """
    gatk GenomicsDBImport \
    --genomicsdb-workspace-path "${meta.study_id}_full_cohort_dbimport_data" \
    --batch-size 50 \
    ${chroms.collect{ f -> "-L ${f}" }.join(' ')} \
    --bypass-feature-reader true \
    --max-num-intervals-to-import-in-parallel 50 \
    --genomicsdb-shared-posixfs-optimizations true \
    --sample-name-map "${SAMPLE_MAP}" \
    --reader-threads 6
    """

    stub:
    """
    mkdir -p demo_full_cohort_dbimport_data
    echo stub > demo_full_cohort_dbimport_data/test.txt
    """

}

process CREATE_DICT {
    container "broadinstitute/gatk:4.2.6.1"
    input: 
        path(ref_genome)
    output:
        tuple path(ref_genome), path("*.dict"), path("*.fai"), emit: ref

    script:
    """
    gatk CreateSequenceDictionary -R $ref_genome
    samtools faidx $ref_genome
    """
    stub:
    """
    echo stub > genome.fai
    echo stub > genome.dict
    """
}

process GATK_GVCF_PER_CHROM {
    tag {CHR[0]}
    container "broadinstitute/gatk:4.2.6.1"
    publishDir "${params.outdir}/Raw_joint_call"
    input: 
    tuple val(meta), path(GENDB)
    tuple path(ref_genome_file), path(ref_genome_dict), path(reference_idx)
    each CHR

    output: 
    tuple val(CHR), path("*vcf.gz"), emit: chrom_vcf

    script:
    """
    gatk --java-options "-Xmx4g -Xms4g" GenotypeGVCFs \
       -R "${ref_genome_file}" \
       -V gendb://${GENDB} \
       -L ${CHR[0]} \
       -O ${CHR[0]}.vcf.gz \
       --genomicsdb-shared-posixfs-optimizations true 
    """

    stub:
    """
    echo stub > ${CHR[0]}.vcf.gz
    """

}

process MERGE_COHORT_VCF {
    container "broadinstitute/gatk:4.2.6.1"
    publishDir "${params.outdir}/Raw_joint_call"
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
    publishDir "${params.outdir}/Raw_joint_call"
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
    publishDir "${params.outdir}/Final_joint_call"
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
    publishDir "${params.outdir}/Final_joint_call"
    input: 
    tuple val(meta), path(raw_genotpye_vcf), path(raw_genotpye_index)
    path(baitset)

    output:
    tuple val(meta), path("*marked.vcf.gz"), path("*marked.vcf.gz.tbi"), emit: marked_variants
    
    script:
    def snp_filter  =  '-filter "QD < 2.0" --filter-name "QD2" ' +
                       '-filter "QUAL < 30.0" --filter-name "QUAL30" ' +
                       '-filter "SOR > 3.0" --filter-name "SOR3" ' +
                       '-filter "FS > 60.0" --filter-name "FS60" ' +
                       '-filter "MQ < 40.0" --filter-name "MQ40" ' +
                       '-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" ' +
                       ' -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"' 
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
    -O "${meta.study_id}${meta.suffix}.marked.vcf.gz"
    """
    stub:
    """
    echo stub > TBC.marked.vcf.gz
    echo stub > TBC.marked.vcf.gz.tbi
    """

}


process FILTER_VARIANTS {
    container "quay.io/biocontainers/bcftools:1.20--h8b25389_0"
    publishDir "${params.outdir}/Final_joint_call"
    input: 
    tuple val(meta), path(marked_genotpye_vcf), path(marked_genotpye_index)
    path(baitset)

    output:
    tuple val(meta), path("*.target.pass.vcf.gz"), path("*.target.pass.vcf.gz.tbi"), emit: filtered_variants
    script:
    """
    bcftools view -f PASS ${marked_genotpye_vcf} \
    -Oz -o "${meta.study_id}${meta.suffix}.marked.target.pass.vcf.gz" \
    -T ${baitset} 
    tabix -p vcf "${meta.study_id}${meta.suffix}.marked.target.pass.vcf.gz"
    """
    stub:
    """
    echo stub > demo.marked.target.pass.vcf.gz
    echo stub > demo.marked.target.pass.vcf.gz.tbi
    """

}


process ANNOTATE_VARIANTS {
    publishDir "${params.outdir}/Final_joint_call"
    container "ensemblorg/ensembl-vep:release_103.1"
    input:
    tuple val(meta), path(vcf_file)
    path(vep_config)
    path(gnomadfile)
    path(clinvarfile)
    path(dbsnpfile)
    path(cosmicfile)
    
    output: 
    tuple val(meta), path("*vep.vcf.gz"), emit: vep_annotation
    script: 
    def custom = "--custom $gnomadfile,gnomAD,vcf,exact,0,FLAG,AF --custom $clinvarfile,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT --custom $dbsnpfile,dbSNP,vcf,exact,0, --custom $cosmicfile,COSMIC,vcf,exact,0,CNT"
    def outfname = meta.id + "vep.vcf.gz"
    """
    vep -i ${vcf_file} \
    --config ${vep_config}
    --output_file TBC
    """
    stub: 
    """
    echo stub > item.vep.vcf.gz
    """
    }

process CONVERT_TO_TSV {
    input: 
    tuple val(meta), path(vep_vcf)

    output: 
    tuple val(meta), path("*.marked.target.pass.vep.tsv.gz")
    
    script:
    """
    zcat ${vep_vcf} | gatk_germline_full_vcf2table.v2.pl -> "${meta.study_id}_cohort_snps.marked.target.pass.vep.tsv"
    gzip "${meta.study_id}_cohort_snps.marked.target.pass.vep.tsv"
    """
}


process COMBINED_SUMMARY{
    publishDir "results", mode: 'copy'

    input: 
    tuple val(meta), path(VEP_SNP_TSVGZ), path(VEP_INDEL_TSVGZ)
    path(NIH_GERMLINE_TSV)
    path(CANCER_GENE_CENSUS)
    path(FLAG_GENES)


    output:
    path("*.tsv")
    path("*.xlsx")


    script:
    """
    get_XSLX_andfilt_COSMIC_and_clinvar_PATHVars.R \
    --study_id ${meta.study_id} \
    --input_snv_tsv ${VEP_SNP_TSVGZ} \
    --input_indel_tsv ${VEP_INDEL_TSVGZ} \
    --germ_pred_tsv ${NIH_GERMLINE_TSV} \
    --cgc_tsv ${CANCER_GENE_CENSUS} \
    --flags_tsv ${FLAG_GENES} \
    --outdir "."
    """
}

process CONVERT_TO_MAF {
    publishDir "results", mode: 'copy'

    input:
    tuple val(meta), path(SNP_INDEL_GNAMDF_TSVGZ), path(NIH_GERMLINE_TSV), path(CANCER_GENE_CENSUS)

    output: 
    path "*cohort_snps.indel.marked.target.pass.vep.gnmadF.maf.gz", emit: maf
    path("oncplots/*"), emit: oncoplots
    script: 
    """
    convert_germ_sumtsv_to_MAF_plot.R \
    --study_id ${meta.study_id} \
    --input_germ_tsv ${SNP_INDEL_GNAMDF_TSVGZ} \
    --germ_pred_tsv ${NIH_GERMLINE_TSV} \
    --cgc_tsv ${CANCER_GENE_CENSUS} \
    --outdir ${SUMTABDIR}
    """

}
