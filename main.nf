#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { GENERATE_GENOMICS_DB;GATK_GVCF_PER_CHROM;CREATE_DICT;MERGE_COHORT_VCF;INDEX_COHORT_VCF;SELECT_SNP_VARIANTS;SELECT_INDEL_VARIANTS;
          MARK_SNP_VARIANTS; MARK_INDEL_VARIANTS;FILTER_SNP_VARIANTS;FILTER_INDEL_VARIANTS} from './subworkflows/germline_post.nf'

workflow {
    sample_map = file(params.sample_map, checkIfExists: true)
    baitset = file(params.baitset, checkIfExists: true)
    reference_genome = file(params.reference_genome, checkIfExists: true)
    // reference_idx = file(params.reference_idx, checkIfExists: true)
    
    chroms = Channel.fromPath("$baseDir/assets/grch38_chromosome.txt")
    | splitCsv(sep:"\t")
    
    
    chrom_list = chroms
    | map { it -> "-L " + it[0] + " "}
    | reduce{a, b -> a + b}
    

    chrom_vals = chroms.map{it[0]}
    
    Channel.fromPath(params.geno_vcf)
    .map { file -> 
            index = file + ".tbi"
            tuple(file, index)}
     .collect()
     .map { file_list -> tuple([study_id: params.study_id], file_list)}
     .set{ vcf_ch }
    
    GENERATE_GENOMICS_DB(sample_map, chrom_list, vcf_ch)
    CREATE_DICT(reference_genome)
    GATK_GVCF_PER_CHROM(GENERATE_GENOMICS_DB.out.genomicsdb, 
                        CREATE_DICT.out.ref, 
                        chrom_vals)
    gvcf_chrom_files = GATK_GVCF_PER_CHROM.out.chrom_vcf
    .groupTuple()
    .view()
    // .map{meta, it -> it}
    // .collect()
    // .map { file_list -> tuple([study_id: params.study_id], file_list)}
    // 


    MERGE_COHORT_VCF(gvcf_chrom_files)
    INDEX_COHORT_VCF(MERGE_COHORT_VCF.out.cohort_vcf)
    
    SELECT_SNP_VARIANTS(INDEX_COHORT_VCF.out.indexed_cohort_vcf)
    SELECT_INDEL_VARIANTS(INDEX_COHORT_VCF.out.indexed_cohort_vcf)

    MARK_SNP_VARIANTS(SELECT_SNP_VARIANTS.out.raw_variants,baitset)
    MARK_INDEL_VARIANTS(SELECT_INDEL_VARIANTS.out.raw_variants, baitset)

    // FILTER_SNP_VARIANTS(MARK_SNP_VARIANTS.out.marked_variants)
    // FILTER_INDEL_VARIANTS(MARK_INDEL_VARIANTS.out.marked_variants)
    
    // ANNOTATE_VARIANTS(FILTER_SNP_VARIANTS.out.filtered_variants)
    // //ANNOTATE_VARIANTS(MARK_INDEL_VARIANTS.out.filtered_variants)

    // CONVERT_TO_TSV(ANNOTATE_VARIANTS.out.vep_annotation)
}