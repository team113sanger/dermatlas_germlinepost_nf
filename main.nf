#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { GENERATE_GENOMICS_DB;GATK_GVCF_PER_CHROM;
          CREATE_DICT;MERGE_COHORT_VCF;INDEX_COHORT_VCF;
          SELECT_VARIANTS as SELECT_SNP_VARIANTS;SELECT_VARIANTS as SELECT_INDEL_VARIANTS;
          MARK_SNP_VARIANTS; MARK_INDEL_VARIANTS;
          FILTER_VARIANTS as FILTER_SNP_VARIANTS;FILTER_VARIANTS as FILTER_INDEL_VARIANTS} from './subworkflows/germline_post.nf'

workflow {
    sample_map = file(params.sample_map, checkIfExists: true)
    baitset = file(params.baitset, checkIfExists: true)
    reference_genome = file(params.reference_genome, checkIfExists: true)
    
    chroms = Channel.fromPath("$baseDir/assets/grch38_chromosome.txt")
    | splitCsv(sep:"\t")
    | collect(flat: true)
    chrom_idx = chroms.withIndex()
    
    Channel.fromPath(params.geno_vcf)
    .map { file -> 
            index = file + ".tbi"
            tuple(file, index)}
     .collect()
     .map { file_list -> tuple([study_id: params.study_id], file_list)}
     .set{ vcf_ch }
    
    GENERATE_GENOMICS_DB(sample_map, chroms, vcf_ch)
    CREATE_DICT(reference_genome)
    GATK_GVCF_PER_CHROM(GENERATE_GENOMICS_DB.out.genomicsdb, 
                        CREATE_DICT.out.ref, 
                        chrom_idx)

    
    gvcf_chrom_files = GATK_GVCF_PER_CHROM.out.chrom_vcf
    | toSortedList { item -> item[0][1]}
    | transpose()
    | last()
    | map {it -> tuple([study_id: params.study_id], it)}
    
    MERGE_COHORT_VCF(gvcf_chrom_files)
    INDEX_COHORT_VCF(MERGE_COHORT_VCF.out.cohort_vcf)
    // INDEX_COHORT_VCF.out.

    INDEX_COHORT_VCF.out.indexed_cohort_vcf
    | map {meta, file, index -> 
       tuple(meta + ["variant_type": "SNP", "suffix": "_cohort_raw_snps"], file, index)}
    | set {snp_ch}

    INDEX_COHORT_VCF.out.indexed_cohort_vcf
    | map {meta, file, index -> 
    tuple(meta + ["variant_type": "INDEL", "suffix": "_cohort_indel_raw"], file, index)}
    | set {indel_ch}


    
    SELECT_SNP_VARIANTS(snp_ch)
    SELECT_INDEL_VARIANTS(indel_ch)

    MARK_SNP_VARIANTS(SELECT_SNP_VARIANTS.out.raw_variants, baitset)
    MARK_INDEL_VARIANTS(SELECT_INDEL_VARIANTS.out.raw_variants, baitset)

    FILTER_SNP_VARIANTS(MARK_SNP_VARIANTS.out.marked_variants, baitset)
    FILTER_INDEL_VARIANTS(MARK_INDEL_VARIANTS.out.marked_variants,  baitset)
    
    // ANNOTATE_VARIANTS(FILTER_SNP_VARIANTS.out.filtered_variants)
    // //ANNOTATE_VARIANTS(MARK_INDEL_VARIANTS.out.filtered_variants)

    // CONVERT_TO_TSV(ANNOTATE_VARIANTS.out.vep_annotation)
}