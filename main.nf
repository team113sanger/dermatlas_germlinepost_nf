#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { GENERATE_GENOMICS_DB;GATK_GVCF_PER_CHROM;
          CREATE_DICT;MERGE_COHORT_VCF;INDEX_COHORT_VCF;
          SELECT_VARIANTS as SELECT_SNP_VARIANTS;SELECT_VARIANTS as SELECT_INDEL_VARIANTS;
          MARK_VARIANTS as MARK_SNP_VARIANTS; MARK_VARIANTS as  MARK_INDEL_VARIANTS;
          FILTER_VARIANTS as FILTER_SNP_VARIANTS;FILTER_VARIANTS as FILTER_INDEL_VARIANTS;
          ANNOTATE_VARIANTS as ANNOTATE_SNP_VARIANTS;ANNOTATE_VARIANTS as ANNOTATE_INDEL_VARIANTS;
          CONVERT_TO_TSV as CONVERT_SNP_TO_TSV;CONVERT_TO_TSV as CONVERT_INDEL_TO_TSV;
          COMBINED_SUMMARY; CONVERT_TO_MAF} from './subworkflows/germline_post.nf'

workflow {
    sample_map = file(params.sample_map, checkIfExists: true)
    baitset = file(params.baitset, checkIfExists: true)
    reference_genome = file(params.reference_genome, checkIfExists: true)
    vep_cache = file(params.vep_cache, checkIfExists: true)
    vep_config = file(params.vep_config, checkIfExists: true)
    gnomad_file = file(params.gnomad_file, checkIfExists: true)
    dbsnp_file = file(params.dbsnp_file, checkIfExists: true)
    clinvar_file = file(params.clinvar_file, checkIfExists: true)
    cosmic_file = file(params.cosmic_file, checkIfExists: true)
    nih_germline_resource =   file(params.nih_germline_resource, checkIfExists: true)
    cancer_gene_census_resoruce = file(params.cancer_gene_census_resoruce, checkIfExists: true)
    flag_genes =  file(params.flag_genes, checkIfExists: true)
    
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
    
    ANNOTATE_SNP_VARIANTS(FILTER_SNP_VARIANTS.out.filtered_variants, 
                      vep_cache, 
                      vep_config,
                      gnomad_file,
                      dbsnp_file,
                      clinvar_file,
                      cosmic_file,
                      reference_genome)
    ANNOTATE_INDEL_VARIANTS(FILTER_INDEL_VARIANTS.out.filtered_variants, 
                      vep_cache, 
                      vep_config,
                      gnomad_file,
                      dbsnp_file,
                      clinvar_file,
                      cosmic_file,
                      reference_genome)

    CONVERT_SNP_TO_TSV(ANNOTATE_SNP_VARIANTS.out.vep_annotation)
    CONVERT_INDEL_TO_TSV(ANNOTATE_INDEL_VARIANTS.out.vep_annotation)
    COMBINED_SUMMARY(CONVERT_SNP_TO_TSV.out.tsv_file,
                     CONVERT_INDEL_TO_TSV.out.tsv_file,
                    nih_germline_resource,
                    cancer_gene_census_resoruce,
                    flag_genes)
    CONVERT_TO_MAF(COMBINED_SUMMARY.out.outfile,
                        nih_germline_resource,
                        cancer_gene_census_resoruce)
}