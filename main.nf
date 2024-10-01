#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GENERATE_GENOMICS_DB } from "./modules/generate_genomicsdb.nf"
include { CREATE_DICT } from "./modules/generate_gatk_ref.nf"
include { GATK_GVCF_PER_CHROM;MERGE_COHORT_VCF;INDEX_COHORT_VCF} from "./modules/gatk_variant_handling.nf"
include { PROCESS_VARIANT_SET as PROCESS_SNPS; PROCESS_VARIANT_SET as PROCESS_INDELS } from "./subworkflows/process_variant_type.nf"
include { NF_DEEPVARIANT } from "./subworkflows/hgi_nfdeepvariant.nf"
include { GERMLINE } from "./subworkflows/germline.nf"
include { POSTPROCESS_ONLY } from "./subworkflows/postprocess_only.nf"

include { COMBINED_SUMMARY;CONVERT_TO_MAF } from "./modules/summarise_results.nf"

workflow {

    baitset = file(params.baitset, checkIfExists: true)
    reference_genome = file(params.reference_genome, checkIfExists: true)
    vep_cache = file(params.vep_cache, checkIfExists: true)
    vep_config = file(params.vep_config, checkIfExists: true)
    gnomad_file = file(params.gnomad_file, checkIfExists: true)
    dbsnp_file = file(params.dbsnp_file, checkIfExists: true)
    clinvar_file = file(params.clinvar_file, checkIfExists: true)
    cosmic_file = file(params.cosmic_file, checkIfExists: true)
    nih_germline_resource =   file(params.nih_germline_resource, checkIfExists: true)
    cancer_gene_census_resource = file(params.cancer_gene_census_resource, checkIfExists: true)
    flag_genes =  file(params.flag_genes, checkIfExists: true)
    
    chroms = Channel.fromPath("$baseDir/assets/grch38_chromosome.txt")
    | splitCsv(sep:"\t")
    | collect(flat: true)
    chrom_idx = chroms.withIndex()

    CREATE_DICT(reference_genome)
    
    if (params.post_process_only){
        POSTPROCESS_ONLY()
        sample_map = POSTPROCESS_ONLY.out.sample_map
        db_ch = GENERATE_GENOMICS_DB(sample_map, chroms, POSTPROCESS_ONLY.out.vcf_ch)
    } else {
        GERMLINE()
        sample_map = GERMLINE.out.sample_map
        db_ch = GENERATE_GENOMICS_DB(sample_map, chroms, GERMLINE.out.vcf_ch)
    }

    GATK_GVCF_PER_CHROM(db_ch, 
                    CREATE_DICT.out.ref, 
                    chrom_idx)
    gvcf_chrom_files = GATK_GVCF_PER_CHROM.out.chrom_vcf
    | toSortedList { item -> item[0][1]}
    | transpose()
    | last()
    | map {it -> tuple([study_id: params.study_id], it)}
    
    MERGE_COHORT_VCF(gvcf_chrom_files)
    INDEX_COHORT_VCF(MERGE_COHORT_VCF.out.cohort_vcf)


    INDEX_COHORT_VCF.out.indexed_cohort_vcf
    | multiMap {meta, file, index -> 
       snp_ch:   tuple(meta + ["variant_type": "SNP",   "suffix": "_cohort_raw_snps"], file, index)
       indel_ch: tuple(meta + ["variant_type": "INDEL", "suffix": "_cohort_indel_raw"], file, index)}
    | set {variant_multi}

    PROCESS_SNPS(variant_multi.snp_ch,
                 baitset,
                 vep_cache,
                 vep_config,
                 gnomad_file,
                 dbsnp_file,
                 clinvar_file,
                 cosmic_file,
                 reference_genome)
    
    PROCESS_INDELS(variant_multi.indel_ch,
                   baitset,
                   vep_cache,
                   vep_config,
                   gnomad_file,
                   dbsnp_file,
                   clinvar_file,
                   cosmic_file,
                   reference_genome)
    
    COMBINED_SUMMARY(PROCESS_SNPS.out.publish_vars,
                     PROCESS_INDELS.out.publish_vars,
                    nih_germline_resource,
                    cancer_gene_census_resource,
                    flag_genes)
    CONVERT_TO_MAF(COMBINED_SUMMARY.out.outfile,
                        nih_germline_resource,
                        cancer_gene_census_resource)
}