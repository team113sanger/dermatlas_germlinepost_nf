#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GENERATE_GENOMICS_DB } from "./modules/generate_genomicsdb.nf"
include { CREATE_DICT } from "./modules/generate_gatk_ref.nf"
include { GATK_GVCF_PER_CHROM; MERGE_COHORT_VCF;INDEX_COHORT_VCF} from "./modules/gatk_variant_handling.nf"
include { PROCESS_VARIANT_SET as PROCESS_SNPS; PROCESS_VARIANT_SET as PROCESS_INDELS } from "./subworkflows/process_variant_type.nf"
include { NF_DEEPVARIANT } from "./subworkflows/hgi_nfdeepvariant.nf"
include { GERMLINE } from "./subworkflows/germline.nf"
include { POSTPROCESS_ONLY } from "./subworkflows/postprocess_only.nf"

include { COMBINED_SUMMARY;CONVERT_TO_MAF } from "./modules/summarise_results.nf"

workflow DERMATLAS_GERMLINE {
    main:
    baitset = file(params.baitset, checkIfExists: true)
    reference_genome = file(params.reference_genome, checkIfExists: true)
    vep_cache = file(params.vep_cache, checkIfExists: true)
    
    custom_files = Channel.of(params.custom_files.split(';'))
    .map(it -> file(it, checkIfExists: true))
    .collect()
    custom_args = Channel.of(params.custom_args.split(';'))
    .collect()
    .map { '--custom ' + it.join(' --custom ') }


    
    chroms = Channel.fromPath(params.chrom_list)
    | splitCsv(sep:"\t")
    | collect(flat: true)
    chrom_idx = chroms.flatten().toList().map { it.withIndex() }

    CREATE_DICT(reference_genome)
    
    if (params.post_process_only){
        POSTPROCESS_ONLY()
        sample_map = POSTPROCESS_ONLY.out.sample_map
        db_ch = GENERATE_GENOMICS_DB(sample_map, chroms, POSTPROCESS_ONLY.out.vcf_ch)
    } else {
        GERMLINE(CREATE_DICT.out.ref, baitset)
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
    | map {it -> [[study_id: params.study_id], it]}
    
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
                 custom_files,
                 custom_args,
                 reference_genome,
                 params.species,
                 params.assembly,
                 params.db_version)
    
    PROCESS_INDELS(variant_multi.indel_ch,
                   baitset,
                   vep_cache,
                   custom_files,
                   custom_args,
                   reference_genome, 
                   params.species,
                   params.assembly,
                   params.db_version)
    
    if (params.summarise_results){
    nih_germline_resource = file(params.nih_germline_resource, checkIfExists: true)
    cancer_gene_census_resource = file(params.cancer_gene_census_resource, checkIfExists: true)
    flag_genes =  file(params.flag_genes, checkIfExists: true)
    
    
    PROCESS_SNPS.out.publish_vars.join(PROCESS_INDELS.out.publish_vars)
    .collectFile(storeDir: "${params.outdir}/${params.release_version}"){
        meta, file, index -> 
        new File("${params.outdir}/summary").mkdirs()
        def filename = "sample_list.tsv"
        [filename, "${meta["sample_id"]}\n"]} 
    | set {sample_list}

    VCF_TO_KEEP_MAF(sample_list)
    VCF_TO_PASS_MAF(sample_list)

                    nih_germline_resource,
                    cancer_gene_census_resource,
                    flag_genes)
    CONVERT_TO_MAF(COMBINED_SUMMARY.out.outfile,
                        nih_germline_resource,
                        cancer_gene_census_resource)
    }
    emit:
    indel_file = PROCESS_INDELS.out.annotated_vars
    snp_file = PROCESS_SNPS.out.annotated_vars
    
    summary_files = COMBINED_SUMMARY.out.outfile
    cohort_maf = CONVERT_TO_MAF.out.maf

}

workflow {
    DERMATLAS_GERMLINE()
}