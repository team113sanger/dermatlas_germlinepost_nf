#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GENERATE_GENOMICS_DB } from "./modules/generate_genomicsdb.nf"
include { CREATE_DICT } from "./modules/generate_gatk_ref.nf"
include { GATK_GVCF_PER_CHROM; MERGE_COHORT_VCF;INDEX_COHORT_VCF} from "./modules/gatk_variant_handling.nf"
include { PROCESS_VARIANT_SET as PROCESS_SNPS; PROCESS_VARIANT_SET as PROCESS_INDELS } from "./subworkflows/process_variant_type.nf"
include { GERMLINE } from "./subworkflows/call_variants.nf"
include { SETUP_CALLING_INPUTS } from "./subworkflows/setup_calling_inputs.nf"
include { GERMLINE_COHORT_ANALYSIS } from "./subworkflows/summarise_germline_analysis.nf"



workflow DERMATLAS_GERMLINE {
    main:
    // Setup parameters and variants
    log.info("Checking baitset path: ${params.baitset}")
    baitset = file(params.baitset, checkIfExists: true)
    log.info("Checking refgenome path: ${params.reference_genome}")
    reference_genome = file(params.reference_genome, checkIfExists: true)
    log.info("Checking VEP cache path: ${params.vep_cache}")
    vep_cache = file(params.vep_cache, checkIfExists: true)
    log.info("Checking alt transcript path: ${params.alternative_transcripts}")
    alternative_transcripts = file(params.alternative_transcripts, checkIfExists: true)

    nih_genes = file(params.nih_germline_resource, checkIfExists: true)
    cgc_genes = file(params.cancer_gene_census_resource, checkIfExists: true)
    flags = file(params.flag_genes, checkIfExists: true)


    custom_files = Channel.of(params.custom_files.split(';'))
        .map(it -> file(it, checkIfExists: true))
        .collect()
    log.info("Custom files exist")
    
    custom_args = Channel.of(params.custom_args.split(';'))
        .collect()
        .map { '--custom ' + it.join(' --custom ') }
    
    log.info("Checking chrom list path: ${params.chrom_list}")
    
    // Read chromosomes and create deterministic indexed channel
    chroms_sorted = Channel.fromPath(params.chrom_list, checkIfExists: true)
        .splitCsv(sep:"\t")
        .flatten()
        .toSortedList()
    
    // Create chromosome channel for GENERATE_GENOMICS_DB
    chroms = chroms_sorted.flatten().collect(flat: true)
    
    // Create indexed chromosome channel for GATK_GVCF_PER_CHROM  
    chrom_idx = chroms_sorted.flatten()

    CREATE_DICT(reference_genome)
    
    if (params.post_process_only){
        SETUP_CALLING_INPUTS()
        sample_map = SETUP_CALLING_INPUTS.out.sample_map
        db_ch = GENERATE_GENOMICS_DB(sample_map, chroms, SETUP_CALLING_INPUTS.out.vcf_ch)
    } else {
        GERMLINE(CREATE_DICT.out.ref, baitset)
        sample_map = GERMLINE.out.sample_map
        db_ch = GENERATE_GENOMICS_DB(sample_map, chroms, GERMLINE.out.vcf_ch)
    }

    log.info("Setups complete")
    GATK_GVCF_PER_CHROM(db_ch, 
                    CREATE_DICT.out.ref, 
                    chrom_idx)

    // Deterministic collection and sorting of VCF files in genomic order  
    gvcf_chrom_files = GATK_GVCF_PER_CHROM.out.chrom_vcf
        .collect()
        .map { vcf_tuples -> 
            // Read chromosome order from file
            def chrom_order = file(params.chrom_list).readLines()
            
            def sorted_vcfs = vcf_tuples.sort { a, b -> 
                chrom_order.indexOf(a[0]) <=> chrom_order.indexOf(b[0])
            }
            def vcf_files = sorted_vcfs.collect { it[1] }
            return [[study_id: params.study_id], vcf_files]
        }
    
    MERGE_COHORT_VCF(gvcf_chrom_files)
    INDEX_COHORT_VCF(MERGE_COHORT_VCF.out.cohort_vcf)


    INDEX_COHORT_VCF.out.indexed_cohort_vcf
        .multiMap { meta, file, index -> 
            snp_ch:   tuple(meta + ["variant_type": "SNP",   "suffix": "_cohort_raw_snps"], file, index)
            indel_ch: tuple(meta + ["variant_type": "INDEL", "suffix": "_cohort_indel_raw"], file, index)
        }
        .set { variant_multi }

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
        snp_conversion_ch = PROCESS_SNPS.out.annotated_vars
        indel_conversion_ch = PROCESS_INDELS.out.annotated_vars
        
        GERMLINE_COHORT_ANALYSIS(snp_conversion_ch, 
                                 indel_conversion_ch,
                                 sample_map,
                                 params.assembly,
                                 params.filter_col,
                                 nih_genes,
                                 cgc_genes,
                                 flags,
                                 alternative_transcripts)
    }
    
    emit: 
        indel_file = PROCESS_INDELS.out.annotated_vars
        snp_file = PROCESS_SNPS.out.annotated_vars
        results = GERMLINE_COHORT_ANALYSIS.out.results
}

workflow {
    DERMATLAS_GERMLINE()
}