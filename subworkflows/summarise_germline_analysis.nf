include { VCF_TO_MAF as VCF_TO_KEEP_MAF} from "../modules/vcf_to_maf.nf"
include { VCF_TO_MAF as VCF_TO_PASS_MAF} from "../modules/vcf_to_maf.nf"
include { FILTER_AND_ONCOPLOT } from "../modules/summarise_results.nf"

workflow GERMLINE_COHORT_ANALYSIS {
    take: 
    snp_conversion_ch
    indel_conversion_ch
    sample_map
    assembly
    filter_col
    nih_germline_resource
    cancer_gene_census_resource
    flag_genes
    alternative_transcripts
    
    main:
    snp_conversion_ch
    .concat(indel_conversion_ch)
    .map { meta, file, index -> file.baseName + ".gz" }
    .collectFile(name: 'vcf_list.txt', newLine: true)
    .set {vcf_list}

    Channel.fromPath(sample_map)
    .splitCsv(sep:"\t")
    .map { row -> row[0]}
    .collectFile(name: 'sample_list.txt', newLine: true)
    .set {sample_list}


    snp_conversion_ch
    .map { meta, file, index -> tuple(meta + ['variant_set': 'keep'], file, index) }
    .set { keep_vars }

    snp_conversion_ch
    .map { meta, file, index -> tuple(meta + ['variant_set': 'pass'], file, index) }
    .set { pass_vars }

    
    VCF_TO_KEEP_MAF(keep_vars, 
                    indel_conversion_ch,
                    vcf_list,
                    assembly,
                    filter_col,
                    alternative_transcripts)

    VCF_TO_PASS_MAF(pass_vars, 
                    indel_conversion_ch,
                    vcf_list,
                    assembly,
                    filter_col,
                    alternative_transcripts)

    FILTER_AND_ONCOPLOT(VCF_TO_KEEP_MAF.out.maf, 
                        nih_germline_resource, 
                        cancer_gene_census_resource, 
                        flag_genes, 
                        sample_list)
    emit: 
        FILTER_AND_ONCOPLOT.out
}