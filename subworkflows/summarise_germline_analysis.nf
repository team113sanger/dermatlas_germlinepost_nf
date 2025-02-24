include { VCF_TO_KEEP_MAF; VCF_TO_PASS_MAF } from "../modules/vcf_to_maf.nf"
include { FILTER_AND_ONCOPLOT } from "../modules/summarise_results.nf"

workflow GERMLINE_COHORT_ANALYSIS {
    take: 
    snp_conversion_ch
    indel_conversion_ch
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

    VCF_TO_KEEP_MAF(snp_conversion_ch, 
                    indel_conversion_ch,
                    vcf_list,
                    assembly,
                    filter_col,
                    alternative_transcripts)

    VCF_TO_PASS_MAF(snp_conversion_ch, 
                    indel_conversion_ch,
                    vcf_list,
                    assembly,
                    filter_col,
                    alternative_transcripts)

    FILTER_AND_ONCOPLOT(VCF_TO_KEEP_MAF.out.maf, 
                        vcf_list, 
                        nih_germline_resource, 
                        cancer_gene_census_resource, 
                        flag_genes)
}