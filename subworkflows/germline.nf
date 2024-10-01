include { NF_DEEPVARIANT } from "../subworkflows/hgi_nfdeepvariant.nf"
workflow GERMLINE {
    
    channel_inputs_bams = Channel.fromPath(params.tsv_file)
    .splitCsv(header: true, sep: '\t')
    .map{row->tuple(row.sample, row.object, row.object_index)}
    .take(params.samples_to_process)

    NF_DEEPVARIANT(channel_inputs_bams, params.run_mode)
    NF_DEEPVARIANT.out.gatk_haplotypecaller_out.collectFile(
        name: "tmp_sample_map.txt")
        {
        meta, vcf, index ->
    ["tmp_sample_map.txt", "${meta}\t${vcf.baseName}\n"]}
    .set{ sample_map } 
    
    
    emit: 
    vcf_ch = NF_DEEPVARIANT.out.gatk_haplotypecaller_out
    sample_map = sample_map

}