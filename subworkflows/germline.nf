include { NF_DEEPVARIANT } from "../subworkflows/hgi_nfdeepvariant.nf"
workflow GERMLINE {
    take: 
    ref_genome
    baitset

    main:
    channel_inputs_bams = Channel.fromPath(params.tsv_file)
    .splitCsv(header: true, sep: '\t')
    .map{row->tuple(row.sample, row.object, row.object_index)}
    .take(params.samples_to_process)

    NF_DEEPVARIANT(channel_inputs_bams, ref_genome, baitset)
    NF_DEEPVARIANT.out.haplos.collectFile(name: "tmp_sample_map.txt"){
        meta, vcf, index ->
    ["tmp_sample_map.txt", "${meta}\t${vcf.baseName}.gz\n"]
    }
    .set{ sample_map } 
    NF_DEEPVARIANT.out.haplos
    .map{ sample_id, file, index -> tuple(file, index)}
    .collect()
    .map { file_list -> tuple([study_id: params.study_id], file_list)}
    .set{ vcf_ch }

    
    emit: 
    vcf_ch = vcf_ch
    sample_map = sample_map

}