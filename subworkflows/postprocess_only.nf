workflow POSTPROCESS_ONLY {

    main:
    Channel.fromPath(params.geno_vcf)
    .map { file -> 
            index = file + ".tbi"
            tuple(file, index)}
     .collect()
     .map { file_list -> tuple([study_id: params.study_id], file_list)}
     .set{ vcf_ch }


    sample_map = file(params.sample_map, checkIfExists: true)

    Channel.fromPath(sample_map)
    .splitCsv(sep:"\t", header: ["sample", "file"])
    .map{ meta  ->
    def basename = new File(meta.file).getName()
    [ meta.sample, basename ]}
    .collectFile{ meta -> ["sample_base.txt", "${meta[0]}\t${meta[1]}\n"]}
    | set {base_map}
    
    Channel.fromPath(params.geno_vcf)
    .map { file -> 
            index = file + ".tbi"
            tuple(file, index)}
     .collect()
     .map { file_list -> tuple([study_id: params.study_id], file_list)}
     .set{ vcf_ch }


    
    emit: 
    vcf_ch = vcf_ch
    sample_map = base_map

}