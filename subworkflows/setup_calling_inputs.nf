workflow SETUP_CALLING_INPUTS {

    main:
    log.info("Checking Sample map path: ${params.sample_map}")
    sample_map = file(params.sample_map, checkIfExists: true)

    Channel.fromPath(sample_map)
    .splitCsv(sep:"\t", header: ["sample", "file"])
    .map{ meta  ->
    def basename = new File(meta.file).getName()
    [ meta.sample, basename ]}
    .collectFile{ meta -> ["sample_base.txt", "${meta[0]}\t${meta[1]}\n"]}
    | set {base_map}
    log.info("Checking GVCF path: ${params.geno_vcf}")
    Channel.fromPath(params.geno_vcf)
    .map { file -> 
            def index = file + ".tbi"
            tuple(file, index)}
     .collect()
     .map { file_list -> tuple([study_id: params.study_id], file_list)}
     .set{ vcf_ch }


    
    emit: 
    vcf_ch = vcf_ch
    sample_map = base_map

}