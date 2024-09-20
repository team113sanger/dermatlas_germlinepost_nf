process GENERATE_GENOMICS_DB {
    container "broadinstitute/gatk:4.2.6.1"
    publishDir "${params.outdir}", mode: "copy"
    
    input: 
    path(SAMPLE_MAP)
    val(chroms)
    tuple val(meta), path(genotype_vcfs)

    output: 
    tuple val(meta), path("*data"), emit: genomicsdb

    script:
    """
    gatk GenomicsDBImport \
    --genomicsdb-workspace-path "${meta.study_id}_full_cohort_dbimport_data" \
    --batch-size 50 \
    ${chroms.collect{ f -> "-L ${f}" }.join(' ')} \
    --bypass-feature-reader true \
    --max-num-intervals-to-import-in-parallel 50 \
    --genomicsdb-shared-posixfs-optimizations true \
    --sample-name-map "${SAMPLE_MAP}" \
    --reader-threads 6
    """

    stub:
    """
    mkdir -p demo_full_cohort_dbimport_data
    echo stub > demo_full_cohort_dbimport_data/test.txt
    """

}
