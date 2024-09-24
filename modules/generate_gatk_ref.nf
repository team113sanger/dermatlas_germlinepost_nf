
process CREATE_DICT {
    container "broadinstitute/gatk:4.2.6.1"
    label "gatk_steps"
    input: 
        path(ref_genome)
    output:
        tuple path(ref_genome), path("*.dict"), path("*.fai"), emit: ref

    script:
    """
    gatk CreateSequenceDictionary -R $ref_genome
    samtools faidx $ref_genome
    """
    stub:
    """
    echo stub > genome.fai
    echo stub > genome.dict
    """
}