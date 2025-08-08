include { SORT_CRAM } from "../modules/sort_cram.nf"
include { MARK_DUPLICATES } from "../modules/markDuplicates.nf"
include { COORD_SORT_CRAM } from "../modules/coord_sort_cram.nf"
include { BAM_TO_CRAM } from "../modules/bam_to_cram.nf"
include { DEEPVARIANT } from "../modules/deepvariant.nf"
include { GATK_HAPLOTYPECALLER } from "../modules/gatk_haplotypecaller.nf"

workflow NF_DEEPVARIANT {
    take:
    channel_inputs_bams
    ref_genome
    baitset

    main:
        SORT_CRAM(channel_inputs_bams)
        MARK_DUPLICATES(SORT_CRAM.out.sorted_sample_cram, ref_genome)
        COORD_SORT_CRAM(MARK_DUPLICATES.out.markdup_sample_cram)
        BAM_TO_CRAM(COORD_SORT_CRAM.out.markdup_sample_cram_crai, ref_genome)
        GATK_HAPLOTYPECALLER(COORD_SORT_CRAM.out.markdup_sample_cram_crai, ref_genome, baitset)

    
    emit:
    haplos = GATK_HAPLOTYPECALLER.out
}
