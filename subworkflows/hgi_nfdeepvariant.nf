include { sort_cram } from "../modules/sort_cram.nf"
include { markDuplicates } from "../modules/markDuplicates.nf"
include { coord_sort_cram } from "../modules/coord_sort_cram.nf"
include { bam_to_cram } from "../modules/bam_to_cram.nf"
include { deepvariant } from "../modules/deepvariant.nf"
include { gatk_haplotypecaller } from "../modules/gatk_haplotypecaller.nf"

workflow NF_DEEPVARIANT {
    take:
    channel_inputs_bams
    ref_genome
    baitset

    main:
        sort_cram(channel_inputs_bams)
        markDuplicates(sort_cram.out.sorted_sample_cram, ref_genome)
        coord_sort_cram(markDuplicates.out.markdup_sample_cram)
        bam_to_cram(coord_sort_cram.out.markdup_sample_cram_crai, ref_genome)
        gatk_haplotypecaller(coord_sort_cram.out.markdup_sample_cram_crai, ref_genome, baitset)

    
    emit:
    haplos = gatk_haplotypecaller.out
}
