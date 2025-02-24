

include { SELECT_VARIANTS;MARK_VARIANTS;FILTER_VARIANTS} from "../modules/gatk_variant_handling.nf"
include { ANNOTATE_VARIANTS } from "../modules/vep_annotation.nf"
include { CONVERT_TO_TSV } from "../modules/summarise_results.nf"

workflow PROCESS_VARIANT_SET {
    take:
    variant_ch
    baitset
    vep_cache
    custom_files
    custom_args
    reference_genome
    species
    assembly
    db_version

    main:
    SELECT_VARIANTS(variant_ch)
    MARK_VARIANTS(SELECT_VARIANTS.out.raw_variants, baitset)
    FILTER_VARIANTS(MARK_VARIANTS.out.marked_variants, baitset)
    ANNOTATE_VARIANTS(FILTER_VARIANTS.out.filtered_variants, 
                      vep_cache, 
                      custom_files,
                      custom_args,
                      reference_genome,
                      species,
                      assembly,
                      db_version)

    CONVERT_TO_TSV(ANNOTATE_VARIANTS.out.vep_annotation)

    emit:
    publish_vars = CONVERT_TO_TSV.out.tsv_file 
    annotated_vars = ANNOTATE_VARIANTS.out.vep_annotation
}