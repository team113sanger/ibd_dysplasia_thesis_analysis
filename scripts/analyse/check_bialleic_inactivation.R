library(readr)

maf <- read_tsv("data/variants/7100_3235-filtered_mutations_matched_allTum_keepPA.maf")

meta <- read_tsv("metadata/final_metadata_qc_pass.tsv") |>
    select(sanger_dna_id, study_id, group)
cn_loh <- read_tsv("data/copy_number/cn_LOH_check.tsv") |>
    left_join(meta, by = c(sample = "sanger_dna_id"))


apc_loh <- cn_loh |>
    filter(gene == "APC") |>
    filter(cn %in% c("cn-LOH", "LOH-del", "LOH-amp")) |>
    select(sample, study_id, gene, cn, group) |>
    write_tsv("results/tables/apc_cn_loh.tsv")

apc_maf <- maf |>
    filter(Hugo_Symbol == "APC") |>
    select(Tumor_Sample_Barcode, Main_consequence_VEP, Variant_Type, HGVSp_Short, VAF_tum) |>
    rename(sample = Tumor_Sample_Barcode)

combine_apc <- apc_loh |>
    left_join(tp53_maf) |>
    arrange(Main_consequence_VEP)

tp53_loh <- cn_loh |>
    filter(gene == "TP53") |>
    filter(cn %in% c("cn-LOH", "LOH-del", "LOH-amp")) |>
    select(sample, study_id, gene, cn, group) |>
    write_tsv("results/tables/tp53_cn_loh.tsv")

tp53_maf <- maf |>
    filter(Hugo_Symbol == "TP53") |>
    select(Tumor_Sample_Barcode, Main_consequence_VEP, Variant_Type, HGVSp_Short, VAF_tum) |>
    rename(sample = Tumor_Sample_Barcode)

combine_tp53 <- tp53_loh |>
    left_join(tp53_maf) |>
    arrange(Main_consequence_VEP)

combined_df <- bind_rows(combine_apc, combine_tp53) |>
    write_tsv("results/tables/bialleic_inactivation.tsv")
