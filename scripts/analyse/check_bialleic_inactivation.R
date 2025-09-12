library(readr)

cn_loh <- read_tsv("data/copy_number/cn_LOH_check.tsv")
maf <- read_tsv("data/variants/7100_3235-filtered_mutations_matched_allTum_keepPA.maf")

apc_loh <- cn_loh |>
    filter(gene == "APC") |>
    filter(cn %in% c("cn-LOH", "LOH-del", "LOH-amp")) |>
    select(sample, gene, cn) |>
    write_tsv("results/tables/apc_cn_loh.tsv")

apc_maf <- maf |>
    filter(Hugo_Symbol == "APC") |>
    select(Tumor_Sample_Barcode, Main_consequence_VEP, Variant_Type, HGVSp_Short, VAF_tum) |>
    rename(sample = Tumor_Sample_Barcode)

combine <- apc_loh |>
    left_join(apc_maf) |>
    arrange(Main_consequence_VEP) |>
    write_tsv("results/tables/bialleic_inactivation_apc.tsv")