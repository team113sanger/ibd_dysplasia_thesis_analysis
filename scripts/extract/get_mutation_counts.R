library(readr)
library(dplyr)

maf <- read_tsv("data/variants/7100_3235-filtered_mutations_matched_allTum_keep.maf")
mafPA <- read_tsv("data/variants/7100_3235-filtered_mutations_matched_allTum_keepPA.maf")

samples <- read_lines("metadata/sample_lists/all_one_ppat.list")
length(samples)

maf_split <- maf |>
    filter(Tumor_Sample_Barcode %in% samples)

mafPA_split <- mafPA |>
    filter(Tumor_Sample_Barcode %in% samples)

muts_count <- nrow(maf_split)
muts_count

PA_muts_count <- nrow(mafPA_split)
PA_muts_count