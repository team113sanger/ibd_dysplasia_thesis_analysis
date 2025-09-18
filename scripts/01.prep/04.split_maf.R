library(readr)
library(dplyr)

# For running dndscv

maf <- read_tsv("data/variants/7100_3235-filtered_mutations_matched_allTum_keep.maf")

prog_pre_list <- read_lines("metadata/sample_lists/progressor_precursor_samples_ppat.tsv")
n_prog_pre_list <- read_lines("metadata/sample_lists/non_progressor_precursor_samples_ppat.tsv")
prog_fol_list <- read_lines("metadata/sample_lists/progressor_follow_up_samples_ppat.tsv")
n_prog_fol_list <- read_lines("metadata/sample_lists/non_progressor_follow_up_samples_ppat.tsv")

prog_pre_maf <- maf |>
  filter(Tumor_Sample_Barcode %in% prog_pre_list)
write_tsv(prog_pre_maf, "data//variants/split_mafs/progressors_precursors.maf")

n_prog_pre_maf <- maf |>
  filter(Tumor_Sample_Barcode %in% n_prog_pre_list)
write_tsv(n_prog_pre_maf, "data/variants/split_mafs/non_progressors_precursors.maf")


# TP53
maf <- read_tsv("data/variants/7100_3235-filtered_mutations_matched_allTum_keepPA.maf")

prog_pre_maf <- maf |>
  filter(Tumor_Sample_Barcode %in% prog_pre_list)
write_tsv(prog_pre_maf, "data/variants/split_mafs/prog_pre_keepPA.maf")

n_prog_pre_maf <- maf |>
  filter(Tumor_Sample_Barcode %in% n_prog_pre_list)
write_tsv(n_prog_pre_maf, "data/variants/split_mafs/n_prog_pre_keepPA.maf")

prog_fol_maf <- maf |>
  filter(Tumor_Sample_Barcode %in% prog_fol_list)
write_tsv(prog_fol_maf, "data/variants/split_mafs/prog_fol_keepPA.maf")

n_prog_fol_maf <- maf |>
  filter(Tumor_Sample_Barcode %in% n_prog_fol_list)
write_tsv(n_prog_fol_maf, "data/variants/split_mafs/n_prog_fol_keepPA.maf")



