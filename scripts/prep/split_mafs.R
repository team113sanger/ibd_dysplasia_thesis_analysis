library(readr)
library(dplyr)

maf <- read_tsv("data/variants/7100_3235-filtered_mutations_matched_indepTum_keep.maf")

prog_pre_list <- read_lines("metadata/sample_lists/progressor_precursor_samples_ppat.tsv")
n_prog_pre_list <- read_lines("metadata/sample_lists/non_progressor_precursor_samples_ppat.tsv")

prog_pre_maf <- maf|>
    filter(Tumor_Sample_Barcode %in% prog_pre_list)
write_tsv(prog_pre_maf, "data/mafs/progressors_precursors.maf")

n_prog_pre_maf <- maf|>
    filter(Tumor_Sample_Barcode %in% n_prog_pre_list)
write_tsv(n_prog_pre_maf, "data/mafs/non_progressors_precursors.maf")