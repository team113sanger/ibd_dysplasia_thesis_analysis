library(dplyr)
library(readr)

maf <- read_tsv("data/variants/7100_3235-filtered_mutations_all_indepTum_keepPA.maf")

problems(maf)

genes <- maf |>
  filter(Tumor_Sample_Barcode == "PD62027c") |>
  pull(Hugo_Symbol)

write_lines(genes, "data/genes.txt")
