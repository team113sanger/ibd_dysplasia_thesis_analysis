library(readr)
library(dplyr)

# Read in data
meta <- read_tsv("metadata/final_metadata_qc_pass.tsv")
maf <- read_tsv("data/7100_3235-filtered_mutations_all_indepTum_keepPA.maf")

# Get sample lists
p53_samples <- maf |>
    filter(Hugo_Symbol == "TP53") |>
    pull(Tumor_Sample_Barcode)

progressors <- meta |>
    filter(group == "Progressor") |>
    pull(sanger_dna_id)

non_p53_progressors <- progressors[!progressors %in% p53_samples]

write_lines(non_p53_progressors, "metadata/sample_lists/no_p53_progressors.tsv")