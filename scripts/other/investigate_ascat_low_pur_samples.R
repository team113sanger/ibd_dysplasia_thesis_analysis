library(readr)
library(dplyr)

low_pur <- read_tsv("data/qc/ascat_low_purity.txt", col_names = "sanger_dna_id")
meta <- read_tsv("metadata/final_metadata_qc_pass.tsv")

compare <- low_pur |>
    left_join(meta) |>
    select(sanger_dna_id, study_id, group, precursor_or_follow_up, grade_of_dysplasia, diagnosis)

write_tsv(compare, "results/ascat_low_purity_cases.tsv")