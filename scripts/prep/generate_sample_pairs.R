library(readr)
library(dplyr)

metadata <- read_tsv("metadata/final_metadata_qc_pass.tsv")

metadata <- metadata |>
  mutate(patient_id = sub("([A-Za-z]+\\d+)[a-zA-Z]$", "\\1", sanger_dna_id))

filtered_metadata <- metadata |>
  group_by(patient_id) |>
  filter(all(c("Precursor", "Follow up") %in% precursor_or_follow_up)) |>
  ungroup() |>
  filter(!sanger_dna_id %in% c("PD62031d", "PD62037a", "PD62038a",
                               "PD62064c", "PD62065e", "PD62075e"))

write_tsv(filtered_metadata, "metadata/sample_pairs.tsv")

