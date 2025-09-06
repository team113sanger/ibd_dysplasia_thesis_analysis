library(readr)
library(dplyr)

metadata <- read_tsv("metadata/final_metadata_qc_pass.tsv")

metadata <- metadata |>
  mutate(patient_id = sub("([A-Za-z]+\\d+)[a-zA-Z]$", "\\1", sanger_dna_id))

sample_pairs <- metadata |>
  group_by(patient_id) |>
  filter(all(c("Precursor", "Follow up") %in% precursor_or_follow_up)) |>
  ungroup() |>
  # Remove extra cases from the same patient
  filter(!sanger_dna_id %in% c(
    "PD62030c", "PD62031d", "PD62033c", "PD62041c", "PD62047g",
    "PD62064c", "PD62065a", "PD62065e", "PD62068a", "PD62068c",
    "PD62068d", "PD62068e", "PD62075c", "PD62075e"
  )) |>
  # Separate independent lesions from the same patient
  mutate(patient_id = case_when(
    sanger_dna_id %in% c("PD62037c", "PD62037d") ~ paste0(patient_id, "x"),
    sanger_dna_id %in% c("PD62038c", "PD62038d") ~ paste0(patient_id, "x"),
    TRUE ~ patient_id
  ))

write_tsv(sample_pairs, "metadata/sample_pairs.tsv")

n_prog_pre <- sample_pairs |>
  filter(group == "Non-progressor",
        precursor_or_follow_up == "Precursor") |>
        pull(sanger_dna_id)
write_lines(n_prog_pre, "metadata/sample_pairs/n_prog_pre.tsv")


n_prog_fol <- sample_pairs |>
  filter(group == "Non-progressor",
        precursor_or_follow_up == "Follow up") |>
        pull(sanger_dna_id)
write_lines(n_prog_fol, "metadata/sample_pairs/n_prog_fol.tsv")

prog_pre <- sample_pairs |>
  filter(group == "Progressor",
        precursor_or_follow_up == "Precursor") |>
        pull(sanger_dna_id)
write_lines(prog_pre, "metadata/sample_pairs/prog_pre.tsv")


prog_fol <- sample_pairs |>
  filter(group == "Progressor",
        precursor_or_follow_up == "Follow up") |>
        pull(sanger_dna_id)
write_lines(prog_fol, "metadata/sample_pairs/prog_fol.tsv")
