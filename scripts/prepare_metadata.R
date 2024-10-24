library(readr)
library(janitor)
library(stringr)
library(dplyr)

# Load metadata
progressors <- read_tsv("metadata/progressors_raw.tsv") |>
  clean_names() |>
  mutate(group = "Progressor") |>
  select(-macrodissection_coring)

non_progressors <- read_tsv("metadata/non-progressors_raw.tsv") |>
  clean_names() |>
  mutate(group = "Non-progressor") |>
  select(-notes)

meta <- read_tsv("metadata/sanger_metadata.tsv") |>
  clean_names() |>
  select(
    case_id, sanger_dna_id, phenotype, diagnosis,
    sex, age
  ) |>
  rename(study_id = case_id)

qc_pass_samples <- read_tsv("metadata/sample_lists/qc_pass_sample_list.tsv", col_names = F) |>
  pull(X1)

# Combine metadata
meta_temp <- rbind(progressors, non_progressors)
meta_combined <- left_join(meta, meta_temp) |>
  select(
    study_id, sanger_dna_id, phenotype, precursor_or_follow_up,
    grade_of_dysplasia, diagnosis, ibd_diagnosis,
    disease_duration_to_first_dyaplastic_lesion,
    time_between_lesions, site, psc, family_history, pancolitis,
    sex, age, additional_info_of_note, group
  )

# Filter for only samples that passed QC
meta_pass <- meta_combined |>
  filter(sanger_dna_id %in% qc_pass_samples)

# Tidy columns
meta_tidy <- meta_pass |>
  mutate(precursor_or_follow_up = case_when(
    str_detect(precursor_or_follow_up, regex("precursor", ignore_case = TRUE)) ~ "Precursor",
    str_detect(precursor_or_follow_up, regex("follow up", ignore_case = TRUE)) ~ "Follow up",
    TRUE ~ precursor_or_follow_up
  )) |>
  mutate(grade_of_dysplasia = case_when(
    str_detect(grade_of_dysplasia, regex("low grade", ignore_case = TRUE)) ~ "Low grade",
    str_detect(grade_of_dysplasia, regex("LG", ignore_case = TRUE)) ~ "Low grade",
    str_detect(grade_of_dysplasia, regex("NOS", ignore_case = TRUE)) ~ "NOS",
    str_detect(grade_of_dysplasia, regex("high grade", ignore_case = TRUE)) ~ "High grade",
    str_detect(grade_of_dysplasia, regex("adenocarcinoma", ignore_case = TRUE)) ~ "Adenocarcinoma",
    TRUE ~ grade_of_dysplasia
  )) |>
  mutate(ibd_diagnosis = case_when(
    str_detect(ibd_diagnosis, regex("crohn'?s", ignore_case = TRUE)) ~ "Crohn's",
    TRUE ~ ibd_diagnosis
  )) |>
  mutate(
    psc = str_to_title(psc),
    family_history = str_to_title(family_history),
    pancolitis = str_to_title(pancolitis)
  )

write_tsv(meta_tidy, "metadata/final_metadata_qc_pass.tsv")
