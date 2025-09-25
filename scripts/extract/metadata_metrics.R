library(readr)
library(stringr)
library(dplyr)

sample_list <- read_lines("metadata/sample_lists/all_one_ppat.list")
meta <- read_tsv("metadata/final_metadata_qc_pass.tsv") |>
            mutate(patient_id = str_remove(sanger_dna_id, "[a-z]$")) |> #keep to remove multiple lesions per patient
            filter(sanger_dna_id %in% sample_list)

patient_coint <- meta |>
    distinct(patient_id) |>
    count()
patient_coint

lesions_per_patient <- meta |>
  group_by(patient_id) |>
  summarise(n_lesions = n(), .groups = "drop")
lesions_per_patient

ibd_forms_count <- meta |>
    distinct(patient_id, .keep_all = TRUE) |>
    group_by(ibd_diagnosis) |>
    count()
ibd_forms_count

site_count <- meta |>
    distinct(patient_id, .keep_all = TRUE) |>
    count(site_general)
site_count

sex_count <-  meta |>
    distinct(patient_id, .keep_all = TRUE) |>
    group_by(sex) |>
    count()
sex_count

sex_count_by_group <-  meta |>
    distinct(patient_id, .keep_all = TRUE) |>
    group_by(sex, group) |>
    count()
sex_count_by_group

age_metrics <- meta |>
  group_by(group) |>
  summarise(
    median_age = median(age, na.rm = TRUE),
    min_age = min(age, na.rm = TRUE),
    max_age = max(age, na.rm = TRUE),
    .groups = "drop"
  )

age_metrics