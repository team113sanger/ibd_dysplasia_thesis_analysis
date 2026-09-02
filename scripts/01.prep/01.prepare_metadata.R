library(readr)
library(janitor)
library(stringr)
library(dplyr)

# Load metadata
progressors <- read_tsv("metadata/progressors_raw.tsv") |>
  clean_names() |>
  mutate(group = "Progressor") |>
  select(-macrodissection_coring)
  #mutate(study_id = str_replace(study_id, "^16CPC-1$", "16CPC")) #Check later is this is correct (16CPC in sanger metadata, 16CPC123 in other)

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
meta_combined <- left_join(meta, meta_temp, by = "study_id") |>
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
    str_detect(grade_of_dysplasia, regex("adenocaricnoma", ignore_case = TRUE)) ~ "Adenocarcinoma",
    str_detect(grade_of_dysplasia, regex("moderate", ignore_case = TRUE)) ~ "Low grade", # Check later is okay 
    str_detect(grade_of_dysplasia, regex("not specified", ignore_case = TRUE)) ~ "NOS", #Check this is correct
    TRUE ~ grade_of_dysplasia
  )) |>
  mutate(ibd_diagnosis = case_when(
    str_detect(ibd_diagnosis, regex("crohn'?s", ignore_case = TRUE)) ~ "Crohn's",
    TRUE ~ ibd_diagnosis
  )) |>
  mutate(
    psc = str_to_title(psc),
    family_history = str_to_title(family_history),
    pancolitis = str_to_title(pancolitis),
    site = str_to_title(site)
  ) |>
  mutate(site_general = case_when(
        site %in% c("Caecum", "Ascending", "Proximal Ascending", "Distal Ascending") ~ "Right colon",
        site %in% c("Hepatic Flexure", "Transverse", "Splenic Flexure") ~ "Transverse colon",
        site %in% c("Descending", "Sigmoid", "Rectosigmoid") ~ "Left colon",
        site %in% c("Rectum") ~ "Rectum",
        TRUE ~ "Other"
      ))

# Remove extra PD62028a/3CPA case in 'non-progressor' group
# Remove N-prog lesion from patient with lesions in both N-prog & prog groups
# 37CNPA/PD62039d & 37CNPB/PD62039e
# Remove samples that are now LGD in precurosr (update from Shahida on 19/06/2026)
# 5CNP, 11CNP, 12CNP, 18CNPB-1, 18CNPB-2, 27CNPA, 30CNPA, 30CNPB, 34CNPA, 34CNPB
# 5CPA, 5CPB, 5CPC, 7CPA, 26CPA, 26CPB
meta_filtered <- meta_tidy |>
  filter(!(sanger_dna_id == "PD62028a" & group == "Non-progressor")) |>
  filter(!(sanger_dna_id %in% c("PD62039d", "PD62039e"))) |>
  filter(!(study_id %in% c("18CNPB-1", "18CNPB-2", "30CNPA", "30CNPB", "34CNPA", "34CNPB"))) |>
  filter(!(study_id %in% c("5CPA", "5CPB", "5CPC")))

# Fix samples with missing metadata (not in prog/non-prog sheets)
# PD62045c, PD62041d 
meta_edit <- meta_filtered |>
  mutate(
    precursor_or_follow_up = if_else(
      sanger_dna_id %in% c("PD62045c", "PD62041d"), "Follow up", precursor_or_follow_up
    ),
    grade_of_dysplasia = if_else(
      sanger_dna_id %in% c("PD62045c", "PD62041d"), "High grade", grade_of_dysplasia # Need to follow up if PD62041d really is high grade
    ),
    group = if_else(
      sanger_dna_id %in% c("PD62045c", "PD62041d"), "Progressor", group
    )
  ) |>
  # Fix incorrect labels in metadata
  mutate(precursor_or_follow_up = case_when(
    sanger_dna_id == "PD62064d" ~ "Follow up",
    sanger_dna_id == "PD62064c" ~ "Precursor",
    TRUE ~ precursor_or_follow_up
  )) |>
  # Move high grade sample in NP to P 
  mutate(group = case_when(
    sanger_dna_id == "PD62082c" ~ "Progressor", # check later why this was in NP
    TRUE ~ group
  ))

# Add patient ID column
meta_fin <- meta_edit |>
  mutate(patient_id = sub("([A-Za-z]+\\d+)[a-zA-Z]$", "\\1", sanger_dna_id)) |>
  # Separate independent lesions from the same patient
  mutate(patient_id = case_when(
    sanger_dna_id %in% c("PD62037c", "PD62037d") ~ paste0(patient_id, "x"),
    sanger_dna_id %in% c("PD62038c", "PD62038d") ~ paste0(patient_id, "x"),
    TRUE ~ patient_id
  ))

write_tsv(meta_fin, "metadata/final_metadata_qc_pass.tsv")
