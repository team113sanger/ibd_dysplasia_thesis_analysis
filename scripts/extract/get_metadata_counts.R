library(readr)
library(dplyr)

sample_pairs <- read_tsv("metadata/sample_pairs.tsv")

# Compute male-to-female ratio per group
sex_ratio_table <- sample_pairs %>%
  distinct(patient_id, group, sex) %>%  # Count each patient only once
  group_by(group, sex) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = sex, values_from = count, values_fill = 0) %>%
  mutate(male_to_female_ratio = Male / Female)

# View result
print(sex_ratio_table)

metadata <- read_tsv("metadata/final_metadata_qc_pass.tsv")

library(dplyr)
library(tidyr)
library(stringr)

sex_ratio_precursor <- metadata %>%
  # Updated: Extract patient ID from study_id (handles both CP and CNP)
  mutate(patient_id = str_extract(study_id, "^[0-9]+CNP?|^[0-9]+CP")) %>%

  # Filter for precursor lesions only
  filter(precursor_or_follow_up == "Precursor") %>%

  # Keep one record per patient
  distinct(patient_id, group, sex) %>%

  # Count sexes per group
  group_by(group, sex) %>%
  summarise(count = n(), .groups = "drop") %>%

  # Pivot and calculate ratio
  pivot_wider(names_from = sex, values_from = count, values_fill = 0) %>%
  mutate(male_to_female_ratio = Male / Female)

# View result
print(sex_ratio_precursor)
