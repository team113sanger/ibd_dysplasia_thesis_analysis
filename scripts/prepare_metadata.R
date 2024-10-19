library(readr)
library(janitor)
library(stringr)
library(dplyr)

# Load metadata
progressors <- read_tsv("metadata/Progressors.tsv") |>
    clean_names() |>
    mutate(group = "Progressor") |>
    select(-macrodissection_coring)

non_progressors <- read_tsv("metadata/Non-progressors.tsv") |>
    clean_names()|>
    mutate(group = "Non-progressor")|>
    select(-notes)

meta <- read_tsv("metadata/7100_metadata.tsv") |>
    clean_names() |>
    select(case_id, country, diagnosis, phenotype, tissue, tumour_type,
            sex, age, location_of_normal, sanger_dna_id) |>
            rename(study_id = case_id)

# Combine metadata
meta_temp <- rbind(progressors, non_progressors)
meta <- left_join(meta, meta_temp)

# Tidy
# meta_tidy <- meta |>
#   mutate(grade_of_dysplasia = case_when(
#     str_detect(grade_of_dysplasia, regex("low grade", ignore_case = TRUE)) ~ "Low grade dysplasia",
#     str_detect(grade_of_dysplasia, regex("high grade", ignore_case = TRUE)) ~ "High grade dysplasia",
#     str_detect(grade_of_dysplasia, regex("adenocarcinoma", ignore_case = TRUE)) ~ "Adenocarcinoma",
#     TRUE ~ grade_of_dysplasia
#   ))

write_tsv(meta, "metadata/processed_metadata.tsv")




