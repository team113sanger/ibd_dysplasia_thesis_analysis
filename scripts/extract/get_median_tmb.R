library(readr)
library(dplyr)

tmb <- read_tsv("data/variants/mutations_per_Mb.tsv", col_names = c("sample", "tmb")) |>
  left_join(read_tsv("metadata/final_metadata_qc_pass.tsv"), by = c("sample" = "sanger_dna_id")) |>
  filter(sample %in% read_lines("metadata/sample_lists/all_one_ppat.list"))

# Median TMB for LGD Non-Progressors
median_LGD_nonprog <- tmb |>
  filter(grade_of_dysplasia == "Low grade", group == "Non-progressor") |>
  summarise(median_tmb = median(tmb)) |>
  pull(median_tmb)

# Median TMB for LGD Progressors
median_LGD_prog <- tmb |>
  filter(grade_of_dysplasia == "Low grade", group == "Progressor") |>
  summarise(median_tmb = median(tmb)) |>
  pull(median_tmb)

# Median TMB for HGD Progressors
median_HGD_prog <- tmb |>
  filter(grade_of_dysplasia == "High grade", group == "Progressor") |>
  summarise(median_tmb = median(tmb)) |>
  pull(median_tmb)

# Median TMB for AC Progressors
median_AC_prog <- tmb |>
  filter(grade_of_dysplasia == "Adenocarcinoma", group == "Progressor") |>
  summarise(median_tmb = median(tmb)) |>
  pull(median_tmb)

# Print results
median_LGD_nonprog
median_LGD_prog
median_HGD_prog
median_AC_prog
