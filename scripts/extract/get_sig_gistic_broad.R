library(readr)
library(dplyr)
library(tidyr)

broad_gistic_pre <- read_tsv("data/copy_number/gistic/prog_pre/MIN_0/broad_significance_results.txt")
broad_gistic_fol <- read_tsv("data/copy_number/gistic/prog_fol/MIN_0/broad_significance_results.txt")

sig_table <- broad_gistic_pre |>
  pivot_longer(
    cols = c(`Amp q-value`, `Del q-value`),
    names_to = "cn_type",
    values_to = "q_value"
  ) |>
  mutate(cn_type = if_else(cn_type == "Amp q-value", "Amplification", "Deletion")) |>
  filter(q_value < 0.1) |>
  select(Arm, cn_type, q_value)

sig_table


sig_table <- broad_gistic_fol |>
  pivot_longer(
    cols = c(`Amp q-value`, `Del q-value`),
    names_to = "cn_type",
    values_to = "q_value"
  ) |>
  mutate(cn_type = if_else(cn_type == "Amp q-value", "Amplification", "Deletion")) |>
  filter(q_value < 0.1) |>
  select(Arm, cn_type, q_value)

sig_table
