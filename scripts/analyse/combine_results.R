library(tidyr)
library(readr)
library(dplyr)

# Read in data
cn_clusters <- read_tsv("results/copy_number/precursors/cn_clusters.tsv") |>
    select(sample, `2`) |>
    rename(cn_cluster = `2`)
cn_prop <- read_tsv("data/copy_number/proportions/all_cn_props.tsv") |>
    select(Sample, proportion)
tp53_status <- read_tsv("results/p53_mutations/p53_status.tsv")
meta <- read_tsv("metadata/final_metadata_qc_pass.tsv")
pre_samples <- c(
    read_lines("metadata/sample_lists/non_progressor_precursor_samples_ppat.tsv"),
    read_lines("metadata/sample_lists/progressor_precursor_samples_ppat.tsv"))
fol_samples <- c(
    read_lines("metadata/sample_lists/non_progressor_follow_up_samples_ppat.tsv"),
    read_lines("metadata/sample_lists/progressor_follow_up_samples_ppat.tsv"))

# Combine data 
df_pre <- meta |>
  select(sanger_dna_id, group) |>
  filter(sanger_dna_id %in% pre_samples) |>
  mutate(
    TP53_status = if_else(sanger_dna_id %in% tp53_status$Tumor_Sample_Barcode, "Mut", "WT")
  ) |>
  left_join(
    cn_clusters |> rename(sanger_dna_id = sample),
    by = "sanger_dna_id"
  ) |>
  left_join(
    cn_prop |> rename(sanger_dna_id = Sample),
    by = "sanger_dna_id"
  ) |>
  filter(!sanger_dna_id %in% c("PD62028a", "PD62077a")) # Remove samples with no CN results

write_tsv(df_pre, "results/precursor_combined_results.tsv")

  df_fol <- meta |>
  select(sanger_dna_id, group) |>
  filter(sanger_dna_id %in% fol_samples) |>
  mutate(
    TP53_status = if_else(sanger_dna_id %in% tp53_status$Tumor_Sample_Barcode, "Mut", "WT")
  ) |>
  left_join(
    cn_clusters |> rename(sanger_dna_id = sample),
    by = "sanger_dna_id"
  ) |>
  left_join(
    cn_prop |> rename(sanger_dna_id = Sample),
    by = "sanger_dna_id"
  )

write_tsv(df_fol, "results/follow_up_combined_results.tsv")