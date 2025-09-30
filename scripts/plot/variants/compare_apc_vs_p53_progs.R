library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)


maf <- read_tsv("data/variants/7100_3235-filtered_mutations_matched_allTum_keepPA.maf")
meta <- read_tsv("metadata/final_metadata_qc_pass.tsv")
sample_list <- read_lines("metadata/sample_lists/all_one_ppat.list")

ibd_crc_ids <- c("5CPA", "14CPA", "16CPA", "17CPA", "18CPA", "21CPA", "25CPB", "27CPA-1", "3CPA")
sporadic_crc_ids <- c("4CPB", "8CPA", "9CPA", "22CPA")

meta_groups <- meta %>%
  mutate(sample_group = case_when(
    study_id %in% ibd_crc_ids ~ "IBD-CRC",
    study_id %in% sporadic_crc_ids ~ "Sporadic-CRC",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(sample_group)) %>%
  # clean disease_duration_to_first_dyaplastic_lesion
  mutate(disease_duration_num = disease_duration_to_first_dyaplastic_lesion %>%
           # remove "years" and "year"
           gsub("years?", "", .) %>%
           # remove "At least"
           gsub("At least", "", .) %>%
           # trim whitespace
           trimws() %>%
           # convert to numeric
           as.numeric())

p1 <- ggplot(meta_groups, aes(x = sample_group, y = disease_duration_num, fill = sample_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, width = 0.6) +
  geom_jitter(width = 0.1, alpha = 0.7, size = 1) +
  geom_text_repel(aes(label = study_id), size = 2.5, max.overlaps = Inf) +  # add labels
  scale_fill_brewer(palette = "Set3") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    x = NULL,
    y = "Disease duration"
  )

# Subset APC mutations
maf_apc <- maf %>%
  filter(Hugo_Symbol == "APC")

# Join to metadata with groups
maf_apc_meta <- maf_apc %>%
  inner_join(meta_groups, by = c("Tumor_Sample_Barcode" = "sanger_dna_id"))

# APC VAFs plot
p2 <- ggplot(maf_apc_meta, aes(x = sample_group, y = VAF_tum, fill = sample_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, width = 0.6) +
  geom_jitter(width = 0.1, alpha = 0.7, size = 1) +
  geom_text_repel(aes(label = study_id), size = 2.5, max.overlaps = Inf) +  # add labels
  scale_fill_brewer(palette = "Accent") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    x = NULL,
    y = "APC VAF"
  )

combined <- p1 + p2

ggsave("plots/variants/ibd_crc_vs_sporadic.png",
       combined, dpi = 300, width = 6, height = 3)
