library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)

maf <- read_tsv("data/variants/7100_3235-filtered_mutations_matched_allTum_keepPA.maf")
meta <- read_tsv("metadata/final_metadata_qc_pass.tsv")
sample_list <- read_lines("metadata/sample_lists/all_one_ppat.list")

# Subset MAF to APC, KRAS, TP53
genes_of_interest <- c("APC", "KRAS", "TP53")

maf_subset <- maf %>%
  filter(Hugo_Symbol %in% genes_of_interest)

# Join with metadata
maf_meta <- maf_subset %>%
  inner_join(meta, by = c("Tumor_Sample_Barcode" = "sanger_dna_id"))

# ---- Precursor samples ----
precursors <- maf_meta %>%
  filter(precursor_or_follow_up == "Precursor") %>%
  mutate(progression_status = case_when(
    group == "Progressor" ~ "Prog",
    group == "Non-progressor" ~ "N-Prog",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(progression_status))

p_pre <- ggplot(precursors, aes(x = progression_status, y = VAF_tum, fill = progression_status)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.6, width = 0.6) +
  facet_wrap(~ Hugo_Symbol, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
          strip.text = element_text(face ="italic")) +
  labs(
    x = NULL,
    y = "Variant Allele Frequency",
    title = "Precursor samples"
  )

# ---- Follow-up samples ----
followups <- maf_meta %>%
  filter(precursor_or_follow_up == "Follow up") %>%
  mutate(progression_status = case_when(
    group == "Progressor" ~ "Prog",
    group == "Non-progressor" ~ "N-Prog",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(progression_status))

p_follow <- ggplot(followups, aes(x = progression_status, y = VAF_tum, fill = progression_status)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.6, width = 0.6) +
  facet_wrap(~ Hugo_Symbol, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        strip.text = element_text(face ="italic")) +
  labs(
    x = NULL,
    y = "Variant Allele Frequency",
    title = "Follow-up samples"
  )

# ---- Combine vertically ----
combined <- p_pre / p_follow

ggsave("plots/variants/compare_gene_vafs.png",
       combined, dpi = 300, width = 8, height = 6)