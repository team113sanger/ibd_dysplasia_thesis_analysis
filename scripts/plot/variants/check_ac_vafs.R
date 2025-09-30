library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)

maf <- read_tsv("data/variants/7100_3235-filtered_mutations_matched_allTum_keep.maf")
#maf <- read_tsv("data/variants/7100_3235-filtered_mutations_matched_allTum_keepPA.maf")
meta <- read_tsv("metadata/final_metadata_qc_pass.tsv")
sample_list <- read_lines("metadata/sample_lists/all_one_ppat.list")

# 1. Clean metadata
meta_clean <- meta %>%
  filter(!grade_of_dysplasia %in% c("NOS")) %>%                # remove NOS
  filter(!(group == "Non-progressor" & precursor_or_follow_up == "Follow up")) %>% # remove non-progressor follow-ups
  group_by(patient_id) %>%
  # keep only precursor sample if multiple LGD from same patient
  filter(!(grade_of_dysplasia == "Low grade" & precursor_or_follow_up == "Follow up")) %>%
  ungroup()

# 2. Join maf and metadata
maf_meta <- maf %>%
  inner_join(meta_clean, by = c("Tumor_Sample_Barcode" = "sanger_dna_id")) |>
    mutate(
    grade_of_dysplasia = case_when(
      grade_of_dysplasia == "Low grade" ~ "LGD",
      grade_of_dysplasia == "High grade" ~ "HGD",
      grade_of_dysplasia == "Adenocarcinoma" ~ "AC",
      grade_of_dysplasia == "NOS" ~ "NOS",
      TRUE ~ NA_character_
    )) |>
  mutate(grade_of_dysplasia = factor(grade_of_dysplasia, levels = c("LGD", "HGD", "AC")))

# 3. Boxplot
p <- ggplot(maf_meta, aes(x = grade_of_dysplasia, y = VAF_tum, fill = grade_of_dysplasia)) +
        geom_boxplot(outlier.shape = 21, alpha = 0.6, width = 0.6) +
        #geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
        scale_fill_brewer(palette = "Set2") +
        theme_bw() +
        theme(legend.position = "none") +
        labs(
            x = "Grade of Dysplasia",
            y = "Variant Allele Frequency"
        )
ggsave("plots/variants/compare_ac_vafs.png", p, dpi = 300, width = 4, height = 4)

# 4. Adenocarcinoma samples separately
p2 <- maf_meta %>%
  filter(grade_of_dysplasia == "AC") %>%
  ggplot(aes(x = study_id, y = VAF_tum, fill = Tumor_Sample_Barcode)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.6, width = 0.6) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6)
  ) +
  labs(
    x = "Adenocarcinoma Samples",
    y = "Variant Allele Frequency"
  )

# Combine with patchwork, giving p2 more space
combined <- p + p2 + plot_layout(widths = c(1, 1.8))  # p narrower, p2 wider

# Save combined plot
ggsave("plots/variants/compare_ac_vafs_per_sample.png",
       combined, dpi = 300, width = 8, height = 3.5)
