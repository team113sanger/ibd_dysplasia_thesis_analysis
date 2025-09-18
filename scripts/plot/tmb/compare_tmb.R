library(maftools)
library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)

# Load metadata
metadata <- read_tsv("metadata/final_metadata_qc_pass.tsv") |>
  select(sanger_dna_id, group, precursor_or_follow_up)

# Load sample list
samples <- read_tsv("metadata/sample_lists/all_one_ppat.list", col_names = "Sample") |>
  left_join(metadata, by = c(Sample = "sanger_dna_id"))

# Load TMB
tmb <- read_tsv("data/variants/mutations_per_Mb.tsv", col_names = c("Sample", "TMB")) |>
  right_join(samples, by = c(Sample = "Sample")) |>
  mutate(
    facet_group = case_when(
      group == "Non-progressor" & precursor_or_follow_up == "Precursor" ~ "N-Prog Pre",
      group == "Non-progressor" & precursor_or_follow_up == "Follow up" ~ "N-Prog Fol",
      group == "Progressor"     & precursor_or_follow_up == "Precursor" ~ "Prog Pre",
      group == "Progressor"     & precursor_or_follow_up == "Follow up" ~ "Prog Fol",
      TRUE ~ NA_character_
    ),
    facet_group = factor(facet_group, label = c("N-Prog Pre", "N-Prog Fol", "Prog Pre", "Prog Fol")),
    group = case_when(
        group =="Non-progressor" ~ "NP",
        group == "Progressor" ~ "P",
        TRUE ~ group
    )
  )


p <- ggplot(tmb, aes(x = group, y = TMB, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
  facet_wrap(~precursor_or_follow_up) +
  theme_bw(base_size = 10) +
  scale_y_log10() +
  #scale_fill_brewer(palette = "Accent") +
  scale_fill_manual(labels = c("P" = "Progressor", "NP" = "Non-progressor"),
                    values = c("P" = "chocolate", "NP" = "darkseagreen")) +
  labs(
    x = NULL,
    y = "TMB (mut/Mb)",
    fill = NULL
  ) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.border = element_blank()
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("NP", "P")),
    label = "p.signif",
    hide.ns = FALSE,
    size = 3
  )

ggsave("plots/tmb/compare_group_tmb.png", p, width = 4, height = 3.5, dpi = 300)
