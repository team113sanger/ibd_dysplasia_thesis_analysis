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
    group = case_when(
        group =="Non-progressor" ~ "NP",
        group == "Progressor" ~ "P",
        TRUE ~ group
    ),
    precursor_or_follow_up = factor(precursor_or_follow_up, label = c("Precursor", "Follow up")),

  )

n_counts <- tmb |>
  group_by(group, precursor_or_follow_up) |>
  summarise(n = n(), .groups = "drop")

p <- ggplot(tmb, aes(x = group, y = TMB, fill = group)) +
        geom_violin(alpha = 0.3, color = NA) +  
        geom_boxplot(alpha = 0.6, width = 0.5) +  
        geom_jitter(width = 0.15, size = 0.5, alpha = 0.7) +
        facet_wrap(~precursor_or_follow_up) +
        geom_text(
          data = n_counts,
          aes(x = group, y = 30, label = paste0("n=", n)),
          inherit.aes = FALSE,
          vjust = 1.5,
          color = "black",
          size = 3
        ) +
        theme_bw(base_size = 10) +
        scale_y_log10() +
        scale_fill_manual(
          labels = c("P" = "Progressor", "NP" = "Non-progressor"),
          values = c("P" = "chocolate", "NP" = "darkseagreen")
        ) +
        labs(
          x = NULL,
          y = "TMB (Mut/Mb)",
          fill = NULL
        ) +
        theme(
          legend.position = "bottom",
          legend.text = element_text(size = 8),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold")
        ) +
        stat_compare_means(
          method = "wilcox.test",
          comparisons = list(c("NP", "P")),
          hide.ns = FALSE,
          size = 3
        )

ggsave("plots/tmb/tmb_by_group.png", p, width = 4, height = 4, dpi = 300)
