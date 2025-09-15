library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)

# Read data
p53_status <- read_tsv("results/p53_mutations/p53_status.tsv")
cn_props <- read_tsv("data/copy_number/proportions/all_cn_props.tsv") |>
  select(Sample, proportion) |>
  rename(sanger_dna_id = Sample)

# Sample lists
nprog_pre <- read_lines("metadata/sample_lists/non_progressor_precursor_samples_ppat.tsv")
prog_pre  <- read_lines("metadata/sample_lists/progressor_precursor_samples_ppat.tsv")
nprog_fol <- read_lines("metadata/sample_lists/non_progressor_follow_up_samples_ppat.tsv")
prog_fol  <- read_lines("metadata/sample_lists/progressor_follow_up_samples_ppat.tsv")

# Functions
subset_groups <- function(df, set1, set2, set1_label, set2_label, timepoint) {
  df %>%
    filter(sanger_dna_id %in% c(set1, set2)) %>%
    mutate(group = ifelse(sanger_dna_id %in% set1, set1_label, set2_label),
           timepoint = timepoint)
}

add_TP53_status <- function(df, p53_df) {
  df %>%
    mutate(p53 = ifelse(sanger_dna_id %in% p53_df$Tumor_Sample_Barcode, "Mut", "WT"))
}

# Build combined dataset
df_all <- bind_rows(
  subset_groups(cn_props, nprog_pre, prog_pre, "Non-Progressor", "Progressor", "Precursor"),
  subset_groups(cn_props, nprog_fol, prog_fol, "Non-Progressor", "Progressor", "Follow-up")
) %>%
  add_TP53_status(p53_status)

df_all[["timepoint"]] <- factor(
  df_all[["timepoint"]],
  levels = c("Precursor", "Follow-up")
)
# Wilcoxon per facet (by group now)
wilcox_labels <- df_all %>%
  group_by(timepoint) %>%
  summarise(
    p = wilcox.test(proportion ~ group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = paste0("Wilcox p = ", signif(p, 3)))

# Boxplot/violin
p <- ggplot(df_all, aes(x = group, y = proportion, fill = group)) +
    geom_violin(width = 1, alpha = 0.8, linewidth = 0.3) +
    geom_boxplot(color = "#383838ff", alpha = 0.15, linewidth = 0.3,
                 outlier.shape = NA, width = 0.5) +
    geom_jitter(fill = "#585858", size = 1, stroke = 0, alpha = 0.4,
                width = 0.1, height = 0.1) +
    labs(x = NULL, y = "CNA Proportion") +
    scale_fill_manual(values = c("Progressor" = "darkorange",
                                 "Non-Progressor" = "darkseagreen")) +
    theme_classic(base_size = 10) +
    facet_wrap(~timepoint) +
    geom_text(
      data = wilcox_labels,
      aes(x = 1.5, y = 1.05 * max(df_all$proportion), label = label),
      inherit.aes = FALSE,
      size = 3,
      fontface = "italic") +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 20, hjust = 1),
      axis.ticks.x = element_blank())

# Counts plot
p2 <- ggplot(df_all, aes(x = group, fill = p53)) +
        geom_bar(position = "fill", width =0.5) +   
        labs(x = "Group", y = "Proportion", fill = "TP53 Status") +
        scale_fill_viridis_d(end = 0.8) +
        theme_classic(base_size = 10) +
        facet_wrap(~timepoint) +
        theme(
          legend.position = "right",
          strip.text = element_blank())

# Arrange 
p_combined <- p / p2 + plot_layout(heights = c(3, 1))

ggsave("plots/copy_number/proportions/group_cn_props.png",
       p_combined, width = 6, height = 5)

df_progs <- df_all |>
    filter(group == "Progressor") |>
    arrange(desc(proportion)) |>
    arrange(timepoint) |>
    write_tsv("results/tables/progs_cn_status.tsv")