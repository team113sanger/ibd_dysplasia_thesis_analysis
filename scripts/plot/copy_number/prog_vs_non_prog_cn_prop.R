library(readr)
library(dplyr)
library(ggplot2)

# Read data
cn_props <- read_tsv("data/copy_number/proportions/all_cn_props.tsv") %>%
  rename(sanger_dna_id = Sample)

# Sample lists
nprog_pre <- read_lines("metadata/sample_lists/non_progressor_precursor_samples_ppat.tsv")
prog_pre  <- read_lines("metadata/sample_lists/progressor_precursor_samples_ppat.tsv")
nprog_fol <- read_lines("metadata/sample_lists/non_progressor_follow_up_samples_ppat.tsv")
prog_fol  <- read_lines("metadata/sample_lists/progressor_follow_up_samples_ppat.tsv")

# Functions
subset_groups <- function(df, set1, set2, set1_label = "Non-Progressor", set2_label = "Progressor") {
  df %>%
    filter(sanger_dna_id %in% c(set1, set2)) %>%
    mutate(group = ifelse(sanger_dna_id %in% set1, set1_label, set2_label))
}

perform_wilcox <- function(df){
  wilcox_result <- wilcox.test(proportion ~ group, data = df)
  p_label <- paste0("p = ", signif(wilcox_result$p.value, 3))
  return(p_label)
}

plot_boxplot <- function(df, p_label){
  ggplot(df, aes(x = group, y = proportion, fill = group)) +
    geom_violin(width = 1, alpha = 0.8, linewidth = 0.3) +
    geom_boxplot(color = "#585858", alpha = 0.15, linewidth = 0.3, outlier.shape = NA) +
    geom_jitter(fill = "#585858", size = 1, stroke = 0, alpha = 0.4, width = 0.1, height = 0.1) +
    labs(y = "CNA Proportion", x = NULL) +
    scale_fill_manual(values = c("Progressor" = "darkorange", "Non-Progressor" = "darkseagreen")) +
    theme_classic() +
    guides(fill = "none") +
    annotate("text",
             x = 2, y = max(df$proportion, na.rm = TRUE) * 1.05,
             label = paste0("Wilcox Test: ", p_label),
             hjust = 1, vjust = 0, size = 3, fontface = "italic", colour = "black"
    )
}

plot_density <- function(df, p_label){
  ggplot(df, aes(x = proportion, fill = group)) +
    geom_density(alpha = 0.5) +
    labs(x = "CNA Proportion", y = "Density") +
    scale_fill_manual(values = c("Progressor" = "darkorange", "Non-Progressor" = "darkseagreen")) +
    theme_classic() +
    annotate("text",
             x = max(df$proportion, na.rm = TRUE) * 0.9,
             y = max(density(df$proportion)$y) * 0.9,
             label = paste0("Wilcox Test: ", p_label),
             hjust = 1, vjust = 1, size = 3, fontface = "italic", colour = "black"
    )
}

# Compare precursors
df_pre <- subset_groups(cn_props, nprog_pre, prog_pre)
p_label_pre <- perform_wilcox(df_pre)
p_boxplot_pre <- plot_boxplot(df_pre, p_label_pre)
p_density_pre <- plot_density(df_pre, p_label_pre)

ggsave("plots/copy_number/proportions/nprog_vs_prog_precursor_boxplot.png", p_boxplot_pre,
       width = 3, height = 2.5)
ggsave("plots/copy_number/proportions/nprog_vs_prog_precursor_density.png", p_density_pre,
       width = 5, height = 4)

# Compare follow ups
df_fol <- subset_groups(cn_props, nprog_fol, prog_fol)
p_label_fol <- perform_wilcox(df_fol)
p_boxplot_fol <- plot_boxplot(df_fol, p_label_fol)
p_density_fol <- plot_density(df_fol, p_label_fol)

ggsave("plots/copy_number/proportions/nprog_vs_prog_followup_boxplot.png", p_boxplot_fol,
       width = 3, height = 2.5)
ggsave("plots/copy_number/proportions/nprog_vs_prog_followup_density.png", p_density_fol,
       width = 5, height = 4)