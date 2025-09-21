library(dplyr)
library(readr)
library(ggplot2)
library(purrr)
library(tidyr)

# --- Load data ---
sample_list <- read_lines("metadata/sample_lists/all_one_ppat.list")
maf <- read_tsv("data/variants/7100_3235-filtered_mutations_matched_allTum_keepPA.maf")
meta <- read_tsv("metadata/final_metadata_qc_pass.tsv") |>
  filter(sanger_dna_id %in% sample_list) |>
  mutate(group = case_when(
    group == "Non-progressor" ~ "N-Prog",
    group == "Progressor" ~ "Prog",
    TRUE ~ group
  ))

# --- Helper function to run whole pipeline ---
make_gene_plot <- function(status_filter) {
  
  # filter metadata by precursor/follow-up
  meta_sub <- meta %>% filter(precursor_or_follow_up == status_filter)
  
  # filter maf
  filt_maf <- maf %>%
    filter(Hugo_Symbol %in% c("APC", "KRAS", "TP53")) %>%
    filter(Tumor_Sample_Barcode %in% meta_sub$sanger_dna_id) %>%
    distinct(Tumor_Sample_Barcode, Hugo_Symbol)
  
  # mutation status per gene/sample
  gene_status <- meta_sub %>%
    select(Tumor_Sample_Barcode = sanger_dna_id, group) %>%
    crossing(Hugo_Symbol = c("APC", "KRAS", "TP53")) %>%
    mutate(has_mutation = ifelse(
      paste(Tumor_Sample_Barcode, Hugo_Symbol) %in%
        paste(filt_maf$Tumor_Sample_Barcode, filt_maf$Hugo_Symbol),
      1, 0
    ))
  
  # counts
  gene_counts <- gene_status %>%
    group_by(group, Hugo_Symbol) %>%
    summarise(
      mutated = sum(has_mutation),
      non_mutated = n() - mutated,
      .groups = "drop"
    )
  
  # fisher test
  fisher_results <- gene_counts %>%
    group_by(Hugo_Symbol) %>%
    summarise(
      p_value = fisher.test(
        matrix(c(mutated[group == "Prog"],
                 non_mutated[group == "Prog"],
                 mutated[group == "N-Prog"],
                 non_mutated[group == "N-Prog"]),
               nrow = 2, byrow = TRUE)
      )$p.value
    )
  
  # proportions
  gene_props <- gene_counts %>%
    mutate(total = mutated + non_mutated,
           prop_mutated = mutated / total) %>%
    left_join(fisher_results, by = "Hugo_Symbol")
  
  # plot
  p <- ggplot(gene_props, aes(x = group, y = prop_mutated, fill = group)) +
    geom_col(position = "dodge") +
    facet_wrap(~ Hugo_Symbol, ncol = 1) +
    labs(y = "Proportion mutated", x = NULL) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0,1)) +
    scale_fill_manual(values = c("Prog" = "darkorange", "N-Prog" = "darkseagreen")) +
    theme_classic(base_size = 10) +
    theme(legend.position = "none",
          strip.text = element_text(face = "bold"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)) +
    geom_text(
      data = gene_props %>% distinct(Hugo_Symbol, p_value),
      aes(x = 1.5, y = 0.95, label = paste0("Fisher p=", signif(p_value, 2))),
      inherit.aes = FALSE, size = 3
    )
  
  return(p)
}

# --- Run for both precursor and follow-up ---
plots <- map(c("Precursor", "Follow up"), make_gene_plot)

# Save
ggsave("plots/variants/mut_freq_pres.png", plots[[1]], width = 2, height = 5, dpi = 300)
ggsave("plots/variants/mut_freq_fols.png", plots[[2]], width = 2, height = 5, dpi = 300)
