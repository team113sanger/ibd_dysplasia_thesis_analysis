library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

# Load data
sample_list <- read_lines("metadata/sample_lists/all_one_ppat.list")
maf <- read_tsv("data/variants/7100_3235-filtered_mutations_matched_allTum_keepPA.maf")
meta <- read_tsv("metadata/final_metadata_qc_pass.tsv") |>
  filter(sanger_dna_id %in% sample_list) |>
  filter(!(precursor_or_follow_up == "Follow up" & group == "Non-progressor")) |>
  filter(grade_of_dysplasia != "NOS") |>
  mutate(grade_of_dysplasia = case_when(
    grade_of_dysplasia == "Low grade" ~ "LGD",
    grade_of_dysplasia == "High grade" ~ "HGD",
    grade_of_dysplasia == "Adenocarcinoma" ~ "AC",
    TRUE ~ grade_of_dysplasia
  ))

# Filter MAF for target genes
filt_maf <- maf |>
  filter(Hugo_Symbol %in% c("APC", "KRAS", "TP53")) |>
  filter(Tumor_Sample_Barcode %in% sample_list) |>
  distinct(Tumor_Sample_Barcode, Hugo_Symbol)

# Mutation status per sample × gene
gene_status <- meta |>
  select(Tumor_Sample_Barcode = sanger_dna_id, grade_of_dysplasia) |>
  crossing(Hugo_Symbol = c("APC", "KRAS", "TP53")) |>
  mutate(has_mutation = ifelse(
    paste(Tumor_Sample_Barcode, Hugo_Symbol) %in%
      paste(filt_maf$Tumor_Sample_Barcode, filt_maf$Hugo_Symbol),
    1, 0
  ))

# Summarise counts
gene_counts <- gene_status |>
  group_by(grade_of_dysplasia, Hugo_Symbol) |>
  summarise(
    mutated = sum(has_mutation),
    non_mutated = n() - mutated,
    .groups = "drop"
  )

# Chi-sq test per gene (global across 3 grades)
chisq_results <- gene_counts |>
  group_by(Hugo_Symbol) |>
  summarise(
    p_value = chisq.test(
      matrix(c(mutated, non_mutated), ncol = 2, byrow = FALSE)
    )$p.value,
    .groups = "drop"
  )

# Calculate proportions
gene_props <- gene_counts |>
  mutate(total = mutated + non_mutated,
         prop_mutated = mutated / total) |>
  left_join(chisq_results, by = "Hugo_Symbol") |>
  mutate(grade_of_dysplasia = factor(grade_of_dysplasia, levels = c("LGD", "HGD", "AC")))

# Plot
p <- ggplot(gene_props, aes(x = grade_of_dysplasia, y = prop_mutated, fill = grade_of_dysplasia)) +
  geom_col(position = "dodge") +
  facet_wrap(~ Hugo_Symbol, ncol = 1) +
  labs(y = "Proportion mutated", x = "Grade of dysplasia") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  ) +
  geom_text(
    data = gene_props |> distinct(Hugo_Symbol, p_value),
    aes(x = 2, y = 0.95, label = paste0("Chi-sq p=", signif(p_value, 2))),
    inherit.aes = FALSE, size = 3
  )

ggsave("plots/variants/mut_freqs_by_grade.png", p, width = 2, height = 5, dpi = 300)
