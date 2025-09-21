library(dplyr)
library(readr)
library(ggplot2)
library(broom)
library(purrr)
library(tidyr)

# --- Load files ---
sample_list <- read_lines("metadata/sample_lists/all_one_ppat.list")
maf <- read_tsv("data/variants/7100_3235-filtered_mutations_matched_allTum_keepPA.maf")
meta <- read_tsv("metadata/final_metadata_qc_pass.tsv") |>
  filter(sanger_dna_id %in% sample_list)

# --- Identify TP53 mutated samples ---
tp53_maf <- maf |>
  filter(Hugo_Symbol == "TP53") |>
  filter(Tumor_Sample_Barcode %in% sample_list) |>
  distinct(Tumor_Sample_Barcode)

# --- Annotate mutation status ---
tp53_status <- meta |>
  select(Tumor_Sample_Barcode = sanger_dna_id, group, precursor_or_follow_up) |>
  mutate(has_mutation = ifelse(Tumor_Sample_Barcode %in% tp53_maf$Tumor_Sample_Barcode, 1, 0))

tp53_status$precursor_or_follow_up <- factor(tp53_status$precursor_or_follow_up,
                                            levels = c("Precursor", "Follow up"))

# --- Summarise counts ---
tp53_counts <- tp53_status |>
  group_by(group, precursor_or_follow_up) |>
  summarise(
    mutated = sum(has_mutation),
    non_mutated = n() - mutated,
    .groups = "drop"
  )

# --- Fisher’s tests per timepoint ---
tp53_pvals <- tp53_counts |>
  group_by(precursor_or_follow_up) |>
  summarise(
    test = list(
      fisher.test(
        matrix(c(
          mutated[group == "Non-progressor"], non_mutated[group == "Non-progressor"],
          mutated[group == "Progressor"],     non_mutated[group == "Progressor"]
        ), nrow = 2, byrow = TRUE)
      )
    ),
    .groups = "drop"
  ) |>
  mutate(test = map(test, tidy)) |>
  unnest(test) |>
  mutate(label = paste0("Fisher's Test\np = ", signif(p.value, 3)))

# --- Compute proportions for plotting ---
tp53_props <- tp53_counts |>
  mutate(
    prop_mutated = mutated / (mutated + non_mutated),
    group = recode(group, "Non-progressor" = "N-Prog", "Progressor" = "Prog")
  )

# --- Plot proportions with annotated p-values ---
p <- ggplot(tp53_props, aes(x = group, y = prop_mutated, fill = group)) +
        geom_col(position = "dodge") +
        facet_wrap(~precursor_or_follow_up, scales = "fixed") +
        scale_fill_manual(values = c("Prog" = "darkorange", "N-Prog" = "darkseagreen")) +
        scale_y_continuous(labels = scales::percent_format(), limits = c(0,1)) +
        labs(
            x = NULL,
            y = "TP53 mutated samples (%)"  ) +
        theme_classic(base_size = 10) +
        theme(legend.position = "none",
                strip.text = element_text(face = "bold"),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
        ) +
        # Add Fisher’s p-values manually above bars
        geom_text(
            data = tp53_pvals,
            aes(x = 1, y = 0.9, label = label),
            inherit.aes = FALSE,
            size = 3
        )

ggsave("plots/variants/tp53_fishers.png", p, width = 4, height = 2.2, dpi = 300)
ggsave("plots/variants/tp53_fishers.pdf", p, width = 4, height = 2.2)
