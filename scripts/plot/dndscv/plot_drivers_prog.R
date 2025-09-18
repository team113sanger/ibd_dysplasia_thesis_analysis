library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

#### Prep ####
# Read in data
prog <- read_tsv("data/variants/dndscv/dndscv_genes_prog_pre.tsv")
prog_ci <- read_tsv("data/variants/dndscv/prog_pre_gene_ci.tsv") |>
    filter(!gene == "RBM10")

# Filter and reshape data
muts_df <- prog |>
  select(gene_name, n_syn, n_mis, n_non, n_spl, n_ind, qglobal_cv) |>
  filter(qglobal_cv < 0.05) |>
  pivot_longer(
    cols = starts_with("n_"),
    names_to = "mutation_type",
    values_to = "count"
  ) |>
  mutate(
    mutation_type = recode(mutation_type,
      n_syn = "Synonymous", n_mis = "Missense", n_non = "Nonsense",
      n_spl = "Splice", n_ind = "Indel"
    )
  )

ratios_df <- prog |>
  filter(qglobal_cv < 0.05) |>
  select(gene_name, wmis_cv, wnon_cv, wspl_cv, wind_cv) |>
  pivot_longer(
    cols = starts_with("w"),
    names_to = "mutation_type",
    values_to = "w_dnds"
  ) |>
  mutate(
    mutation_type = recode(mutation_type,
      wmis_cv = "Missense",
      wnon_cv = "Nonsense",
      wspl_cv = "Splice",
      wind_cv = "Indel"
    ),
    gene_name = factor(gene_name, levels = unique(gene_name)) # preserve order
  )

ci_df <- prog_ci |>
  rename(gene_name = gene) |>
  tidyr::pivot_longer(cols = starts_with(c("mis", "tru")),
                      names_to = c("type", ".value"),
                      names_pattern = "(mis|tru)_(.*)") |>
  mutate(type = dplyr::recode(type,
                         "mis" = "Missense",
                         "tru" = "Truncating")
        )

# Set order of genes
genes_order <- muts_df |>
    arrange(qglobal_cv) |>
    pull(gene_name) |>
    unique()

muts_df <- muts_df |>
    mutate(gene_name = factor(gene_name, labels = genes_order))

ratios_df <- ratios_df |>
    mutate(gene_name = factor(gene_name, labels = genes_order))


ci_df <- ci_df |>
    mutate(gene_name = factor(gene_name, labels = genes_order))


# Define colour palette
mutation_colours <- c(
  "Missense"   = "#8cb87aff",  
  "Truncating" = "#7c679eff",
  "Indel" = "#77a9baff",
  "Nonsense" = "#7c679eff",
  "Splice" = "#d5ad5eff",
  "Synonymous" ="#696a6bff"
)

# Counts barplot
counts_barplot <- ggplot(muts_df, aes(x = gene_name, y = count, fill = mutation_type)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = mutation_colours) +
  labs(
    x = NULL,
    y = "Total mutations",
    fill = NULL,
    title = "Progressor"
  ) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.text  = element_text(size = 6),
    legend.key.size = unit(0.3, "cm")  
  )

# ratios barplot
ratios_barplot <- ggplot(ratios_df, aes(x = gene_name, y = w_dnds, fill = mutation_type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  labs(
    x = NULL,
    y = "dN/dS ratio",
    fill = NULL
  ) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.text  = element_text(size = 6),
    legend.key.size = unit(0.3, "cm")
  ) +
  scale_fill_manual(values = mutation_colours) 

# dnds ratios barplot
# Plot
ci_barplot <- ggplot(ci_df, aes(x = gene_name, y = mle, fill = type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = low, ymax = high),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30") +
  scale_fill_manual(values = mutation_colours) +
  scale_y_log10(trans = "pseudo_log", breaks = c(1, 10, 100, 1000, 10000)) +
  labs(
    x = "Genes under positive selection",
    y = "dN/dS ratio",
    fill = NULL
  ) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(face = "bold.italic"),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.text  = element_text(size = 6),
    legend.key.size = unit(0.3, "cm")
  )

# Combine plots
combined_plot <- counts_barplot / ratios_barplot / ci_barplot + plot_layout(heights = c(1, 1, 1))

ggsave("plots/dndscv/precursors/dndscv_progressors_plot.png",
        plot = combined_plot, height = 3, width = 3.4, dpi = 300)
