library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(patchwork)

#### Prep ####
# Read in data
prog <- read_tsv("data/variants/dndscv/dndscv_genes_prog_pre.tsv") |>
  mutate(group = "Progressor")
non_prog <- read_tsv("data/variants/dndscv/dndscv_genes_non_prog_pre.tsv") |>
  mutate(group = "Non-Progressor")
prog_ci <- read_tsv("data/variants/dndscv/prog_gene_ci.tsv") |>
  mutate(group = "Progressor")
non_prog_ci <- read_tsv("data/variants/dndscv/non_prog_gene_ci.tsv") |>
  mutate(group = "Non-Progressor")

combined_df <- bind_rows(prog, non_prog)

# Reshape and filter data
muts_df <- combined_df |>
  select(gene_name, qglobal_cv, group, n_syn, n_mis, n_non, n_spl, n_ind) |>
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

sig_muts_sums <- muts_df |>
  filter(qglobal_cv < 0.05) |>
  group_by(gene_name, group, mutation_type) |>
  summarise(total = sum(count), .groups = "drop")

##### Plot ##### 
### Barplot ###
# Set gene order
ordered_genes <- sig_muts_sums |>
  group_by(gene_name) |>
  summarise(TotalMut = sum(total)) |>
  arrange(desc(TotalMut)) |>
  pull(gene_name)

plot_df <- sig_muts_sums |>
  mutate(gene_name = factor(gene_name, levels = ordered_genes))

# Define colours
mutation_colours <- c(
  "Synonymous" = "azure4",
  "Missense" = "forestgreen",
  "Nonsense" = "firebrick",
  "Splice" = "orange2",
  "Indel" = "sienna4"
)

# Plot
counts_barplot <- ggplot(sig_muts_sums, aes(x = gene_name, y = total, fill = mutation_type)) +
  geom_col() +
  facet_wrap(~group, scales = "free_x") +
  scale_fill_manual(values = mutation_colours) +
  labs(
    x = "Gene",
    y = "Mutation Count",
    fill = "Mutation Type"
  ) +
  theme_classic(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.line.y = element_line(colour = "black"),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9)
  )

ggsave("plots/dndscv/precursors/dndscv_barplot.png", plot = counts_barplot, height = 4, width = 6, dpi = 300)

### TP53 Plot ###
tp53_df <- muts_df |>
  filter(gene_name == "TP53") |>
  group_by(gene_name, group, mutation_type) |>
  summarise(total = sum(count), .groups = "drop")

barplot <- ggplot(tp53_df, aes(x = gene_name, y = total, fill = mutation_type)) +
              geom_col(position = "stack") +
              facet_wrap(~group, scales = "free_x") +
              scale_fill_manual(values = mutation_colours) +
              labs(
                x = "Gene",
                y = "Mutation Count",
                fill = "Mutation Type"
              ) +
              theme_classic(base_size = 13) +
              theme(
                panel.grid = element_blank(),
                axis.line.y = element_line(colour = "black"),
                axis.line.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.text.x = element_text(angle = 45, hjust = 1),
                strip.background = element_blank(),
                strip.text = element_text(face = "bold"),
                legend.position = "right",
                legend.key.size = unit(0.4, "cm"),
                legend.text = element_text(size = 9),
                legend.title = element_text(size = 9)
              )

# Prep
ci_df <- bind_rows(prog_ci, non_prog_ci) |>
  rename(gene_name = gene) |>
  filter(gene_name %in% ordered_genes) |>
  mutate(gene_name = factor(gene_name, levels = ordered_genes))

# Dot plot of truncating dN/dS with error bars
dotplot <- ggplot(ci_df, aes(x = gene_name, y = tru_mle, colour = group)) +
  geom_point(size = 2, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = tru_low, ymax = tru_high),
    width = 0.3,
    position = position_dodge(width = 0.6)
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30") +
  scale_y_log10() +
  facet_wrap(~group, scales = "free_x") +
  theme_classic(base_size = 13) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position = "none"
  ) +
  labs(y = "Truncating dN/dS")

# Combine
combined_plot <- dotplot / barplot + plot_layout(heights = c(1, 2))

ggsave("plots/dndscv/precursors/dndscv_tp53_plot.png", plot = combined_plot, height = 6, width = 5, dpi = 300)
