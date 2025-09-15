library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

#### Prep ####
# Read in data
prog <- read_tsv("data/variants/dndscv/dndscv_genes_prog_pre.tsv") |>
  mutate(group = "Progressor")
non_prog <- read_tsv("data/variants/dndscv/dndscv_genes_non_prog_pre.tsv") |>
  mutate(group = "Non-Progressor")
prog_ci <- read_tsv("data/variants/dndscv/prog_pre_gene_ci.tsv") |>
  mutate(group = "Progressor")
non_prog_ci <- read_tsv("data/variants/dndscv/non_prog_pre_gene_ci.tsv") |>
  mutate(group = "Non-Progressor")

# combine dfs
combined_muts <- bind_rows(prog, non_prog)
combine_ci <- bind_rows(prog_ci, non_prog_ci) |>
  rename(gene_name = gene) 

# n in each group
n_prog <- 19
n_nprog <- 21

# Reshape and filter data
muts_df <- combined_muts |>
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

ci_df <- combine_ci |>
  tidyr::pivot_longer(cols = starts_with(c("mis", "tru")),
                      names_to = c("type", ".value"),
                      names_pattern = "(mis|tru)_(.*)") |>
  mutate(type = dplyr::recode(type,
                         "mis" = "Missense",
                         "tru" = "Truncating")
        ) |>
  rename(dnds = mle)

# Only significant genes
sig_muts_sums <- muts_df |>
  filter(qglobal_cv < 0.05) |>
  group_by(gene_name, group, mutation_type) |>
  summarise(total = sum(count), .groups = "drop")

# Get genes for plot
genes_to_plot <- sig_muts_sums |>
  group_by(gene_name) |>
  arrange(desc(total)) |>
  pull(gene_name) |>
  unique()

# Get percentage mutated
gene_pct <- sig_muts_sums |>
  group_by(group, gene_name) |>
  summarise(total_mut = sum(total), .groups = "drop") |>
  mutate(percent_samples = round(100 * total_mut / ifelse(group == "Progressor", n_prog, n_nprog), 1))

##### Plot ##### 
# Define colours
mutation_colours <- c(
  "Synonymous" = "#8DA0CB",   # slate blue-grey
  "Missense"   = "#68a334ff",   # muted green
  "Nonsense"   = "#9f2727ff",   # softer, warm red
  "Splice"     = "#E6AB02",   # goldenrod
  "Indel"      = "#7d6a56ff"    # teal
)

# Mutation count barplot
counts_barplot <- ggplot(sig_muts_sums, aes(x = reorder(gene_name, -total), y = total, fill = mutation_type)) +
  geom_col() +
  geom_text(
    data = gene_pct,
    aes(x = gene_name, y = total_mut + 2.5 , label = paste0(percent_samples, "%")),
    inherit.aes = FALSE, 
    vjust = 1, size = 3, angle = 45
  ) +
  facet_wrap(~group, scales = "free_x", ncol =1, strip.position = "top") +
  scale_fill_manual(values = mutation_colours) +
  labs(
    x = NULL,
    y = "Total Mutations",
    fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
   # axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1.2, face = "italic"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9)
  ) +
  scale_y_continuous(limits = c(0,20))

ggsave(
  "plots/dndscv/precursors/dndscv_counts_barplot.png",
  plot = counts_barplot, height = 5.5, width = 3.5, dpi = 300
)

# dnds ratios barplot
# dnds_plot <- ggplot(ci_df, aes(x = gene_name, y = dnds, fill = type)) +
#   geom_col(position = position_dodge(width = 0.9)) +
#   geom_errorbar(
#     aes(ymin = low, ymax = high),
#     width = 0.3,
#     position = position_dodge(width = 0.9)
#   ) +
#   facet_wrap(~group, ncol = 1) +
#   scale_y_log10() +
#   scale_fill_brewer(palette = "Set2") +
#   labs(
#     x = "Gene",
#     y = "dN/dS (log scale)",
#     fill = NULL
#   ) +
#   theme_classic(base_size = 12) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
#     strip.text = element_text(face = "bold"),
#     legend.position = "right"
#   )

# ggsave("plots/dndscv/precursors/dnds_ratios_barplot.png",
#        plot = dnds_plot, height = 5, width = 6, dpi = 300)
