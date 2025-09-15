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

n_prog <- 19
n_nprog <- 21

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

gene_pct <- sig_muts_sums |>
  group_by(group, gene_name) |>
  summarise(total_mut = sum(total), .groups = "drop") |>
  mutate(percent_samples = round(100 * total_mut / ifelse(group == "prog", n_prog, n_nprog), 1))

##### Plot ##### 
### Barplot ###
# Define colours
mutation_colours <- c(
  "Synonymous" = "#8DA0CB",   # slate blue-grey
  "Missense"   = "#68a334ff",   # muted green
  "Nonsense"   = "#9f2727ff",   # softer, warm red
  "Splice"     = "#E6AB02",   # goldenrod
  "Indel"      = "#7d6a56ff"    # teal
)

# Plot
counts_barplot <- ggplot(sig_muts_sums, aes(x = reorder(gene_name, -total), y = total, fill = mutation_type)) +
  geom_col() +
  geom_text(
    data = gene_pct,
    aes(x = gene_name, y = total_mut + 1.5 , label = paste0(percent_samples, "%")),
    inherit.aes = FALSE, 
    vjust = 1, size = 3, angle = 45
  ) +
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
    axis.text.x = element_text(angle = 90, hjust = 1.2, face = "italic"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9)
  ) +
  scale_y_continuous(limits = c(0,18))

ggsave("plots/dndscv/precursors/dndscv_barplot.png", plot = counts_barplot, height = 4, width = 4.7, dpi = 300)

### TP53 Plot ###
tp53_df <- muts_df |>
  filter(gene_name == "TP53") |>
  group_by(gene_name, group, mutation_type) |>
  summarise(total = sum(count), .groups = "drop")

barplot <- ggplot(tp53_df, aes(x = gene_name, y = total, fill = mutation_type)) +
              geom_col(position = "stack", width = 0.7) +
              facet_wrap(~group, scales = "free_x") +
              scale_fill_manual(values = mutation_colours) +
              labs(
                x = "TP53",
                y = "Mutation Count",
                fill = "Mutation Type"
              ) +
              theme_classic(base_size = 13) +
              theme(
                panel.grid = element_blank(),
                axis.line.y = element_line(colour = "black"),
                axis.line.x = element_blank(),
                axis.title.x = element_text(face = "italic"),
                axis.ticks.x = element_blank(),
                axis.text.x = element_blank(),
                strip.background = element_blank(),
                strip.text = element_blank(),
                legend.position = "right",
                legend.key.size = unit(0.4, "cm"),
                legend.text = element_text(size = 9),
                legend.title = element_text(size = 9)
              )

# Prep
ci_df <- bind_rows(prog_ci, non_prog_ci) |>
  rename(gene_name = gene) |>
  tidyr::pivot_longer(cols = starts_with(c("mis", "tru")),
                      names_to = c("type", ".value"),
                      names_pattern = "(mis|tru)_(.*)") |>
  mutate(type = dplyr::recode(type,
                         "mis" = "Missense",
                         "tru" = "Truncating")
        )

# Dot plot of truncating dN/dS with error bars
dotplot <- ggplot(ci_df, aes(x = gene_name, y = mle, colour = type)) +
  geom_point(size = 2, position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = low, ymax = high),
    width = 0.3,
    position = position_dodge(width = 0.7)
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30") +
  scale_color_manual(values = c("#68a334ff", "#9f2727ff")) +
  scale_y_log10(trans = "pseudo_log", breaks = c(1, 10, 100, 1000, 10000)) +
  facet_wrap(~group, scales = "free_x") +
  theme_classic(base_size = 13) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.position = "right"
  ) +
  labs(y = "dN/dS", colour = "Mutation Type")

# Combine
combined_plot <- dotplot / barplot + plot_layout(heights = c(1, 1.5))

ggsave("plots/dndscv/precursors/dndscv_tp53_plot.png", plot = combined_plot, height = 5, width = 5.2, dpi = 300)
