library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)

non_prog <- read_tsv("/lustre/scratch125/casm/team113da/projects/IBD_Associated_Dysplasia/7100_3235_IBD-associated_dysplasia/analysis/dndscv/release_v1/non_progressor_precursor/dndscv_genes_non_prog_pre.tsv")
prog <- read_tsv("/lustre/scratch125/casm/team113da/projects/IBD_Associated_Dysplasia/7100_3235_IBD-associated_dysplasia/analysis/dndscv/release_v1/progressor_precursor/dndscv_genes_prog_pre.tsv")

# Filter genes with qglobal_cv < 0.01
filtered_data <- non_prog %>%
  filter(qglobal_cv < 0.05)

# Pivot mutation columns to long format
long_data <- filtered_data %>%
  select(gene_name, qglobal_cv, n_syn, n_mis, n_non, n_spl, n_ind) %>%
  pivot_longer(
    cols = c(n_syn, n_mis, n_non, n_spl, n_ind),
    names_to = "mutation_type",
    values_to = "count"
  )

# Create the barplot
mutation_plot <- ggplot(long_data, aes(x = gene_name, y = count, fill = mutation_type)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    aes(label = ifelse(mutation_type == "n_ind", paste0("q=", round(qglobal_cv, 3)), "")),
    position = position_stack(vjust = 1.25), # Slightly higher placement
    size = 2.5
  ) +
  scale_fill_brewer(palette = "Dark2") +
  labs(
    x = NULL, # Remove x-axis title
    y = "Mutation Count",
    fill = "Mutation Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, size = 10), # Horizontal x-axis text
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    panel.grid.major.x = element_blank(), # Remove vertical gridlines for clarity
    panel.grid.minor.x = element_blank()
  )

# Save the plot
ggsave("plots/dndscv/mutation_numbers_qglobal_cv.png", plot = mutation_plot, height = 4, width = 4, dpi = 300)

#### Progressors ####
# Filter genes with qglobal_cv < 0.01
filtered_data <- prog %>%
  filter(qglobal_cv < 0.05)

# Pivot mutation columns to long format
long_data <- filtered_data %>%
  select(gene_name, qglobal_cv, n_syn, n_mis, n_non, n_spl, n_ind) %>%
  pivot_longer(
    cols = c(n_syn, n_mis, n_non, n_spl, n_ind),
    names_to = "mutation_type",
    values_to = "count"
  )

# Create the barplot
mutation_plot <- ggplot(long_data, aes(x = gene_name, y = count, fill = mutation_type)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    aes(label = ifelse(mutation_type == "n_ind", paste0("q=", round(qglobal_cv, 3)), "")),
    position = position_stack(vjust = 1.25), # Slightly higher placement
    size = 2.5
  ) +
  scale_fill_brewer(palette = "Dark2") +
  labs(
    x = NULL, # Remove x-axis title
    y = "Mutation Count",
    fill = "Mutation Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, size = 10), # Horizontal x-axis text
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    panel.grid.major.x = element_blank(), # Remove vertical gridlines for clarity
    panel.grid.minor.x = element_blank()
  )

# Save the plot
ggsave("plots/dndscv/prog_mutation_numbers_qglobal_cv.png", plot = mutation_plot, height = 4, width = 3, dpi = 300)






library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(patchwork)

# --- STEP 1: Filter significant genes and label groups ---
non_prog_filtered <- non_prog %>%
  filter(qglobal_cv < 0.05) %>%
  mutate(Group = "Non-Progressor")

prog_filtered <- prog %>%
  filter(qglobal_cv < 0.05) %>%
  mutate(Group = "Progressor")

# --- STEP 2: Combine and reshape mutation counts ---
combined <- bind_rows(non_prog_filtered, prog_filtered)

mutation_long <- combined %>%
  select(gene_name, Group, n_syn, n_mis, n_non, n_spl, n_ind) %>%
  pivot_longer(
    cols = starts_with("n_"),
    names_to = "Mutation_Type",
    values_to = "Count"
  ) %>%
  mutate(
    Mutation_Type = recode(Mutation_Type,
      n_syn = "Synonymous",
      n_mis = "Missense",
      n_non = "Nonsense",
      n_spl = "Splice",
      n_ind = "Indel"
    )
  )

mutation_summed <- mutation_long %>%
  group_by(gene_name, Group, Mutation_Type) %>%
  summarise(Total = sum(Count), .groups = "drop")

# Set consistent gene order
ordered_genes <- mutation_summed %>%
  group_by(gene_name) %>%
  summarise(TotalMut = sum(Total)) %>%
  arrange(desc(TotalMut)) %>%
  pull(gene_name)

mutation_summed <- mutation_summed %>%
  mutate(gene_name = factor(gene_name, levels = ordered_genes))

# --- STEP 3: Barplot ---
mutation_colours <- c(
  "Synonymous" = "azure4",
  "Missense" = "forestgreen",
  "Nonsense" = "firebrick",
  "Splice" = "orange2",
  "Indel" = "sienna4"
)

barplot <- ggplot(mutation_summed, aes(x = gene_name, y = Total, fill = Mutation_Type)) +
  geom_col() +
  facet_wrap(~Group, scales = "free_x") +
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

ggsave("plots/dndscv/dndscv_plot.png", plot = barplot, height = 5, width = 7, dpi = 300)






mutation_summed_p53 <- mutation_summed |>
  filter(gene_name == "TP53")

barplot <- ggplot(mutation_summed_p53, aes(x = gene_name, y = Total, fill = Mutation_Type)) +
  geom_col() +
  facet_wrap(~Group, scales = "free_x") +
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

# --- STEP 4: Prepare truncating dN/dS confidence intervals ---
prog_ci <- prog_ci %>% mutate(Group = "Progressor")
non_prog_ci <- non_prog_ci %>% mutate(Group = "Non-Progressor")

combined_ci <- bind_rows(prog_ci, non_prog_ci) %>%
  rename(gene_name = gene) %>%
  filter(gene_name %in% ordered_genes) %>%
  mutate(gene_name = factor(gene_name, levels = ordered_genes))

# --- STEP 5: Dot plot of truncating dN/dS with error bars ---
dotplot <- ggplot(combined_ci, aes(x = gene_name, y = tru_mle, colour = Group)) +
  geom_point(size = 2, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = tru_low, ymax = tru_high),
    width = 0.3,
    position = position_dodge(width = 0.6)
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30") +
  scale_y_log10() +
  facet_wrap(~Group, scales = "free_x") +
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

# --- STEP 6: Combine plots ---
combined_plot <- dotplot / barplot + plot_layout(heights = c(1, 2))

ggsave("plots/dndscv/dndscv_combined_plot.png", plot = combined_plot, height = 6, width = 5, dpi = 300)
