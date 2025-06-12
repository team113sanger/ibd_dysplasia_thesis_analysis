library(tidyr)
library(readr)
library(dplyr)
library(ggplot2)

# Load data
maf <- read_tsv("data/variants/7100_3235-filtered_mutations_all_indepTum_keepPA.maf")
meta <- read_tsv("metadata/sample_pairs.tsv")

# Reduce to the genes of interest
genes_of_interest <- c("TP53", "APC", "KRAS")

# Step 1: Get actual mutations
patient_mutations <- maf %>%
  filter(Hugo_Symbol %in% genes_of_interest) %>%
  inner_join(meta, by = c("Tumor_Sample_Barcode" = "sanger_dna_id")) %>%
  distinct(patient_id, precursor_or_follow_up, Hugo_Symbol) %>%
  mutate(mutation = 1)

# Step 2: Get all expected combinations — only across the required dimensions
all_combinations <- meta %>%
  distinct(patient_id, precursor_or_follow_up) %>%
  tidyr::crossing(Hugo_Symbol = genes_of_interest)

# Step 3: Join + fill missing mutation = 0
patient_mutations_complete <- all_combinations %>%
  left_join(patient_mutations, by = c("patient_id", "precursor_or_follow_up", "Hugo_Symbol")) %>%
  mutate(mutation = replace_na(mutation, 0)) %>%
  left_join(meta %>% distinct(patient_id, group), by = "patient_id")

# Step 4: Create patient index (sorted by group then patient_id)
patient_ids <- patient_mutations_complete %>%
  arrange(group) |>
  distinct(patient_id) %>%
  mutate(patient_index = row_number())

# Add index to main data
plot_data <- patient_mutations_complete %>%
  left_join(patient_ids, by = c("patient_id"))

plot_data[["precursor_or_follow_up"]] <- factor(
  plot_data[["precursor_or_follow_up"]],
  levels = c("Precursor", "Follow up")
)

p <- ggplot(plot_data, aes(
  x = precursor_or_follow_up,
  y = patient_index,
  group = patient_id
)) +
  geom_line(aes(colour = group), linewidth = 0.6, alpha = 0.5) +
  geom_point(
    aes(shape = factor(mutation), fill = factor(mutation)),
    size = 3, colour = "black"
  ) +
  facet_wrap(~ Hugo_Symbol, scales = "free_y") +
  scale_fill_manual(
    values = c("0" = "white", "1" = "black"),
    name = "Mutation"
  ) +
  scale_shape_manual(
    values = c("0" = 21, "1" = 24),
    name = "Mutation"
  ) +
  scale_colour_manual(
    values = c("Progressor" = "#E64B35", "Non-progressor" = "#3CB371")
  ) +
  labs(x = NULL, y = "Patient", colour = "Group") +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "italic"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave("plots/variants/patient_mutation_dotplot.png", p, width = 8, height = 5.5)
