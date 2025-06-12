library(readr)
library(dplyr)
library(ggplot2)

maf <- read_tsv("data/variants/7100_3235-filtered_mutations_all_indepTum_keepPA.maf")
meta <-read_tsv("metadata/sample_pairs.tsv")

### Prep ### 
# Filter MAF
genes_of_interest <- c("APC", "KRAS", "TP53")

filtered_maf <- maf %>%
  filter(Hugo_Symbol %in% genes_of_interest)

# Add metadata
maf_annotated <- filtered_maf %>%
  inner_join(meta, by = c("Tumor_Sample_Barcode" = "sanger_dna_id"))

mutation_presence <- maf_annotated %>%
  distinct(patient_id, precursor_or_follow_up, Hugo_Symbol, group) %>%
  mutate(mutation = 1) %>%
  tidyr::complete(patient_id, precursor_or_follow_up, Hugo_Symbol, group, fill = list(mutation = 0))

mutation_presence[["precursor_or_follow_up"]] <- factor(
  mutation_presence[["precursor_or_follow_up"]],
  levels = c("Precursor", "Follow up")
)

### Barplot ###
bar_data <- mutation_presence %>%
  group_by(Hugo_Symbol, precursor_or_follow_up, group) %>%
  summarise(freq = mean(mutation), .groups = "drop")

p <- ggplot(bar_data, aes(x = freq, y = Hugo_Symbol, fill = group)) +
      geom_col(position = "dodge") +
      scale_x_continuous(labels = scales::percent_format()) +
      theme_classic(base_size = 13) +
      labs(x = "Mutation Frequency", fill = "Sample") +
      facet_wrap(~precursor_or_follow_up) +
      theme(axis.text.y = element_text(face = "italic"),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold"),
            axis.title.y = element_blank(),
            legend.position = "right") +
      scale_fill_manual(values = c("Progressor" = "darkorange", "Non-progressor" = "darkseagreen"))

ggsave("plots/variants/mutation_freq_barplot.png", plot = p, width = 7, height = 3.5)