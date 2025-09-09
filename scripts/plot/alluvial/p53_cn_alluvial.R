library(tidyr)
library(readr)
library(dplyr)
library(ggalluvial)

# Read in data
cn_clusters <- read_tsv("results/copy_number/cn_clusters.tsv") |>
    select(sample, `4`) |>
    rename(cn_cluster = `4`)
tp53_status <- read_tsv("results/p53_mutations/p53_status.tsv")
meta <- read_tsv("metadata/final_metadata_qc_pass.tsv")
samples <- c(
    read_lines("metadata/sample_lists/non_progressor_precursor_samples_ppat.tsv"),
    read_lines("metadata/sample_lists/progressor_precursor_samples_ppat.tsv"))

# Combine data 
df_all <- meta |>
  select(sanger_dna_id, group) |>
  filter(sanger_dna_id %in% samples) |>
  mutate(
    p53 = if_else(sanger_dna_id %in% tp53_status$Tumor_Sample_Barcode, "Mut", "WT")
  ) |>
  left_join(
    cn_clusters |> rename(sanger_dna_id = sample),
    by = "sanger_dna_id"
  )

df_alluvial <- df_all |>
  group_by(p53, cn_cluster, group) |>
  summarise(Freq = n(), .groups = "drop")

# Plot
axis_labels <- data.frame(
  x = c(1, 2, 3),                  
  y = 42,
  label = c("Group", "TP53 Status", "CN Cluster")
)

p <- ggplot(df_alluvial,
        aes(axis1 = group, axis2 = p53, axis3 = cn_cluster, y = Freq, fill = group)) +
        geom_alluvium(width = 0.5, alpha = 0.8) +
        geom_stratum(width = 0.5, color = "black") +
        geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
        geom_text(data = axis_labels, aes(x = x, y = y, label = label),
                  inherit.aes = FALSE, fontface = "bold", size = 4) +
        scale_fill_manual(values = c("Non-progressor" = "#8FB996",
                             "Progressor" = "#E6A272"), na.value = "white") +
        theme_classic() +
        theme(
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            legend.position = "none"
        )

ggsave("plots/alluvial/p53_cn_alluvial.png", p, height = 5)




