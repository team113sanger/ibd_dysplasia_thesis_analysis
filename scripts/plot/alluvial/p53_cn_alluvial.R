library(tidyr)
library(readr)
library(dplyr)
library(ggalluvial)

# Read in data
results_df <- read_tsv("results/precursor_combined_results.tsv")

df_alluvial <- results_df |>
  group_by(p53, cn_cluster, group) |>
  summarise(Freq = n(), .groups = "drop")

# Plot
axis_labels <- data.frame(
  x = c(1, 2, 3),                  
  y = 40,
  label = c("Group", "TP53 Status", "CN Cluster")
)

p <- ggplot(df_alluvial,
        aes(axis1 = group, axis2 = p53, axis3 = cn_cluster, y = Freq, fill = group)) +
        geom_alluvium(width = 0.5, alpha = 0.7) +
        geom_stratum(width = 0.5, color = "black", alpha = 0.9) +
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

ggsave("plots/alluvial/group_cn_p53_alluvial.png", p, height = 5)

axis_labels <- data.frame(
  x = c(1, 2, 3),                  
  y = 40,
  label = c("TP53 Status", "CN Cluster", "Group")
)

p2 <- ggplot(df_alluvial,
        aes(axis1 = p53, axis2 = cn_cluster, axis3 = group, y = Freq, fill = p53)) +
        geom_alluvium(width = 0.5, alpha = 0.7) +
        geom_stratum(width = 0.5, color = "black", alpha = 0.9) +
        geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
        geom_text(data = axis_labels, aes(x = x, y = y, label = label),
                  inherit.aes = FALSE, fontface = "bold", size = 4) +
        scale_fill_manual(values = c("#41B6C4","#edf8b1"), na.value = "white") +
        theme_classic() +
        theme(
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            legend.position = "none"
        )

ggsave("plots/alluvial/p53_cn_group_alluvial.png", p2, height = 5)