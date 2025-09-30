library(tidyr)
library(readr)
library(dplyr)
library(ggalluvial)

### Functions ###
prep_alluvial_data <- function(df) {
  df_alluvial <- df |>
    select(-cn_proportion) |>
    group_by(group, cn_cluster, TP53_status) |>
    summarise(Freq = n(), .groups = "drop")
  return(df_alluvial)
}

plot_alluvial <- function(plot_df, axes_order, axis_labels, colours) {
  axis_label_df <- data.frame(
    x = seq_along(axes_order),
    y = 43,
    label = axis_labels
  )
  
  ggplot(plot_df,
         aes_string(axis1 = axes_order[1], axis2 = axes_order[2], axis3 = axes_order[3],
                    y = "Freq", fill = axes_order[1])) +
    geom_alluvium(width = 0.5, alpha = 0.7) +
    geom_stratum(width = 0.5, color = "black", alpha = 0.9) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    geom_text(data = axis_label_df, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, fontface = "bold", size = 4) +
    scale_fill_manual(values = colours, na.value = "white") +
    theme_classic() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.position = "none"
    )
}

### Precursors ###
pre_df <- read_tsv("results/precursor_combined_results.tsv") |>
     mutate(group = recode(group,"Non-progressor" = "N-Prog", "Progressor" = "Prog"))


pre_alluvial <- prep_alluvial_data(pre_df)

# By group
p1 <- plot_alluvial(
  pre_alluvial,
  axes_order = c("group", "TP53_status", "cn_cluster"),
  axis_labels = c("Group", "TP53 Status", "CN Cluster"),
  colours = c("N-Prog" = "#8FB996", "Prog" = "#E6A272")
)
ggsave("plots/alluvial/precursors/group_cn_TP53.png", p1, height = 3, width = 5)

# By TP53
p2 <- plot_alluvial(
  pre_alluvial,
  axes_order = c("TP53_status", "cn_cluster", "group"),
  axis_labels = c("TP53 Status", "CN Cluster", "Group"),
  colours = c("WT" = "#79AF97", "Mut" = "#6A6599")
)
ggsave("plots/alluvial/precursors/TP53_cn_group.png", p2, height = 3, width = 5)

### Follow Ups ###
fol_df <- read_tsv("results/follow_up_combined_results.tsv") |>
     mutate(group = recode(group,"Non-progressor" = "N-Prog", "Progressor" = "Prog"))

fol_alluvial <- prep_alluvial_data(fol_df)

# By group
f1 <- plot_alluvial(
  fol_alluvial,
  axes = c("group", "TP53_status", "cn_cluster"),
  axis_labels = c("Group", "TP53 Status", "CN Cluster"),
  colours = c("N-Prog" = "#8FB996", "Prog" = "#E6A272")
)
ggsave("plots/alluvial/follow_ups/group_cn_TP53.png", f1, height = 3, width = 5)

# By TP53
f2 <- plot_alluvial(
  fol_alluvial,
  axes = c("TP53_status", "cn_cluster", "group"),
  axis_labels = c("TP53 Status", "CN Cluster", "Group"),
  colours = c("WT" = "#79AF97", "Mut" = "#6A6599")
)
ggsave("plots/alluvial/follow_ups/TP53_cn_group.png", f2, height = 3, width = 5)