library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)

cn_props <- read_tsv("data/copy_number/proportions/all_cn_props.tsv") |>
    rename(sanger_dna_id = Sample)

sample_pairs <- read_tsv("metadata/sample_pairs.tsv") |>
  left_join(props)

sample_pairs[["precursor_or_follow_up"]] <- factor(
  sample_pairs[["precursor_or_follow_up"]],
  levels = c("Precursor", "Follow up")
)

### Line Plot ###
p_pairs <- ggplot(sample_pairs, aes(x = precursor_or_follow_up, y = proportion, group = patient_id)) +
  geom_line(aes(colour = group), alpha = 0.6, linewidth = 0.8) +
  geom_point(aes(fill = group), shape = 21, size = 3, stroke = 0.2) +
  geom_text_repel(aes(label = patient_id, colour = group), 
                  size = 2, max.overlaps = 20, show.legend = FALSE) +
  scale_y_log10() +
  scale_fill_manual(values = c("Progressor" = "darkorange", "Non-progressor" = "darkseagreen")) +
  scale_colour_manual(values = c("Progressor" = "darkorange", "Non-progressor" = "darkseagreen")) +
  facet_wrap(~group) +
  labs(
    x = NULL,
    y = "TMB (Mutations/Mb)"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

ggsave("plots/copy_number/proportions/paired_cn_prop_lineplot.png", p_pairs, width = 6, height = 4)

### Wilcox Test ###
props_pre <- sample_pairs |>
  filter(precursor_or_follow_up == "Precursor")
wilcox_result <- wilcox.test(proportion ~ group, data = props_pre)
print(wilcox_result)
p_label <- paste0("p = ", signif(wilcox_result$p.value, 3))

props_fol <- sample_pairs |>
  filter(precursor_or_follow_up == "Follow up")
wilcox_result <- wilcox.test(proportion ~ group, data = props_fol)
print(wilcox_result)
p_label_2 <- paste0("p = ", signif(wilcox_result$p.value, 3))

pvals_df <- tibble(
  precursor_or_follow_up = c("Precursor", "Follow up"),
  p_label = c(
    paste0("Wilcox Test: ", p_label),
    paste0("Wilcox Test: ", p_label_2)
  ),
  x = 1,  # Adjust X and Y as needed for plot positioning
  y = 0.9
) |>
  mutate(
    precursor_or_follow_up = factor(
      precursor_or_follow_up,
      levels = c("Precursor", "Follow up")
    )
  )

### Boxplot ###
sample_pairs <- sample_pairs |>
    mutate(group = recode(group,"Non-progressor" = "N-Prog", "Progressor" = "Prog")) |>
    mutate(
      precursor_or_follow_up = factor(
        precursor_or_follow_up,
        levels = c("Precursor", "Follow up") 
      )
    )

p_boxplot <- ggplot(sample_pairs, aes(x = group, y = proportion, fill = group)) +
  geom_boxplot() +
  labs(y = "CNA Proportion") +
  scale_fill_manual(values = c("Prog" = "darkorange", "N-Prog" = "darkseagreen")) +
  theme_classic() +
  guides(fill = "none") +
  facet_wrap(~precursor_or_follow_up, ncol = 1) +
  geom_text(data = pvals_df, aes(x = x, y = y, label = p_label), 
            inherit.aes = FALSE, size = 2.5, fontface = "italic", hjust = 0)

ggsave("plots/copy_number/proportions/paired_cn_prop_boxplot.png", p_boxplot, width = 2.5, height = 4)
#ggsave("plots/copy_number/proportions/paired_cn_prop_boxplot.pdf", p_boxplot, width = 2.5, height = 4)
