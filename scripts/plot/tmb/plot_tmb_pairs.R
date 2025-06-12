library(readr)
library(dplyr)
library(ggplot2)

tmb <- read_tsv("data/variants/mutations_per_Mb.tsv", col_names = c("sanger_dna_id", "tmb"))
  
sample_pairs <- read_tsv("metadata/sample_pairs.tsv") |>
  left_join(tmb)

sample_pairs[["precursor_or_follow_up"]] <- factor(
    sample_pairs[["precursor_or_follow_up"]], 
    levels = c("Precursor", "Follow up")
    )

### Boxplot ###
p_boxplot <- ggplot(sample_pairs, aes(x = precursor_or_follow_up, y = tmb, fill = precursor_or_follow_up)) +
        geom_boxplot(alpha = 0.8) +
        labs(
            y = "TMB (Mutations/Mb)",
            fill = "Sample Status"
        ) +
        theme_bw(base_size = 14) +
        theme(legend.position = "none",
              axis.title.x = element_blank()) +
        scale_fill_brewer(palette = "Dark2") +
        scale_y_log10() +
        facet_wrap(~ group)

ggsave("plots/tmb/samples_pairs_tmb_boxplot.png", p_boxplot, width = 5.6, height = 4)

### Wilcox Test ###
wilcox_prog <- wilcox.test(
  tmb ~ precursor_or_follow_up,
  data = filter(sample_pairs, group == "Progressor")
)

wilcox_nprog <- wilcox.test(
  tmb ~ precursor_or_follow_up,
  data = filter(sample_pairs, group == "Non-progressor")
)

### Density Plot ###
p_density <- ggplot(sample_pairs, aes(x = tmb, fill = precursor_or_follow_up)) +
    geom_density(alpha = 0.6) +
    scale_x_log10() +
    labs(
        x = "TMB (Mutations/Mb)",
        y = "Density",
        fill = "Sample Status"
    ) +
    facet_wrap(~ group) + 
    theme_bw(base_size = 14) +
    theme(legend.position = "top") +
    scale_fill_brewer(palette = "Dark2")

ggsave("plots/tmb/sample_pairs_tmb_density.png", p_density, width = 7, height = 5)


### Individual Pairs ###

p_pairs <- ggplot(sample_pairs, aes(x = precursor_or_follow_up, y = tmb, group = patient_id)) +
    geom_line(aes(colour = group), alpha = 0.6, linewidth = 0.8) +
    geom_point(aes(fill = group), shape = 21, size = 3, stroke = 0.2) +
    scale_y_log10() +
    scale_fill_manual(values = c("Progressor" = "darkorange", "Non-progressor" = "darkseagreen")) +
    scale_colour_manual(values = c("Progressor" = "darkorange", "Non-progressor" = "darkseagreen")) +
    facet_wrap(~ group) +
    labs(
        x = NULL,
        y = "TMB (Mutations/Mb)"
    ) +
    theme_bw(base_size = 14) +
    theme(legend.position = "none")

ggsave("plots/tmb/paired_tmb_lineplot.png", p_pairs, width = 6, height = 4)