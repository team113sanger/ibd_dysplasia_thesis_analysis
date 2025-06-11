library(readr)
library(dplyr)
library(ggplot2)

prog_props <- read_tsv("data/copy_number/proportions/progressor_precursors_cn_props.tsv") |>
    mutate(group = "Progressor")
nprog_props <- read_tsv("data/copy_number/proportions/non-progressors_cn_props.tsv") |>
    mutate(group = "Non-Progressor")

props <- bind_rows(prog_props, nprog_props)


### Boxplot ###
p_boxplot <- ggplot(props, aes(x = group, y = proportion, fill = group)) +
        geom_boxplot() +
        labs(y = "CNA Proportion") +
        scale_fill_manual(values = c("Progressor" = "darkseagreen", "Non-Progressor" = "darkorange")) +
        theme_classic() +
        guides(fill="none")

ggsave("plots/copy_number/prog_vs_nprog_cn_prop_boxplot.png", p_boxplot,
        width = 4, height = 4)


### Wilcox Test ### 
wilcox_result <- wilcox.test(proportion ~ group, data = props)
print(wilcox_result)

p_label <- paste0("p = ", signif(wilcox_result$p.value, 3))

### Density Plot ### 
p_density <- ggplot(props, aes(x = proportion, fill = group)) +
    geom_density(alpha = 0.5) + 
    labs(x = "CNA Proportion", y = "Density") +
    scale_fill_manual(values = c("Progressor" = "darkseagreen", "Non-Progressor" = "darkorange")) +
    theme_classic() +
    annotate("text", x = Inf, y = Inf, label = paste0("Wilcox Test: ", p_label), 
             hjust = 1.1, vjust = 1.5, size = 3, fontface = "italic", colour = "black")

ggsave("plots/copy_number/prog_vs_nprog_cn_prop_density.png", p_density,
       width = 5, height = 4)
