library(dplyr)
library(ggplot2)
library(readr)

#TODO/Scores don't look right and need to colour on those with -log10 q value greater than 1

gistic <- read_tsv("data/copy_number/gistic/scores.gistic.tsv")

# order chromosomes 1–22, X, Y
chrom_order <- c(1:22, "X", "Y")

# compute chromosome lengths from your tibble
chrom_lengths <- gistic %>%
  group_by(Chromosome) %>%
  summarise(chr_len = max(End), .groups = "drop") %>%
  arrange(factor(Chromosome, levels = chrom_order))

# cumulative offset for plotting
chrom_lengths <- chrom_lengths %>%
  mutate(chr_start = lag(cumsum(chr_len), default = 0),
         chr_mid   = chr_start + chr_len / 2)

gistic_plot <- gistic %>%
  inner_join(chrom_lengths, by = "Chromosome") %>%
  mutate(
    Chromosome = factor(Chromosome, levels = chrom_order),
    start_cum = Start + chr_start,
    end_cum   = End + chr_start,
    score     = ifelse(Type == "Del", -`G-score`, `G-score`)
  )

p <- ggplot(gistic_plot) +
        geom_rect(aes(
            xmin = start_cum, xmax = end_cum,
            ymin = 0, ymax = score,
            fill = Type
        )) +
        scale_fill_manual(values = c("Amp" = "forestgreen", "Del" = "steelblue")) +
        scale_x_continuous(
            breaks = chrom_lengths$chr_mid,
            labels = chrom_lengths$Chromosome,
            expand = c(0,0)
        ) +
        labs(x = "Chromosome", y = "GISTIC score") +
        theme_classic(base_size = 10) +
        theme(
            axis.line.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        )


ggsave("plots/copy_number/gistic_plot.png", p, height = 3, width = 4)