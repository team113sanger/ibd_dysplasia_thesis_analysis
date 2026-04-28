library(readr)
library(dplyr)
library(ggplot2)
# read in tg data
tg_calls <- read_tsv("data/copy_number/trevor_graham/tg_discovery_calls.tsv") |>
    tibble::column_to_rownames("Sample") |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("Bin")
tg_bin_sizes <- read_tsv("data/copy_number/trevor_graham/tg_bin_sizes.tsv")
tg_samples <- read_csv("data/copy_number/trevor_graham/tg_discovery_samples.csv") |>
    rename(sample = Sample)

# read in our data
calls <- read_tsv("results/copy_number/cnv_all_calls.txt")

# Plot freq
# Single Plot
tidy_tg_calls <- reshape2::melt(tg_calls)
colnames(tidy_tg_calls) <- c("bin", "sample", "cn_state")

tidy_tg_calls$cn_state <- as.character(tidy_tg_calls$cn_state)

tidy_tg_calls$bin <- sub("^V", "", tidy_tg_calls$bin)
tidy_tg_calls$sample <- sub(".*\\.", "", tidy_tg_calls$sample)

tidy_tg_calls <- tidy_tg_calls |>
  mutate(bin = as.numeric(bin)) |>
  left_join(tg_samples |> select(sample, Progression)) |>
   filter(!is.na(Progression))

sample_list <- unique(tidy_tg_calls$sample)

# heatmap
p <- ggplot(data = tidy_tg_calls) +
  geom_tile(aes(x = bin, y = sample, fill = cn_state)) +
  geom_hline(yintercept = ((1:length(sample_list)) - 0.5), color = "gray80", size = 0.3) +
  geom_vline(xintercept = 4401, color = "gray80", size = 0.3) +
  scale_fill_manual(values = c("0" = "white", "1" = "#B04968", "2" = "#ac3659", "-1" = "#4776A2", "-2" = "#306ba2"),
                    labels = c("0" = "No change", "1" = "Gain", "-1" = "Loss"),
                    guide = guide_legend(
                    override.aes = list(color = "black", size = 0.5))) +
  theme_bw(base_size = 10) +
  scale_x_continuous(expand = c(0, 0), breaks = tg_bin_sizes$cum_bins, labels = c(1:22)) +
  scale_y_discrete(expand = c(0, 0)) +
  facet_wrap(~Progression, ncol = 1, scales = "free") +
  labs(x = "Chromosomes", y = NULL, fill = NULL) +
  theme(axis.ticks = element_blank(), axis.text.x = element_text(size = 5, hjust = 1),
        axis.text.y = element_blank(), axis.title.y = element_text(size = 9),
        legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_blank(), strip.text = element_text(face = "bold"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 9),
        plot.title = element_text(size = 9, face = "bold"),
        legend.key.size = unit(0.4, "cm"))

  ggsave("plots/copy_number/trevor_graham/tg_cn_heatmap.pdf", p, width = 10, heigh = 8)

# frequency plot
freq_df <- tidy_tg_calls |>
  mutate(gain = cn_state %in% c("1", "2"), loss = cn_state %in% c("-1", "-2")) |>
  group_by(Progression, bin) |>
  summarise(gain_freq = mean(gain), loss_freq = mean(loss))

p <- ggplot(freq_df) +
  geom_area(aes(x = bin, y = gain_freq), fill = "#B04968", alpha = 0.6) +
  geom_area(aes(x = bin, y = -loss_freq), fill = "#4776A2", alpha = 0.6) +
        scale_x_continuous(
            expand = c(0,0),
            breaks = tg_bin_sizes$cum_bins,
            labels = 1:22) +
  facet_wrap(~Progression, ncol =1) +
  theme_bw() +
  labs(y = "Frequency", x = "Chromosome") +
  theme(strip.background = element_blank())
ggsave("plots/copy_number/test_tg.pdf", p, width = 12, height = 5)