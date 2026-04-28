library(readr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(patchwork)

# Load data
tg_calls <- read_tsv("data/copy_number/trevor_graham/tg_discovery_calls.tsv") |>
  tibble::column_to_rownames("Sample") |>
  t() |>
  as.data.frame() |>
  tibble::rownames_to_column("Bin")
tg_bin_sizes <- read_tsv("data/copy_number/trevor_graham/tg_bin_sizes.tsv")
tg_samples <- read_csv("data/copy_number/trevor_graham/tg_discovery_samples.csv") |>
  rename(sample = Sample)

cnv_calls <- read_tsv("results/copy_number/cnv_all_calls.txt")
metadata  <- read_tsv("metadata/final_metadata_qc_pass.tsv")
our_bin_sizes <- read_tsv("metadata/rescources/our_bin_sizes.tsv")
nprog_pre <- read_tsv("metadata/sample_lists/non_progressor_precursor_samples_ppat.tsv",
                      col_names = "sample") |> mutate(Progression = "Non-Progressor")
prog_pre  <- read_tsv("metadata/sample_lists/progressor_precursor_samples_ppat.tsv",
                      col_names = "sample") |> mutate(Progression = "Progressor")
our_sample_groups <- bind_rows(nprog_pre, prog_pre)

# Tidy helpers
make_freq_df <- function(tidy_calls, gain_states, loss_states) {
  tidy_calls |>
    mutate(
      gain = cn_state %in% gain_states,
      loss = cn_state %in% loss_states
    ) |>
    group_by(Progression, bin) |>
    summarise(gain_freq = mean(gain), loss_freq = mean(loss), .groups = "drop")
}

freq_plot <- function(freq_df, bin_sizes, title) {
  ggplot(freq_df) +
    geom_area(aes(x = bin, y =  gain_freq), fill = "#B04968", alpha = 0.6) +
    geom_area(aes(x = bin, y = -loss_freq), fill = "#4776A2", alpha = 0.6) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
    scale_x_continuous(
      expand = c(0, 0),
      breaks = bin_sizes$cum_bins,
      labels = 1:22
    ) +
    scale_y_continuous(
      limits = c(-1, 1),
      labels = abs
    ) +
    facet_wrap(~Progression, ncol = 1) +
    theme_bw() +
    labs(title = title, y = "Frequency", x = "Chromosome") +
    theme(
      strip.background = element_blank(),
      plot.title       = element_text(face = "bold", size = 11)
    )
}

# TG tidy
tidy_tg_calls <- melt(tg_calls, id.vars = "Bin") |>
  rename(bin = Bin, sample = variable, cn_state = value) |>
  mutate(
    cn_state = as.character(cn_state),
    bin      = as.numeric(sub("^V", "", bin)),
    sample   = sub(".*\\.", "", sample)
  ) |>
  left_join(tg_samples |> select(sample, Progression), by = "sample") |>
  filter(!is.na(Progression))

# Our data tidy
tidy_our_calls <- melt(cnv_calls, id.vars = "bin") |>
  rename(sample = variable, cn_state = value) |>
  mutate(
    cn_state = as.character(cn_state),
    bin      = as.numeric(bin)
  ) |>
  inner_join(our_sample_groups, by = "sample")  # filters to prog/non-prog pre only

# Freq data
freq_tg  <- make_freq_df(tidy_tg_calls,  gain_states = c("1", "2"),  loss_states = c("-1", "-2"))
freq_our <- make_freq_df(tidy_our_calls, gain_states = c("1"),        loss_states = c("2"))

# Plots
p_tg  <- freq_plot(freq_tg,  tg_bin_sizes,  "TG Discovery Cohort")
p_our <- freq_plot(freq_our, our_bin_sizes, "Our Cohort")

p_combined <- p_tg / p_our +
  plot_layout(heights = c(2, 2))  # adjust if one cohort has more groups

ggsave("plots/copy_number/cnv_freq_tg_comparison.pdf", p_combined, width = 10, height = 8)