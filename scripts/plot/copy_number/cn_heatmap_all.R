library(ggplot2)
library(cowplot)
library(readr)
library(dplyr)
library(tibble)

# load in chromosome sizes
chrom_sizes <- read.table("metadata/rescources/GRCh38_chrom_sizes.csv", header = 0, stringsAsFactors = F, sep = ",")
chrom_sizes[, 3] <- cumsum(as.numeric(chrom_sizes[, 2]))
chrom_sizes$V1 <- gsub("^chr", "", chrom_sizes$V1)

# Read in metadata
meta <- read_tsv("metadata/final_metadata_qc_pass.tsv") |>
  select(sanger_dna_id, study_id, group, precursor_or_follow_up) |>
  mutate(facet_group = paste(group, precursor_or_follow_up, sep = " "))

samples <- read_lines("metadata/sample_lists/all_one_ppat.list")

# Load in CN segments
segments <- read_tsv("data/copy_number/segments/all_segments.tsv")

loh_segments <- read_tsv("data/copy_number/segments/all_cn-loh_segments.tsv") |>
  mutate(CN = "cn-LOH")

# Create a sample list (vector) of unique sample IDs
sample_list <- intersect(samples, unique(segments$Sample))

segments <- segments |>
  rows_update(loh_segments, by = c("Sample", "chr", "startpos", "endpos")) |>
  filter(Sample %in% sample_list)

# Prep
bin_size <- 100000
filter_size <- 1000000

cnv_calls <- array(data = 0, dim = c((chrom_sizes[24, 3] / bin_size) + 1, length(sample_list)))
colnames(cnv_calls) <- unique(sample_list)

segments$startcumpos <- segments$startpos
segments$endcumpos <- segments$endpos

# Call segments
for (i in 1:nrow(segments)) {
  if (segments$chr[i] != "1") {
    segments$startcumpos[i] <- segments$startpos[i] + chrom_sizes$V3[which(chrom_sizes$V1 == segments$chr[i]) - 1]
    segments$endcumpos[i] <- segments$endpos[i] + chrom_sizes$V3[which(chrom_sizes$V1 == segments$chr[i]) - 1]
  }
}

for (i in 1:length(sample_list)) {
  # Initialize vectors for gains and losses
  gain_temp <- vector()
  loss_temp <- vector()
  loh_temp <- vector()

  # Filter the segments for the current sample
  sample_segments <- segments[segments$Sample == sample_list[i], ]

  # Loop through each row of the filtered sample segments
  for (j in 1:nrow(sample_segments)) {
    # Check the CN status for the current segment
    if (sample_segments$CN[j] == "gain") {
      # Add bins for gain
      gain_temp <- c(gain_temp, floor(sample_segments$startcumpos[j] / bin_size):(floor(sample_segments$endcumpos[j] / bin_size) + 1))
    } else if (sample_segments$CN[j] == "loss") {
      # Add bins for loss
      loss_temp <- c(loss_temp, floor(sample_segments$startcumpos[j] / bin_size):(floor(sample_segments$endcumpos[j] / bin_size) + 1))
    } else if (sample_segments$CN[j] == "cn-LOH") {
      loh_temp <- c(loh_temp, floor(sample_segments$startcumpos[j] / bin_size):(floor(sample_segments$endcumpos[j] / bin_size) + 1))
    }
  }

  # Update the cnv_calls matrix for the current sample
  cnv_calls[unique(gain_temp), sample_list[i]] <- 1
  cnv_calls[unique(loss_temp), sample_list[i]] <- 2
  cnv_calls[unique(loh_temp), sample_list[i]] <- 3
}

cnv_calls_1 <- cnv_calls[1:28751, ]

# Write out cnv
cnv_calls_1_df <- as.data.frame(cnv_calls_1) |>
  rownames_to_column("bin")

#write_tsv(cnv_calls_1_df, "results/copy_number/cnv_all_calls.txt")

# Single Plot
tidy_cnv_calls <- reshape2::melt(cnv_calls_1)
colnames(tidy_cnv_calls) <- c("bin", "sample", "cn_state")

tidy_cnv_calls$cn_state <- as.character(tidy_cnv_calls$cn_state)
tidy_cnv_calls <- tidy_cnv_calls |>
  left_join(meta, by = c("sample" = "sanger_dna_id"))

precursor_calls <- tidy_cnv_calls |>
  filter(precursor_or_follow_up == "Precursor")

follow_up_calls <- tidy_cnv_calls |>
  filter(precursor_or_follow_up == "Follow up")

pp <- ggplot(data = precursor_calls) +
  geom_tile(aes(x = bin, y = study_id, fill = cn_state)) +
  geom_hline(yintercept = ((1:length(sample_list)) - 0.5), color = "gray80", size = 0.3) +
  geom_vline(xintercept = chrom_sizes$V3[1:22] / bin_size, color = "gray80", size = 0.3) +
  facet_wrap(~group, scales = "free_y", ncol = 1) +
  scale_fill_manual(
    values = c("0" = "white", "1" = "#B04968", "2" = "#4776A2", "3" = "#69C2A7"),
    labels = c("0" = "No change", "2" = "Loss", "1" = "Gain", "3" = "cn-LOH"),
    guide = guide_legend(
      override.aes = list(color = "black", size = 0.5),
    )
  ) +
  theme_bw(base_size = 10) +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = c(chrom_sizes$V3[1:22] - (chrom_sizes$V2[1:22]) / 2) / bin_size,
    labels = c(1:22)
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Chromosomes", y = NULL, fill = NULL, title = "Precursor") +
  theme(
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 6.5),
    legend.position = "none",
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(size = 9, face = "bold"),
    legend.key.size = unit(0.4, "cm"),
    axis.title.y = element_text(size = 9) 
  )

# ggsave("plots/copy_number/heatmap/cn_heatmap_all.pdf", plot = p, width = 10, height = length(sample_list) / 6)
ggsave("plots/copy_number/heatmap/precursors/cn_heatmap_pre.png", plot = p, width = 5.5, height = 5, dpi = 300)

pf <- ggplot(data = follow_up_calls) +
  geom_tile(aes(x = bin, y = study_id, fill = cn_state)) +
  geom_hline(yintercept = ((1:length(sample_list)) - 0.5), color = "gray80", size = 0.3) +
  geom_vline(xintercept = chrom_sizes$V3[1:22] / bin_size, color = "gray80", size = 0.3) +
  facet_wrap(~group, scales = "free_y", ncol = 1) +
  scale_fill_manual(
    values = c("0" = "white", "1" = "#B04968", "2" = "#4776A2", "3" = "#69C2A7"),
    labels = c("0" = "No change", "2" = "Loss", "1" = "Gain", "3" = "cn-LOH"),
    guide = guide_legend(
      override.aes = list(color = "black", size = 0.5),
    )
  ) +
  theme_bw(base_size = 10) +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = c(chrom_sizes$V3[1:22] - (chrom_sizes$V2[1:22]) / 2) / bin_size,
    labels = c(1:22)
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Chromosomes", y = NULL, fill = "Copy number", title = "Follow Up") +
  theme(
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 6.5),
    legend.position = "bottom",
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    plot.title = element_text(size = 9, face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.key.size = unit(0.4, "cm"),
    axis.title.y = element_text(size = 9) 
  )

# ggsave("plots/copy_number/heatmap/cn_heatmap_all.pdf", plot = p, width = 10, height = length(sample_list) / 6)
ggsave("plots/copy_number/heatmap/follow_ups/cn_heatmap_fol.png", plot = p, width = 5.5, height = 5, dpi = 300)

library(patchwork)
combined <- pp / pf  

ggsave("plots/copy_number/heatmap/cn_heatmap_combined.png",
       plot = combined, width = 6.5, height = 9.2, dpi = 300)