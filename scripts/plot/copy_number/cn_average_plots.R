library(ggplot2)
library(cowplot)
library(readr)
library(dplyr)

# load in chromosome sizes
chrom_sizes <- read.table("data/other/GRCh38_chrom_sizes.csv", header = 0, stringsAsFactors = F, sep = ",")
chrom_sizes[, 3] <- cumsum(as.numeric(chrom_sizes[, 2]))
chrom_sizes$V1 <- gsub("^chr", "", chrom_sizes$V1)

# Load in CN segments
non_prog_segments <- read_tsv("/lustre/scratch125/casm/team113da/projects/IBD_Associated_Dysplasia/7100_3235_IBD-associated_dysplasia/analysis/ASCAT/release_v1/non_progressor_precursors/PLOTS_ONE_PER_PATIENT/7100_3235_segments.tsv")
prog_segments <- read_tsv("/lustre/scratch125/casm/team113da/projects/IBD_Associated_Dysplasia/7100_3235_IBD-associated_dysplasia/analysis/ASCAT/release_v1/progressor_precursors/PLOTS_ONE_PER_PATIENT/7100_3235_segments.tsv")

non_prog_cn_loh_segments <- read_tsv("/lustre/scratch125/casm/team113da/projects/IBD_Associated_Dysplasia/7100_3235_IBD-associated_dysplasia/analysis/ASCAT/release_v1/non_progressor_precursors/PLOTS_ONE_PER_PATIENT/7100_3235_cn-loh_segments.tsv") |>
  mutate(CN = "cn-LOH")
prog_cn_loh_segments <- read_tsv("/lustre/scratch125/casm/team113da/projects/IBD_Associated_Dysplasia/7100_3235_IBD-associated_dysplasia/analysis/ASCAT/release_v1/progressor_precursors/PLOTS_ONE_PER_PATIENT/7100_3235_cn-loh_segments.tsv") |>
  mutate(CN = "cn-LOH")

non_prog_segments <- non_prog_segments |>
  rows_update(non_prog_cn_loh_segments, by = c("Sample", "chr", "startpos", "endpos"))
prog_segments <- prog_segments |>
  rows_update(prog_cn_loh_segments, by = c("Sample", "chr", "startpos", "endpos"))

combined_segments <- bind_rows(
  non_prog_segments %>% mutate(Group = "Non-Progressor"),
  prog_segments %>% mutate(Group = "Progressor")
)

# Create a sample list (vector) of unique sample IDs
metadata <- read_tsv("metadata/final_metadata_qc_pass.tsv")

precursors <- metadata %>%
  filter(precursor_or_follow_up == "Precursor") %>%
  pull(sanger_dna_id)

# Sample list 2: Follow-up samples
follow_ups <- metadata %>%
  filter(precursor_or_follow_up == "Follow up") %>%
  pull(sanger_dna_id)

sample_list <- intersect(precursors, unique(combined_segments$Sample))
sample_list <- setdiff(sample_list, c("PD62031d", "PD62038c"))

# Prep
bin_size <- 100000
filter_size <- 1000000

cnv_calls <- array(data = 0, dim = c((chrom_sizes[24, 3] / bin_size) + 1, length(sample_list)))
colnames(cnv_calls) <- unique(sample_list)

combined_segments$startcumpos <- combined_segments$startpos
combined_segments$endcumpos <- combined_segments$endpos

for (i in 1:nrow(combined_segments)) {
  if (combined_segments$chr[i] != "1") {
    combined_segments$startcumpos[i] <- combined_segments$startpos[i] + chrom_sizes$V3[which(chrom_sizes$V1 == combined_segments$chr[i]) - 1]
    combined_segments$endcumpos[i] <- combined_segments$endpos[i] + chrom_sizes$V3[which(chrom_sizes$V1 == combined_segments$chr[i]) - 1]
  }
}

for (i in 1:length(sample_list)) {
  # Initialize vectors for gains and losses
  gain_temp <- vector()
  loss_temp <- vector()
  loh_temp <- vector()

  # Filter the combined_segments for the current sample
  sample_segments <- combined_segments[combined_segments$Sample == sample_list[i], ]

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

cnv_calls_1 <- cnv_calls[1:28776, ]

# Single Plot
tidy_cnv_calls <- reshape2::melt(cnv_calls_1)
colnames(tidy_cnv_calls) <- c("bin", "sample", "cn_state")

tidy_cnv_calls$cn_state <- as.character(tidy_cnv_calls$cn_state)
tidy_cnv_calls <- tidy_cnv_calls %>%
  left_join(combined_segments %>%
    select(Sample, Group) %>%
    distinct(), by = c("sample" = "Sample"))

p <- ggplot(data = tidy_cnv_calls) +
  geom_tile(aes(x = bin, y = sample, fill = cn_state)) +
  geom_hline(yintercept = ((1:length(sample_list)) - 0.5), color = "gray90") +
  geom_vline(xintercept = chrom_sizes$V3[1:22] / bin_size) +
  facet_wrap(~Group, nrow = 2, scales = "free_y") +
  scale_fill_manual(
    values = c("0" = "white", "1" = "lightsalmon2", "2" = "skyblue2", "3" = "mediumaquamarine"),
    labels = c("0" = "No change", "2" = "Loss", "1" = "Gain", "3" = "cn-LOH"),
    guide = guide_legend(
      override.aes = list(color = "black", size = 0.5) # Add outline to legend squares
    )
  ) +
  theme_bw() +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = c(chrom_sizes$V3[1:22] - (chrom_sizes$V2[1:22]) / 2) / bin_size,
    labels = c(1:22)
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "", y = "", fill = "Copy number") +
  theme(
    axis.ticks = element_blank(),
    legend.position = "top"
  )

ggsave("plots/CN/cn_heatmap.pdf", plot = p, width = 10, height = length(sample_list) / 5)
ggsave("plots/CN/cn_heatmap.png", plot = p, width = 8, height = length(sample_list) / 5)


# # By variable
# groups <- unique(combined_segments$Group)

# matrix <- array(data = 0, dim = c((chrom_sizes[24,3] / bin_size) + 1, length(groups)))
# colnames(matrix) <- groups

# per_variable_gains <- matrix
# per_variable_gains_summary <- matrix
# per_variable_losses <- matrix
# per_variable_losses_summary <- matrix

# for(i in 1:length(groups)){
#   per_variable_sample_ids <- unique(combined_segments$Sample[combined_segments$Group == groups[i]])

#   cnv_calls_variable_inc <- cnv_calls[,per_variable_sample_ids, drop = F]

#   gains_variable_inc <- cnv_calls_variable_inc
#   gains_variable_inc[gains_variable_inc == 1] <- 1
#   gains_variable_inc[gains_variable_inc != 1] <- 0

#   per_variable_gains[,i] <- rowSums(gains_variable_inc)
#   per_variable_gains_summary[,i] <- rowSums(gains_variable_inc) / length(per_variable_sample_ids)

#   losses_variable_inc <- cnv_calls_variable_inc
#   losses_variable_inc[losses_variable_inc != 2] <- 0
#   losses_variable_inc[losses_variable_inc == 2] <- 1

#   per_variable_losses[,i] <- rowSums(losses_variable_inc)
#   per_variable_losses_summary[,i] <- rowSums(losses_variable_inc) / length(per_variable_sample_ids)

# }

# per_variable_gains <-  per_variable_gains[1:28776,]
# per_variable_losses <-  per_variable_losses[1:28776,]

# per_variable_gains_summary <-  per_variable_gains_summary[1:28776,]
# per_variable_losses_summary <-  per_variable_losses_summary[1:28776,]


# tidy_variable_gains_summary <- reshape2::melt(per_variable_gains_summary)
# tidy_variable_losses_summary <- reshape2::melt(per_variable_losses_summary)


# colnames(tidy_variable_gains_summary) <- c("Chromosome", "Variable", "Gains")
# colnames(tidy_variable_losses_summary) <- c("Chromosome", "Variable", "Losses")

# #Plotting
# #gains only
#  p1 <-  ggplot(data = tidy_variable_gains_summary) +
#         geom_tile(aes(x = Chromosome, y = Variable, fill = Gains)) +
#         scale_fill_viridis_c(option = "rocket", direction = -1, limits = c(0, 1)) +
#         geom_hline(yintercept = ((1:length(groups)) - 0.5), color = "white") +
#         geom_vline(xintercept = chrom_sizes$V3[1:22] / bin_size) +
#         theme_bw() +
#         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#         scale_x_continuous(expand = c(0,0), breaks = c(chrom_sizes$V3[1:22] - (chrom_sizes$V2[1:22])/2) / bin_size, labels = c(1:22)) +
#         scale_y_discrete(expand = c(0,0)) +
#         labs(fill = "Copy number", title = "Gains") +
#         ylab(var_title) +
#         theme(axis.ticks = element_blank(), text = element_text(size = 15))
# #ggsave(filename = paste0("rcc_all_countries_only_gains_battenberg_no_dup.pdf"), width = 15, height = length(countries)/4)
# #ggsave(filename = paste0("rcc_all_countries_only_gains_battenberg_no_dup.png"), width = 15, height = length(countries)/4)

# #losses only
#   p2 <- ggplot(data = tidy_variable_losses_summary) +
#         geom_tile(aes(x = Chromosome, y = Variable, fill = Losses)) +
#         #facet_grid(cols = vars(Chromosome)) +
#         scale_fill_viridis_c(option = "mako", direction = -1, limits = c(0, 1)) +
#         geom_hline(yintercept = ((1:length(groups)) - 0.5), color = "white") +
#         geom_vline(xintercept = chrom_sizes$V3[1:22] / bin_size) +
#         theme_bw() +
#         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#         scale_x_continuous(expand = c(0,0), breaks = c(chrom_sizes$V3[1:22] - (chrom_sizes$V2[1:22])/2) / bin_size, labels = c(1:22)) +
#         scale_y_discrete(expand = c(0,0)) +
#         labs(fill = "Copy number", title = "Losses") +
#         ylab(var_title) +
#         theme(axis.ticks = element_blank(), text = element_text(size = 15))
# #ggsave(filename = paste0("rcc_all_countries_only_losses_battenberg_no_dup.pdf"), width = 15, height = 10)
# #ggsave(filename = paste0("test_rcc_all_countries_only_losses_battenberg_no_dup.png"), width = 15, height = length(countries)/4)


# plot_grid(p1, p2, labels = c('A', 'B'), nrow = 2, label_size = 12)
# ggsave(filename = paste0("plots/CN/copy_number_gains_and_losses.pdf"), width = 15, height = 6)
# #ggsave(filename = paste0("rcc_country_copy_number_gains_and_losses.png"), width = 15, height = length(countries))
