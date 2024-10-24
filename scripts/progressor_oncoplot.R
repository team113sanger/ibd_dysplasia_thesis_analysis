source("/lustre/scratch125/casm/team113da/projects/dermatlas_analysis_methods/somatic-variant-plots/src/oncoplot_utils.R")

library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(cowplot)
library(ggpubr)
library(tibble)

# Load data 
MAF_PATH <- "data/7100_3235-filtered_mutations_all_indepTum_keepPA.maf"
METADATA_PATH <- "metadata/processed_metadata.tsv"
TMB_PATH <- "data/mutations_per_Mb.tsv"
SAMPLE_LIST <- "metadata/sample_lists/qc_pass_samples_list.txt"

maf <- read_tsv(MAF_PATH)
metadata <- read_tsv(METADATA_PATH, col_types = cols(
    Sex = col_character()
  ))
tmb <- read_tsv(TMB_PATH, col_names = FALSE)
sample_list <- read_tsv(SAMPLE_LIST, col_names = F)

# Get sample lists
progressors <- metadata |>
    filter(group == "Progressor") |>
    filter(sanger_dna_id %in% sample_list[["X1"]]) |>
    pull(sanger_dna_id)

######## Prepare Variants Data ########
# Rename consequences in Main_consequence_VEP
variants_short <- shorten_consequence(maf)

# Select samples, only oncoKB genes and determine min hits required
variants_plot <- variants_short |>
  filter(Tumor_Sample_Barcode %in% progressors) |>
  # filter(Hugo_Symbol %in% oncoKB) |>
  group_by(Hugo_Symbol) |>
  filter(n_distinct(Tumor_Sample_Barcode) >= 5) |>
  ungroup() |>
  select(Hugo_Symbol, Tumor_Sample_Barcode, Main_consequence_VEP)

# Get order for plotting genes
genes_order <- get_order(variants_plot, Hugo_Symbol)

# Remove duplicates
variants_unique <- remove_duplicates(variants_plot)

# Add missing samples
variants_unique_full <- add_missing_samples(progressors, variants_unique, "No mutation")
samples_order <- c(samples_order, missing_samples)

# Set levels
variants_unique_full[["Hugo_Symbol"]] <-
  factor(variants_unique_full[["Hugo_Symbol"]], levels = genes_order)

# Expand variants for plotting
variants_expanded <- expand_dataframe(variants_unique_full, "No mutation")

# Add group info
meta_select <- metadata |>
  select(sanger_dna_id, precursor_or_follow_up)

variants_plot_data <- left_join(variants_expanded, meta_select,
                                by = c("Tumor_Sample_Barcode" = "sanger_dna_id"))

variants_plot_data$precursor_or_follow_up <- factor(
  variants_plot_data$precursor_or_follow_up,
  levels = c("Precursor", "Follow up")
)

# Get sample order for plotting
samples_order_1 <- variants_plot_data |>
  filter(precursor_or_follow_up == "Precursor") |>
  group_by(Hugo_Symbol) |>
  mutate(n = n()) |>
  arrange(desc(n)) |>
  pull(Tumor_Sample_Barcode)

samples_order_1 <- unique(samples_order_1)

# Get sample order for plotting
samples_order_2 <- variants_plot_data |>
  filter(precursor_or_follow_up == "Follow up") |>
  group_by(Hugo_Symbol) |>
  mutate(n = n()) |>
  arrange(desc(n)) |>
  pull(Tumor_Sample_Barcode)
  
samples_order_2 <- unique(samples_order_2)

samples_order <- c(samples_order_1, samples_order_2)

variants_plot_data[["Tumor_Sample_Barcode"]] <-
  factor(variants_plot_data[["Tumor_Sample_Barcode"]], levels = samples_order)

######## Variant Plot #########
consequence_palette <- c(
  "Missense" = "forestgreen",
  "Nonsense" = "firebrick",
  "Splice site" = "orange2",
  "Frameshift" = "mediumpurple3",
  "Inframe indel" = "sienna4",
  "Multi hit" = "#f7b4ae",
  "Start codon lost" = "dodgerblue4",
  "No mutation" = "snow3"
)
legend_order <- names(consequence_palette)
# pclo_y <- which(levels(variants_expanded$Hugo_Symbol) == "PCLO")

variants_plot <- ggplot(data = variants_plot_data) +
  geom_tile(aes(
    x = Tumor_Sample_Barcode, y = Hugo_Symbol,
    fill = Main_consequence_VEP
  ), color = "white", lwd = 0.3) +
  scale_fill_manual(values = consequence_palette, breaks = legend_order) +
  scale_color_identity(guide = "none") +
  theme(
    legend.position = "right", axis.ticks = element_blank(),
    legend.justification = c("left", "top"),
    axis.text.y = element_text(face = "italic", size = 8),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    panel.background = element_rect(fill = "snow3"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  ) +
  scale_x_discrete(expand = c(0, 0), guide = guide_axis(angle = 90)) +
  scale_y_discrete(expand = c(0, 0)) +
  facet_wrap(~ precursor_or_follow_up, scales = "free_x") +
  # geom_hline(yintercept = pclo_y + 0.5, color = "grey30",
  #           linewidth = 0.4, linetype = "dashed") +
  labs(fill = "Consequence")
# coord_fixed()

variants_plot
ggsave('plots/test.pdf')
variants_plot_leg <- as_ggplot(get_legend(variants_plot + theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9))))

######## Prepare Metadata ########
meta <- metadata |>
  filter(sanger_dna_id %in% progressors) |>
  rename(
    sample = sanger_dna_id
  ) |>
  select(sample, diagnosis, grade_of_dysplasia, ibd_diagnosis, precursor_or_follow_up)

meta[["sample"]] <-
  factor(meta[["sample"]], levels = samples_order)

metadata_expanded <- meta |>
  pivot_longer(cols = c(grade_of_dysplasia,
                        ibd_diagnosis, precursor_or_follow_up),
                names_to = "Category", values_to = "Value")

######## Metadata Plots ########
metadata_palette <- c(
  "Male" = "lightskyblue3",
  "Female" = "bisque3",
  "Extremities" = "lightgoldenrod",
  "Head and neck" = "rosybrown2",
  "Trunk" = "palegreen4",
  "Papillary Eccrine Adenoma" = "darkseagreen",
  "Apocrine Tubular Adenoma" = "mediumpurple3"
)
legend_order <- names(metadata_palette)

metadata_plot <- ggplot(data = metadata_expanded) +
  geom_tile(aes(
    x = sample, y = Category,
    fill = Value
  ), color = "white", lwd = 0.3) +
  # scale_fill_manual(values = metadata_palette, breaks = legend_order) +
  scale_fill_viridis_d() +
  theme(
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 8),
    legend.position = "bottom", legend.box = "horizontal",
    legend.justification = c("left", "top"),
    legend.text = element_text(size = 8), legend.title = element_text(size = 9),
    panel.background = element_rect(fill = "snow3"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  ) +
  scale_x_discrete(expand = c(0, 0), guide = guide_axis(angle = 90)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(fill = "Metadata")
# coord_fixed()

metadata_plot
metadata_plot_leg <- as_ggplot(
  get_legend(
    metadata_plot + theme(
      legend.position = "right", legend.text = element_text(size = 8),
      legend.title = element_text(size = 9)
    )
  )
)

######## TMB Plot #########
# missing_sample <- c("PD59552a", 0.00)
# tmb <- rbind(tmb, missing_sample)
tmb[["X2"]] <- as.numeric(tmb[["X2"]])
tmb_data <- tmb |>
    filter(X1 %in% progressors)

tmb_data[["X1"]] <- factor(tmb_data[["X1"]], levels = samples_order)

tmb_plot <- ggplot(tmb_data, aes(x = X1, y = X2)) +
  geom_bar(stat = "identity", fill = "grey40") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank()
  ) +
  scale_x_discrete(expand = c(0, 0), guide = guide_axis(angle = 90)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "TMB")

tmb_plot

######## Combine Plots ########
all_legends <- list(
  variants_plot_leg, metadata_plot_leg
)

all_plots <- list(
  tmb_plot,
  variants_plot + theme(legend.position = "none"),
  metadata_plot + theme(legend.position = "none")
)

p <- plot_grid(plotlist = all_plots, ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 5, 1.5))
legend_p1 <- plot_grid(plotlist = all_legends, ncol = 1, align = "vh", axis = "l")
#legend_p2 <- plot_grid(legend_p1, , ncol = 2, rel_heights = c(1, 0.75))
plot_grid(p, legend_p1, ncol = 2, rel_widths = c(1, 0.65))

ggsave("plots/progressors_oncoplot.pdf",
  height = 5.5, width = 7
)