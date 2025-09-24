source("scripts/plot/oncoplots/oncoplot_utils.R")

library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(cowplot)
library(ggpubr)
library(tibble)

MAF_PATH <- "data/variants/7100_3235-filtered_mutations_matched_allTum_keepPA.maf"
METADATA_PATH <- "metadata/final_metadata_qc_pass.tsv"
TMB_PATH <- "data/variants/mutations_per_Mb.tsv"
SAMPLE_LIST_PATH <- "metadata/sample_lists/non_progressor_precursor_samples_ppat.tsv"
TP53_LOH <- "results/tables/tp53_cn_loh.tsv"
APC_LOH <- "results/tables/apc_cn_loh.tsv"
CN_PROP <- "data/copy_number/proportions/all_cn_props.tsv"


######## Read in data ########
maf <- read_tsv(MAF_PATH)
sample_list <- read_lines(SAMPLE_LIST_PATH) 
metadata <- read_tsv(METADATA_PATH)
tmb <- read_tsv(TMB_PATH, col_names = FALSE) |>
  filter(X1 %in% sample_list)
tp53_loh <- read_tsv(TP53_LOH)
apc_loh <- read_tsv(APC_LOH)
cn_props <- read_tsv(CN_PROP)


# Define genes to plot
plot_genes <- c("TP53", "APC", "KRAS", "RNF43", "RBM10", "LRP1B", "FBXW7", "PIK3CA", "SMAD4", "IDH1", "MSH3", "POLD1", "MLH3")

######## Prepare Variants Data ########
# Rename consequences in Main_consequence_VEP
variants_short <- shorten_consequence(maf)

# Select samples, only oncoKB genes and determine min hits required
variants_filt <- variants_short |>
  filter(Tumor_Sample_Barcode %in% sample_list) |>
  filter(Hugo_Symbol %in% plot_genes) |>
  group_by(Hugo_Symbol) |>
  filter(n_distinct(Tumor_Sample_Barcode) >= 1) |>
  ungroup() |>
  select(Hugo_Symbol, Tumor_Sample_Barcode, Main_consequence_VEP)

# Remove duplicates
variants_unique <- remove_duplicates(variants_filt)

# Get order for plotting genes
genes_order <- get_order(variants_unique, Hugo_Symbol)

# Get sample order for plotting
samples_order <- variants_unique |>
  group_by(Hugo_Symbol) |>
  mutate(n = n()) |>
  arrange(desc(n)) |>
  pull(Tumor_Sample_Barcode)

samples_order <- unique(samples_order)

# Add missing samples
variants_unique_full <- add_missing_samples(sample_list, variants_unique, "No mutation")
samples_order <- c(samples_order, missing_samples)

# Set levels
variants_unique_full[["Tumor_Sample_Barcode"]] <-
  factor(variants_unique_full[["Tumor_Sample_Barcode"]], levels = samples_order)

variants_unique_full[["Hugo_Symbol"]] <-
  factor(variants_unique_full[["Hugo_Symbol"]], levels = genes_order)

# Expand variants for plotting
variants_expanded <- expand_dataframe(variants_unique_full, "No mutation")

######## Variant Plot #########
consequence_palette <- c(
  "Missense" = "#68a334ff",
  "Nonsense" = "#9f2727ff",
  "Splice site" = "#E6AB02",
  "Frameshift" = "mediumpurple4",
  "Inframe indel" = "#7d6a56ff",
  "Multi hit" = "#EF9A9A",
  "Start codon lost" = "dodgerblue4",
  "No mutation" = "snow3"
)
legend_order <- names(consequence_palette)
# pclo_y <- which(levels(variants_expanded$Hugo_Symbol) == "PCLO")

variants_plot <- ggplot(data = variants_expanded) +
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
  labs(fill = "Consequence") 
  # facet_grid(~Pipeline, scales = "free_x")
# coord_fixed()

variants_plot
variants_plot_leg <- as_ggplot(get_legend(variants_plot + theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9))))

######## Prepare CN Props data ########
cn_props <- cn_props |>
  select(Sample, proportion) |>
  filter(Sample %in% sample_list)

cnprops_long <- cn_props |>
  pivot_longer(cols = c(proportion), names_to = "Category", values_to = "Value") |>
  mutate(Category = "CN Proportion") |>
  bind_rows(
    tibble(
      Sample   = "PD62077a",
      Category = "CN Proportion",
      Value    = NA_real_  
    )
  )

cnprops_long[["Sample"]] <-
  factor(cnprops_long[["Sample"]], levels = samples_order)

######## CN Props Plot ########
cnprops_plot <- ggplot(data = cnprops_long) +
  geom_tile(aes(
    x = Sample, y = Category,
    fill = Value
  ), color = "white", lwd = 0.3) +
  scale_fill_viridis_c(direction = -1, limits = c(0, 1), na.value = "snow2") +
  scale_color_identity(guide = "none") +
  theme(
    axis.ticks = element_blank(),
    axis.title.x = element_blank(), axis.text.x = element_blank(),
    legend.justification = c("left", "top"),
    axis.title.y = element_blank(),
    legend.position = "right", legend.box = "horizontal",
    panel.background = element_rect(fill = "snow3"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  ) +
  scale_x_discrete(expand = c(0, 0), guide = guide_axis(angle = 90)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(fill = "CN Proportion Scale")
# coord_fixed()

cnprops_plot
cnprops_plot_leg <- as_ggplot(get_legend(cnprops_plot + theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9))))

######## Prepare LOH data ########
## TP53 ##
loh_data <- tp53_loh |>
  rbind(apc_loh) |>
  select(gene, sample, cn) |>
  filter(sample %in% sample_list)

loh_complete <- add_missing_samples(sample_list, loh_data, "No LOH")
loh_expand <- expand_dataframe(loh_complete, "No LOH") |>
  mutate(Category = gene) |>
  select(sample, Category, cn)

loh_expand[["sample"]] <-
  factor(loh_expand[["sample"]], levels = samples_order)

######## CN LOH Plot ########
loh_palette <- c(
  # "cn-LOH" = "palegreen3",
  "LOH-del" = "#925E9F",
  "cn-LOH" = "#79AF97",
  "LOH-amp" = "rosybrown1",
  "No LOH" = "snow3",
  "No data" = "snow2"
)
legend_order <- names(loh_palette)

loh_plot <- ggplot(data = loh_expand) +
  geom_tile(aes(
    x = sample, y = Category,
    fill = cn
  ), color = "white", lwd = 0.3) +
  scale_fill_manual(values = loh_palette, breaks = legend_order) +
  scale_color_identity(guide = "none") +
  theme(
    axis.ticks = element_blank(),
    axis.title.x = element_blank(), axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(face = "italic", size = 8),
    legend.position = "right", legend.box = "horizontal",
    legend.justification = c("left", "top"),
    panel.background = element_rect(fill = "snow3"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  ) +
  scale_x_discrete(expand = c(0, 0), guide = guide_axis(angle = 90)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(fill = "LOH Status")
# coord_fixed()

loh_plot
loh_plot_leg <- as_ggplot(get_legend(loh_plot + theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9))))


######## Prepare Metadata ########
meta <- metadata |>
  filter(sanger_dna_id %in% sample_list) |>
  dplyr::rename(
    Sample = sanger_dna_id
  ) |>
  select(Sample, study_id, ibd_diagnosis, site_general, sex, psc, pancolitis)

meta <- meta |>
  mutate(
    study_id = factor(
      study_id,
      levels = study_id[match(samples_order, Sample)]
    )
  )

metadata_expanded <- meta |>
  pivot_longer(
    cols = -c(Sample, study_id), 
    names_to = "Category",
    values_to = "Value"
  ) |>
  mutate(
    Category = str_replace_all(Category, "_", " "),  
    Category = str_to_title(Category)               
  )

######## Metadata Plots ########
metadata_palette <- c(
  # Sex
  "Male"   = "#AED581",
  "Female" = "#A1887F",
  
  # PSC (Yes/No)
  "Yes" = "#9575CD",
  "No"  = "gray80",
  
  # Diagnosis
  "UC"     = "#83A8CC",
  "Crohn's" = "#DF8F44",
  "IBDU"   = "gray60",
  
  # Site
  "Left colon"       = "#B1746F",
  "Right colon"      = "#A5D6A7",
  "Transverse colon" = "#EFC000",
  "Rectum"           = "#7986CB"
)

legend_order <- names(metadata_palette)

metadata_plot <- ggplot(data = metadata_expanded) +
  geom_tile(aes(
    x = study_id, y = Category,
    fill = Value
  ), color = "white", lwd = 0.3, height = 0.7) +
  scale_fill_manual(values = metadata_palette, breaks = legend_order) +
    theme(
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 8),
    legend.position = "bottom", legend.box = "horizontal",
    legend.justification = c("left", "top"),
    legend.text = element_text(size = 8), legend.title = element_text(size = 9),
    panel.background = element_rect(fill = "white"),
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

tmb[["X1"]] <- factor(tmb[["X1"]], levels = samples_order)

tmb_plot <- ggplot(tmb, aes(x = X1, y = X2)) +
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
  cnprops_plot_leg, loh_plot_leg, variants_plot_leg
)

all_plots <- list(
  tmb_plot,
  cnprops_plot + theme(legend.position = "none"),
  loh_plot + theme(legend.position = "none"),
  variants_plot + theme(legend.position = "none"),
  metadata_plot + theme(legend.position = "none")
)

p <- plot_grid(plotlist = all_plots, ncol = 1, align = "v", axis = "lr", rel_heights = c(1.2, 0.5, 0.8, 4, 2.5))
legend_p1 <- plot_grid(plotlist = all_legends, ncol = 1, align = "vh", axis = "l", rel_heights = c(1.2, 0.8, 2))
legend_p2 <- plot_grid(legend_p1, metadata_plot_leg, ncol = 2, rel_heights = c(1, 0.75))
plot_grid(p, legend_p2, ncol = 2, rel_widths = c(1, 0.7))

ggsave("plots/oncoplots/combined_pre_npro.pdf",
  height = 5, width = 6.5
)

