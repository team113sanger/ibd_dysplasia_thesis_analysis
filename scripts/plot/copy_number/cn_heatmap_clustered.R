# Load libraries
library(readr)
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)

# ================================
# 1. Read data
# ================================

cnv_mat <- read_tsv("results/copy_number/precursors/cnv_calls.txt")
clusters <- read_tsv("results/copy_number/precursors/cn_clusters.tsv")
metadata <- read_tsv("metadata/final_metadata_qc_pass.tsv")

# ================================
# 2. Prepare matrices
# ================================
# Remove bin column
num_mat <- as.matrix(cnv_mat[ , -1])
rownames(num_mat) <- cnv_mat$bin
colnames(num_mat) <- colnames(cnv_mat)[-1]

# Merge clusters with metadata
annot_data <- clusters |>
  left_join(metadata, by = c("sample" = "sanger_dna_id"))

# ================================
# 3. Prepare chromosomes
# ================================
chrom_sizes <- read_csv("metadata/rescources/GRCh38_chrom_sizes.csv") |>
  rename(chr = `#Chr`) |>
  filter(!chr %in% c("chrX", "chrY")) |>
  mutate(cum_size = cumsum(Size))
bin_size <- 100000

chrom_sizes$n_bins <- floor(chrom_sizes$Size / bin_size)

chr_vector <- rep(chrom_sizes$chr, times = round(chrom_sizes$Size / bin_size))

# Make same size... correct later
chr_vector <- c(chr_vector, "chr22", "chr22")
chr_vector <- factor(chr_vector, levels = paste0("chr", 1:22))

# Check
length(chr_vector)
table(chr_vector)

# ================================
# 4. Define annotations
# ================================
# Top annotation: cluster labels
top_annotation <- HeatmapAnnotation(
  "CN cluster" = annot_data$`2`,  # <- or use `3` or `8` depending on cut
  col = list(
    "CN cluster" = c(
      "1" = "#225ea8",
      "2" = "#41b6c4",
      "3" = "#a1dab4",
      "4" = "#ffffcc" 
    )
  ),
  annotation_name_side = "left"
)

# Bottom annotation: metadata
bottom_annotation <- HeatmapAnnotation(
  "Group" = annot_data$group,
  "Grade" = annot_data$grade_of_dysplasia,
  "Sex" = annot_data$sex,
  col = list(
    "Group" = c(
      "Non-progressor" = "#6aa84f",  
      "Progressor"     = "#e69138"   
    ),
    "Grade" = c(
      "Low grade"       = "#edf8b1", 
      "Moderate"        = "#7fcdbb", 
      "High grade"      = "#225ea8", 
      "NOS"             = "#cccccc",  
      "Not specified"   = "#cccccc"   
    ),
    "Sex" = c(
      "Male"   = "#7866c2ff",  
      "Female" = "#a9cb8dff"   
    )
  ),
  annotation_name_side = "left"
)
# ================================
# 5. Build heatmap
# ================================
cnv_colors <- c(
  "0" = "#FFFFFF",      # No change (white)
  "1" = "#BB3737",      # Gain (red)
  "2" = "#4176A1",      # Loss (blue)
  "3" = "#66c2b7ff"       # cn-LOH (orange)
)

HM <- Heatmap(
  num_mat,
  name = "CNV",
  col = cnv_colors,
  cluster_rows = FALSE,        # keep bins in genomic order
  cluster_columns = FALSE,     # do not cluster samples
  show_column_names = FALSE,
  show_row_names = FALSE,
  top_annotation = top_annotation,
  bottom_annotation = bottom_annotation,
  border = TRUE,
  column_split = annot_data$`2`,  # split samples by cluster assignment
  column_title = NULL,     
  row_split = chr_vector,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 10),
  row_gap = unit(0, "points")
)

# ================================
# 6. Save plot
# ================================
pdf("plots/copy_number/heatmap/precursors/cn_clustered_heatmap_2.pdf", width = 5.5, height = 7)
draw(
    HM,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legend = TRUE)
dev.off()

