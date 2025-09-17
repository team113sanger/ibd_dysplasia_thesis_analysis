source("scripts/plot/oncoplots/oncoplot_utils.R")

library(dplyr)
library(tidyr)
library(readr)
library(ComplexHeatmap)
library(grid)

# Load in data
oncoKB <- read_tsv("metadata/rescources/cancerGeneList.tsv") |>
  pull(`Hugo Symbol`)

metadata <- read_tsv("metadata/final_metadata_qc_pass.tsv") |>
  filter(group == "Progressor")

metadata$precursor_or_follow_up <- factor(
  metadata$precursor_or_follow_up,
  levels = c("Precursor", "Follow up")
)

progressors <- metadata |>
  pull(sanger_dna_id)

tmb_raw <- read_tsv("data/variants/mutations_per_Mb.tsv", col_names = c("Sample", "TMB")) |>
  filter(Sample %in% progressors)

maf <- read_tsv("data/variants/7100_3235-filtered_mutations_all_indepTum_keepPA.maf") |>
  filter(Hugo_Symbol %in% oncoKB) |>
  filter(Tumor_Sample_Barcode %in% progressors) |>
  group_by(Hugo_Symbol) |>
  filter(n_distinct(Tumor_Sample_Barcode) >= 3) |>
  ungroup() |>
  select(Hugo_Symbol, Tumor_Sample_Barcode, Main_consequence_VEP) |>
  shorten_consequence() |>
  left_join(metadata |> select(sanger_dna_id, study_id),
    by = c("Tumor_Sample_Barcode" = "sanger_dna_id")
  ) |>
  select(Hugo_Symbol, study_id, Main_consequence_VEP) |>
  arrange(study_id)

# Transform data
mat_df <- maf |>
  group_by(Hugo_Symbol, study_id) |>
  summarise(Main_consequence_VEP = paste(unique(Main_consequence_VEP), collapse = ";"), .groups = "drop") |>
  pivot_wider(names_from = study_id, values_from = Main_consequence_VEP, values_fill = "")

mat <- as.matrix(mat_df[, -1])
rownames(mat) <- mat_df$Hugo_Symbol

mat <- mat[, order(colnames(mat))]

# Define plot parameters
col <- c(
  "Missense" = "forestgreen",
  "Nonsense" = "firebrick",
  "Splice site" = "orange2",
  "Frameshift" = "mediumpurple3",
  "Inframe indel" = "sienna4",
  "Start codon lost" = "dodgerblue4",
  "Stop codon lost" = "pink"
)

alter_fun <- list(
  `Missense` = function(x, y, w, h) grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = col["Missense"], col = NA)),
  `Nonsense` = function(x, y, w, h) grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = col["Nonsense"], col = NA)),
  `Splice site` = function(x, y, w, h) grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = col["Splice site"], col = NA)),
  `Frameshift` = function(x, y, w, h) grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = col["Frameshift"], col = NA)),
  `Inframe indel` = function(x, y, w, h) grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = col["Inframe indel"], col = NA)),
  `Start codon lost` = function(x, y, w, h) grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = col["Start codon lost"], col = NA)),
  `Stop codon lost` = function(x, y, w, h) grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = col["Stop codon lost"], col = NA))
)

# Add TMB
tmb <- tmb_raw |>
  left_join(metadata |> select(sanger_dna_id, study_id),
    by = c("Sample" = "sanger_dna_id")
  ) |>
  select(study_id, TMB) |>
  filter(study_id %in% colnames(mat)) |>
  arrange(study_id)

top_anno <- HeatmapAnnotation(
  "TMB\n(per Mb)" = anno_barplot(tmb$TMB, baseline = 0),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10),
  border = T,
  height = unit(1, "cm")
)

# Add metadata
metadata <- metadata |>
  filter(study_id %in% colnames(mat)) |>
  arrange(study_id)

bottom_anno <- HeatmapAnnotation(
  "Grade" = as.vector(metadata$grade_of_dysplasia),
  "IBD" = as.vector(metadata$ibd_diagnosis),
  show_legend = c(
    "Grade" = T,
    "IBD" = T
  ),
  col = list(
    "Grade" = c("Low grade" = "#ffffcc", "High grade" = "#66c2a5", "Adenocarcinoma" = "#225ea8"),
    "IBD" = c("Crohn's" = "#ffffcc", "IBDU" = "#66c2a5", "UC" = "#225ea8")
  ),
  gap = unit(c(0, 0, 1), "mm"),
  border = T, na_col = "#e6e6e6",
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10),
  simple_anno_size = unit(0.4, "cm")
)

# Create plot
p <- oncoPrint(
  mat,
  alter_fun = alter_fun,
  show_pct = FALSE,
  row_names_side = "left",
  pct_digits = 1,
  show_column_names = TRUE,
  column_split = factor(metadata$precursor_or_follow_up),
  column_gap = unit(1, "mm"),
  column_names_gp = gpar(fontsize = 5),
  top_annotation = top_anno,
  bottom_annotation = bottom_anno,
  width = unit(10, "cm"),
  height = unit(8, "cm"),
  # border = T,
  row_names_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 12)
)

pdf(
  file = paste0("plots/oncoplots/progressors_all.pdf"),
  width = 6, height = 5
)
draw(p,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = TRUE,
)
dev.off()
