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
  filter(sanger_dna_id %in% read_lines("metadata/sample_lists/all_one_ppat.list")) |>
  filter(precursor_or_follow_up == "Precursor") |>
  filter(!grade_of_dysplasia == "High grade")

metadata$group <- factor(
  metadata$group,
  levels = c("Non-progressor", "Progressor")
)

samples <- metadata |>
  pull(sanger_dna_id)

tmb_raw <- read_tsv("data/variants/mutations_per_Mb.tsv", col_names = c("Sample", "TMB")) |>
  filter(Sample %in% samples)

maf <- read_tsv("data/variants/7100_3235-filtered_mutations_matched_allTum_keepPA.maf") |>
  filter(Hugo_Symbol %in% oncoKB) |>
  filter(Tumor_Sample_Barcode %in% samples) |>
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
  "Missense" = "#20854E",
  "Nonsense" = "#9f2727ff",
  "Splice site" = "#E6AB02",
  "Frameshift" = "mediumpurple3",
  "Inframe indel" = "#7d6a56ff",
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
  "TMB" = anno_barplot(tmb$TMB, baseline = 0),
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
  "Group" = as.vector(metadata$group),
  "IBD" = as.vector(metadata$ibd_diagnosis),
  "Site" = as.vector(metadata$site_general),
  show_legend = c(
    "Grade" = T,
    "IBD" = T,
    "Site" = T,
    "Group" =T
  ),
  col = list(
    "Grade" = c(
      "Low grade" = "#ffffcc",
      "High grade" = "#66c2a5",
      "Adenocarcinoma" = "#225ea8",
      "NOS" = "lightgrey",
      "Moderate" = "#F1E6C8",
      "Not specified" = 'grey'
    ),
    "IBD" = c(
        "Crohn's" = "#EE4C97",
        "IBDU" = "#FFDC91",
        "UC" = "#5B8FA8"),
    "Site" = c(
      "Left colon"  = "#7375B5",
      "Transverse colon" = "#919C4C",
      "Right colon" = "#F5C04A",
      "Rectum" = "#E68C7C"
    ),
    "Group" = c(
      "Progressor" = "#dd882dff",
      "Non-progressor" = "#388643ff"
    )
  ),
  gap = unit(c(1, 1, 1, 1), "mm"),
  border = T, na_col = "#e6e6e6",
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10),
  simple_anno_size = unit(0.3, "cm")
)

# Create plot
p <- oncoPrint(
  mat,
  alter_fun = alter_fun,
  show_pct = FALSE,
  row_names_side = "left",
  pct_digits = 1,
  show_column_names = TRUE,
  column_split = factor(metadata$group),
  column_gap = unit(1, "mm"),
  column_names_gp = gpar(fontsize = 5),
  top_annotation = top_anno,
  bottom_annotation = bottom_anno,
  width = unit(12, "cm"),
  height = unit(7, "cm"),
  # border = T,
  row_names_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 12)
)

pdf(
  file = paste0("plots/oncoplots/precursors_oncoplot.pdf"),
  width = 7, height = 6
)
draw(p,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legend = TRUE,
)
dev.off()
