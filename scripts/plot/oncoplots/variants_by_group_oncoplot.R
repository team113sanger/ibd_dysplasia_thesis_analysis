source("scripts/plot/oncoplots/oncoplot_utils.R")

library(dplyr)
library(tidyr)
library(readr)
library(ComplexHeatmap)
library(grid)
library(maftools)

# Load in data
oncoKB <- read_tsv("metadata/rescources/cancerGeneList.tsv") |>
  pull(`Hugo Symbol`)

tmb_raw <- read_tsv("data/variants/mutations_per_Mb.tsv", col_names = c("Sample", "TMB"))
samples <- read_lines("metadata/sample_lists/all_one_ppat.list")

metadata <- read_tsv("metadata/final_metadata_qc_pass.tsv") |>
  mutate(group = paste(group, " ", precursor_or_follow_up, sep = "")) |>
  filter(sanger_dna_id %in% samples)

# Add levels
metadata[["group"]] <- factor(
  metadata[["group"]],
  levels = c(
    "Non-progressor Precursor",  "Non-progressor Follow up",
    "Progressor Precursor", "Progressor Follow up"
  )
)

# metadata[["grade_of_dysplasia"]] <- factor(
#   metadata[["grade_of_dysplasia"]],
#   levels = c("NOS", "Low grade", "High grade", "Adenocarcinoma")
# )

maf <- read_tsv("data/variants/7100_3235-filtered_mutations_matched_allTum_keepPA.maf") |>
  filter(Tumor_Sample_Barcode %in% samples) |>
  filter(!Hugo_Symbol %in% maftools:::flags(top = 20)) |>
  filter(Hugo_Symbol %in% oncoKB) |>
  group_by(Hugo_Symbol) |>
  filter(n_distinct(Tumor_Sample_Barcode) >= 4) |>
  ungroup() |>
  select(Hugo_Symbol, Tumor_Sample_Barcode, Main_consequence_VEP) |>
  shorten_consequence() |>
  left_join(metadata |> select(sanger_dna_id),
    by = c("Tumor_Sample_Barcode" = "sanger_dna_id")
  ) |>
  arrange(Tumor_Sample_Barcode)

# Transform data
mat_df <- maf |>
  group_by(Hugo_Symbol, Tumor_Sample_Barcode) |>
  summarise(Main_consequence_VEP = paste(unique(Main_consequence_VEP), collapse = ";"), .groups = "drop") |>
  pivot_wider(names_from = Tumor_Sample_Barcode, values_from = Main_consequence_VEP, values_fill = "")

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
  filter(Sample %in% colnames(mat)) |>
  arrange(Sample)

top_anno <- HeatmapAnnotation(
  "TMB\n(per Mb)" = anno_barplot(tmb$TMB, baseline = 0),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10),
  border = T,
  height = unit(1, "cm")
)

# Add metadata
metadata <- metadata |>
  filter(sanger_dna_id %in% colnames(mat)) |>
  arrange(sanger_dna_id)

bottom_anno <- HeatmapAnnotation(
  "Grade" = as.vector(metadata$grade_of_dysplasia),
  "IBD" = as.vector(metadata$ibd_diagnosis),
  "Site" = as.vector(metadata$general_site),
  show_legend = c(
    "Grade" = T,
    "IBD" = T,
    "Site" = TRUE
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
        "Crohn's" = "lightpink",
        "IBDU" = "darkorchid",
        "UC" = "darkseagreen"),
    "Site" = c(
      "Left colon"  = "#d0e3e8",
      "Transverse colon" = "#A4C8E1",
      "Right colon" = "#7FB1CC",
      "Rectum" = "#5B9BB1"
    )
  ),
  gap = unit(c(0, 0, 1), "mm"),
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
  width = unit(17, "cm"),
  height = unit(7, "cm"),
  remove_empty_columns = FALSE,
  # border = T,
  row_names_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 9)
)

pdf(
  file = paste0("plots/oncoplots/variants_by_group.pdf"),
  width = 8, height = 6
)
draw(p,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legend = TRUE,
)
dev.off()
