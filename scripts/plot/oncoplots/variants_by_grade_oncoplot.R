source("scripts/plotting/oncoplots/oncoplot_utils.R")

library(dplyr)
library(tidyr)
library(readr)
library(ComplexHeatmap)
library(grid)

# Load in data
oncoKB <- read_tsv("metadata/rescources/cancerGeneList.tsv") |>
  pull(`Hugo Symbol`)

tmb_raw <- read_tsv("data/variants/mutations_per_Mb.tsv", col_names = c("Sample", "TMB"))

metadata <- read_tsv("metadata/final_metadata_qc_pass.tsv")

metadata$grade_of_dysplasia <- factor(
  metadata$grade_of_dysplasia,
  levels = c("NOS", "Low grade", "High grade", "Adenocarcinoma")
)

samples_to_plot <- metadata %>%
  filter(
    (grade_of_dysplasia == "Low grade" & precursor_or_follow_up == "Precursor") |
    grade_of_dysplasia == "High grade" |
    grade_of_dysplasia == "Adenocarcinoma"
  ) %>%
  pull(sanger_dna_id)

#plot_genes <- c("TP53", "KRAS", "APC", "RNF43", "RBM10", "MSH3", "POLD1", "APOBEC3A", "PIK3CA")

maf <- read_tsv("data/variants/7100_3235-filtered_mutations_all_indepTum_keepPA.maf") |>
  filter(Hugo_Symbol %in% oncoKB) |>
  filter(Tumor_Sample_Barcode %in% samples_to_plot) |>
  group_by(Hugo_Symbol) |>
  filter(n_distinct(Tumor_Sample_Barcode) >= 4) |>
  ungroup() |>
  select(Hugo_Symbol, Tumor_Sample_Barcode, Main_consequence_VEP) |>
  shorten_consequence() |>
  left_join(metadata |> select(sanger_dna_id, study_id), 
            by = c("Tumor_Sample_Barcode" = "sanger_dna_id")) |>
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
            by = c("Sample" = "sanger_dna_id")) |>
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
  "Status" = as.vector(metadata$precursor_or_follow_up),
  "Group" = as.vector(metadata$group),
  "IBD" = as.vector(metadata$ibd_diagnosis),
  "Site" = as.vector(metadata$site),
  show_legend = c(
    "Status" = TRUE,
    "Group" = TRUE,
    "IBD" = TRUE,
    "Site" = TRUE
  ),
  col = list(
    "Status" = c("Precursor" = "#ffffcc", "Follow up" = "#225ea8"),
    "Group" = c("Progressor" = "#225ea8",
                 "Non-progressor" = "#ffffcc"), 
    "IBD" = c("Crohn's" = "#ffffcc", 
               "IBDU" = "#66c2a5", 
               "UC" = "#225ea8"),
    "Site" = c(
      "Sigmoid" = "#d0e3e8",  # Light Blue
      "Transverse" = "#A4C8E1",  # Medium Blue
      "Rectum" = "#7FB1CC",  # Blue
      "Ascending" = "#5B9BB1",  # Darker Blue
      "Distal Ascending" = "#009688",  # Teal
      "Proximal Ascending" = "#66C2A5",  # Light Green
      "Splenic Flexure" = "#A0DAB3",  # Mint Green
      "Hepatic Flexure" = "#F1E6C8",  # Light Yellow
      "Caecum" = "#FFEA78",  # Pale Yellow
      "Descending" = "#FFD700",  # Gold
      "Rectosigmoid" = "#F4A582"
    )
  ),
  gap = unit(c(0, 1, 1), "mm"),
  border = TRUE, na_col = "#e6e6e6",
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
  column_split = factor(metadata$grade_of_dysplasia),
  column_gap = unit(1, "mm"),
  column_names_gp = gpar(fontsize = 5),
  top_annotation = top_anno,
  bottom_annotation = bottom_anno,
  remove_empty_columns = FALSE,
  width = unit(17, "cm"),
  height = unit(8, "cm"),
  # border = T,
  row_names_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 9)
)

pdf(
  file = paste0("plots/oncoplots/variants_by_grade.pdf"),
  width = 9, height = 6
)
draw(p,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = TRUE,
)
dev.off()