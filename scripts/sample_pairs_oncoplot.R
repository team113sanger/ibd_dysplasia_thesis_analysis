source("/lustre/scratch125/casm/team113da/projects/dermatlas_analysis_methods/somatic-variant-plots/src/oncoplot_utils.R")

library(dplyr)
library(tidyr)
library(readr)
library(ComplexHeatmap)
library(grid)

#### Load in data ####
oncoKB <- read_tsv("/lustre/scratch124/casm/team113/secure-lustre/resources/dermatlas/oncoKB/cancerGeneList.tsv") |>
  pull(`Hugo Symbol`)

tmb_raw <- read_tsv("data/mutations_per_Mb.tsv", col_names = c("Sample", "TMB"))

metadata <- read_tsv("metadata/sample_pairs.tsv")
metadata$precursor_or_follow_up <- factor(
  metadata$precursor_or_follow_up,
  levels = c("Precursor", "Follow up")
)
pair_ids <- metadata[["sanger_dna_id"]]
progressors <- metadata |>
    filter(group == "Progressor") |>
    pull(sanger_dna_id)
non_progressors <- metadata |>
    filter(group == "Non-progressor") |>
    pull(sanger_dna_id)

#### Progressors ####
maf <- read_tsv("data/7100_3235-filtered_mutations_all_indepTum_keepPA.maf") |>
  # filter(Hugo_Symbol %in% oncoKB) |>
  filter(Tumor_Sample_Barcode %in% pair_ids) |>
  filter(Tumor_Sample_Barcode %in% progressors) |>
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
  "Grade" = as.vector(metadata$grade_of_dysplasia),
  "IBD" = as.vector(metadata$ibd_diagnosis),
  "Site" = as.vector(metadata$site),
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
      "NOS" = "lightgrey"
    ),
    "IBD" = c("Crohn's" = "#ffffcc", "IBDU" = "#66c2a5", "UC" = "#225ea8"),
    "Site" = c(
      "Sigmoid" = "#d0e3e8",  
      "Transverse" = "#A4C8E1", 
      "Rectum" = "#7FB1CC",  
      "Ascending" = "#5B9BB1", 
      "Distal Ascending" = "#009688",  
      "Proximal Ascending" = "#66C2A5",  
      "Splenic Flexure" = "#A0DAB3", 
      "Hepatic Flexure" = "#F1E6C8",  
      "Caecum" = "#FFEA78", 
      "Descending" = "#FFD700",
      "Rectosigmoid" = "#F4A582"
    )
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
  column_split = factor(metadata$patient_id),
  column_gap = unit(1, "mm"),
  column_names_gp = gpar(fontsize = 5),
  top_annotation = top_anno,
  bottom_annotation = bottom_anno,
  width = unit(12, "cm"),
  height = unit(10, "cm"),
  remove_empty_columns = FALSE,
  # border = T,
  row_names_gp = gpar(fontsize = 8),
  column_title = c("Sample Pairs"),
  column_title_gp = gpar(fontsize = 12)
)

pdf(
  file = paste0("plots/sample_pairs_oncoplot.pdf"),
  width = 8, height = 6
)
draw(p,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = TRUE,
)
dev.off()

#### Non-progressors ####
# maf <- read_tsv("data/7100_3235-filtered_mutations_all_indepTum_keepPA.maf") |>
#   # filter(Hugo_Symbol %in% oncoKB) |>
#   filter(Tumor_Sample_Barcode %in% pair_ids) |>
#   filter(Tumor_Sample_Barcode %in% non_progressors) |>
#   group_by(Hugo_Symbol) |>
#   filter(n_distinct(Tumor_Sample_Barcode) >= 4) |>
#   ungroup() |>
#   select(Hugo_Symbol, Tumor_Sample_Barcode, Main_consequence_VEP) |>
#   shorten_consequence() |>
#   left_join(metadata |> select(sanger_dna_id, study_id), 
#             by = c("Tumor_Sample_Barcode" = "sanger_dna_id")) |>
#   select(Hugo_Symbol, study_id, Main_consequence_VEP) |>
#   arrange(study_id)

# # Transform data
# mat_df <- maf |>
#   group_by(Hugo_Symbol, study_id) |>
#   summarise(Main_consequence_VEP = paste(unique(Main_consequence_VEP), collapse = ";"), .groups = "drop") |>
#   pivot_wider(names_from = study_id, values_from = Main_consequence_VEP, values_fill = "")

# mat <- as.matrix(mat_df[, -1])
# rownames(mat) <- mat_df$Hugo_Symbol

# mat <- mat[, order(colnames(mat))]

# # Define plot parameters
# col <- c(
#   "Missense" = "forestgreen",
#   "Nonsense" = "firebrick",
#   "Splice site" = "orange2",
#   "Frameshift" = "mediumpurple3",
#   "Inframe indel" = "sienna4",
#   "Start codon lost" = "dodgerblue4",
#   "Stop codon lost" = "pink"
# )

# alter_fun <- list(
#   `Missense` = function(x, y, w, h) grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = col["Missense"], col = NA)),
#   `Nonsense` = function(x, y, w, h) grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = col["Nonsense"], col = NA)),
#   `Splice site` = function(x, y, w, h) grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = col["Splice site"], col = NA)),
#   `Frameshift` = function(x, y, w, h) grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = col["Frameshift"], col = NA)),
#   `Inframe indel` = function(x, y, w, h) grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = col["Inframe indel"], col = NA)),
#   `Start codon lost` = function(x, y, w, h) grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = col["Start codon lost"], col = NA)),
#   `Stop codon lost` = function(x, y, w, h) grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = col["Stop codon lost"], col = NA))
# )

# # Add TMB
# tmb <- tmb_raw |>
#   left_join(metadata |> select(sanger_dna_id, study_id), 
#             by = c("Sample" = "sanger_dna_id")) |>
#   select(study_id, TMB) |>
#   filter(study_id %in% colnames(mat)) |>
#   arrange(study_id)

# top_anno <- HeatmapAnnotation(
#   "TMB\n(per Mb)" = anno_barplot(tmb$TMB, baseline = 0),
#   annotation_name_side = "left",
#   annotation_name_gp = gpar(fontsize = 10),
#   border = T,
#   height = unit(1, "cm")
# )

# # Add metadata
# metadata <- metadata |>
#   filter(study_id %in% colnames(mat)) |>
#   arrange(study_id)

# bottom_anno <- HeatmapAnnotation(
#   "Grade" = as.vector(metadata$grade_of_dysplasia),
#   "IBD" = as.vector(metadata$ibd_diagnosis),
#   "Site" = as.vector(metadata$site),
#   show_legend = c(
#     "Grade" = T,
#     "IBD" = T,
#     "Site" = TRUE
#   ),
#   col = list(
#     "Grade" = c(
#       "Low grade" = "#ffffcc",
#       "High grade" = "#66c2a5",
#       "Adenocarcinoma" = "#225ea8",
#       "NOS" = "lightgrey"
#     ),
#     "IBD" = c("Crohn's" = "#ffffcc", "IBDU" = "#66c2a5", "UC" = "#225ea8"),
#     "Site" = c(
#       "Sigmoid" = "#d0e3e8",  
#       "Transverse" = "#A4C8E1", 
#       "Rectum" = "#7FB1CC",  
#       "Ascending" = "#5B9BB1", 
#       "Distal Ascending" = "#009688",  
#       "Proximal Ascending" = "#66C2A5",  
#       "Splenic Flexure" = "#A0DAB3", 
#       "Hepatic Flexure" = "#F1E6C8",  
#       "Caecum" = "#FFEA78", 
#       "Descending" = "#FFD700",
#       "Rectosigmoid" = "#F4A582"
#     )
#   ),
#   gap = unit(c(0, 0, 1), "mm"),
#   border = T, na_col = "#e6e6e6",
#   annotation_name_side = "left",
#   annotation_name_gp = gpar(fontsize = 10),
#   simple_anno_size = unit(0.4, "cm")
# )

# # Create plot
# p <- oncoPrint(
#   mat,
#   alter_fun = alter_fun,
#   show_pct = FALSE,
#   row_names_side = "left",
#   pct_digits = 1,
#   show_column_names = TRUE,
#   column_split = factor(metadata$patient_id),
#   column_gap = unit(1, "mm"),
#   column_names_gp = gpar(fontsize = 5),
#   top_annotation = top_anno,
#   bottom_annotation = bottom_anno,
#   width = unit(12, "cm"),
#   height = unit(10, "cm"),
#   remove_empty_columns = FALSE,
#   # border = T,
#   row_names_gp = gpar(fontsize = 8),
#   column_title_gp = gpar(fontsize = 9)
# )

# pdf(
#   file = paste0("plots/np_sample_pairs_oncoplot.pdf"),
#   width = 8, height = 6
# )
# draw(p,
#   heatmap_legend_side = "right",
#   annotation_legend_side = "right",
#   merge_legend = TRUE,
# )
# dev.off()