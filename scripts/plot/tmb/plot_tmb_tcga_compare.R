library(maftools)
library(readr)
library(dplyr)

maf <- read_tsv("data/variants/7100_3235-filtered_mutations_matched_indepTum_keep.maf")
metadata <- read_tsv("metadata/final_metadata_qc_pass.tsv")

#### By Group ###
progressors <- metadata |>
  filter(group == "Progressor") |>
  pull(sanger_dna_id)

progressor_maf <- maf |>
  filter(Tumor_Sample_Barcode %in% progressors)

non_progressors <- metadata |>
  filter(group == "Non-progressor") |>
  pull(sanger_dna_id)

non_progressor_maf <- maf |>
  filter(Tumor_Sample_Barcode %in% non_progressors)

maf_1 <- read.maf(progressor_maf, verbose = FALSE)
maf_2 <- read.maf(non_progressor_maf, verbose = FALSE)

pdf("plots/TMB/tmb_by group.pdf", width = 4, height = 6)

tcgaCompare(
  maf = c(maf_1, maf_2), cohortName = c("PRG", "N-PRG"), logscale = TRUE,
  capture_size = 50, tcga_cohorts = "COAD"
)

dev.off()

#### By Grade ####
low_grade <- metadata |>
  filter(grade_of_dysplasia == "Low grade") |>
  pull(sanger_dna_id)

low_grade_maf <- maf |>
  filter(Tumor_Sample_Barcode %in% low_grade)

high_grade <- metadata |>
  filter(grade_of_dysplasia == "High grade") |>
  pull(sanger_dna_id)

high_grade_maf <- maf |>
  filter(Tumor_Sample_Barcode %in% high_grade)

adenocarcinoma <- metadata |>
  filter(grade_of_dysplasia == "Adenocarcinoma") |>
  pull(sanger_dna_id)

adenocarcinoma_maf <- maf |>
  filter(Tumor_Sample_Barcode %in% adenocarcinoma)

maf_1 <- read.maf(low_grade_maf, verbose = FALSE)
maf_2 <- read.maf(high_grade_maf, verbose = FALSE)
maf_3 <- read.maf(adenocarcinoma_maf, verbose = FALSE)

pdf("plots/TMB/tmb_by_grade.pdf", width = 4, height = 6)

tcgaCompare(
  maf = c(maf_1, maf_2, maf_3), cohortName = c("LGD", "HGD", "IBD-CRC"), logscale = TRUE,
  capture_size = 50, tcga_cohorts = "COAD"
)

dev.off()
