library(maftools)
library(readr)
library(dplyr)

# Load metadata
metadata <- read_tsv("metadata/final_metadata_qc_pass.tsv") |>
  select(sanger_dna_id, group, precursor_or_follow_up)

samples <- read_tsv("metadata/sample_lists/all_one_ppat.list", col_names = "Sample") |>
  left_join(metadata, by = c(Sample = "sanger_dna_id"))

# Load maf
maf_all <- read.maf(maf = "data/variants/7100_3235-filtered_mutations_matched_allTum_keep.maf", verbose = FALSE)

# Sample group
prog_pre_ids <- samples |> filter(group == "Progressor", precursor_or_follow_up == "Precursor") |> pull(Sample)
prog_fol_ids <- samples |> filter(group == "Progressor", precursor_or_follow_up == "Follow up") |> pull(Sample)
nprog_pre_ids <- samples |> filter(group == "Non-progressor", precursor_or_follow_up == "Precursor") |> pull(Sample)
nprog_fol_ids <- samples |> filter(group == "Non-progressor", precursor_or_follow_up == "Follow up") |> pull(Sample)

# Split mafs
prog_pre_maf  <- subsetMaf(maf_all, tsb = prog_pre_ids)
prog_fol_maf  <- subsetMaf(maf_all, tsb = prog_fol_ids)
nprog_pre_maf <- subsetMaf(maf_all, tsb = nprog_pre_ids)
nprog_fol_maf <- subsetMaf(maf_all, tsb = nprog_fol_ids)

pdf("plots/tmb/tmb_by group.pdf", width = 4, height = 6)
tcgaCompare(
  maf = c(nprog_pre_maf, prog_pre_maf), cohortName = c("N-PRG", "PRG"), logscale = TRUE,
  capture_size = 50, tcga_cohorts = "COAD"
)
dev.off()

 #TODO: rest of script still needs updating
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
