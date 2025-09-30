library(readr)
library(dplyr)
library(stringr)
library(maftools)

### Prog Pre ###
maf <- read.maf("data/variants/split_mafs/prog_pre_keepPA.maf")

pdf("plots/variants/prog_pre_mut_vafs.pdf",
    width = 5, height =4)
plotVaf(maf = maf, vafCol = 'VAF_tum')
dev.off()

pdf("plots/variants/tp53_lollipop.pdf",
    width = 5, height =4)
lollipopPlot(
  maf = maf,
  gene = 'TP53',
  defaultYaxis = TRUE
)
dev.off()

### Prog Fol ###
maf <- read.maf("data/variants/split_mafs/prog_fol_keepPA.maf")

pdf("plots/variants/prog_fol_mut_vafs.pdf",
    width = 5, height =4)
plotVaf(maf = maf, vafCol = 'VAF_tum')
dev.off()

pdf("plots/variants/tp53_fol_lollipop.pdf",
    width = 5, height =4)
lollipopPlot(
  maf = maf,
  gene = 'TP53',
  defaultYaxis = TRUE
)
dev.off()

### N-Prog Pre ###
maf <- read.maf("data/variants/split_mafs/n_prog_pre_keepPA.maf")

pdf("plots/variants/n_prog_pre_mut_vafs.pdf",
    width = 5, height =4)
plotVaf(maf = maf, vafCol = 'VAF_tum')
dev.off()

### N-Prog Fol ###
maf <- read.maf("data/variants/split_mafs/n_prog_fol_keepPA.maf")

pdf("plots/variants/n_prog_fol_mut_vafs.pdf",
    width = 5, height =4)
plotVaf(maf = maf, vafCol = 'VAF_tum')
dev.off()

#### TP53 lollipop plot ####

# Load data
maf <- read_tsv("data/variants/7100_3235-filtered_mutations_matched_allTum_keepPA.maf")
tp53_muts <- read_tsv("results/p53_mutations/p53_status.tsv")
meta <- read_tsv("metadata/final_metadata_qc_pass.tsv")

# Filter for TP53 mutations and join metadata
tp53_maf <- maf |>
  filter(Hugo_Symbol == "TP53") |>
  left_join(meta, by = c("Tumor_Sample_Barcode" = "sanger_dna_id")) |>
  distinct(patient_id, HGVSp_Short, .keep_all = TRUE)

# --- (1) Non-progressor TP53 mutations ---
tp53_maf_nonprog <- tp53_maf |> filter(group == "Non-progressor")

tp53_maf_nonprog_maf <- read.maf(tp53_maf_nonprog, verbose = FALSE)

pdf("plots/variants/tp53_lollipop_nonprogressor.pdf", width = 5, height = 4)
lollipopPlot(
  maf = tp53_maf_nonprog_maf,
  gene = "TP53",
  defaultYaxis = TRUE
)
dev.off()

# --- (2) Progressor TP53 mutations ---
tp53_maf_prog <- tp53_maf |> filter(group == "Progressor")

tp53_maf_prog_maf <- read.maf(tp53_maf_prog, verbose = FALSE)

pdf("plots/variants/tp53_lollipop_progressor.pdf", width = 5, height = 4)
lollipopPlot(
  maf = tp53_maf_prog_maf,
  gene = "TP53",
  defaultYaxis = TRUE
)
dev.off()