library(readr)
library(dplyr)
library(GenomicRanges)
library(IRanges)
library(trackViewer)

# Load data ---------------------------------------------------------------

tp53_muts <- read_tsv("results/p53_mutations/p53_muts.txt")
mut_status <- read_tsv("results/p53_mutations/p53_status.tsv")

samples <- read_lines("metadata/sample_lists/all_one_ppat.list")

# Merge tables on sample + protein change --------------------------------

muts <- tp53_muts %>%
  inner_join(mut_status, 
             by = c("Sample_ID" = "Tumor_Sample_Barcode",
                    "HGVSp_Short" = "HGVSp_Short")) |>
 filter(!is.na(group),
        Sample_ID %in% samples,
        precursor_or_follow_up == "Precursor") 

# Count how many times each mutation appears ------------------------------

mut_summary <- muts %>%
  group_by(Chromosome, Start_Position, HGVSp_Short, Main_consequence_VEP, group) %>%
  summarise(count = n(), .groups = "drop")

# Define Colours ------------------------------

mutation_colors <- c(
  missense_variant = "#9BCD9B",
  stop_gained = "#e27396",
  frameshift_variant = "#7EC0EE",
  splice_donor_variant = "#7077FF",
  splice_acceptor_variant = "#7077FF"
)

# Build GRanges for mutations --------------------------------------------

plot.snp <- GRanges(
  seqnames = mut_summary$Chromosome,
  ranges   = IRanges(
    start = mut_summary$Start_Position,
    width = 1,
    names = mut_summary$HGVSp_Short
  ),
  score    = mut_summary$count
)

plot.snp$SNPsideID <- ifelse(
  mut_summary$group == "Progressor", 
  "top",
  "bottom"
)

plot.snp$color <- mutation_colors[ mut_summary$Main_consequence_VEP ]
plot.snp$label.parameter.rot <- 45
plot.snp$label.parameter.gp <- gpar(cex = 0.8)


# TP53 gene structure (hg38) ----------------------------------------------
# Exon coordinates are hard-coded—no external database needed.

gene.features <- GRanges(
  seqnames = "chr17",
  ranges = IRanges(
    start = c(7676521, 7676382, 7675994, 7675053, 7674859, 7674181, 7673701, 7673535, 7670609, 7668421),
    end   = c(7676622, 7676403, 7676272, 7675236, 7674971, 7674290, 7673837, 7673608, 7670715, 7669690)
  ),
  names = paste("Exon", 1:10),
  fill  = "#CDC0B0",
  height = 0.03
)

yaxis <- c(0, 1, 2, 3, 4, 5)

# Plot --------------------------------------------------------------------

pdf("plots/variants/tp53_lollipop.pdf")
lolliplot(
  plot.snp,
  gene.features,
  lollipop_style_switch_limit = 10, 
  yaxis = yaxis,
  lolliplot_style = "caterpillar",
  ylab = ""
)

grid.text(label = "Mutation Frequency", 
      x = unit(0.02, "npc"),
      y = unit(0.275, "npc"), 
      rot = 90,
      gp = gpar(fontsize = 12, col = "black"))


dev.off()