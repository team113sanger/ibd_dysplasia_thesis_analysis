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
