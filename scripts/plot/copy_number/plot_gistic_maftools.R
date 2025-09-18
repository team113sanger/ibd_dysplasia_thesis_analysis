library(maftools)

pre.gistic = readGistic(
    gisticAllLesionsFile = "data/copy_number/gistic/prog_pre/MIN_0/all_lesions.conf_95.txt",
    gisticAmpGenesFile = "data/copy_number/gistic/prog_pre/MIN_0/amp_genes.conf_95.txt",
    gisticDelGenesFile = "data/copy_number/gistic/prog_pre/MIN_0/del_genes.conf_95.txt",
    gisticScoresFile = "data/copy_number/gistic/prog_pre/MIN_0/scores.gistic")

pdf("plots/copy_number/gistic/pre_gistic_plot.pdf", width = 7, height = 4)
gisticChromPlot(gistic = pre.gistic, markBands = "all", ref.build = "hg38")
dev.off()

pdf("plots/copy_number/gistic/pre_gistic_plot_2.pdf")
gisticBubblePlot(gistic = pre.gistic)
dev.off()

fol.gistic = readGistic(
    gisticAllLesionsFile = "data/copy_number/gistic/prog_fol/MIN_0/all_lesions.conf_95.txt",
    gisticAmpGenesFile = "data/copy_number/gistic/prog_fol/MIN_0/amp_genes.conf_95.txt",
    gisticDelGenesFile = "data/copy_number/gistic/prog_fol/MIN_0/del_genes.conf_95.txt",
    gisticScoresFile = "data/copy_number/gistic/prog_fol/MIN_0/scores.gistic")

pdf("plots/copy_number/gistic/fol_gistic_plot.pdf", width = 7, height = 4)
gisticChromPlot(gistic = fol.gistic, markBands = "all", ref.build = "hg38")
dev.off()

pdf("plots/copy_number/gistic/fol_gistic_plot_2.pdf")
gisticBubblePlot(gistic = fol.gistic)
dev.off()

pdf("plots/copy_number/gistic/gistic_plot.pdf", width = 6, height = 5)
par(mfrow=c(2,1))
gisticChromPlot(gistic = pre.gistic, markBands = "all", ref.build = "hg38")
gisticChromPlot(gistic = fol.gistic, markBands = "all", ref.build = "hg38")
dev.off()