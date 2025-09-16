library(maftools)

laml.gistic = readGistic(
    gisticAllLesionsFile = "data/copy_number/gistic/MIN_0/all_lesions.conf_95.txt",
    gisticAmpGenesFile = "data/copy_number/gistic/MIN_0/amp_genes.conf_95.txt",
    gisticDelGenesFile = "data/copy_number/gistic/MIN_0/del_genes.conf_95.txt",
    gisticScoresFile = "data/copy_number/gistic/MIN_0/scores.gistic")

pdf("plots/copy_number/gistic/gistic_plot.pdf")
gisticChromPlot(gistic = laml.gistic, markBands = "all", ref.build = "hg38")
dev.off()

pdf("plots/copy_number/gistic/gistic_plot_2.pdf")
gisticBubblePlot(gistic = laml.gistic)
dev.off()
