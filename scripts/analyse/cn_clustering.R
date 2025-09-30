
library(dplyr)
library(tibble)

# Read in cnv calls
cnv_mat <- read_tsv("results/copy_number/cnv_calls.txt")

# Convert to matrix and drop bins column
num_mat <- as.matrix(cnv_mat[ , -1])

# Compute pairwise distances between samples
dist_mat <- dist(t(num_mat), method = "euclidean")

# Hierarchical clustering (build dendrogram)
hc_mat <- hclust(dist_mat, method = "ward.D")

# Cut into different cluster numbers (2, 3 and 4 clusters)
CN_clusters <- as.data.frame(cutree(hc_mat, k = c(2, 3, 4))) |>
  tibble::rownames_to_column("sample")

# Save results
write_tsv(CN_clusters, "results/copy_number/precursors/cn_clusters.tsv")

# Plot dendrogram
pdf("results/copy_number/precursors/clusters_dendrogram.pdf")
plot(hc_mat, labels = colnames(num_mat), main = "CNV clustering")
dev.off()