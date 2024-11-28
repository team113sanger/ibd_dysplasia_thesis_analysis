library(dplyr)
library(readr)

maf <- read_tsv("/lustre/scratch125/casm/team113da/projects/IBD_Associated_Dysplasia/ibd_dysplasia_analysis/data/7100_3235-filtered_mutations_all_indepTum_keepPA.maf")

problems(maf)

genes <- maf |> 
    filter(Tumor_Sample_Barcode == "PD62027c") |>
    pull(Hugo_Symbol)

write_lines(genes, "data/test_genes.txt")