library(readr)
library(dplyr)

maf <- read_tsv("data/variants/7100_3235-filtered_mutations_matched_allTum_keepPA.maf")
metadata <- read_tsv("metadata/final_metadata_qc_pass.tsv")

p53_muts <- maf |>
  filter(Hugo_Symbol == "TP53") |>
  select(Tumor_Sample_Barcode, Main_consequence_VEP, Variant_Type, HGVSp_Short, VAF_tum) |>
  left_join(metadata |> select(sanger_dna_id, group, precursor_or_follow_up, grade_of_dysplasia), by = c("Tumor_Sample_Barcode" = "sanger_dna_id"))

write_tsv(p53_muts, "results/p53_mutations/p53_status.tsv")


maf_new <- maf |>
  filter(Hugo_Symbol == "TP53") |>
  select(
    Sample_ID = Tumor_Sample_Barcode, Chromosome, Start_Position, End_Position,
    Reference_Allele, Variant_Allele = Tumor_Seq_Allele2, Hugo_Symbol
  ) |>
  arrange(Start_Position)

write_tsv(maf_new, "results/p53_mutations/p53_muts.txt")
