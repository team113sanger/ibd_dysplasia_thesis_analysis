library(readr)
library(dplyr)

# Read in data
samples <- read_lines("metadata/sample_lists/all_one_ppat.list")
meta <- read_tsv("metadata/final_metadata_qc_pass.tsv") |>
    select(sanger_dna_id, group, precursor_or_follow_up, grade_of_dysplasia) |>
    filter(sanger_dna_id %in% samples) |>
    mutate(group_status = paste(group, precursor_or_follow_up, sep = "_"))


maf <- read_tsv("data/variants/7100_3235-filtered_mutations_matched_allTum_keep.maf")
mafPA <- read_tsv("data/variants/7100_3235-filtered_mutations_matched_allTum_keepPA.maf") |>
    left_join(meta, by = c(Tumor_Sample_Barcode = "sanger_dna_id"))

# Number of samples
length(samples)
table(meta$grade_of_dysplasia)
table(meta$group_status)

# Total counts
maf_split <- maf |>
    filter(Tumor_Sample_Barcode %in% samples)

mafPA_split <- mafPA |>
    filter(Tumor_Sample_Barcode %in% samples)

muts_count <- nrow(maf_split)
muts_count

PA_muts_count <- nrow(mafPA_split)
PA_muts_count

# Gene counts by group
all_gene_counts <- mafPA_split %>%
  filter(Hugo_Symbol %in% c("APC", "KRAS", "TP53", "RNF43")) %>%
  count(group_status, Hugo_Symbol) %>%
  arrange(Hugo_Symbol, group_status)

all_gene_counts

distinct_gene_counts <- mafPA_split %>%
  filter(Hugo_Symbol %in% c("APC", "KRAS", "TP53", "RNF43")) %>%
  mutate(group_status = paste(group, precursor_or_follow_up, sep = "_")) %>%
  distinct(Tumor_Sample_Barcode, Hugo_Symbol, group_status) %>% 
  count(group_status, Hugo_Symbol) %>%
  arrange(Hugo_Symbol, group_status)

distinct_gene_counts

# Gene counts by grade
all_gene_counts <- mafPA_split %>%
  filter(Hugo_Symbol %in% c("APC", "KRAS", "TP53")) %>%
  count(grade_of_dysplasia, Hugo_Symbol) %>%
  arrange(Hugo_Symbol, grade_of_dysplasia)

all_gene_counts

distinct_gene_counts <- mafPA_split %>%
  filter(Hugo_Symbol %in% c("APC", "KRAS", "TP53")) %>%
  distinct(Tumor_Sample_Barcode, Hugo_Symbol, grade_of_dysplasia) %>% 
  count(grade_of_dysplasia, Hugo_Symbol) %>%
  arrange(Hugo_Symbol, grade_of_dysplasia)

distinct_gene_counts