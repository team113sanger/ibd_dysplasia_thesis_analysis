library(readr)
library(dplyr)
library(ggplot2)
library(stringr)

sample_list <- read_lines("metadata/sample_lists/all_one_ppat.list")
meta <- read_tsv("metadata/final_metadata_qc_pass.tsv") |>
            mutate(patient_id = str_remove(sanger_dna_id, "[a-z]$")) |>
            filter(sanger_dna_id %in% sample_list) |>
            mutate(disease_duration_years = as.numeric(str_extract(disease_duration_to_first_dyaplastic_lesion, "\\d+")))


duration_stats <- meta |>
  group_by(group) |>
  summarise(
    n = n(),
    mean_duration = mean(disease_duration_years, na.rm = TRUE),
    median_duration = median(disease_duration_years, na.rm = TRUE),
    sd_duration = sd(disease_duration_years, na.rm = TRUE),
    min_duration = min(disease_duration_years, na.rm = TRUE),
    max_duration = max(disease_duration_years, na.rm = TRUE)
  )
duration_stats

ggplot(meta, aes(x = group, y = disease_duration_years, fill = group)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.6) +
  scale_fill_manual(values = c("Non-progressor" = "darkseagreen", 
                                "Progressor" = "darkorange")) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs(
    x = NULL,
    y = "Disease Duration (years)"
  ) +
  theme_classic() +
  theme(legend.position = "none")

ggsave("plots/metadata/disease_duration.png", width =3, height =3)