library(readr)
library(dplyr)
library(ggplot2)
library(stringr)

sample_list <- read_lines("metadata/sample_lists/all_one_ppat.list")
meta <- read_tsv("metadata/final_metadata_qc_pass.tsv") |>
            mutate(patient_id = str_remove(sanger_dna_id, "[a-z]$")) |> #keep to remove multiple lesions per patient
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

p1 <- ggplot(meta, aes(x = group, y = disease_duration_years, fill = group)) +
        geom_boxplot(outlier.shape = 21, alpha = 0.6, outlier.fill = "black") +
        scale_fill_manual(values = c("Non-progressor" = "darkseagreen", 
                                      "Progressor" = "chocolate")) +
        # geom_jitter(width = 0.2, alpha = 0.6) +
        labs(
          x = NULL,
          y = "Disease Duration (years)"
        ) +
        theme_classic() +
        theme(legend.position = "none")

ggsave("plots/metadata/disease_duration.png", p1, width =3, height =3)

prog <- meta |>
  filter(group =="Progressor")
hist(prog$disease_duration_years)

nprog <- meta |>
  filter(group =="Non-progressor")
hist(nprog$disease_duration_years)

p2 <- ggplot(meta, aes(x = disease_duration_years, color = group)) +
        stat_ecdf(size = 1) +
        theme_classic() +
        labs(colour = NULL, x = "Disease duration to dysplasia") +
        theme(legend.position = c(0.75, 0.3))
ggsave("plots/metadata/disease_duration_cdf.png", p2, width = 3.5, height =3.5)

p3 <- ggplot(meta, aes(x = disease_duration_years, fill = group)) +
        geom_density(alpha = 0.4) +
        labs(y = "Density", fill = NULL, x = "Disease duration to dysplasia") +
        scale_fill_manual(values = c("Non-progressor" = "darkseagreen", 
                               "Progressor" = "chocolate")) +
        theme_classic() +
        theme(legend.position = c(0.75, 0.7))
ggsave("plots/metadata/disease_duration_density.png", p3, width = 4, height = 4)

