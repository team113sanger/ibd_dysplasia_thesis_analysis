library(readr)
library(dplyr)
library(ggplot2)
library(stringr)

sample_list <- read_lines("metadata/sample_lists/all_one_ppat.list")
meta <- read_tsv("metadata/final_metadata_qc_pass.tsv") |>
            filter(sanger_dna_id %in% sample_list)

meta_time <- meta |>
  mutate(
    time_clean = case_when(
      str_detect(time_between_lesions, "^[0-9]+.*(year|month)") ~ time_between_lesions,
      TRUE ~ NA_character_
    ),
    time_clean = str_to_lower(time_clean),
    years = as.numeric(str_extract(time_clean, "(?<=\\b)\\d+(?=\\s*years?)")),
    months = as.numeric(str_extract(time_clean, "(?<=\\b)\\d+(?=\\s*months?)")),
    # convert to total months
    total_months = coalesce(years, 0) * 12 + coalesce(months, 0)
    )

time_stats <- meta_time |>
  group_by(group) |>
  summarise(
    n = n(),
    mean_months = mean(total_months, na.rm = TRUE),
    median_months = median(total_months, na.rm = TRUE),
    sd_months = sd(total_months, na.rm = TRUE),
    min_months = min(total_months, na.rm = TRUE),
    max_months = max(total_months, na.rm = TRUE),
    IQR = IQR(total_months, na.rm = TRUE)
  )
time_stats

p1 <- ggplot(meta_time, aes(x = group, y = total_months, fill = group)) +
        geom_boxplot(alpha = 0.6) +
        # geom_jitter(width = 0.2, alpha = 0.7) +
        scale_fill_manual(values = c("Non-progressor" = "darkseagreen", 
                            "Progressor" = "chocolate")) +
        labs(
            x = NULL,
            y = "Time Between Lesions (months)"
        ) +
        theme_classic() +
        theme(legend.position = "none")
ggsave("plots/metadata/time_between_boxplot.png", p1, width = 3, height = 3)

p2 <- ggplot(meta_time, aes(x = total_months, fill = group)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = c("Non-progressor" = "darkseagreen", 
                               "Progressor" = "chocolate")) +
  theme_classic() +
  labs(x = "Time Between Lesions (months)", y = "Density", fill = NULL) +
  theme(
    legend.position = c(0.75, 0.7)
  )
ggsave("plots/metadata/time_between_hist.png", p2, width =4, height = 4)


nonprgs <- meta_time %>%
  filter(group == "Non-progressor")

summary(nonprgs$total_months)
hist(nonprgs$total_months, breaks = 10, main = "Non-progressors: Time Between Lesions", xlab = "Months")

threshold <- quantile(nonprgs$total_months, 0.1, na.rm = TRUE)  # 10th percentile

tail_patients <- nonprgs %>%
  filter(total_months <= threshold)

p <- ggplot(nonprgs, aes(x = total_months)) +
        geom_histogram(binwidth = 5, fill = "darkseagreen", color = "black", alpha = 0.6) +
        geom_vline(xintercept = threshold, color = "red", linetype = "dashed") +
        labs(x = "Months", y = "Count") +
        theme_classic()
ggsave("plots/metadata/time_between_nprog.png", p, width =4, height = 4)
