library(readr)
library(ggplot2)
library(dplyr)

cov <- read_tsv("data/qc/7100-coverage.txt")

adenocarcinomas <- cov |>
  filter(sample %in% c(
    "PD62036a", "PD62030d", "PD62038e", "PD62045a", "PD62038d",
    "PD62039c", "PD62031c", "PD62051c", "PD62027c"
  ))

write_tsv(adenocarcinomas, "results/adenocarcinomas_cov.tsv")


# plot barplot of pass/failed samples
# Replace values in the 'pass' column using mutate
cov_plot <- cov |>
  mutate(pass = ifelse(pass == TRUE, "PASS", "FAIL"))

# Create the barplot
pass_fail_plot <- ggplot(cov_plot, aes(x = pass, fill = pass)) +
  geom_bar() +
  scale_fill_manual(values = c("PASS" = "#56B4E9", "FAIL" = "#E69F00")) + # Custom colours for pass/fail
  theme_minimal() +
  labs(
    y = "Sample Count" # Remove x-axis title, keep y-axis title
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.x = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16)
  )

# Save the plot
ggsave("plots/cov_pass_fail_barplot.png", plot = pass_fail_plot, height = 4, width = 3, dpi = 300)
