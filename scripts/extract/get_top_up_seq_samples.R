library(readr)
library(dplyr)

cov <- read_tsv("data/qc/7100-coverage.txt")

rescue_samples <- cov |>
  filter(pass == FALSE) |>
  arrange(desc(`21+`)) |>
  filter(`21+` > 70) |>
  select(sample, `21+`)

write_tsv(rescue_samples, "metadata/rescue_samples.tsv")