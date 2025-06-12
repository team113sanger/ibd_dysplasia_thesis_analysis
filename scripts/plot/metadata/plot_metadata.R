library(ggplot2)
library(patchwork)
library(stringr)

# Read metadata
metadata_raw <- read_tsv("metadata/final_metadata_qc_pass.tsv")

metadata_raw <- metadata_raw %>%
  mutate(
    precursor_or_follow_up = factor(
      precursor_or_follow_up,
      levels = c("Precursor", "Follow up")
    )
  )

# Extract patient ID from sanger_dna_id by removing the trailing letter (e.g. PD62033a -> PD62033)
metadata_raw <- metadata_raw %>%
  mutate(patient_id = str_remove(sanger_dna_id, "[a-zA-Z]$"))

# Separate into Progressor and Non-progressor groups
prog_meta <- metadata_raw %>%
  filter(group == "Progressor") %>%
  mutate(Group = "Progressor")

nonprog_meta <- metadata_raw %>%
  filter(group == "Non-progressor") %>%
  mutate(Group = "Non-progressor")

# Helper function for basic barplots
make_barplot <- function(data, var, title) {
  ggplot(data, aes_string(x = var)) +
    geom_bar(fill = "steelblue") +
    labs(title = title, y = "Sample Count", x = NULL) +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
}

# Barplot coloured by sex for "Status" (precursor or follow-up)
make_status_plot <- function(data, title) {
  ggplot(data, aes(x = precursor_or_follow_up, fill = sex)) +
    geom_bar(position = "stack") +
    labs(title = title, y = "Sample Count", x = NULL, fill = "Sex") +
    scale_fill_manual(values = c("mediumpurple4", "darkolivegreen4")) +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
}

# IBD diagnosis plot using unique patients only
make_ibd_unique_plot <- function(data, title) {
  data_unique <- data %>%
    distinct(patient_id, .keep_all = TRUE)

  ggplot(data_unique, aes(x = ibd_diagnosis)) +
    geom_bar(fill = "steelblue") +
    labs(title = title, y = "Patient Count", x = NULL) +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
}

# Create all plots

# Status by sex
status_plot_nonprog <- make_status_plot(nonprog_meta, "Non-Progressors\nStatus by Sex")
status_plot_prog <- make_status_plot(prog_meta, "Progressors\nStatus by Sex")

# IBD Diagnosis (unique patients)
ibd_plot_nonprog <- make_ibd_unique_plot(nonprog_meta, "IBD Diagnosis")
ibd_plot_prog <- make_ibd_unique_plot(prog_meta, "IBD Diagnosis")

# Combine and save
p <- (status_plot_nonprog | ibd_plot_nonprog) /
  (status_plot_prog | ibd_plot_prog)

ggsave("plots/metadata/metadata_plots.png", p, dpi = 300, width = 6, height = 7)
