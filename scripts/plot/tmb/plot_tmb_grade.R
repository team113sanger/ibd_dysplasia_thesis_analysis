library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)

samples <- read_lines("metadata/sample_lists/all_one_ppat.list")
metadata <- read_tsv("metadata/final_metadata_qc_pass.tsv") |>
  filter(sanger_dna_id %in% samples)

tmb <- read_tsv("data/variants/mutations_per_Mb.tsv", col_names = c("sanger_dna_id", "tmb")) |>
  left_join(metadata) |>
  filter(!sanger_dna_id == "PD62082c") |>
  mutate(
    grade_of_dysplasia = case_when(
      grade_of_dysplasia == "Low grade" ~ "LGD",
      grade_of_dysplasia == "High grade" ~ "HGD",
      grade_of_dysplasia == "Adenocarcinoma" ~ "AC",
      grade_of_dysplasia == "NOS" ~ "NOS",
      TRUE ~ NA_character_
    )
  ) |>
  filter(
    (group == "Non-progressor" & precursor_or_follow_up == "Precursor") |
      (group == "Progressor" & (
        grade_of_dysplasia != "LGD" |
          (grade_of_dysplasia == "LGD" & precursor_or_follow_up == "Precursor")
      ))
  ) |>
  filter(!grade_of_dysplasia == "NOS") |>
  select(sanger_dna_id, tmb, grade_of_dysplasia, group)

remap_tmb <- read_tsv("data/variants/ibd_crc_2018_remap/mutations_per_Mb.tsv",
  col_names = c("sanger_dna_id", "tmb")
) |>
  mutate(
    grade_of_dysplasia = "AC_2018",
    group = "Progressor"
  )

tmb <- rbind(tmb, remap_tmb)

tmb[["grade_of_dysplasia"]] <- factor(
  tmb[["grade_of_dysplasia"]],
  levels = c("LGD", "HGD", "AC", "AC_2018")
)

n_counts <- tmb |>
  group_by(group, grade_of_dysplasia) |>
  summarise(n = n(), .groups = "drop")

# Calculate Kruskal–Wallis p-values
tmb_sub1 <- tmb |> filter(grade_of_dysplasia %in% c("LGD", "HGD", "AC"))
tmb_sub2 <- tmb |> filter(grade_of_dysplasia %in% c("LGD", "HGD", "AC", "AC_2018"))

pval1 <- kruskal.test(tmb ~ grade_of_dysplasia, data = tmb_sub1)$p.value
pval2 <- kruskal.test(tmb ~ grade_of_dysplasia, data = tmb_sub2)$p.value

# Format nicely
pval1_label <- paste0("Kruskal-Wallis (LGD–HGD–AC): p = ", signif(pval1, 2))
pval2_label <- paste0("Kruskal-Wallis (LGD–HGD–AC–AC_2018): p = ", signif(pval2, 2))

# Count samples per group
n_counts <- tmb |>
  group_by(grade_of_dysplasia) |>
  summarise(n = n(), .groups = "drop")

# Plot
p <- ggplot(tmb, aes(x = grade_of_dysplasia, y = tmb, fill = grade_of_dysplasia)) +
  geom_violin(alpha = 0.3, color = NA) +
  geom_boxplot(alpha = 0.8, width = 0.5, outlier.shape = NA) +
  geom_jitter(size = 0.5, alpha = 0.8, width = 0.15) +
  geom_text(
    data = n_counts,
    aes(x = grade_of_dysplasia, y = 30, label = paste0("n=", n)),
    inherit.aes = FALSE,
    vjust = 1.5,
    color = "black",
    size = 3
  ) +
  labs(
    y = "TMB (Mut/Mb)",
    fill = NULL
  ) +
  theme_bw(base_size = 10) +
  scale_y_log10() +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  # Add the two Kruskal test labels manually
  annotate("text", x = 2, y = 400, label = pval1_label, size = 3, fontface = "italic") +
  annotate("text", x = 2.5, y = 250, label = pval2_label, size = 3, fontface = "italic") +
  scale_fill_brewer(palette = "Dark2") 

ggsave("plots/tmb/tmb_by_grade_with_group.png", p, width = 4, height = 4)


#### Just by grade ####
tmb <- tmb |> filter(grade_of_dysplasia != "AC_2018")

n_counts <- tmb |>
  group_by(grade_of_dysplasia) |>
  summarise(n = n(), .groups = "drop")

p <- ggplot(tmb, aes(x = grade_of_dysplasia, y = tmb, fill = grade_of_dysplasia)) +
  geom_violin(alpha = 0.3, color = NA) +   # transparent violin behind
  geom_boxplot(alpha = 0.8, width = 0.5, outlier.shape = NA) +  # narrower boxes
  geom_jitter(size = 0.5, alpha = 0.8, width = 0.15) +  # optional points
  geom_text(
    data = n_counts,
    aes(x = grade_of_dysplasia, y = 30, label = paste0("n=", n)),
    inherit.aes = FALSE,
    vjust = 1.5,
    color = "black",
    size = 3
  ) +
  labs(
    y = "TMB (Mut/Mb)", fill = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  scale_fill_brewer(palette = "Dark2") +
  scale_y_log10() +
  stat_compare_means(
    method = "wilcox.test", 
    comparisons = list(
      c("LGD", "HGD"),
      c("LGD", "AC"),
      c("HGD", "AC")
    ),
    size = 2.5
  ) +
  stat_compare_means(method = "kruskal.test", label.y = 3, size = 3) 

ggsave("plots/tmb/tmb_by_grade.png", p, width = 4, height = 4)