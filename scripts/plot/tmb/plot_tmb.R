library(readr)
library(dplyr)
library(ggplot2)

metadata <- read_tsv("metadata/final_metadata_qc_pass.tsv")
tmb <- read_tsv("data/variants/mutations_per_Mb.tsv", col_names = c("sanger_dna_id", "tmb")) |>
  left_join(metadata) |>
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
    (group == "Non-progressor" & precursor_or_follow_up == "Follow up") |
    (group == "Progressor" & (
      grade_of_dysplasia != "LGD" |
      (grade_of_dysplasia == "LGD" & precursor_or_follow_up == "Precursor")
    ))
  )

tmb[["grade_of_dysplasia"]] <- factor(
    tmb[["grade_of_dysplasia"]], 
    levels = c("LGD", "HGD", "AC")
    )

n_counts <- tmb %>%
    group_by(group, grade_of_dysplasia) %>%
    summarise(n = n(), .groups = "drop")

p <- ggplot(tmb, aes(x = grade_of_dysplasia, y = tmb, fill = grade_of_dysplasia)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.8) +
        # geom_jitter(size = 1, alpha = 0.8) +
       #             position = position_dodge(width = 0)) +
        facet_grid(~ group, scales = "free_x", space = "free") +
        geom_text(
            data = n_counts,
            aes(x = grade_of_dysplasia, y = 70, label = paste0("n=", n)),
            inherit.aes = FALSE,
            vjust = 1.5,
            color = "black",
            size = 4
        ) +
        labs(
            y = "TMB (Mutations/Mb)"
        ) +
        theme_bw(base_size = 14) +
        theme(legend.position = "none",
              axis.title.x = element_blank()) +
        scale_fill_brewer(palette = "Dark2") +
        scale_y_log10()

ggsave("plots/tmb/tmb_per_group.png", p, width = 5.6, height = 6)

#### Add 2018 remapped adenocarcinomas #### 
remap_tmb <- read_tsv("/lustre/scratch125/casm/team113da/projects/IBD_Associated_Dysplasia/3361_3511_IBD_CRC_ReMap/analysis/variants_combined/release_v1/all_tumours/matched_samples/mutations_per_Mb.tsv",
                        col_names = c("sanger_dna_id", "tmb")) |>
                mutate(
                    grade_of_dysplasia = "AC_2018",
                    group = "Progressor"
                )

tmb <- tmb |>
        select("sanger_dna_id", "tmb", "grade_of_dysplasia", "group")

combined_tmb <- bind_rows(tmb, remap_tmb)

combined_tmb[["grade_of_dysplasia"]] <- factor(
    combined_tmb[["grade_of_dysplasia"]], 
    levels = c("LGD", "HGD", "AC", "AC_2018")
    )

n_counts <- combined_tmb %>%
    group_by(group, grade_of_dysplasia) %>%
    summarise(n = n(), .groups = "drop")

p2 <- ggplot(combined_tmb, aes(x = grade_of_dysplasia, y = tmb, fill = grade_of_dysplasia)) +
            geom_boxplot( alpha = 0.8) +
            # geom_dotplot(binaxis='y', position=position_jitterdodge(jitter.width=0, dodge.width = 0.3), dotsize=0.5) +
        #    geom_jitter(size = 1, alpha = 0.8) +
        #             position = position_dodge(width = 0)) +
            facet_grid(~ group, scales = "free_x", space = "free") +
            geom_text(
                data = n_counts,
                aes(x = grade_of_dysplasia, y = 200, label = paste0("n=", n)),
                inherit.aes = FALSE,
                color = "black",
                size = 4
            ) +
            labs(
                y = "TMB (Mutations/Mb)"
            ) +
            theme_bw(base_size = 14) +
            theme(legend.position = "none",
                axis.title.x = element_blank(),
                strip.text = element_text(size = 10)) +
            scale_fill_brewer(palette = "Dark2") +
            scale_y_log10()

ggsave("plots/TMB/tmb_per_group_with_remap.png", p2, width = 6, height = 5)

