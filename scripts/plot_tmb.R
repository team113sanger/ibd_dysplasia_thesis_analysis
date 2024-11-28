library(readr)
library(dplyr)
library(ggplot2)

metadata <- read_tsv("metadata/final_metadata_qc_pass.tsv")
tmb <- read_tsv("data/mutations_per_Mb.tsv", col_names = c("sanger_dna_id", "tmb")) |>
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
    filter(!grade_of_dysplasia == "NOS")

tmb[["grade_of_dysplasia"]] <- factor(
    tmb[["grade_of_dysplasia"]], 
    levels = c("LGD", "HGD", "AC")
    )

p <- ggplot(tmb, aes(x = grade_of_dysplasia, y = tmb, fill = grade_of_dysplasia)) +
        geom_violin(trim = TRUE, alpha = 0.6) +
        geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.3, fill = "white", color = "black") +
        geom_jitter(size = 2, alpha = 0.8, color = "grey25",
                    position = position_dodge(width = 0)) +
        facet_grid(~ group, scales = "free_x", space = "free") +
        labs(
            y = "TMB (Mutations/Mb)"
        ) +
        theme_bw(base_size = 14) +
        theme(legend.position = "none",
              axis.title.x = element_blank()) +
        scale_fill_brewer(palette = "Dark2") +
        scale_y_log10()

ggsave("plots/TMB/test.png", p)