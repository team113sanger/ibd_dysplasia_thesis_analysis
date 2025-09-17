library(readr)
library(ggplot2)

sample_list <- read_lines("metadata/sample_lists/all_one_ppat.list")

meta <- read_tsv("metadata/final_metadata_qc_pass.tsv") |>
            filter(sanger_dna_id %in% sample_list) |>
            mutate(patient_id = str_remove(sanger_dna_id, "[a-z]$")) |> #keep to remove multiple lesions per patient
            distinct(patient_id, .keep_all = TRUE) |>
            mutate(site_general = case_when(
                site %in% c("Caecum", "Ascending", "Proximal Ascending", "Distal Ascending") ~ "Right colon",
                site %in% c("Hepatic Flexure", "Transverse", "Splenic Flexure") ~ "Transverse colon",
                site %in% c("Descending", "Sigmoid", "Rectosigmoid") ~ "Left colon",
                site %in% c("Rectum") ~ "Rectum",
                TRUE ~ "Other"
            ))

p <- ggplot(meta, aes(x = patient_id, y = site_general, colour = ibd_diagnosis, shape = sex)) +
        geom_point(size = 2, alpha = 0.99) +
        facet_wrap(~ group, scales = "free_x") +
        theme_bw() +
        scale_colour_brewer(palette = "Dark2") +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6),
            panel.grid.major = element_line(colour = "grey90"),
            panel.grid.minor = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold"),
            text = element_text(family = "serif", colour = "black")
        ) +
        labs(
            x = NULL,
            y = NULL,
            colour = "Sex",
            shape = "IBD Diagnosis"
        )
ggsave("plots/metadata/patient_overview.png", p, width = 6.8, height = 3)

