library(readr)
library(dplyr)
library(ggplot2)

p53_status <- read_tsv("results/p53_status.tsv")
props <- read_tsv("data/progressor_precursors_cn_props.tsv") |>
  mutate(p53 = if_else(Sample %in% p53_status$Tumor_Sample_Barcode, "Mut", "WT")) |>
  left_join(p53_status %>% select(Tumor_Sample_Barcode, VAF_tum), by = c("Sample" = "Tumor_Sample_Barcode"))


# p <- ggplot(props, aes(x = p53, y = proportion, fill = p53)) +
#         geom_bar(stat = "identity", position = "dodge") +
#         labs(x = "TP53 Status", y = "Proportion") +
#         scale_fill_manual(values = c("WT" = "blue", "Mut" = "red")) +
#         theme_classic()
#  ggsave("plots/cn_props/progressor_precursors_p53_cn_prop.png", p)

p2 <- ggplot(props, aes(x = p53, y = proportion, fill = p53)) +
        geom_violin(trim = TRUE) +
        labs(x = "p53 Status", y = "CNA Proportion") +
        scale_fill_manual(values = c("WT" = "darkseagreen", "Mut" = "darkorange")) +
        theme_classic()
ggsave("plots/cn_props/progressor_precursors_p53_cn_prop_viol.png", p2,
        width = 4, height = 4)

write_tsv(props, "results/p53_cn_props.tsv")

# Follow-ups
p53_status <- read_tsv("results/p53_status.tsv")
props <- read_tsv("data/progressors_followups_cn_props.tsv") |>
  mutate(p53 = if_else(Sample %in% p53_status$Tumor_Sample_Barcode, "Mut", "WT")) |>
  left_join(p53_status %>% select(Tumor_Sample_Barcode, VAF_tum), by = c("Sample" = "Tumor_Sample_Barcode"))


p2 <- ggplot(props, aes(x = p53, y = proportion, fill = p53)) +
        geom_violin(trim = TRUE) +
        labs(x = "p53 Status", y = "CNA Proportion") +
        scale_fill_manual(values = c("WT" = "darkseagreen", "Mut" = "darkorange")) +
        theme_classic()
ggsave("plots/cn_props/progressor_followups_p53_cn_prop_viol.png", p2,
        width = 4, height = 4)

write_tsv(props, "results/p53_cn_props_fups.tsv")

