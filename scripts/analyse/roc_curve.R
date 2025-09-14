library(pROC)
library(dplyr)
library(readr)
library(ggplot2)

# Read in data
df <- read_tsv("results/precursor_combined_results.tsv") |>
    select(sanger_dna_id, group, proportion)

# Convert to binary outcome (1 = Progressor, 0 = Non-progressor)
df <- df |>
  mutate(group_outcome = ifelse(group == "Progressor", 1, 0)) 

# ROC analysis
roc_obj <- roc(df$group_outcome, df$proportion,
               levels = c(0, 1), direction = "<")

# See results
# Area under curve (AUC)
auc_val <- auc(roc_obj)
auc_val

# See all results
all <- coords(roc_obj, ret = "all", transpose = FALSE) 

# Best threshold (Youden index)
best_thresh <- coords(
    roc_obj, "best",
    ret = c("threshold", "sensitivity", "1-specificity", "specificity", "precision", "recall"),
    transpose = FALSE
    )
best_thresh

# Plot ROC
# pdf("results/regression_analysis/roc.pdf")
# plot(roc_obj, print.auc = TRUE, col = "darkorange", lwd = 2, legacy.axes = TRUE)
# dev.off()

p <- ggroc(list(roc_obj), size = 1, legacy.axes = TRUE)+
        theme_classic(base_size = 12) +
        geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey", linetype="dashed") + #random classifier baseline
        theme(axis.text.x = element_text(vjust = 0.5, hjust=1, colour="black", size=10),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_blank(),
              legend.position="none") +
        annotate("text", x=0.7, y=0.2, label= paste0("functional score\n\ AUC: ",round(auc_val,4)),family = "mono", colour="black") +
        scale_colour_manual(values = c("black")) +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
        annotate("point", x=best_thresh$`1-specificity`, y=best_thresh$sensitivity, colour="red", size=2.5) +
        annotate("text", x=0.7, y=(best_thresh$sensitivity-0.05), label= paste0(" threshold: ",round(best_thresh$threshold,5)),family = "mono", colour="black") +
        annotate("text", x=0.7, y=(best_thresh$sensitivity-0.1), label= paste0("  sensitivity: ", round(best_thresh$sensitivity,5)),family = "mono", colour="black") +
        annotate("text", x=0.7, y=(best_thresh$sensitivity-0.15), label= paste0("1-specificity: ", round(best_thresh$`1-specificity`,5)),family = "mono", colour="black") +
        ylab("Sensitivity (True Positive Rate)") +
        xlab("1-Specificity (False Positve Rate)")

ggsave("results/regression_analysis/roc_curve.png", p, dpi = 300)

# Get positve and negative predictive values 
# Apply threshold to classify samples
df <- df |>
  mutate(pred_progression = ifelse(proportion >= best_thresh$threshold, "Progressor", "Non-progressor"))

# Confusion matrix
conf <- table(Observed = df$group, Predicted = df$pred_progression)
conf

# Extract values
TP <- conf["Progressor", "Progressor"]
FP <- conf["Non-progressor", "Progressor"]
TN <- conf["Non-progressor", "Non-progressor"]
FN <- conf["Progressor", "Non-progressor"]

# Predictive values
PPV <- TP / (TP + FP)
NPV <- TN / (TN + FN)

PPV
NPV

# Write results
# Build results summary
results_summary <- tibble::tibble(
  Metric = c("AUC",
             "Threshold",
             "Sensitivity",
             "1-Specificity",
             "PPV",
             "NPV",
             "TP", "FP", "TN", "FN"),
  Value = c(
    as.numeric(auc_val),
    best_thresh$threshold,
    best_thresh$sensitivity,
    best_thresh$`1-specificity`,
    PPV,
    NPV,
    TP, FP, TN, FN
  )
)

write_tsv(results_summary, "results/regression_analysis/roc_metrics.tsv")

