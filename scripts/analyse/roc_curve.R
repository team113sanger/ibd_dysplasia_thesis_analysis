library(pROC)
library(dplyr)
library(readr)
library(ggplot2)

set.seed(123) # for reproducibility

# Read in data
df <- read_tsv("results/precursor_combined_results.tsv") |>
  select(sanger_dna_id, group, cn_proportion) |>
  mutate(group_outcome = ifelse(group == "Progressor", 1, 0))

# Train-test split (75/25)
train_idx <- sample(seq_len(nrow(df)), size = 0.75 * nrow(df))
train_df <- df[train_idx, ]
test_df  <- df[-train_idx, ]

# ROC analysis on training set
roc_train <- roc(train_df$group_outcome, train_df$cn_proportion,
                 levels = c(0, 1), direction = "<")

auc_train <- auc(roc_train)

# Best threshold from training
best_thresh <- coords(
  roc_train, "best",
  ret = c("threshold", "sensitivity", "specificity"),
  transpose = FALSE
)

# Apply threshold to test set
test_df <- test_df |>
  mutate(pred_progression = ifelse(cn_proportion >= best_thresh$threshold, 1, 0))

# Confusion matrix on test set
conf <- table(Observed = test_df$group_outcome, Predicted = test_df$pred_progression)

TP <- conf["1", "1"]
FP <- conf["0", "1"]
TN <- conf["0", "0"]
FN <- conf["1", "0"]

PPV <- TP / (TP + FP)
NPV <- TN / (TN + FN)

# ROC on test set
roc_test <- roc(test_df$group_outcome, test_df$cn_proportion,
                levels = c(0, 1), direction = "<")
auc_test <- auc(roc_test)

p_both <- ggroc(list(Train = roc_train, Test = roc_test), size = 1, legacy.axes = TRUE) +
  theme_classic(base_size = 12) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color = "grey", linetype = "dashed") +
  scale_color_manual(values = c("Train" = "blue", "Test" = "red")) +
  labs(color = "Dataset") +
  annotate("text", x = 0.7, y = 0.3,
           label = paste0("Train AUC: ", round(auc_train, 3)),
           family = "mono", colour = "blue") +
  annotate("text", x = 0.7, y = 0.2,
           label = paste0("Test AUC: ", round(auc_test, 3)),
           family = "mono", colour = "red") +
  theme(axis.line = element_blank(),
        legend.position = c(0.9, 0.55)) +
  ylab("Sensitivity (TPR)") +
  xlab("1-Specificity (FPR)")

ggsave("results/regression_analysis/roc_curve_train_test.png", p_both, dpi = 300, 
        width = 5.5, height = 5.5)

# Results summary
results_summary <- tibble::tibble(
  Metric = c("Train AUC", "Test AUC", "Threshold (train)",
             "Sensitivity (train)", "Specificity (train)",
             "PPV (test)", "NPV (test)", "TP", "FP", "TN", "FN"),
  Value = c(
    as.numeric(auc_train),
    as.numeric(auc_test),
    best_thresh$threshold,
    best_thresh$sensitivity,
    best_thresh$specificity,
    PPV, NPV, TP, FP, TN, FN
  )
)

write_tsv(results_summary, "results/regression_analysis/roc_metrics_train_test.tsv")