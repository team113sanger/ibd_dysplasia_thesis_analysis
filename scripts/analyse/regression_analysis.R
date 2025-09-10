library(readr)
library(dplyr)
library(ggplot2)

# Read in data
pre_df <- read_tsv("results/precursor_combined_results.tsv") |>
    rename(PGA = proportion)

# Reformat to 0/1 outcomes
df <- pre_df |>
  dplyr::mutate(
    group = ifelse(group == "Progressor", 1, 0),
    TP53_status = ifelse(TP53_status == "Mut", 1, 0)
  )

# Linear regression
lm_model <- lm(PGA ~ group + TP53_status, data = df)

# Model summary
summary(lm_model)

# Diagnostic plots
pdf("results/regression_analysis/lm.pdf")
par(mfrow=c(2,2))
plot(lm_model)
dev.off()

# Logistic regression
glm_model <- glm(group ~ PGA + TP53_status,
                 data = df, family = binomial)

# Model summary
summary(glm_model)

# Odds ratios with confidence intervals
exp(cbind(OR = coef(glm_model), confint(glm_model)))


# Plot 
# PGA by progression status
p1 <- ggplot(df, aes(x = factor(group), y = PGA, fill = factor(group))) +
        geom_boxplot(alpha=0.6) +
        labs(x = "Progression status (0=non, 1=progressor)",
            y = "Proportion genome altered (PGA)") +
        theme_classic()
ggsave("results/regression_analysis/pga_by_group.png", p1)

# PGA by TP53 status
p2 <- ggplot(df, aes(x = factor(TP53_status), y = PGA, fill = factor(TP53_status))) +
        geom_boxplot(alpha=0.6) +
        labs(x = "TP53 status (0=WT, 1=mutated)",
            y = "Proportion genome altered (PGA)") +
        theme_classic()
ggsave("results/regression_analysis/pga_by_tp53.png", p2)
