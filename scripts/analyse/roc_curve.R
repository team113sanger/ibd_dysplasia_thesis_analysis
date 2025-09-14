library(pROC)
library(plyr)
library(readr)
library(dplyr)
library(ggplot2)
library(PRROC)

# Read in data
df <- read_tsv("results/precursor_combined_results.tsv") |>
    select(sanger_dna_id, group, proportion)

result <-roc(df$group, df$proportion, ci=FALSE)

#can use this to get values and filter for which parameters you want
all <- coords(result, ret = "all", transpose = FALSE) #%>% select(precision, recall)
ideal <- coords(result, "best", ret=c("threshold", "sensitivity", "1-specificity", "specificity","precision","recall"))
write.csv(all,"results/regression_analysis/group_cna_output.csv", row.names = FALSE)
write.csv(ideal,"results/regression_analysis//group_cna_output_thresholds.csv", row.names = FALSE)

g <- ggroc(list(result), size = 1, legacy.axes = TRUE)+
        theme_classic()+
        #theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
        theme(axis.text.x = element_text(vjust = 0.5, hjust=1, colour="black", size=10))+
        theme(axis.text.y = element_text(vjust = 0.5, hjust=1, colour="black", size=10))+
        theme(axis.title.x = element_text(size = 12, colour = "black"))+
        theme(axis.title.y = element_text(size = 12, colour = "black"))+
        geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey", linetype="dashed")+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_blank())+
        theme(legend.position="none")+
        annotate("text", x=0.7, y=0.25, label= paste0("functional score\n\ AUC: ",round(result$auc,4)),family = "mono", colour="black") +
        scale_colour_manual(values = c("black"))+
       # scale_y_continuous(breaks = seq(0,0.25,0.5,0.75,1), limits = c(0,1))
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25))+
        annotate("point", x=ideal$`1-specificity`, y=ideal$sensitivity, colour="red", size=2.5)+
        annotate("text", x=0.7, y=(ideal$sensitivity-0.05), label= paste0(" threshold: ",round(ideal$threshold,5)),family = "mono", colour="black") +
        annotate("text", x=0.7, y=(ideal$sensitivity-0.1), label= paste0("  sensitivity: ", round(ideal$sensitivity,5)),family = "mono", colour="black") +
        annotate("text", x=0.7, y=(ideal$sensitivity-0.15), label= paste0("1-specificity: ", round(ideal$`1-specificity`,5)),family = "mono", colour="black") +
        ylab("Sensitivity")+
        xlab("1-Specificity")

ggsave("results/regression_analysis/roc.png", g)