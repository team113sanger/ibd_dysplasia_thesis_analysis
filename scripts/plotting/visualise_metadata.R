library(ggplot2)
library(cowplot)
library(dplyr)
library(readr)
library(grid)
library(gridExtra)
library(rcartocolor)

metadata <- read_tsv("metadata/final_metadata_qc_pass.tsv")

p1 <- ggplot(metadata, aes(x = sex, fill = sex)) +
  geom_bar() +
  scale_fill_carto_d(palette = "Fall") + 
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank(),
        axis.line.y = element_line(color = 'black'), axis.ticks.y = element_line(color = 'black')) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Sex") +
  guides(fill = "none") 

p2 <- ggplot(metadata, aes(x = site, fill = site)) +
  geom_bar() +
  #scale_fill_manual(values = pastel_palette) +  # Use more colors for site categories
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank(),
        axis.line.y = element_line(color = 'black'), axis.ticks.y = element_line(color = 'black')) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Site") +
  guides(fill = "none") 

p3 <- ggplot(metadata, aes(x = grade_of_dysplasia, fill = grade_of_dysplasia)) +
  geom_bar() +
  scale_fill_carto_d(palette = "Fall") + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.line.y = element_line(color = 'black'),
    axis.ticks.y = element_line(color = 'black'),
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Grade") +
  guides(fill = "none") 

p4 <- ggplot(metadata, aes(x = group, fill = precursor_or_follow_up)) +
  geom_bar() +
  scale_fill_carto_d(palette = "Fall") + 
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank(),
        axis.line.y = element_line(color = 'black'), axis.ticks.y = element_line(color = 'black'),
        legend.position = "top") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Group") +
  guides(fill = guide_legend(title = NULL))

# Age histogram
pdf("plots/metadata/age.pdf", width = 4, height = 4)
hist(metadata$age, main = "Age", xlab = "Age", ylab = NULL, col = "peachpuff4", border = "white")
dev.off()

# Arrange all plots together
plot <- plot_grid(p1, p2, p3, p4, nrow = 2, label_size = 12)
y.grob <- textGrob("Frequency", gp=gpar(fontsize=12), rot=90)
plot <- grid.arrange(arrangeGrob(plot, left = y.grob))
ggsave("plots/metadata/metadata_plots.pdf", plot = plot, height = 6, width = 6)