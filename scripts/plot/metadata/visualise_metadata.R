library(ggplot2)
library(cowplot)
library(dplyr)
library(readr)
library(grid)
library(gridExtra)
library(rcartocolor)

metadata_raw <- read_tsv("metadata/final_metadata_qc_pass.tsv")

metadata <- metadata_raw |>
  filter(precursor_or_follow_up == "Precursor") |>
  mutate(
    colon_region = case_when(
      site %in% c("Caecum", "Ascending", "Proximal Ascending", "Distal Ascending", "Transverse") ~ "Right colon",
      site %in% c("Sigmoid", "Splenic Flexure") ~ "Left colon",
      site %in% c("Rectum", "Rectosigmoid") ~ "Rectum",
      TRUE ~ "Other" # In case there are any sites not listed
    )
  )


p1 <- ggplot(metadata, aes(x = sex, fill = sex)) +
  geom_bar() +
  scale_fill_carto_d(palette = "Fall") + 
  facet_wrap(~ group) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank(),
        axis.line.y = element_line(color = 'black'), axis.ticks.y = element_line(color = 'black')) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Sex") +
  guides(fill = "none") 
  ggsave("plots/metadata/metadata_sex_plot.pdf", plot = p1, height = 6, width = 6)

p2 <- ggplot(metadata, aes(x = colon_region, fill = colon_region)) +
  geom_bar() +
  scale_fill_carto_d(palette = "Fall") + 
  theme_minimal() +
  facet_wrap(~ group) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank(),
        axis.line.y = element_line(color = 'black'), axis.ticks.y = element_line(color = 'black')) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Site") +
  guides(fill = "none") 
  ggsave("plots/metadata/metadata_site_plot.pdf", plot = p2, height = 6, width = 6)

  # Arrange all plots together
  plot <- plot_grid(p1, p2, nrow = 2, label_size = 12)
  y.grob <- textGrob("Frequency", gp=gpar(fontsize=12), rot=90)
  plot <- grid.arrange(arrangeGrob(plot, left = y.grob))
  ggsave("plots/metadata/metadata_plots.pdf", plot = plot, height = 5, width = 5)

# p3 <- ggplot(metadata, aes(x = grade_of_dysplasia, fill = grade_of_dysplasia)) +
#   geom_bar() +
#   scale_fill_carto_d(palette = "Fall") + 
#   theme_minimal() +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.title.y = element_blank(),
#     axis.title.x = element_blank(),
#     axis.line.y = element_line(color = 'black'),
#     axis.ticks.y = element_line(color = 'black'),
#     axis.text.x = element_text(angle = 30, hjust = 1)
#   ) +
#   scale_y_continuous(expand = c(0, 0)) +
#   labs(title = "Grade") +
#   guides(fill = "none") 


p4 <- ggplot(metadata_raw, aes(x = group, fill = precursor_or_follow_up)) +
  geom_bar() +
  scale_fill_carto_d(palette = "Fall") + 
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.line.y = element_line(color = 'black'), axis.ticks.y = element_line(color = 'black'),
        legend.position = "top") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Group") +
  guides(fill = guide_legend(title = NULL))
ggsave("plots/metadata/metadata_group_plot.png", plot = p4, height = 4, width = 3, dpi = 300)

p5 <- ggplot(metadata, aes(x = group, y = age, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.8) +
  labs(y = "Age") +
  scale_fill_carto_d(palette = "Fall") + 
  theme_minimal() +
  theme(
    legend.position = "none", 
  #  axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("plots/metadata/metadata_age_plot.png", plot = p5, height = 4, width = 3, dpi = 300)



# Age histogram
# pdf("plots/metadata/age.pdf", width = 4, height = 4)
# hist(metadata$age, main = "Age", xlab = "Age", ylab = NULL, col = "peachpuff4", border = "white")
# dev.off()


# library(ggmosaic)

# p <- ggplot(metadata) +
#   geom_mosaic(aes(x = product(sex, colon_region), fill = colon_region)) +
#   facet_wrap(~ group) + # Facet by group
#   scale_fill_carto_d(palette = "Fall") +
#   theme_classic() +
#   theme(
#     panel.grid.major = element_blank(), 
#     panel.grid.minor = element_blank(),
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.title.y = element_blank(),
#     axis.title.x = element_blank(),
#     legend.position = "top"
#   ) +
#   labs(title = "Mosaic Plot of Sex and Colon Region by Group") +
#   guides(fill = guide_legend(title = "Colon Region"))

# # Save the mosaic plot as a PNG file
# ggsave("plots/metadata/metadata_combined_mosaic_plot.png", plot = p, height = 6, width = 8, dpi = 300)


# Stacked barplot for Sex
p1 <- ggplot(metadata, aes(x = group, fill = sex)) +
  geom_bar(position = "fill") +  # Use "fill" to stack proportionally
  scale_fill_carto_d(palette = "Fall") + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.line.y = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    legend.position = "top",  # Move legend to the top
    legend.box.margin = margin(t = 10, r = 20, b = 10, l = 10)  # Add margin around the legend
  ) +
  guides(fill = guide_legend(title = "Sex"))

# Stacked barplot for Colon Region (Site)
p2 <- ggplot(metadata, aes(x = group, fill = colon_region)) +
  geom_bar(position = "fill") +  # Use "fill" for proportional stacking
  scale_fill_carto_d(palette = "Fall") + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.line.y = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    legend.position = "top",  # Move legend to the top
    legend.box.margin = margin(t = 10, r = 20, b = 10, l = 10)  # Add margin around the legend
  ) +
  guides(fill = guide_legend(title = "Site"))

# Arrange the two plots together
plot <- plot_grid(p1, p2, ncol = 2, align = "v", label_size = 12)

# Add a shared y-axis label
y.grob <- textGrob("Proportion", gp = gpar(fontsize = 12), rot = 90)
plot <- grid.arrange(arrangeGrob(plot, left = y.grob))

# Save the combined plot
ggsave("plots/metadata/metadata_stacked_plots.png", plot = plot, height = 4, width = 6, dpi = 300)
