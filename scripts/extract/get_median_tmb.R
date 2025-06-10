# Median TMB for LGD Progressors
median_LGD_prog <- tmb %>%
  filter(grade_of_dysplasia == "LGD", group == "Progressor") %>%
  summarise(median_tmb = median(tmb)) %>%
  pull(median_tmb)

# Median TMB for LGD Non-Progressors
median_LGD_nonprog <- tmb %>%
  filter(grade_of_dysplasia == "LGD", group == "Non-progressor") %>%
  summarise(median_tmb = median(tmb)) %>%
  pull(median_tmb)

# Median TMB for HGD Progressors
median_HGD_prog <- tmb %>%
  filter(grade_of_dysplasia == "HGD", group == "Progressor") %>%
  summarise(median_tmb = median(tmb)) %>%
  pull(median_tmb)

# Median TMB for AC Progressors
median_AC_prog <- tmb %>%
  filter(grade_of_dysplasia == "AC", group == "Progressor") %>%
  summarise(median_tmb = median(tmb)) %>%
  pull(median_tmb)

# Print results
median_LGD_prog
median_LGD_nonprog
median_HGD_prog
median_AC_prog
