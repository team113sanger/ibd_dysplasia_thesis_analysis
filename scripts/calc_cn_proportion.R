library(dplyr)

# Sample  chr     startpos        endpos  nMajor  nMinor  Sex     Purity  Ploidy  CN      Size
# PD56526a        1       826681  248918428       1       1       M       0.76    2.07625432667259        neutral 248091748
# PD56526a        2       41509   11344147        1       1       M       0.76    2.07625432667259        neutral 11302639
# PD56526a        2       11446553        28850976        1       0       M       0.76    2.07625432667259        loss    17404424

args <- commandArgs(trailing = TRUE)

file <- args[1]
outfile <- args[2]

if (is.na(file) | is.na(outfile)) {
  stop("No inputs provided.\n\nUsage: Rscript calc_cn_proportion.R  input_file output_file\n\n")
}


# Read in ASCAT annotated file

df <- read.table(file, header = T, sep = "\t", stringsAsFactors = F)

# Get the total size of all segments by sample and CN type (gain/loss/neutral)

df_grouped <- df %>%
  group_by(Sample, CN) %>%
  summarise(tot = sum(Size))

# Get non-neutral CN segments by samples and get the total segment size

gainloss <- ungroup(df_grouped) %>%
  filter(CN != "neutral") %>%
  group_by(Sample) %>%
  summarise(gain_loss = sum(tot))

# Get the neutral CN segments and sum total of these segments

neutral <- df_grouped %>%
  filter(CN == "neutral") %>%
  select(-CN) %>%
  rename("neutral" = "tot")

joined <- gainloss %>%
  full_join(., neutral, by = "Sample") %>%
  replace(is.na(.), 0)

# Calculate proportion
joined$proportion <- joined$gain_loss / (joined$gain_loss + joined$neutral)

# Write output

write.table(data.frame(joined), file = outfile, sep = "\t", row.names = F, quote = F)
