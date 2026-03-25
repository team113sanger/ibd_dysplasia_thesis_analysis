library(readr)
library(dplyr)

# Read in metadata
meta <- read_tsv("metadata/final_metadata_qc_pass.tsv")

get_samples <- function(metadata_df, filter_variable, filter_value,
                        filter_variable2 = NULL, filter_value2 = NULL) {
  filtered_df <- metadata_df |>
    filter({{ filter_variable }} == filter_value)

  if (!is.null(filter_value2)) {
    filtered_df <- filtered_df |>
      filter({{ filter_variable2 }} == filter_value2)
  }

  sample_list <- filtered_df |>
    select(sanger_dna_id, study_id)

  return(sample_list)
}

outdir <- "metadata/sample_lists"

extra_samples <- c("PD62065a", "PD62065e", "PD62068d", "PD62068e", "PD62075c", "PD62075e", "PD62083d", #non-prog fol
                   "PD62064c", "PD62068a", "PD62068c", "PD62081c", #non-prog pre
                   "PD62030c", "PD62036a", "PD62047g", "PD62048c", #prog fol
                   "PD62031d", "PD62033c", "PD62041c") # prog pre
#### Sample lists by Group ####
## Progressors ##
# All
progressors <- get_samples(meta, group, "Progressor")
write_tsv(progressors, file.path(outdir, "progressor_samples.tsv"))

progressors_ppat <- get_samples(meta, group, "Progressor") |>
  filter(!sanger_dna_id %in% extra_samples)
write_tsv(progressors_ppat, file.path(outdir, "progressor_samples_ppat.tsv"))

# Precursors
progressor_precursors <- get_samples(
  meta, group, "Progressor",
  precursor_or_follow_up, "Precursor"
)
write_tsv(progressor_precursors, file.path(outdir, "progressor_precursor_samples.tsv"))

progressor_precursors_ppat <- progressor_precursors |> 
  filter(!sanger_dna_id %in% extra_samples)
write_tsv(progressor_precursors_ppat, file.path(outdir, "progressor_precursor_samples_ppat.tsv"))

# Follow Ups
progressor_follow_ups <- get_samples(
  meta, group, "Progressor",
  precursor_or_follow_up, "Follow up"
)
write_tsv(progressor_follow_ups, file.path(outdir, "progressor_follow_up_samples.tsv"))

progressor_follow_ups_ppat <- progressor_follow_ups |> 
  filter(!sanger_dna_id %in% extra_samples)
write_tsv(progressor_follow_ups_ppat, file.path(outdir, "progressor_follow_up_samples_ppat.tsv"))

## Non-progressors ##
# All
non_progressors <- get_samples(meta, group, "Non-progressor")
write_tsv(non_progressors, file.path(outdir, "non_progressor_samples.tsv"))

non_progressors_ppat <- non_progressors |> 
  filter(!sanger_dna_id %in% extra_samples)
write_tsv(non_progressors_ppat, file.path(outdir, "non_progressor_samples_ppat.tsv"))

# Precursors
non_progressor_precursors <- get_samples(
  meta, group, "Non-progressor",
  precursor_or_follow_up, "Precursor"
)
write_tsv(non_progressor_precursors, file.path(outdir, "non_progressor_precursor_samples.tsv"))

non_progressor_precursors_ppat <- non_progressor_precursors |> 
  filter(!sanger_dna_id %in% extra_samples)
write_tsv(non_progressor_precursors_ppat, file.path(outdir, "non_progressor_precursor_samples_ppat.tsv"))

# Follow Ups
non_progressor_follow_ups <- get_samples(
  meta, group, "Non-progressor",
  precursor_or_follow_up, "Follow up"
)
write_tsv(non_progressor_follow_ups, file.path(outdir, "non_progressor_follow_up_samples.tsv"))

non_progressor_follow_ups_ppat <- non_progressor_follow_ups |> 
  filter(!sanger_dna_id %in% extra_samples)
write_tsv(non_progressor_follow_ups_ppat, file.path(outdir, "non_progressor_follow_up_samples_ppat.tsv"))

#### Sample lists by grade ####
low_grade <- get_samples(meta, grade_of_dysplasia, "Low grade")
write_tsv(low_grade, file.path(outdir, "low_grade_samples.tsv"))

low_grade_ppat <- low_grade |> 
  filter(!sanger_dna_id %in% extra_samples)
write_tsv(low_grade_ppat, file.path(outdir, "low_grade_samples_ppat.tsv"))

high_grade <- get_samples(meta, grade_of_dysplasia, "High grade")
write_tsv(high_grade, file.path(outdir, "high_grade_samples.tsv"))

high_grade_ppat <- high_grade |> 
  filter(!sanger_dna_id %in% extra_samples)
write_tsv(high_grade_ppat, file.path(outdir, "high_grade_samples_ppat.tsv"))

carinoma <- get_samples(meta, grade_of_dysplasia, "Adenocarcinoma")
write_tsv(carinoma, file.path(outdir, "adenocarcinoma_samples.tsv"))

carinoma_ppat <- carinoma |> 
  filter(!sanger_dna_id %in% extra_samples)
write_tsv(carinoma_ppat, file.path(outdir, "adenocarcinoma_samples_ppat.tsv"))

all_ppat <- rbind(non_progressors_ppat, progressors_ppat) |>
  write_tsv(file.path(outdir, "all_one_ppat.list"))