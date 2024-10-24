library(readr)
library(dplyr)

# Read in metadata
meta <- read_tsv("metadata/final_metadata_qc_pass.tsv")

get_samples <- function(metadata_df, filter_variable, filter_value,
                        filter_variable2 = NULL, filter_value2 = NULL) {
    filtered_df <- metadata_df |>
        filter({{filter_variable}} == filter_value)
    
    if (!is.null(filter_value2)) {
        filtered_df <- filtered_df |>
            filter({{filter_variable2}} == filter_value2)
    }

    sample_list <- filtered_df |>
        pull(sanger_dna_id)
    
    return(sample_list)
}

outdir <- "metadata/sample_lists"

#### Sample lists by Group ####
## Progressors ##
# All
progressors <- get_samples(meta, group, "Progressor")
write_lines(progressors, file.path(outdir, "progressor_samples.tsv"))

# Precursors
progressor_precursors <- get_samples(meta, group, "Progressor",
                                    precursor_or_follow_up, "Precursor")
write_lines(progressor_precursors, file.path(outdir, "progressor_precursor_samples.tsv"))

# Follow Ups
progressor_follow_ups <- get_samples(meta, group, "Progressor",
                                    precursor_or_follow_up, "Follow up")
write_lines(progressor_follow_ups, file.path(outdir, "progressor_follow_up_samples.tsv"))

## Non-progressors ##
# All
non_progressors <- get_samples(meta, group, "Non-progressor")
write_lines(non_progressors, file.path(outdir, "non_progressor_samples.tsv"))

# Precursors
non_progressor_precursors <- get_samples(meta, group, "Non-progressor",
                                          precursor_or_follow_up, "Precursor")
write_lines(non_progressor_precursors, file.path(outdir, "non_progressor_precursor_samples.tsv"))

# Follow Ups
non_progressor_follow_ups <- get_samples(meta, group, "Non-progressor",
                                        precursor_or_follow_up, "Follow up")
write_lines(non_progressor_follow_ups, file.path(outdir, "non_progressor_follow_up_samples.tsv"))

#### Sample lists by grade ####
low_grade <- get_samples(meta, grade_of_dysplasia, "Low grade")
write_lines(low_grade, file.path(outdir, "low_grade_samples.tsv"))

high_grade <- get_samples(meta, grade_of_dysplasia, "High grade")
write_lines(high_grade, file.path(outdir, "high_grade_samples.tsv"))

carinoma <- get_samples(meta, grade_of_dysplasia, "Adenocarcinoma")
write_lines(carinoma, file.path(outdir, "adenocarcinoma_samples.tsv"))