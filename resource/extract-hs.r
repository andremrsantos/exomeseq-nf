library(tidyverse)
## Load arguments
args = commandArgs(trailingOnly=TRUE)
## Validate input
if (length(args) == 0) {
	stop("You must provide at least 1 arguments <ProjectDir>.")
} else if (length(args) == 1) {
	args[2] = "exomeseq.hs_metrics.tsv"
}

## Load data
report_dir   <- file.path(args[1], "reports")
report_files <- list.files(report_dir, pattern = ".hs_metrics")
file.path(report_dir, report_files) %>%
	set_names(report_files) %>%
	map_dfr(read_tsv, comment = "#", skip = 4, .id = "sample") %>%
	mutate(sample = gsub(".hs_metrics", "", sample)) %>%
	write_tsv(args[2]) 