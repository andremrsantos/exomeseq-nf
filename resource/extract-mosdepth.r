library(tidyverse)
## Load arguments
args = commandArgs(trailingOnly=TRUE)
## Validate input
if (length(args) == 0) {
	stop("You must provide at least 1 arguments <ProjectDir>.")
} else if (length(args) == 1) {
	args[2] = "exomeseq.mosdepth.dist.tsv"
}

## Load data
report_dir   <- file.path(args[1], "reports", "cov")
report_files <- list.files(report_dir, pattern = ".mosdepth.dist.txt")
file.path(report_dir, report_files) %>%
	set_names(report_files) %>%
	map_dfr(read_tsv, 
		col_types = "cin",
		col_names = c("region", "coverage", "proportion"), 
		.id = "sample") %>%
	filter(region == "total") %>%
	mutate(sample = gsub(".mosdepth.dist.txt", "", sample)) %>%
	write_tsv(args[2]) 