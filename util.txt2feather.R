#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 2) q()

in.data.file <- argv[1]
out.data.file <- argv[2]

library(readr)
library(feather)

in.data <- read_tsv(in.data.file, progress = TRUE, col_names = FALSE)
write_feather(in.data, path = out.data.file)
