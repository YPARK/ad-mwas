#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

library(dplyr)

sample.file <- argv[1] # e.g., 'data/raw/samples.txt'
masterid.file <- argv[2] # e.g., 'data/raw/masterid.csv'
out.file <- argv[3]

## Match sample names
options(stringsAsFactors = FALSE)

masterid <- read.table(masterid.file, sep = ',', header = TRUE)
samples <- read.table(sample.file)
samples <- cbind(samples, meth.pos = 1:NROW(samples))

colnames(samples) <- c('subjectid', 'meth.pos')

ret <- samples %>% left_join(masterid, by = 'subjectid') %>%
    select(IID, subjectid, meth.pos, gwas1709, methylation)

write.table(ret, file = gzfile(out.file), row.names = FALSE, col.names = FALSE, quote = FALSE,
            sep = '\t')
