#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

library(dplyr)
options(stringsAsFactors = FALSE)

sample.file <- argv[1] # e.g., sample.file = 'data/raw/samples.txt'
masterid.file <- argv[2] # e.g., masterid.file = 'data/raw/masterid.csv'
out.file <- argv[3]

## Match sample names
masterid <- read.table(masterid.file, sep = ',', header = TRUE)
samples <- read.table(sample.file)
colnames(samples) <- c('subjectid', 'meth.pos', 'data.pos')

ret <- samples %>% left_join(masterid, by = 'subjectid') %>%
    select(IID, subjectid, meth.pos, gwas1709, methylation)

write.table(ret, file = gzfile(out.file), row.names = FALSE, col.names = FALSE, quote = FALSE,
            sep = '\t')
