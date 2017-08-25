#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 3) { q() }

data.probes.file <- argv[1] # e.g., 'data/raw/all-probes.txt.gz'
chr.probe.file <- argv[2] # e.g., 'data/probes/chr21-probes.txt.gz'
out.file <- argv[3]

library(dplyr)
options(stringsAsFactors = FALSE)

data.probes <- read.table(data.probes.file, sep = '\t')
chr.probes <- read.table(chr.probe.file, sep = '\t')

colnames(data.probes) <- c('cg', 'data.pos')
colnames(chr.probes) <- c('cg', 'chr', 'loc')

ret <- chr.probes %>% left_join(data.probes, by = 'cg')

write.table(ret, file = gzfile(out.file), quote = FALSE, sep = '\t', col.names = FALSE,
            row.names = FALSE)
