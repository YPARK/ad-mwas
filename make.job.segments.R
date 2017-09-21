#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)
cpg.file <- argv[1] # e.g., cpg.file <- 'data/probes/chr21-probes.txt.gz'

options(stringsAsFactors = FALSE)
library(dplyr)

block.size0 <- 1e6
max.size <- 20

cpg.cols <- c('cg', 'chr', 'cg.loc')

## 1. assign cpgs into blocks
cpgs <- read.table(cpg.file, col.names = cpg.cols) %>%
    mutate(idx = 1:n()) %>%
        mutate(block.size = block.size0, block = round(cg.loc/block.size)) %>%
            arrange(cg.loc)

## 2. break large blocks into smaller chunks
while(TRUE) {
    big.blocks <- cpgs %>% group_by(block) %>%
        summarize(block.size = n()) %>%
            filter(block.size > max.size) %>% select(block)

    if(nrow(big.blocks) == 0) break

    ret1 <- cpgs %>% filter(! block %in% big.blocks$block)

    ret2 <-
        cpgs %>% right_join(big.blocks, by = 'block') %>%
            mutate(block.size = block.size / 2) %>%
                mutate(block = round(cg.loc / block.size * 4) / 4)

    cpgs <- rbind(ret1, ret2) %>% arrange(cg.loc)
}

out <- cpgs %>% group_by(block) %>%
    summarize(idx.vec = paste('c(',paste(idx, collapse = ','),')', sep=''))

cat(paste(out$idx.vec, collapse = '\n'), file = stdout())
cat('\n', file = stdout())
