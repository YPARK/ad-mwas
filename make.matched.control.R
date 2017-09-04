#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 3) q()

source('util.R')
library(dplyr)
library(reshape2)
library(feather)
options(stringsAsFactors = FALSE)

matched.sample.file <- argv[1] # e.g., matched.sample.file = 'data/raw/matched.samples.txt.gz'
ctrl.file <- argv[2] # e.g., ctrl.file = 'data/raw/control-probes.txt.gz'
out.file <- argv[3] # e.g., out.ft

ctrl.raw.tab <- read.table(ctrl.file, header = TRUE, sep = '\t', check.names = FALSE)
samples.tab <- read.table(matched.sample.file, header = TRUE, sep = '\t', check.names = FALSE)
colnames(samples.tab) <- c('iid', 'meth.id', 'meth.pos', 'has.geno', 'has.meth')

.tag <- function(s) {
    ret <- strsplit(s, split = '[.]')[[1]]
    ret[length(ret)]
}

.rm.tag <- function(s) {
    ret <- strsplit(s, split = '[.]')[[1]]
    paste(ret[-length(ret)], sep = '')
}

name.tab <- data.frame(col.name = names(ctrl.raw.tab), col.pos = 1:ncol(ctrl.raw.tab)) %>%
    mutate(col.tag = sapply(col.name, .tag)) %>%
        mutate(col.name = sapply(col.name, .rm.tag))

.take.mat <- function(xx, tag.name) {
    ret <- ctrl.raw.tab %c% (name.tab$col.tag == tag.name) %>% t()
    rownames(ret) <- sapply(rownames(ret), .rm.tag)
    colnames(ret) <- xx$ProbeID
    ret
}

ctrl.red <- .take.mat(ctrl.raw.tab, 'Signal_Red')
ctrl.green <- .take.mat(ctrl.raw.tab, 'Signal_Grn')
ctrl.p <- .take.mat(ctrl.raw.tab, 'Detection Pval')

ctrl.logit <- (log(ctrl.red) - log(ctrl.green)) %>% as.matrix()
ctrl.logit[ ctrl.p > 0.05 ] <- NA

ctrl.type <- ctrl.raw.tab[, c(3,2)]

ctrl.melt <- ctrl.logit %>%
    melt(varname = c('meth.id', 'ProbeID'), stringsAsFactors = FALSE) %>%
        left_join(ctrl.type, by = 'ProbeID') %>%
            group_by(meth.id, TargetID) %>%
                summarize(ctrl.type = mean(value, na.rm = TRUE),
                          ctrl.type.se = sd(value, na.rm = TRUE)) %>%
                              as.data.frame()

uniq.target <- unique(ctrl.melt$TargetID)
uniq.meth.id <- unique(ctrl.melt$meth.id)

ret <- matrix(NA, nrow = length(uniq.meth.id), ncol = length(uniq.target))
ret[cbind(match(ctrl.melt$meth.id, uniq.meth.id), match(ctrl.melt$TargetID, uniq.target))] <-
    ctrl.melt$ctrl.type
colnames(ret) <- uniq.target
ret <- data.frame(meth.id = as.character(uniq.meth.id), ret, check.names = FALSE)
ret <- samples.tab %>% left_join(ret, by = 'meth.id') %>%
    as.data.frame()
    
write_feather(ret, path = out.file)
