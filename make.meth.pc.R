#!/usr/bin/env Rscript

source('util.R')
library(feather)
library(dplyr)

sample.file <- 'data/raw/matched.samples.txt.gz'
raw.data.files <- 'data/raw/chr' %&&% 1:22 %&&% '-beta.txt.gz'
data.files <- 'data/meth/chr' %&&% 1:22 %&&% '-logit.ft'

samples <- read.table(sample.file, sep='\t',
                      col.names=c('iid', 'meth.id', 'meth.pos', 'has.geno', 'has.meth'))

take.sd <- function(ff) {
    in.data <- scan(ff, sep = '\t')
    n.r <- nrow(samples)
    n.c <- length(in.data) / n.r
    ret <- apply(matrix(data = in.data, nrow = n.r, ncol = n.c, byrow = TRUE),
                 2, sd, na.rm = TRUE)
    rm(in.data)
    log.msg('Read %s\n', ff)
    return(ret)
}

Y.sd <- lapply(raw.data.files, take.sd)
gc()

cutoff <- sort(do.call(c, Y.sd), decreasing = TRUE)[1800]

Y.probes <- lapply(Y.sd, function(yy) which(yy > cutoff))
Y.probes <- lapply(Y.probes, function(yy) 'V' %&&% yy)
names(Y.probes) <- data.files


Y <- lapply(data.files, function(ff) read_feather(ff, columns=Y.probes[[ff]]))
Y <- do.call(cbind, Y) %>% scale() %>% as.matrix()
Y <- Y %c% (apply(is.na(Y), 2, sum) == 0)

Y.svd <- svd(Y)

write.tab(Y.svd$u, file = gzfile('data/meth/PC.txt.gz'))
