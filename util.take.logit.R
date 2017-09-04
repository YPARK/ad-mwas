#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 4) q()

probe.file <- argv[1]
data.file <- argv[2]
sample.file <- argv[3]
out.file <- argv[4]

library(readr)
library(feather)

log.msg <- function(...) {
    cat(sprintf(...), file = stderr())
}

take.logit <- function(y) {
    ret <- log(y) - log(1 - y)
}

logit.std <- function(x, lb = -4, ub = 4) {
    ret <- scale(take.logit(x))
    ret <- pmax(pmin(ret, ub), lb)
}


probes <- read.table(probe.file, sep='\t', col.names=c('cg', 'chr', 'loc', 'meth.pos'))
samples <- read.table(sample.file, sep='\t', col.names=c('subj.id', 'meth.pos', 'data.pos'))
in.data <- scan(data.file, sep = '\t')
stopifnot((nrow(probes) * nrow(samples)) == length(in.data))

Y <- matrix(data = in.data, ncol = nrow(probes), byrow = TRUE)
rm(in.data)

out <- logit.std(Y, lb = -4, ub = 4)
write_feather(as.data.frame(out), path = out.file)
log.msg('Read Y: %d x %d\n', nrow(Y), ncol(Y))
