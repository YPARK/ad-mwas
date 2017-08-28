#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 2) q()

in.file <- argv[1]
out.file <- argv[2]

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

Y <- read_feather(in.file)
out <- logit.std(Y, lb = -4, ub = 4)
write_feather(as.data.frame(out), path = out.file)
log.msg('Read Y: %d x %d\n', nrow(Y), ncol(Y))

