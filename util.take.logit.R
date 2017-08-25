#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 2) q()

in.file <- argv[1]
out.file <- argv[2]

library(feather)

log.msg <- function(...) {
    cat(sprintf(...), file = stderr())
}

take.logit <- function(y, TOL = 1e-4) {
    logit.y <- y
    logit.y[y < TOL] <- TOL
    logit.y[y > (1 - TOL)] <- 1 - TOL
    logit.y <- log(logit.y) - log(1 - logit.y)

    log.msg('finished logit transformation: %d x %d\n',
            nrow(y), ncol(y))

    return(logit.y)
}

`%qc%` <- function(X, q.cutoff) {

    ret.X <- apply(X, 2, function(x) {
        lb <- quantile(x, probs = q.cutoff * 0.5, na.rm = TRUE)
        ub <- quantile(x, probs = 1 - q.cutoff * 0.5, na.rm = TRUE)
        pmin(pmax(x, lb), ub)
    })

    log.msg('finished quantile-based truncation: %d x %d\n',
            nrow(X), ncol(X))
    return(ret.X)
}

logit.std <- function(x, q.cutoff = 0.01) {
    ret <- scale(take.logit(x) %qc% q.cutoff)
    log.msg('finished standardization: %d x %d\n',
            nrow(x), ncol(x))
    return(ret)
}

Y <- read_feather(in.file)
log.msg('Read Y: %d x %d\n', nrow(Y), ncol(Y))
out <- logit.std(Y, q.cutoff = 0.01)
write_feather(as.data.frame(out), path = out.file)
