#!/usr/bin/env Rscript

## Generate QTL data

library(fqtl)
library(feather)
library(dplyr)
options(stringsAsFactors = FALSE)

meth.pos <- 1:50
chr <- 21
n.top.y0 <- 2

log.msg <- function(...) {
    cat(sprintf(...), file = stderr())
}

`%&&%` <- function(a, b) paste(a, b, sep = '')

`%c%` <- function(a, b) a[, b, drop = FALSE]

`%r%` <- function(a, b) a[b, , drop = FALSE]

read.meth.chr <- function(chr, ...) {
    in.file <- 'data/meth/chr' %&&% chr %&&% '-logit.ft'
    ret <- read_feather(in.file, ...)
}

.NA <- function(nrow, ncol) {
    matrix(NA, nrow, ncol)
}

fast.cov <- function(x, y) {
    n.obs <- crossprod(!is.na(x), !is.na(y))
    ret <- crossprod(replace(x, is.na(x), 0),
                     replace(y, is.na(y), 0)) / n.obs
}

fast.cor <- function(x, y) {
    x.sd <- apply(x, 2, sd, na.rm = TRUE)
    y.sd <- apply(y, 2, sd, na.rm = TRUE)
    ret <- fast.cov(x, y)
    ret <- sweep(sweep(ret, 1, x.sd, `/`), 2, y.sd, `/`)    
}

################################################################
## Find correlated CpGs in other chromosomes

Y1 <- as.matrix(read.meth.chr(chr, columns = 'X' %&&% meth.pos))

sample.info <- read_feather('data/geno/chr' %&&% chr %&&% '.samples.ft')

take.y0 <- function(.chr, y1, n.top) {
    y0 <- as.matrix(read.meth.chr(.chr))

    log.msg('correlation between y1 and y0 in chr%d\n', .chr)

    abs.cov.mat <- t(abs(fast.cov(y1, y0)))
    y0.idx <- apply(abs.cov.mat, 2, function(x) order(x, decreasing = TRUE)[1:n.top])
    y0.idx <- unique(as.vector(y0.idx))

    ret <- y0[, y0.idx, drop = FALSE]
    rm(y0); rm(abs.cov.mat); rm(y0)
    log.msg('constructed y0 in chr%d\n', .chr)

    return(ret)
}

other.chr <- setdiff(1:22, chr)
y0.list <- lapply(other.chr, take.y0, y1 = Y1, n.top = n.top.y0)
Y0 <- do.call(cbind, y0.list)
gc()

## further refine Y0 -> n x (n.top * n.cpg) at most
abs.cov.y10 <- t(abs(fast.cov(Y1, Y0)))
y0.idx <- apply(abs.cov.y10, 2, function(x) order(x, decreasing=TRUE)[1:n.top.y0])
Y0.ref <- do.call(cbind, lapply(1:ncol(y0.idx), function(j) Y0 %c% y0.idx[, j]))

################################################################
Y1 <- as.matrix(Y1 %r% sample.info$meth.pos)
Y0 <- as.matrix(Y0.ref %r% sample.info$meth.pos)

probes <- read.table('data/raw/chr' %&&% chr %&&% '-probes.txt.gz') %r% meth.pos

geno.bim <- read_feather('data/geno/chr' %&&% chr %&&% '.geno.bim.ft')

colnames(geno.bim) <- c('chr', 'rs', '.', 'snp.pos', 'a1', 'a2')

## Take SNPs in cis-regulatory region
cis.dist <- 1e6

colnames(probes) <- c('cg', 'chr', 'cg.pos', '.')

snp.idx.str <- function(cg.pos) {
    ret <- which(geno.bim$snp.pos > cg.pos - cis.dist &
                     geno.bim$snp.pos < cg.pos + cis.dist)
    paste(ret, collapse = '|')
}


temp <- probes %>% group_by(cg) %>% mutate(geno.rows = snp.idx.str(cg.pos)) %>%
    mutate(geno.rows = sapply(geno.rows, strsplit, split = '[|]'))

geno.rows <- as.integer(unique(unlist(temp$geno.rows)))

x.file <- 'data/geno/chr' %&&% chr %&&% '.geno.mat.ft'
X <- read_feather(x.file, 'V' %&&% geno.rows) %r% 
    sample.info$geno.pos %>%
        as.data.frame()

x.bim <- geno.bim %r% geno.rows %>%
    as.data.frame()

