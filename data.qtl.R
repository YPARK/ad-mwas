#!/usr/bin/env Rscript
## Generate QTL data
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 4) {
    q()
}

chr <- as.integer(argv[1])              # e.g., 21
meth.pos <- eval(parse(text = argv[2])) # e.g., 1:50
n.top.y0 <- as.integer(argv[3])         # e.g., 2
out.hdr <- argv[4]                      # e.g., 'temp'

cis.dist <- 1e6

library(fqtl)
library(feather)
library(dplyr)
source('util.R')
options(stringsAsFactors = FALSE)

dir.create(dirname(out.hdr), recursive = TRUE, showWarnings = FALSE)

################################################################
read.meth.chr <- function(chr, ...) {
    in.file <- 'data/meth/chr' %&&% chr %&&% '-logit.ft'
    ret <- read_feather(in.file, ...)
}

y1.out.file <- out.hdr %&&% '.y1.ft'
y0.out.file <- out.hdr %&&% '.y0.ft'
x.out.file <- out.hdr %&&% '.x.ft'
x.bim.out.file <- out.hdr %&&% '.x.bim.ft'
probe.out.file <- out.hdr %&&% '.y.prb.ft'

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
Y1 <- Y1 %r% sample.info$meth.pos %>% as.data.frame()
Y0 <- Y0.ref %r% sample.info$meth.pos %>% as.data.frame()

probes <- read.table('data/raw/chr' %&&% chr %&&% '-probes.txt.gz') %r% meth.pos

geno.bim <- read_feather('data/geno/chr' %&&% chr %&&% '.geno.bim.ft')

colnames(geno.bim) <- c('chr', 'rs', '.', 'snp.pos', 'a1', 'a2')

## Take SNPs in cis-regulatory region
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

colnames(X) <- x.bim[geno.rows, 'rs']
colnames(Y1) <- probes[, 'cg']

################################################################
## Write them down
write_feather(Y1, path = y1.out.file)
write_feather(Y0, path = y0.out.file)
write_feather(X, path = x.out.file)
write_feather(x.bim, path = x.bim.out.file)
write_feather(probes, path = probe.out.file)
