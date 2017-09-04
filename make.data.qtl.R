#!/usr/bin/env Rscript
## Generate QTL data
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 6) {
    q()
}

meth.file <- argv[1]                     # e.g., meth.file = 'data/meth/chr2-logit.ft'
meth.probe.file <- argv[2]               # e.g., meth.probe.file = 'data/raw/chr2-probes.txt.gz'
meth.cols <- eval(parse(text = argv[3])) # e.g., meth.cols = 1:20
n.top.y0 <- as.integer(argv[4])          # e.g., n.top.y0 = 3
plink.hdr <- argv[5]                     # e.g., plink.hdr = 'rosmap-geno/gen/impute/rosmap1709-chr2'
out.hdr <- argv[6]                       # e.g., out.hdr = 'temp'

cis.dist <- 1e6

library(fqtl)
library(feather)
library(dplyr)
source('util.R')
options(stringsAsFactors = FALSE)

dir.create(dirname(out.hdr), recursive = TRUE, showWarnings = FALSE)
dir.create(out.hdr, recursive = TRUE)
temp.dir <- system('mktemp -d ' %&&% out.hdr %&&% '/temp.XXXX', intern = TRUE, ignore.stderr = TRUE)

sample.file <- 'data/raw/matched.samples.txt.gz'
pheno.file <- 'phenotype/pheno_cov_n3033_032315.csv'
control.probes <- read_feather('data/meth/control.ft')

################################################################
y1.out.file <- out.hdr %&&% '.y1.ft'
y0.out.file <- out.hdr %&&% '.y0.ft'
x.out.file <- out.hdr %&&% '.x.ft'
x.bim.out.file <- out.hdr %&&% '.x.bim.ft'
probe.out.file <- out.hdr %&&% '.y.prb.ft'
cov.out.file <- out.hdr %&&% '.cov.ft'
ctrl.out.file <- out.hdr %&&% '.ctrl.ft'

################################################################
## Find correlated CpGs in other chromosomes

probes <- read.table(meth.probe.file, sep='\t', col.names=c('cg', 'chr', 'loc', 'meth.pos'))
probes <- probes %r% meth.cols
chr <- unique(probes$chr)
samples <- read.table(sample.file, sep='\t', col.names=c('iid', 'meth.id', 'meth.pos', 'has.geno', 'has.meth'))

pheno.tab <- read.table(pheno.file, sep = ',', header = TRUE) %>% rename(fid = FID, iid = IID)
pheno.tab[pheno.tab == -9] <- NA
pheno.tab <- pheno.tab %>% select(-fid)

plink.lb <- max(min(probes$loc) - cis.dist, 0)
plink.ub <- max(probes$loc) + cis.dist

plink.cmd <- sprintf('./bin/plink --bfile %s --make-bed --geno 0.05 --chr %d --from-bp %d --to-bp %d --out %s',
                     plink.hdr, chr, plink.lb, plink.ub, temp.dir %&&% 'plink')
system(plink.cmd)
plink <- read.plink(temp.dir %&&% 'plink')
system('rm -r ' %&&% temp.dir)
colnames(plink$FAM) <- c('fid', 'iid', '.', '..', 'msex', '...')
fam.tab <- data.frame(plink$FAM %>% select(fid, iid, msex), geno.pos = 1:nrow(plink$FAM))

log.msg('probes:\n%s\n\n', paste(probes$cg, collapse=', '))

if(nrow(plink$BIM) < 1) {
    log.msg('No variants in cis\n\n')
    q()
}

Y1 <- as.matrix(read_feather(meth.file, columns = 'V' %&&% meth.cols))
colnames(Y1) <- probes$cg

take.y0 <- function(.chr, y1, n.top) {

    .meth.file <- gsub(meth.file, patter = 'chr' %&&% chr, replacement = 'chr' %&&% .chr)
    y0 <- as.matrix(read_feather(.meth.file))
    log.msg('correlation between y1 and y0 in chr%d\n', .chr)

    abs.cor.mat <- t(abs(fast.cor(y1, y0)))
    y0.idx <- apply(abs.cor.mat, 2, function(x) order(x, decreasing = TRUE)[1:n.top])
    y0.idx <- unique(as.vector(y0.idx))

    ret <- y0[, y0.idx, drop = FALSE]
    rm(y0); rm(abs.cor.mat); rm(y0)
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
## Take individuals with genotypes
sample.info <- samples %>%
    left_join(fam.tab, by = 'iid') %>%
        na.omit() %>%
            left_join(pheno.tab, by = 'iid')

Y1 <- Y1 %r% sample.info$meth.pos %>% as.data.frame()
Y0 <- Y0.ref %r% sample.info$meth.pos %>% as.data.frame()

x.bim <- plink$BIM
X <- plink$BED %r% sample.info$geno.pos %>% as.data.frame()
colnames(X) <- x.bim[, 2]

ctrl.out <- sample.info %>% select(iid) %>%
    left_join(control.probes, by = 'iid') %c% -(1:5) %>%
        scale() %>% as.data.frame()

################################################################
## Write them down
write_feather(sample.info, path = cov.out.file)
write_feather(ctrl.out, path = ctrl.out.file)
write_feather(Y1, path = y1.out.file)
write_feather(Y0, path = y0.out.file)
write_feather(X, path = x.out.file)
write_feather(x.bim, path = x.bim.out.file)
write_feather(probes, path = probe.out.file)
