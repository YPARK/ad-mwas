#!/usr/bin/env Rscript
## Estimate QTL effects
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 2) q()

data.hdr <- argv[1] # e.g., '/broad/hptmp/ypp/AD/mwas/1/141-data'
out.hdr <- argv[2]  # e.g., 'temp'

dir.create(dirname(out.hdr), recursive = TRUE, showWarnings = FALSE)

do.permutation <- FALSE
if(length(argv) > 2) {
    do.permutation <- as.logical(argv[3])
}
    
################################################################
source('util.R')

y1.data.file <- data.hdr %&&% '.y1.ft'
y0.data.file <- data.hdr %&&% '.y0.ft'
x.data.file <- data.hdr %&&% '.x.ft'
x.bim.data.file <- data.hdr %&&% '.x.bim.ft'
probe.data.file <- data.hdr %&&% '.y.prb.ft'

in.files <- c(y1.data.file,
              y0.data.file,
              x.data.file,
              x.bim.data.file,
              probe.data.file)

if(!all(sapply(in.files, file.exists))) {
    log.msg('Insufficient input files\n')
    q()
}

resid.out.file <- out.hdr %&&% '.resid.txt.gz'
probe.out.file <- out.hdr %&&% '.probes.txt.gz'
qtl.resid.out.file <- out.hdr %&&% '.resid.qtl.gz'
qtl.raw.out.file <- out.hdr %&&% '.raw.qtl.gz'
snps.out.file <- out.hdr %&&% '.snps.gz'

out.files <- c(resid.out.file,
               probe.out.file,
               qtl.resid.out.file,
               qtl.raw.out.file,
               snps.out.file)
               
if(all(sapply(out.files, file.exists))) {
    log.msg('All the files exist\n')
    q()
}

################################################################
library(feather)
library(fqtl)
library(dplyr)
source('util.matrixqtl.R')

y1 <- read_feather(y1.data.file) %>% as.matrix()
y0 <- read_feather(y0.data.file) %>% as.matrix()
X <- read_feather(x.data.file) %>% as.matrix()
x.bim <- read_feather(x.bim.data.file) %>% as.data.frame()
probes <- read_feather(probe.data.file) %>% as.data.frame()

if(do.permutation) {
    X <- X %r% sample(nrow(X))
}

################################################################
## Remove Y0 ~ X
vb.opt <- list(vbiter = 2000, gammax = 1e3,
               tol = 1e-8, rate = 1e-2, pi.ub = -2, pi.lb = -4)

out.0 <- fqtl.regress(y = y0, x.mean = X, factored = FALSE, options = vb.opt)

W.proxy <- scale(out.0$resid$theta)

## Remove Y1 ~ W using Beta model
y1.beta <- 1/(1 + exp(-y1))

vb.opt <- list(vbiter = 3000, gammax = 1e3, k = 5, model = 'beta',
               tol = 1e-8, rate = 1e-2, pi.ub = -2, pi.lb = -4)

out.1 <- fqtl.regress(y = y1.beta, x.mean = W.proxy, factored = TRUE, options = vb.opt)

## Save residuals
residual <- out.1$resid$theta %>% as.data.frame()
write.mat(residual, file = gzfile(resid.out.file))

## calculate PVE inside the link function
y1.hat <- W.proxy %*% out.1$mean.left$theta %*% t(out.1$mean.right$theta)
y1.tot <- y1.hat + residual

v.tot <- apply(y1.tot, 2, var, na.rm = TRUE)
v.hat <- apply(y1.hat, 2, var, na.rm = TRUE)
out.probes <- data.frame(probes[, -4], pve = round(data.frame(pve = v.hat / v.tot), 2))

write.tab(out.probes, file = gzfile(probe.out.file))

################################################################
## marginal QTL calculation

rm.gz <- function(x) gsub(pattern='.gz', replacement = '', x)

me.out <- run.matrix.qtl(X, y1,
                         snp.names = x.bim$rs,
                         out.names = probes$cg,
                         me.out.file = rm.gz(qtl.raw.out.file))

system(paste('gzip', rm.gz(qtl.raw.out.file), '-9 -f'))

me.out <- run.matrix.qtl(X, residual,
                         snp.names = x.bim$rs,
                         out.names = probes$cg,
                         me.out.file = rm.gz(qtl.resid.out.file))

system(paste('gzip', rm.gz(qtl.resid.out.file), '-9 -f'))

write.tab(x.bim, file = gzfile(snps.out.file))
