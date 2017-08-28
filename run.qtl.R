#!/usr/bin/env Rscript
## Estimate QTL effects
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 4) q()

data.hdr <- argv[1] # e.g., data.hdr = '/broad/hptmp/ypp/AD/mwas/1/1411-data'
sample.file <- argv[2] # e.g., sample.file = 'data/geno/chr1.samples.ft'
ctrl.probe.file <- argv[3] # e.g., ctrl.probe.file = 'data/meth/control.ft'
out.hdr <- argv[4]  # e.g., out.hdr = 'temp'

dir.create(dirname(out.hdr), recursive = TRUE, showWarnings = FALSE)

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
              probe.data.file,
              sample.file,
              ctrl.probe.file)

if(!all(sapply(in.files, file.exists))) {
    log.msg('Insufficient input files\n')
    q()
}

resid0.out.file <- out.hdr %&&% '.resid0.txt.gz'
resid.out.file <- out.hdr %&&% '.resid.txt.gz'
resid.ctrl.out.file <- out.hdr %&&% '.resid-ctrl.txt.gz'
probe.out.file <- out.hdr %&&% '.probes.txt.gz'

qtl.resid0.out.file <- out.hdr %&&% '.qtl-resid0.gz'
qtl.resid.out.file <- out.hdr %&&% '.qtl-resid.gz'
qtl.raw.out.file <- out.hdr %&&% '.qtl-raw.gz'
qtl.ctrl.out.file <- out.hdr %&&% '.qtl-ctrl.gz'

perm.resid0.out.file <- out.hdr %&&% '.qtl-perm-resid0.gz'
perm.resid.out.file <- out.hdr %&&% '.qtl-perm-resid.gz'
perm.raw.out.file <- out.hdr %&&% '.qtl-perm-raw.gz'
perm.ctrl.out.file <- out.hdr %&&% '.qtl-perm-ctrl.gz'

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

y1 <- read_feather(y1.data.file) %>% as.matrix()
y0 <- read_feather(y0.data.file) %>% as.matrix()
X <- read_feather(x.data.file) %>% as.matrix() %>% scale(scale = FALSE)

x.bim <- read_feather(x.bim.data.file) %>% as.data.frame()
probes <- read_feather(probe.data.file) %>% as.data.frame()
colnames(X) <- x.bim$rs

log.msg('Y1 : %d x %d\n\n', nrow(y1), ncol(y1))
log.msg('Y0 : %d x %d\n\n', nrow(y0), ncol(y0))
log.msg('X : %d x %d\n\n', nrow(X), ncol(X))

take.theta <- function(qtl.out) {
    if('mean.left' %in% names(qtl.out)){
        theta <- qtl.out$mean.left$theta %*% t(qtl.out$mean.right$theta)
    } else {
        theta <- qtl.out$mean$theta
    }
    return(theta)
}

take.pve <- function(qtl.out, xx, cc) {

    theta <- take.theta(qtl.out)
    resid <- qtl.out$resid$theta

    eta.1 <- xx %*% theta
    eta.2 <- cc %*% qtl.out$mean.cov$theta

    eta.tot <- resid + eta.1 + eta.2

    data.frame(v1 = apply(eta.1, 2, var, na.rm = TRUE),
               v2 = apply(eta.2, 2, var, na.rm = TRUE),
               vr = apply(resid, 2, var, na.rm = TRUE))
}

################################################################
samples <- read_feather(sample.file)
ctrl.tab <- read_feather(ctrl.probe.file)

ctrl.tab <- samples %>% select(iid) %>%
    left_join(ctrl.tab, by = 'iid') %>%
        select(-iid, -meth.id, -meth.pos, -has.geno, -has.meth) %>%
            as.matrix()

ctrl.tab <- scale(ctrl.tab, scale = FALSE)

pheno.tab <- samples %>%
    select(np_sqrt, nft_sqrt, cogn_ep_random_slope, msex, age_death, educ) %>%
        mutate(intercept = 1) %>%
            as.matrix()

n.pheno <- ncol(pheno.tab)
pheno.tab[, -n.pheno] <- scale(pheno.tab[, -n.pheno])

################################################################
## Remove Y0 ~ X
vb.opt <- list(vbiter = 1500, gammax = 1e3, out.residual = TRUE,
               tol = 1e-8, rate = 1e-2, pi.ub = -0, pi.lb = -4)

out.W <- fqtl.regress(y = y0, x.mean = X,
                      factored = FALSE,
                      options = vb.opt)
                            
W.proxy <- out.W$resid$theta
rm(out.W)

################################################################
## Remove Y1 ~ Y0 + Pheno + 1
vb.opt <- list(vbiter = 2000, gammax = 1e3, k = 10, out.residual = TRUE,
               tol = 1e-8, rate = 1e-2, pi.ub = -0, pi.lb = -4)

out.0 <- fqtl.regress(y = y1,
                      x.mean = y0,
                      c.mean = pheno.tab,
                      x.var = pheno.tab,
                      factored = TRUE,
                      options = vb.opt)

## Save residuals
residual0 <- out.0$resid$theta # i.e., y1 - y0 %*% take.theta(out.0)
residual0[is.na(y1)] <- NA

write.mat(residual0 %>% as.data.frame(), file = gzfile(resid0.out.file))

################################################################
## Remove Y1 ~ W + Pheno + 1
out.1 <- fqtl.regress(y = y1,
                      x.mean = W.proxy, 
                      c.mean = pheno.tab,
                      x.var = pheno.tab,
                      factored = TRUE,
                      options = vb.opt)

## Save residuals
residual <- out.1$resid$theta # i.e., y1 - W.proxy %*% take.theta(out.1)
residual[is.na(y1)] <- NA

write.mat(residual %>% as.data.frame(), file = gzfile(resid.out.file))

################################################################
## Correction by known control probes
out.2 <- fqtl.regress(y = y1,
                      x.mean = ctrl.tab,
                      c.mean = pheno.tab,
                      x.var = pheno.tab,
                      factored = FALSE,
                      options = vb.opt)

residual.ctrl <- out.2$resid$theta #
residual.ctrl[is.na(y1)] <- NA

write.mat(residual.ctrl %>% as.data.frame(), file = gzfile(resid.ctrl.out.file))

log.msg('Finished confounder correction\n\n')

################################################################
## calculate PVE inside the link function

pve.0 <- take.pve(out.0, y0, pheno.tab)
pve.1 <- take.pve(out.1, W.proxy, pheno.tab)
pve.2 <- take.pve(out.2, ctrl.tab, pheno.tab)

out.probes <- rbind(data.frame(probes[, -4], pve.0, model = 'proxy0'),
                    data.frame(probes[, -4], pve.1, model = 'proxy'),
                    data.frame(probes[, -4], pve.2, model = 'ctrl'))

write.tab(out.probes, file = gzfile(probe.out.file))

log.msg('Finished PVE calculation\n\n')

################################################################
colnames(y1) <- probes$cg
colnames(residual0) <- probes$cg
colnames(residual) <- probes$cg
colnames(residual.ctrl) <- probes$cg

################################################################
## marginal QTL calculation
qtl.raw.out <- get.marginal.qtl(X, y1)
write.tab(qtl.raw.out, file = gzfile(qtl.raw.out.file))
rm(qtl.raw.out)

qtl.resid.out <- get.marginal.qtl(X, residual)
write.tab(qtl.resid.out, file = gzfile(qtl.resid.out.file))
rm(qtl.resid.out)

qtl.resid0.out <- get.marginal.qtl(X, residual0)
write.tab(qtl.resid0.out, file = gzfile(qtl.resid0.out.file))
rm(qtl.resid0.out)

qtl.ctrl.out <- get.marginal.qtl(X, residual.ctrl)
write.tab(qtl.ctrl.out, file = gzfile(qtl.ctrl.out.file))
rm(qtl.ctrl.out)

## permuted QTL calculation
x.perm <- X %r% sample(nrow(X))
colnames(x.perm) <- x.bim$rs

qtl.raw.out <- get.marginal.qtl(x.perm, y1)
write.tab(qtl.raw.out, file = gzfile(perm.raw.out.file))
rm(qtl.raw.out)

qtl.resid.out <- get.marginal.qtl(x.perm, residual)
write.tab(qtl.resid.out, file = gzfile(perm.resid.out.file))
rm(qtl.resid.out)

qtl.resid0.out <- get.marginal.qtl(x.perm, residual0)
write.tab(qtl.resid0.out, file = gzfile(perm.resid0.out.file))
rm(qtl.resid0.out)

qtl.ctrl.out <- get.marginal.qtl(x.perm, residual.ctrl)
write.tab(qtl.ctrl.out, file = gzfile(perm.ctrl.out.file))
rm(qtl.ctrl.out)

log.msg('Finished\n')
