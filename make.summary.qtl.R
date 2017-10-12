#!/usr/bin/env Rscript
## Estimate QTL effects
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 2) q()

data.hdr <- argv[1] # e.g., data.hdr = '/broad/hptmp/ypp/AD/mwas/qtl/1/4-data'
out.hdr <- argv[2]  # e.g., out.hdr = 'temp'

dir.create(dirname(out.hdr), recursive = TRUE, showWarnings = FALSE)

################################################################
source('util.R')

y1.data.file <- data.hdr %&&% '.y1.ft'
y0.data.file <- data.hdr %&&% '.y0.ft'
x.data.file <- data.hdr %&&% '.x.ft'
x.bim.data.file <- data.hdr %&&% '.x.bim.ft'
probe.data.file <- data.hdr %&&% '.y.prb.ft'
cov.data.file <- data.hdr %&&% '.cov.ft'
ctrl.data.file <- data.hdr %&&% '.ctrl.ft'

PC.file <- 'data/meth/PC.txt.gz'

in.files <- c(y1.data.file,
              y0.data.file,
              x.data.file,
              x.bim.data.file,
              probe.data.file,
              cov.data.file,
              ctrl.data.file,
              PC.file)

if(!all(sapply(in.files, file.exists))) {
    log.msg('Insufficient input files:\n%s\n',
            paste(in.files, collapse='\n'))
    q()
}

################################################################
library(dplyr)
library(feather)
library(fqtl)
library(glmnet)
library(methods)
.read.tab <- function(...) read_feather(...) %>% as.data.frame()
.read.mat <- function(...) read_feather(...) %>% as.matrix()
.cor <- function(...) cor(..., use = 'pairwise.complete.obs')

Y1 <- .read.mat(y1.data.file)
Y0 <- .read.mat(y0.data.file)
X <- .read.mat(x.data.file)
x.bim <- .read.tab(x.bim.data.file)
colnames(x.bim) <- c('chr', 'rs', '.', 'snp.loc', 'qtl.a1', 'qtl.a2')

probes <- .read.tab(probe.data.file) %>% select(-cg.pos)
cpg.names <- probes$cg

V <- .read.tab(cov.data.file)
C <- .read.mat(ctrl.data.file)

pheno <- V %>%
    select(amyloid_sqrt, tangles_sqrt, cogn_ep_random_slope, msex.y, studyn, age_death, educ) %>%
        scale() %>% as.matrix()

pheno <- cbind(pheno, 1)

PC <- read.table(PC.file) %>% as.matrix() %r% V$meth.pos %c% 1:30 %>%
    scale()

################################################################
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
    theta.cov <- qtl.out$mean.cov$theta
    resid <- qtl.out$resid$theta

    xx[is.na(xx)] <- 0
    theta[is.na(theta)] <- 0
    theta.cov[is.na(theta.cov)] <- 0
    cc[is.na(cc)] <- 0
    eta.1 <- xx %*% theta
    eta.2 <- cc %*% theta.cov

    data.frame(v1 = apply(eta.1, 2, var, na.rm = TRUE),
               v2 = apply(eta.2, 2, var, na.rm = TRUE),
               vr = apply(resid, 2, var, na.rm = TRUE))
}

run.glmnet <- function(y, x, alpha = 1){
    valid <- !is.na(y)
    xx <- x[valid,,drop=FALSE]
    yy <- as.matrix(y[valid])
    cv.out <- cv.glmnet(x=xx, y=yy, alpha=alpha, nfolds=5)
    bet <- glmnet(x=xx, y=yy, alpha=alpha, lambda=cv.out$lambda.min)$beta
    cat(mean(abs(bet) > 0), '\n')

    resid <- matrix(NA, nrow = length(y), ncol = 1)
    resid[valid,] <- as.matrix(yy - xx %*% bet)
    return(list(beta = as.matrix(bet), resid = resid))
}

estimate.lm <- function(yy, yy.ctrl,
                        out.tag = '',
                        .pheno = pheno,
                        .cpg.names = cpg.names) {

    .pheno[is.na(.pheno)] <- 0
    .nc <- ncol(yy.ctrl)
    .xx <- cbind(yy.ctrl, .pheno)
    .xx[is.na(.xx)] <- 0

    glmnet.list <- apply(yy, 2, run.glmnet, x = .xx, alpha = 0.5)

    out <- list()
    out$resid$theta <- do.call(cbind, lapply(glmnet.list, function(x) x$resid))
    out$mean$theta <- do.call(cbind, lapply(glmnet.list, function(x) x$beta %r% 1:.nc))
    out$mean.cov$theta <- do.call(cbind, lapply(glmnet.list, function(x) x$beta %r% (-(1:.nc))))

    pve <- cbind(.cpg.names, take.pve(out, as.matrix(yy.ctrl), as.matrix(.pheno)))
    colnames(out$resid$theta) <- .cpg.names

    list(R = out$resid$theta, PVE = pve, out.tag = out.tag, lm = glmnet.list)
}

write.confounder <- function(.obj, .out.hdr) {
    resid.out.file <- .out.hdr %&&% '.resid.gz'
    pve.out.file <- .out.hdr %&&% '.pve.gz'
    write.tab(.obj$R, file = gzfile(resid.out.file))
    write.tab(.obj$PVE, file = gzfile(pve.out.file))
}

write.tab.gz <- function(.tab, .out.file) {
    write.tab(.tab, file = gzfile(.out.file))
}

filter.qtl <- function(qtl.tab, cis.dist = 1e6, .probes = probes, .snps = x.bim) {
    qtl.tab %>% mutate(cg = as.character(cg), rs = as.character(rs)) %>%
        left_join(.probes, by = 'cg') %>%
            left_join(.snps, by = 'rs') %>%
                filter(abs(snp.loc - loc) < cis.dist) %>%
                    select(snp.loc, cg, beta, beta.z)
}

conf.1 <- estimate.lm(Y1, Y0, out.tag = 'hs-lm') 

conf.2 <- estimate.lm(Y1, PC, out.tag = 'pc-lm') 

conf.list <- list(conf.1, conf.2)

sapply(conf.list, function(cc) write.confounder(cc, out.hdr %&&% '-' %&&% cc$out.tag))

qtl.y1 <- get.marginal.qtl(X, Y1) %>% filter.qtl()
write.tab.gz(qtl.y1, out.hdr %&&% '.qtl-raw-y1.gz')

qtl.y0 <- get.marginal.qtl(X, Y0) %>% filter.qtl()
write.tab.gz(qtl.y0, out.hdr %&&% '.qtl-raw-y0.gz')

sapply(conf.list, function(cc) {
    write.tab.gz(get.marginal.qtl(X, cc$R) %>% filter.qtl(),
                 out.hdr %&&% '.qtl-' %&&% cc$out.tag %&&% '.gz')
})

write.tab.gz(x.bim, out.hdr %&&% '.snps.gz')
write.tab.gz(probes, out.hdr %&&% '.probes.gz')

log.msg('Finished\n')
