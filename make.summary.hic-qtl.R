#!/usr/bin/env Rscript
## Estimate QTL effects
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 3) q()

data.hdr <- argv[1] # e.g., data.hdr = '/broad/hptmp/ypp/AD/mwas/hic-qtl/1/4-data'
hic.file <- argv[2] # e.g., hic.file = 'hic/CO/chr1.pairs.gz'
out.hdr <- argv[3]  # e.g., out.hdr = 'temp'

dir.create(dirname(out.hdr), recursive = TRUE, showWarnings = FALSE)

hic.resol <- 40000

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
library(readr)

################################################################
hic.tab <- read_tsv(hic.file, col_names = FALSE)
hic.pairs <- hic.tab[, c(2, 4)] / hic.resol
colnames(hic.pairs) <- c('u', 'v')
.pairs1 <- hic.pairs %>% unique() %>% mutate(pp = paste(u, v, sep='-'))
.pairs2 <- hic.pairs %>% unique() %>% mutate(pp = paste(v, u, sep='-'))
hic.pairs <- c(.pairs1$pp, .pairs2$pp)

################################################################
.read.tab <- function(...) read_feather(...) %>% as.data.frame()
.read.mat <- function(...) read_feather(...) %>% as.matrix()
.cor <- function(...) cor(..., use = 'pairwise.complete.obs')

Y1 <- .read.mat(y1.data.file)
Y0 <- .read.mat(y0.data.file)

x.bim <- .read.tab(x.bim.data.file)
colnames(x.bim) <- c('chr', 'rs', '.', 'snp.loc', 'qtl.a1', 'qtl.a2')
X <- .read.mat(x.data.file) %>% scale()
colnames(X) <- x.bim[, 4]

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

take.pve <- function(qtl.out, xx) {

    theta <- take.theta(qtl.out)
    resid <- qtl.out$resid$theta

    xx[is.na(xx)] <- 0
    theta[is.na(theta)] <- 0

    eta.1 <- xx %*% theta
    data.frame(v1 = apply(eta.1, 2, var, na.rm = TRUE),
               vr = apply(resid, 2, var, na.rm = TRUE))
}

run.glmnet <- function(y, x, alpha = 1){
    valid <- !is.na(y)
    xx <- x[valid,,drop=FALSE]
    yy <- as.matrix(y[valid])
    cv.out <- cv.glmnet(x=xx, y=yy, alpha=alpha, nfolds=5)
    bet <- glmnet(x=xx, y=yy, alpha=alpha, lambda=cv.out$lambda.min)$beta
    bet <- as.matrix(bet)

    cat(mean(abs(bet) > 0), '\n')

    if(sum(abs(bet) > 0) == 0) {
        return(list(beta = bet, resid = y))
    }

    ## re-estimate lm on non-zero coefficients
    resid <- matrix(NA, nrow = length(y), ncol = 1)
    xx.sub <- xx %c% which(abs(bet)>0)
    lm.out <- lm(yy ~ xx.sub - 1)

    resid[valid,] <- as.matrix(lm.out$residuals)

    return(list(beta = bet, resid = resid))
}

estimate.lm <- function(yy, yy.ctrl,
                        out.tag = '',
                        .cpg.names = cpg.names) {

    .xx <- yy.ctrl
    .xx[is.na(.xx)] <- 0

    glmnet.list <- apply(yy, 2, run.glmnet, x = .xx, alpha = 0.5)

    out <- list()
    out$resid$theta <- do.call(cbind, lapply(glmnet.list, function(x) x$resid))
    out$mean$theta <- do.call(cbind, lapply(glmnet.list, function(x) x$beta))

    pve <- cbind(.cpg.names, take.pve(out, as.matrix(yy.ctrl)))
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

################################################################
## Hi-C based filtering
hic.filter.qtl <- function(qtl.tab,
                           probes.tab = probes,
                           hic.pairs.tab = hic.pairs,
                           resol = hic.resol) {

    ret <- qtl.tab %>%
        mutate(cg = as.character(cg)) %>%
            left_join(probes.tab, by = 'cg') %>%
                mutate(u = ceiling(loc/resol), v = ceiling(snp/resol))
    
    ret.self <- ret %>% filter(u == v) %>%
        select(snp, cg, beta, beta.z)

    ret.remote <- ret %>%        
        mutate(pp = paste(u, v, sep = '-')) %>%
            filter(pp %in% hic.pairs.tab)
    
    ret.remote <- ret.remote %>%  select(snp, cg, beta, beta.z)

    ret <- rbind(ret.self, ret.remote) %>% arrange(snp, desc(abs(beta.z)))

    return(ret %>% as.data.frame())
}

conf.1 <- estimate.lm(Y1, Y0, out.tag = 'hs-lm', cpg.names) 
conf.2 <- estimate.lm(Y1, PC, out.tag = 'pc-lm', cpg.names) 
conf.list <- list(conf.1, conf.2)

sapply(conf.list, function(cc) write.confounder(cc, out.hdr %&&% '-' %&&% cc$out.tag))

qtl.y1 <- get.marginal.qtl(X, Y1) %>% hic.filter.qtl()
write.tab.gz(qtl.y1, out.hdr %&&% '.qtl-raw-y1.gz')

qtl.y0 <- get.marginal.qtl(X, Y0)
write.tab.gz(qtl.y0, out.hdr %&&% '.qtl-raw-y0.gz')

sapply(conf.list, function(cc) {
    write.tab.gz(get.marginal.qtl(X, cc$R) %>% hic.filter.qtl(),
                 out.hdr %&&% '.qtl-' %&&% cc$out.tag %&&% '.gz')
})

write.tab.gz(x.bim, out.hdr %&&% '.snps.gz')
write.tab.gz(probes, out.hdr %&&% '.probes.gz')

log.msg('Finished\n')
