#!/usr/bin/env Rscript
## Estimate QTL effects
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 2) q()

data.hdr <- argv[1] # e.g., data.hdr = '/broad/hptmp/ypp/AD/mwas/qtl/19/2-data'
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
.read.tab <- function(...) read_feather(...) %>% as.data.frame()
.read.mat <- function(...) read_feather(...) %>% as.matrix()
.cor <- function(...) cor(..., use = 'pairwise.complete.obs')

Y1 <- .read.mat(y1.data.file)
Y0 <- .read.mat(y0.data.file)
X <- .read.mat(x.data.file)
probes <- .read.tab(probe.data.file) %>% select(-meth.pos)

V <- .read.tab(cov.data.file)
C <- .read.mat(ctrl.data.file)

pheno <- V %>%
    select(amyloid_sqrt, tangles_sqrt, cogn_ep_random_slope, msex.y, studyn, age_death, educ) %>%
        scale() %>% as.matrix()

pheno <- cbind(pheno, 1)

PC <- read.table(PC.file) %>% as.matrix() %r% V$meth.pos %c% 1:30 %>%
    scale()

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

##############
## method 1 ##
##############

## 1. W = Y0 - X * theta
## 2. Y1 = f(W + R)
estimate.confounder <- function(yy, yy.ctrl, xx,
                                .pheno = pheno,
                                .probes = probes,
                                .iterations = 5000,
                                clean.confounder = FALSE) {

    vb.opt <- list(vbiter = ceiling(.iterations/2),
                   gammax = 1e4,
                   out.residual = TRUE,
                   tol = 1e-8,
                   rate = 5e-3,
                   pi.ub = -0,
                   pi.lb = -4,
                   model = 'gaussian',
                   k = 10)
    
    if(clean.confounder) {
        out.W <- fqtl.regress(y = yy.ctrl,
                              x.mean = xx,
                              c.mean = .pheno,
                              factored = FALSE,
                              options = vb.opt)
        
        W.proxy <- out.W$resid$theta %>% scale()
    } else {
        W.proxy <- yy.ctrl %>% scale()
    }

    ## Remove Y1 ~ W + .Pheno + 1    
    out.1 <- fqtl.regress(y = yy,
                          x.mean = W.proxy, 
                          c.mean = .pheno,
                          factored = TRUE,
                          options = vb.opt)

    residual <- out.1$resid$theta
    residual[is.na(yy)] <- NA
    colnames(residual) <- .probes[, 'cg']
    pve <- cbind(.probes, take.pve(out.1, W.proxy, pheno))
    list(R = residual, W = W.proxy, PVE = pve,
         mean.left = out.1$mean.left,
         mean.right = out.1$mean.right)
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

conf.1 <- estimate.confounder(Y1, Y0, X, clean.confounder = TRUE)
conf.2 <- estimate.confounder(Y1, Y0, X, clean.confounder = FALSE)
conf.3 <- estimate.confounder(Y1, C, X, clean.confounder = FALSE)
conf.4 <- estimate.confounder(Y1, PC, X, clean.confounder = FALSE)

write.confounder(conf.1, out.hdr %&&% '.conf-y0-clean')
write.confounder(conf.2, out.hdr %&&% '.conf-y0')
write.confounder(conf.3, out.hdr %&&% '.conf-ctrl')
write.confounder(conf.4, out.hdr %&&% '.conf-pc')

qtl.0 <- get.marginal.qtl(X, Y1)
qtl.1 <- get.marginal.qtl(X, conf.1$R)
qtl.2 <- get.marginal.qtl(X, conf.2$R)
qtl.3 <- get.marginal.qtl(X, conf.3$R)
qtl.4 <- get.marginal.qtl(X, conf.4$R)

x.perm <- X %r% sample(nrow(X))

qtl.perm.0 <- get.marginal.qtl(x.perm, Y1)
qtl.perm.1 <- get.marginal.qtl(x.perm, conf.1$R)
qtl.perm.2 <- get.marginal.qtl(x.perm, conf.2$R)
qtl.perm.3 <- get.marginal.qtl(x.perm, conf.3$R)
qtl.perm.4 <- get.marginal.qtl(x.perm, conf.4$R)

write.tab.gz(qtl.0, out.hdr %&&% '.qtl-raw.gz')
write.tab.gz(qtl.1, out.hdr %&&% '.qtl-y0-clean.gz')
write.tab.gz(qtl.2, out.hdr %&&% '.qtl-y0.gz')
write.tab.gz(qtl.3, out.hdr %&&% '.qtl-ctrl.gz')
write.tab.gz(qtl.4, out.hdr %&&% '.qtl-pc.gz')

write.tab.gz(qtl.perm.0, out.hdr %&&% '.qtl.perm-raw.gz')
write.tab.gz(qtl.perm.1, out.hdr %&&% '.qtl.perm-y0-clean.gz')
write.tab.gz(qtl.perm.2, out.hdr %&&% '.qtl.perm-y0.gz')
write.tab.gz(qtl.perm.3, out.hdr %&&% '.qtl.perm-ctrl.gz')
write.tab.gz(qtl.perm.4, out.hdr %&&% '.qtl.perm-pc.gz')

log.msg('Finished\n')
