#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 4) q()

plink.hdr <- argv[1]   # e.g., plink.hdr = 'genotype/step08/BED/chr21'
sample.file <- argv[2] # e.g., sample.file = 'data/raw/matched.samples.txt.gz'
pheno.file <- argv[3] # e.g., pheno.file = 'phenotype/pheno_cov_n3033_032315.csv'
out.hdr <- argv[4]

library(feather)
library(fqtl)
library(dplyr)
options(stringsAsFactors = FALSE)

pheno.tab <- read.table(pheno.file, sep = ',', header = TRUE) %>%
    rename(fid = FID, iid = IID)

pheno.tab[pheno.tab == -9] <- NA

pheno.tab <- pheno.tab %>% select(-fid)

plink <- read.plink(plink.hdr)
samples <- read.table(sample.file)
colnames(samples) <- c('iid', 'meth.id', 'meth.pos', 'has.geno', 'has.meth')

colnames(plink$FAM)[1:2] <- c('fid', 'iid')
fam.tab <- data.frame(plink$FAM[, 1:2], geno.pos = 1:NROW(plink$FAM))

meth.samples <- samples %>%
    left_join(fam.tab, by = 'iid') %>%
        na.omit() %>%
            left_join(pheno.tab, by = 'iid')

geno.mat <- plink$BED[meth.samples$geno.pos, ]
geno.fam <- plink$FAM[meth.samples$geno.pos, ]

## remove SNPs with too many missing values
## and take common variants MAF < 5%
maf <- 0.05

missing.rates <- apply(is.na(geno.mat), 2, mean)

mu <- apply(geno.mat, 2, mean, na.rm = TRUE) * 0.5
vv <- mu * (1 - mu)
var.cutoff <- maf * (1 - maf)

valid.snp.pos <- which(missing.rates < 0.05 & vv > var.cutoff)

geno.mat <- geno.mat[, valid.snp.pos, drop = FALSE]

geno.bim <- plink$BIM[valid.snp.pos, ]

glue <- function(...) paste(..., sep = '')

write_feather(geno.bim, path = glue(out.hdr, '.geno.bim.ft'))
write_feather(geno.fam, path = glue(out.hdr, '.geno.fam.ft'))
write_feather(meth.samples, path = glue(out.hdr, '.samples.ft'))
write_feather(as.data.frame(geno.mat), path = glue(out.hdr, '.geno.mat.ft'))
