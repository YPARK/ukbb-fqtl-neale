#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 7) {
    q()
}

ld.idx <- as.integer(argv[1])
ld.file <- argv[2]            # e.g., 'LD/fourier_ls-all.bed'
data.dir <- argv[3]           # e.g., 'data/'
pheno.file <- argv[4]         # e.g., 'phenotypes/ukbb_pheno_CC.txt'
conf.file <- argv[5]          # e.g., 'result/20180606/confounder/50/1.txt.gz'
K <- as.integer(argv[6])      # e.g., 10
out.hdr <- argv[7]            # e.g., 'temp'

non.neg <- FALSE

if(length(argv) > 7) {
    non.neg <- as.logical(argv[8])
}

options(stringsAsFactors = FALSE)
source('util.R')

dir.create(dirname(out.hdr), recursive = TRUE, showWarnings = FALSE)

snp.factor.file <- out.hdr %&&% '.snp-factor.gz'
trait.factor.file <- out.hdr %&&% '.trait-factor.gz'
zscore.file <- out.hdr %&&% '.zscore.gz'
resid.file <- out.hdr %&&% '.resid.gz'
var.file <- out.hdr %&&% '.var.gz'
conf.out.file <- out.hdr %&&% '.conf.gz'

.files <- c(snp.factor.file, trait.factor.file, zscore.file, resid.file, var.file, conf.out.file)

if(all(sapply(.files, file.exists))) {
    log.msg('Files exists:\n%s\n', paste(.files, collapse = '\n'))
    q()
}

################################################################
source('util.fgwas.R')
library(dplyr)
library(readr)

ld.info <- read.ld.info(ld.idx, ld.file)

temp.dir <- system('mkdir -p /broad/hptmp/ypp/ukbb-fgwas/' %&&% out.hdr %&&%
                   '; mktemp -d /broad/hptmp/ypp/ukbb-fgwas/' %&&% out.hdr %&&%
                   '/temp.XXXXXXXX',
                   intern = TRUE,
                   ignore.stderr = TRUE)
plink.hdr <- 'UK10K_1KG/chr' %&&% ld.info$chr.input
plink <- subset.plink(plink.hdr = plink.hdr,
                      chr = ld.info$chr.input,
                      plink.lb = ld.info$lb.input,
                      plink.ub = ld.info$ub.input,
                      temp.dir = temp.dir,
                      maf = 0.01)
if(dir.exists(temp.dir)) system('rm -r ' %&&% temp.dir)

if(is.null(plink)) {
    log.msg('Failed to read plink\n')
    q()
}

chr.input <- ld.info$chr.input
pheno.tab <- read_tsv(pheno.file)
traits <- pheno.tab$Field.code
gwas.files <- traits %>%
    (function(x) data.dir %&&% '/' %&&% chr.input %&&% '/' %&&% ld.idx %&&% '/' %&&% x %&&% '.stat.gz')

if(!all(file.exists(gwas.files))) {
    missing.files <- gwas.files %>%
        sapply(FUN = function(x) ifelse(file.exists(x), '', x)) %>%
            paste(collapse = '\n')
    log.msg('Missing GWAS files :\n%s\n', missing.files)
    q()
}

gwas.stat.tab <- gwas.files %>% read.gwas.files(plink.obj = plink)
if(is.null(gwas.stat.tab)) {
    log.msg('Failed to read GWAS files\n')
    q()
}

if(nrow(gwas.stat.tab) < 1) {
    write_tsv(data.frame(), path = gzfile(conf.out.file))
    write_tsv(data.frame(), path = gzfile(var.file))
    write_tsv(data.frame(), path = gzfile(snp.factor.file))
    write_tsv(data.frame(), path = gzfile(trait.factor.file))
    write_tsv(data.frame(), path = gzfile(zscore.file))
    write_tsv(data.frame(), path = gzfile(resid.file))
    q()
}

gwas.data <- construct.data.matrix(gwas.stat.tab, plink, traits)
if(is.null(gwas.data)) {
    log.msg('Failed to construct GWAS data\n')
    q()
}

################################################################
## z-scores for potential confounding signals
conf.tab <- read_tsv(conf.file)

n <- nrow(gwas.data$X)
z.conf.0 <- t(gwas.data$X) %*% matrix(1, n, 1) / sqrt(n)

if(ncol(conf.tab) < 2) {
    Z.conf <- z.conf.0
} else {
    C <- plink$FAM %>% select(iid) %>%
        left_join(conf.tab) %>%
            select(-iid) %>%
                as.matrix()

    C <- scale(C)
    
    stopifnot(all(!is.na(C)))
    log.msg('Found %d covaraites\n\n', ncol(C))
    Z.conf <- t(gwas.data$X) %*% C / sqrt(n)
    Z.conf <- cbind(Z.conf, z.conf.0)
}

################################################################
K <- max(min(c(length(traits) - 1, ncol(gwas.data$X) - 1, K)), 1)

vb.opt <- list(pi.ub = -0.1, pi.lb = -3, tau = -5, do.hyper = TRUE,
               do.stdize = TRUE, eigen.tol = 1e-2, gammax = 1e2,
               svd.init = TRUE, do.rescale = TRUE,
               jitter = 0.1, vbiter = 7500, rate = 1e-2, decay = -1e-2,
               tol = 0, right.nn = non.neg, k = K)

z.out <- fit.zqtl(effect = gwas.data$beta, effect.se = gwas.data$se,
                  X = gwas.data$X, C.delta = Z.conf, factored = TRUE,
                  options = vb.opt)

log.msg('Finished fQTL estimation\n\n')

LD.info <- gwas.data$snps %>%
    summarize(chr = min(chr),
              ld.idx,
              LB = ld.info$lb.input,
              UB = ld.info$ub.input)

var.tab <- estimate.variance(z.out, traits, se.mat = gwas.data$se)
var.tab <- data.frame(LD.info, var.tab)

right.tab <- melt.effect(z.out$param.right, traits, 1:K) %>%
    rename(theta.se = theta.var, trait = row, factor = col) %>%
        mutate(theta.se = signif(sqrt(theta.se), 2)) %>%
            as.data.frame()

right.tab <- data.frame(LD.info, right.tab)

left.tab <- melt.effect(z.out$param.left, gwas.data$zz$rs, 1:K) %>%
    rename(theta.se = theta.var, rs = row, factor = col) %>%
        mutate(theta.se = signif(sqrt(theta.se), 2)) %>%
            left_join(gwas.data$zz %>% select(rs, snp.loc)) %>%
                left_join(gwas.data$snps %>% select(-missing, -chr)) %>%
                    as.data.frame()

if(nrow(left.tab) > 0){
    left.tab <- data.frame(LD.info, left.tab)
} else {
    left.tab <- data.frame()
}

n.c <- ncol(Z.conf)

conf.assoc.tab <- melt.effect(z.out$conf.delta, 1:n.c, traits) %>%
    rename(theta.se = theta.var, conf = row, trait = col) %>%
        mutate(theta.se = signif(sqrt(theta.se), 2)) %>%
            as.data.frame()

## estimate residual z-score matrix w/o confounding effects
resid.tab <- gwas.data$zz
resid.hat <- t(z.out$Vt) %*% z.out$Y - Z.conf %*% z.out$conf.delta$theta
resid.tab[, -(1:4)] <- round(resid.hat, 4)

################################################################
write_tsv(conf.assoc.tab, path = gzfile(conf.out.file))
write_tsv(var.tab, path = gzfile(var.file))
write_tsv(left.tab, path = gzfile(snp.factor.file))
write_tsv(right.tab, path = gzfile(trait.factor.file))
write_tsv(gwas.data$zz, path = gzfile(zscore.file))
write_tsv(resid.tab, path = gzfile(resid.file))

log.msg('Successfully finished fQTL\n\n')
