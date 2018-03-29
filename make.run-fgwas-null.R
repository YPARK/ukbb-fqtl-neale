#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) != 5) q()

in.dir <- argv[1]               # e.g., in.dir = 'data/1/1/'
plink.hdr <- argv[2]            # e.g., plink.hdr = '1KG_EUR/chr1'
K.max <- as.integer(argv[3])    # e.g., K.max = 50
NON.NEG <- as.logical(argv[4])  # e.g., NON.NEG = TRUE
out.hdr <- argv[5]              # e.g., out.hdr = 'temp'

options(stringsAsFactors = FALSE)
source('util.R')

dir.create(dirname(out.hdr), recursive = TRUE, showWarnings = FALSE)
snp.factor.file <- out.hdr %&&% '.snp-factor.gz'
trait.factor.file <- out.hdr %&&% '.trait-factor.gz'
zscore.file <- out.hdr %&&% '.zscore.gz'
var.file <- out.hdr %&&% '.var.gz'
conf.out.file <- out.hdr %&&% '.conf.gz'

.files <- c(snp.factor.file, trait.factor.file, zscore.file, var.file, conf.out.file)

if(all(sapply(.files, file.exists))) {
    log.msg('Files exists:\n%s\n', paste(.files, collapse = '\n'))
    q()
}

library(readr)
library(dplyr)
library(tidyr)
library(zqtl)

traits <- read_tsv('phenotypes/ukbb_pheno.txt') %>%
    select(Field.code) %>%
        unlist(use.names = FALSE) %>%
            unique() %>% sort()

stat.files <- in.dir %&&% '/' %&&% traits %&&% '.stat.gz'

if(!all(sapply(stat.files, file.exists))) {
    temp <- paste(stat.files[!sapply(stat.files, file.exists)], collapse = ', ')
    log.msg('Missing stat files : %s\n', temp)
    q()
}

################################################################
## read plink data
plink <- read.plink(plink.hdr)
colnames(plink$BIM) <- c('chr', 'rsid', 'remove', 'snp.loc', 'plink.a1', 'plink.a2')
plink.snps <- plink$BIM %>% mutate(plink.pos = 1:n()) %>%
    select(-remove)
plink.fam <- plink$FAM
colnames(plink.fam) <- c('fid', 'iid', 'father', 'mother', 'sex', 'pheno')

################################################################
## read statistics
read.stat <- function(stat.file) {
    stat.cols <- c('chr', 'rsid', 'snp.loc', 'a1', 'a2', 'beta', 'se')
    stat.types <- 'iciccdd'

    trait <- basename(stat.file) %>%
        gsub(pattern='.stat.gz', replacement = '')

    ret <- read_tsv(stat.file, col_names = stat.cols, col_types = stat.types, skip = 1) %>%
        group_by(chr, rsid, snp.loc, a1, a2) %>%
            summarize(beta = mean(beta), se = mean(se)) %>%
                mutate(trait = trait)
    return(ret)
}

## match directionality between 1kg and ukbb
stat.tab <- bind_rows(lapply(stat.files, read.stat)) %>%
    left_join(plink.snps) %>% na.omit() %>%
        filter((a1 == plink.a2 & a2 == plink.a1) | (a1 == plink.a1 & a2 == plink.a2)) %>%
            mutate(beta = if_else(a1 == plink.a1, -beta, beta)) %>%
                as.data.frame()

if(nrow(stat.tab) < 1) {
    log.msg('Empty stat\n\n')
    write_tsv(data.frame(), path = gzfile(var.file))
    write_tsv(data.frame(), path = gzfile(conf.out.file))
    write_tsv(data.frame(), path = gzfile(snp.factor.file))
    write_tsv(data.frame(), path = gzfile(trait.factor.file))
    write_tsv(data.frame(), path = gzfile(zscore.file))
    q()
}

nsnps.traits <- stat.tab %>% select(trait, snp.loc) %>% unique() %>%
    group_by(trait) %>% summarize(n = n())

## check if there were different number of rsids
max.nsnps <- max(nsnps.traits$n)
discrepancy <- nsnps.traits %>% filter(n < (max.nsnps/2))

if(discrepancy %>% nrow() > 0) {
    library(pander)
    log.msg('Some traits might have preprocessing errors:\n%s\n%s\n\n',
            pandoc.table.return(nsnps.traits, style = 'simple'),
            paste(stat.files, collapse = '\n'))
    q()
}

## construct beta and standard error matrix
beta.tab <- stat.tab %>% select(chr, rsid, snp.loc, a1, a2, plink.pos, beta, trait) %>%
    spread(key = trait, value = beta) %>% arrange(snp.loc, a1) %>% na.omit() %>%
        as.data.frame()

se.tab <- stat.tab %>% select(chr, rsid, snp.loc, a1, a2, plink.pos, se, trait) %>%
    spread(key = trait, value = se) %>% arrange(snp.loc, a1) %>% na.omit() %>%
        as.data.frame()

zscore.tab <- stat.tab %>% select(chr, rsid, snp.loc, a1, a2, plink.pos, beta, se, trait) %>%
    mutate(z = signif(beta / se, 2)) %>% select(-beta, -se) %>%
        spread(key = trait, value = z) %>% arrange(snp.loc, a1) %>% na.omit() %>%
            as.data.frame()

if(nrow(beta.tab) < 2) {
    log.msg('Too few number of rsids\n\n')
    write_tsv(data.frame(), path = gzfile(var.file))
    write_tsv(data.frame(), path = gzfile(conf.out.file))
    write_tsv(data.frame(), path = gzfile(snp.factor.file))
    write_tsv(data.frame(), path = gzfile(trait.factor.file))
    write_tsv(data.frame(), path = gzfile(zscore.file))
    q()
}

rm(stat.tab)
gc()

X <- plink$BED %c% beta.tab$plink.pos %>% scale()
X[is.na(X)] <- 0
plink.snps <- plink.snps %r% beta.tab$plink.pos
rm(plink)
gc()

################################################################
tab2mat <- function(tab) {
    .loc <- match(traits, names(tab))
    ret <- tab %c% .loc  %>% as.matrix()
    return(ret)
}

beta.mat <- plink.snps %>% select(chr, rsid, snp.loc) %>% left_join(beta.tab) %>% na.omit() %>%
    tab2mat()

se.mat <- plink.snps %>% select(chr, rsid, snp.loc) %>% left_join(se.tab) %>% na.omit() %>%
    tab2mat()

rm(se.tab)
rm(beta.tab)
gc()

zscore.tab <- plink.snps %>% select(chr, rsid, snp.loc) %>% left_join(zscore.tab) %>% na.omit()
log.msg('Constructed data\n\n')

################################################################
K <- max(min(c(length(traits) - 1, ncol(X) - 1, K.max)), 1)

################################################################
## generate null data
beta.null <- make.zqtl.null(X = X, beta.mat = beta.mat, se.mat = se.mat, eig.tol = 1e-2, stdize = TRUE)

################################################################
## correct confounder
vb.opt <- list(pi.ub = -1, pi.lb = -4, tau = -5, do.hyper = TRUE,
               do.stdize = TRUE, eigen.tol = 1e-2, gammax = 1e4,
               svd.init = TRUE, jitter = 0.1, vbiter = 5000,
               tol = 1e-8, rate = 1e-2, right.nn = FALSE,
               k = K)

z.fact.out <- fit.zqtl.factorize(beta.null, se.mat, X, options = vb.opt)

cutoff <- log(0.9) - log(0.1)
valid <- which(z.fact.out$param.trait$lodds > cutoff)

if(length(valid) < 1) {
    n <- nrow(X)
    Z.conf <- t(X) %*% matrix(1, n, 1) / sqrt(n)
} else {
    Z.conf <- t(X) %*% (z.fact.out$param.indiv$theta %c% valid) / sqrt(nrow(X))
}

## fit the model
vb.opt <- list(pi.ub = -1, pi.lb = -3, tau = -5, do.hyper = TRUE,
               do.stdize = TRUE, eigen.tol = 1e-2, gammax = 1e4,
               svd.init = TRUE, do.rescale = TRUE,
               jitter = 0.1, vbiter = 5000,
               tol = 1e-8, rate = 1e-2, right.nn = NON.NEG,
               k = K)

z.out <- fit.zqtl(effect = beta.null, effect.se = se.mat,
                  X = X, C.delta = Z.conf, factored = TRUE, options = vb.opt)

log.msg('Finished fQTL estimation\n\n')

LD.info <- plink.snps %>% summarize(chr = min(chr), LB = min(snp.loc), UB = max(snp.loc))

################################################################
## Variance of each factor
theta.snp <- z.out$param.left$theta
theta.trait <- z.out$param.right$theta

take.var.k <- function(k) {

    theta.k <- (theta.snp %c% k) %*% t(theta.trait %c% k)
    theta.k <- theta.k * se.mat ## scale by SE
    eta.k <- sweep(z.out$Vt %*% theta.k, 1, sqrt(z.out$D2), `*`)
    var.k <- apply(eta.k^2, 2, sum)

    return(data.frame(trait = as.character(names(var.k)), factor = k, var = var.k))
}

theta.tot <- (theta.snp %*% t(theta.trait)) * se.mat ## scale by SE
eta.tot <- sweep(z.out$Vt %*% theta.tot, 1, sqrt(z.out$D2), `*`)
var.tot <- data.frame(trait = as.character(traits), factor = 'total', var = apply(eta.tot^2, 2, sum))

var.tab <- bind_rows(lapply(1:K, take.var.k))
rownames(var.tab) <- NULL
var.tab <- rbind(var.tab, var.tot) %>%
    mutate(var = signif(var, 4))

var.tab <- data.frame(LD.info, var.tab)

log.msg('Calculated Variance\n\n')

right.tab <- melt.effect(z.out$param.right, traits, 1:K) %>%
    rename(theta.se = theta.var, trait = row, factor = col) %>%
        mutate(theta.se = signif(sqrt(theta.se), 2)) %>%
            as.data.frame()

right.tab <- data.frame(LD.info, right.tab)

left.tab <- melt.effect(z.out$param.left, zscore.tab$rsid, 1:K) %>%
    rename(theta.se = theta.var, rsid = row, factor = col) %>%
        mutate(theta.se = signif(sqrt(theta.se), 2)) %>%
            left_join(zscore.tab %>% select(rsid, snp.loc, a1, a2)) %>%
                as.data.frame()

if(nrow(left.tab) > 0){
    left.tab <- data.frame(LD.info, left.tab)
} else {
    left.tab <- data.frame()
}

conf.assoc.tab <- melt.effect(z.out$conf.delta, 1:ncol(Z.conf), traits) %>%
    rename(theta.se = theta.var, conf = row, trait = col) %>%
        mutate(theta.se = signif(sqrt(theta.se), 2)) %>%
            as.data.frame()

write_tsv(conf.assoc.tab, path = gzfile(conf.out.file))
write_tsv(var.tab, path = gzfile(var.file))
write_tsv(left.tab, path = gzfile(snp.factor.file))
write_tsv(right.tab, path = gzfile(trait.factor.file))
write_tsv(zscore.tab, path = gzfile(zscore.file))

log.msg('Successfully finished fQTL\n\n')
