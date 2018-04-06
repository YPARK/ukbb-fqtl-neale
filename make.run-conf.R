#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) != 5) q()

in.dir <- argv[1]               # e.g., in.dir = 'data/1/1/'
plink.hdr <- argv[2]            # e.g., plink.hdr = '1KG_EUR/chr1'
re.K <- as.integer(argv[3])     # e.g., re.K = 50
pheno.file <- argv[4]           # e.g., 'phenotypes/ukbb_pheno.txt'
out.hdr <- argv[5]              # e.g., out.hdr = 'temp'

options(stringsAsFactors = FALSE)
source('util.R')

dir.create(dirname(out.hdr), recursive = TRUE, showWarnings = FALSE)
conf.ind.file <- out.hdr %&&% '.conf-ind.gz'
conf.loading.file <- out.hdr %&&% '.conf-trait.gz'
conf.lodds.file <- out.hdr %&&% '.conf-lodds.gz'

.files <- c(conf.ind.file, conf.loading.file, conf.lodds.file)

if(all(sapply(.files, file.exists))) {
    log.msg('Files exists:\n%s\n', paste(.files, collapse = '\n'))
    q()
}

library(readr)
library(dplyr)
library(tidyr)
library(zqtl)

traits <- read_tsv(pheno.file) %>%
    select(Field.code) %>%
        unlist(use.names = FALSE) %>%
            unique() %>% sort()

stat.files <- in.dir %&&% '/' %&&% traits %&&% '.stat.gz'

################################################################
## read plink data
plink <- read.plink(plink.hdr)
colnames(plink$BIM) <- c('chr', 'rsid', 'remove', 'snp.loc', 'plink.a1', 'plink.a2')
plink.snps <- plink$BIM %>% mutate(plink.pos = 1:n()) %>%
    select(-remove)

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
    write_tsv(data.frame(), path = conf.loading.file)
    write_tsv(data.frame(), path = conf.ind.file)
    write_tsv(data.frame(), path = conf.lodds.file)
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

if(nrow(beta.tab) < 2) {
    log.msg('Too few number of rsids\n\n')
    write_tsv(data.frame(), path = conf.loading.file)
    write_tsv(data.frame(), path = conf.ind.file)
    write_tsv(data.frame(), path = conf.lodds.file)
    q()
}

rm(stat.tab)
gc()

fam <- plink$FAM
colnames(fam) <- c('fid', 'iid', 'father', 'mother', 'sex', 'pheno')

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

log.msg('Constructed data\n\n')

################################################################
re.K <- max(min(c(length(traits) - 1, ncol(X) - 1, re.K)), 1)

## Just run factorization to check if there were any confounders
vb.opt <- list(pi.ub = -1, pi.lb = -3, tau = -5, do.hyper = TRUE,
               right.nn = FALSE, do.stdize = TRUE, do.rescale = TRUE, 
               eigen.tol = 1e-1, gammax = 1e3, vbiter = 5000,
               svd.init = TRUE, jitter = 0.1,
               tol = 1e-8, rate = 1e-2, k = re.K)

z.out <- fit.zqtl.factorize(effect = beta.mat, effect.se = se.mat,
                            X = X, options = vb.opt)

log.msg('Finished factorization estimation\n\n')

colnames(z.out$param.indiv$theta) <- 'factor.' %&&% 1:re.K
Z.conf.ind.tab <- data.frame(iid = fam$iid, signif(z.out$param.indiv$theta, 4))

colnames(z.out$param.trait$theta) <- 'factor.' %&&% 1:re.K
Z.conf.loading.tab <- data.frame(traits, signif(z.out$param.trait$theta, 4))

conf.lodds <- signif(as.numeric(z.out$param.trait$lodds), 4)
Z.conf.lodds.tab <- data.frame(conf.factor = 1:re.K, conf.lodds)

write_tsv(Z.conf.loading.tab, path = conf.loading.file)
write_tsv(Z.conf.ind.tab, path = conf.ind.file)
write_tsv(Z.conf.lodds.tab, path = conf.lodds.file)

log.msg('Successfully finished fQTL\n\n')
