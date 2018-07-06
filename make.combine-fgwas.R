#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 5) {
    q()
}

result.dir <- argv[1]               # e.g., result.dir = 'result/20180613/fgwas/50/'
ld.file <- argv[2]                  # e.g., ld.file = 'LD/fourier_ls-all.bed'
CHR <- as.integer(argv[3])          # e.g., CHR = 1
pip.cutoff <- as.numeric(argv[4])   # e.g., pip.cutoff = 0.9
out.file <- argv[5]                 # e.g., out.file = 'temp'

source('util.R')
library(dplyr)
library(readr)
library(tidyr)

dir.create(dirname(out.file), recursive = TRUE, showWarnings = FALSE)

read.summary <- function(ld.idx, result.dir, lodds.cutoff) {

    .collapse <- function(...) paste(..., collapse = '|')

    snp.file <- result.dir %&&% '/' %&&% ld.idx %&&% '.snp-factor.gz'
    trait.file <- result.dir %&&% '/' %&&% ld.idx %&&% '.trait-factor.gz'
    z.file <- result.dir %&&% '/' %&&% ld.idx %&&% '.zscore.gz'
    resid.file <- result.dir %&&% '/' %&&% ld.idx %&&% '.resid.gz'
    var.file <- result.dir %&&% '/' %&&% ld.idx %&&% '.var.gz'
    snp.tab <- suppressMessages(read_tsv(snp.file))

    if(nrow(snp.tab) < 1) {
        return(NULL)
    }

    valid.factors <- snp.tab %>%
        group_by(factor) %>%
            slice(which.max(lodds)) %>%
                filter(lodds > lodds.cutoff) %>%
                    select(factor) %>% unlist()

    if(length(valid.factors) > 0) {

        trait.tab <- suppressMessages(read_tsv(trait.file))
        traits <- trait.tab %>% select(trait) %>% unique() %>% unlist()
        z.tab <- suppressMessages(read_tsv(z.file))
        resid.tab <- suppressMessages(read_tsv(resid.file))

        var.tab <- suppressMessages(read_tsv(var.file, col_types = 'iiiiccd'))

        var.tot <- var.tab %>% filter(factor == 'total') %>% select(trait, var) %>%
            rename(var.gen = var)
        var.conf <- var.tab %>% filter(factor == 'conf') %>% select(trait, var) %>%
            rename(var.conf = var)
        var.trait <- var.tab %>% filter(factor != 'total', factor != 'conf') %>%
            mutate(factor = as.integer(factor)) %>%
                filter(factor %in% valid.factors) %>%
                    rename(var.factor = var)

        var.summary <-
            var.trait %>%
                select(factor, trait, var.factor) %>%
                    left_join(var.conf) %>%
                        left_join(var.tot) %>%
                            mutate(var.factor = signif(var.factor, 2),
                                   var.conf = signif(var.conf, 2),
                                   var.gen = signif(var.gen, 2))

        z.summary <- snp.tab %>%
            filter(lodds > lodds.cutoff) %>%
                select(snp.loc, factor) %>%
                    left_join(z.tab %>% select(-chr,-rs,-plink.pos)) %>%
                        gather_(key_col= 'snp.best.trait',
                                value_col = 'snp.best.z',
                                gather_col = traits)

        z.summary <- z.summary %>%
            group_by(snp.loc, factor) %>%
                slice(which.max(abs(snp.best.z))) %>%
                    select(starts_with('snp'), factor) %>%
                        as.data.frame()

        ## number of GWAS SNPs within this LD block
        .gwas <- z.tab %>%
            gather_(key_col = 'trait', value_col = 'z', gather_cols = traits) %>%
                mutate(p = 2 * pnorm(abs(z), lower.tail = FALSE))


        .gwas.trait <- .gwas %>%
            group_by(trait) %>%
                summarize(gwas.4 = sum(p < 1e-4),
                          gwas.6 = sum(p < 1e-6),
                          gwas.8 = sum(p < 1e-8))

        .gwas.tot <- .gwas %>%
            group_by(snp.loc) %>%
                slice(which.min(p)) %>%
                    as.data.frame() %>%
                        summarize(gwas.4 = sum(p < 1e-4),
                                  gwas.6 = sum(p < 1e-6),
                                  gwas.8 = sum(p < 1e-8))

        ## number of GWAS SNPs after correcting out
        .gwas.resid <- resid.tab %>%
            gather_(key_col = 'trait', value_col = 'z', gather_cols = traits) %>%
                mutate(p = 2 * pnorm(abs(z), lower.tail = FALSE))

        .gwas.trait.resid <- .gwas.resid %>%
            group_by(trait) %>%
                summarize(resid.4 = sum(p < 1e-4),
                          resid.6 = sum(p < 1e-6),
                          resid.8 = sum(p < 1e-8))

        trait.summary <- trait.tab %>%
            filter(factor %in% valid.factors) %>%
                left_join(.gwas.trait) %>%
                    left_join(.gwas.trait.resid) %>%
                        left_join(var.summary) %>%
                            group_by(chr, ld.idx, LB, UB, factor) %>%
                                summarize(trait = .collapse(trait),
                                          trait.theta = .collapse(theta),
                                          trait.theta.se = .collapse(theta.se),
                                          trait.lodds = .collapse(lodds),
                                          var.factor = .collapse(var.factor),
                                          var.conf = .collapse(var.conf),
                                          var.gen = .collapse(var.gen),
                                          gwas.4 = .collapse(gwas.4),
                                          gwas.6 = .collapse(gwas.6),
                                          gwas.8 = .collapse(gwas.8),
                                          resid.4 = .collapse(resid.4),
                                          resid.6 = .collapse(resid.6),
                                          resid.8 = .collapse(resid.8))

        snp.summary <- snp.tab %>% filter(lodds > lodds.cutoff) %>%
            left_join(z.summary) %>%
                group_by(chr, ld.idx, LB, UB, factor) %>%
                    summarize(snp = .collapse(rs),
                              snp.loc = .collapse(snp.loc),
                              snp.theta = .collapse(theta),
                              snp.theta.se = .collapse(theta.se),
                              snp.lodds = .collapse(lodds),
                              snp.best.trait = .collapse(snp.best.trait),
                              snp.best.z = .collapse(snp.best.z))

        tot.summary <- snp.summary %>%
            left_join(trait.summary)

        return(tot.summary %>% as.data.frame())
    } else {
        return(NULL)
    }
}

ld.info <- read.ld.info(ld.idx = NULL, ld.file) %>%
    filter(chr.input == CHR)

out.tab <-
    ld.info$ld.idx %>%
        lapply(FUN = read.summary,
               result.dir = result.dir,
               lodds.cutoff = log(pip.cutoff) - log(1 - pip.cutoff)) %>%
                   bind_rows()

write_tsv(out.tab, out.file)
