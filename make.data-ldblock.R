#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 2) q()

stat.file <- argv[1] # e.g., stat.file = 'ukbiobank_summary/1677.assoc.tsv.gz'
out.dir <- argv[2]   # e.g., out.dir = '/broad/hptmp/ypp/ukbb/tempdata'

LD.FILE = 'ldblocks/EUR/fourier_ls-all.bed'

options(stringsAsFactors = FALSE)
library(readr)
library(dplyr)
library(tidyr)
source('util.R')

trait <-
    basename(stat.file) %>%
        gsub(pattern = '.assoc.tsv.gz', replacement = '')

stat.tab <- read_tsv(stat.file) %>%
    select(variant, rsid, beta, se) %>%
        separate(variant, c('chr', 'snp.loc', 'a1', 'a2')) %>%
            filter(chr %in% as.character(1:22)) %>%
                mutate(chr = as.integer(chr),
                       snp.loc = as.integer(snp.loc))

gc()

################################################################
## select variants existing in the 1KG data
ref.cols <- c('chr', 'rsid', 'miss', 'snp.loc', 'ref.a1', 'ref.a2')
ref.types <- 'iciicc'
ref.var.tab <- bind_rows(lapply('1KG_EUR/chr' %&&% 1:22 %&&% '.bim', read_tsv,
                                col_names = ref.cols, col_types = ref.types))

valid.stat.tab <- ref.var.tab %>%
    left_join(stat.tab) %>%
        na.omit() %>%
            mutate(beta = if_else(a1 == ref.a2, -beta, beta))

rm(stat.tab)
gc()

################################################################
## break down LD by LD
ld.tab <- read_tsv(LD.FILE) %>%
    rename(CHR = chr, START = start, END = stop) %>%
    mutate(CHR = gsub(CHR, pattern = 'chr', replacement = ''))

for(rr in 1:nrow(ld.tab)) {

    CHR <- ld.tab$CHR[rr]
    SS <- ld.tab$START[rr]
    EE <- ld.tab$END[rr]

    out.file <- out.dir %&&% '/' %&&% CHR %&&% '/' %&&% rr %&&% '/' %&&% trait %&&% '.stat.gz'
    dir.create(dirname(out.file), recursive = TRUE, showWarnings = FALSE)

    if(!file.exists(out.file)) {
        .temp <- valid.stat.tab %>%
            filter(chr == CHR, snp.loc >= SS, snp.loc < EE) %>%
                select(chr, rsid, snp.loc, a1, a2, beta, se) %>%
                    mutate(beta = signif(beta, 4), se = signif(se, 4))

        write_tsv(.temp, path = gzfile(out.file))
        rm(.temp)
        gc()
        log.msg('Wrote %s file\n', out.file)
    }
}

log.msg('Successfully broke down data\n\n')
