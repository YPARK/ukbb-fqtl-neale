#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
source('util.R')

library(dplyr)
library(readr)
library(tidyr)

.files <- 'FQTL-v6/chr' %&&% 1:22 %&&% '/50/combined.txt.gz'

gtex.cols <- c('ensg', 'chr', 'tss',
               'tis.idx', 'tis.names', 'tis.theta', 'tis.se', 'tis.lodds',
               'snp.names', 'snp.theta', 'snp.se', 'snp.lodds',
               'factor', 'pip.cutoff')

tot.tab <- bind_rows(lapply(.files, read_tsv, col_names = gtex.cols)) %>%
    filter(pip.cutoff == 0.5)

.split.bar <- function(s) strsplit(s, split = '[|]')

tis.tab <- tot.tab %>%
    select(ensg, factor, starts_with('tis')) %>%
    unnest(idx = .split.bar(tis.idx),
           name = .split.bar(tis.names),
           theta = .split.bar(tis.theta),
           se = .split.bar(tis.se),
           lodds = .split.bar(tis.lodds)) %>%
    select(-starts_with('tis')) %>%
    separate(ensg, c('ensg', 'remove'), sep = '[.]') %>%
    select(-remove)

write_tsv(tis.tab, path = 'gtex-v6-fqtl-tis.txt.gz')
    

    
