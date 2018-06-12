#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 3) q()

result.dir <- argv[1]               # e.g., result.dir = 'result/fgwas_nn/100/Q/'
lodds.cutoff <- as.numeric(argv[2]) # e.g., lodds.cutoff = 0
out.hdr <- argv[3]                  # e.g., out.hdr = 'temp'

source('util.R')
library(dplyr)
library(readr)
library(tidyr)

dir.create(dirname(out.hdr), recursive = TRUE, showWarnings = FALSE)

trait.file <- out.hdr %&&% '_traits.txt.gz'
trait.full.file <- out.hdr %&&% '_full_traits.txt.gz'
var.file <- out.hdr %&&% '_var.txt.gz'
snp.file <- out.hdr %&&% '_snps.txt.gz'

.take.trait.tab <- function(.blk, lodds.cutoff = 0) {

    .trait.file <- .blk %&&% '.trait-factor.gz'
    .blk.idx <- .blk %>% basename() %>% as.integer()

    if(!file.exists(.trait.file)) return(NULL)

    if(lodds.cutoff <= -Inf) {
        .trait.tab <- suppressMessages(read_tsv(.trait.file, col_names = TRUE))
        return(.trait.tab %>% as.data.frame())
    }

    .snp.file <- .blk %&&% '.snp-factor.gz'

    if(!file.exists(.snp.file)) return(NULL)

    .snp.tab <- suppressMessages(read_tsv(.snp.file, col_names= TRUE) %>%
                                     dplyr::filter(lodds > lodds.cutoff))

    if(nrow(.snp.tab) > 0) {
        .factors <- .snp.tab %>% dplyr::select(factor) %>% unique()
        .trait.tab <- suppressMessages(read_tsv(.trait.file, col_names = TRUE) %>%
                                           dplyr::filter(factor %in% .factors$factor))
        .trait.tab <- .trait.tab %>% mutate(ld.idx = .blk.idx)
        return(.trait.tab %>% as.data.frame())
    }
    return(NULL)
}

.take.var.tab <- function(.blk) {
    .var.file <- .blk %&&% '.var.gz'
    .blk.idx <- .blk %>% basename() %>% as.integer()

    if(!file.exists(.var.file)) return(NULL)
    .var.tab <- suppressMessages(read_tsv(.var.file, col_names = TRUE, col_types = 'iiiccd'))
    .var.tab <- .var.tab %>% mutate(ld.idx = .blk.idx)
    return(.var.tab %>% as.data.frame())
}

.take.snp.tab <- function(.blk, lodds.cutoff = 0) {
    .snp.file <- .blk %&&% '.snp-factor.gz'
    .blk.idx <- .blk %>% basename() %>% as.integer()

    if(!file.exists(.snp.file)) return(NULL)

    .snp.tab <- suppressMessages(read_tsv(.snp.file, col_names= TRUE) %>%
                                     dplyr::filter(lodds > lodds.cutoff))
    if(nrow(.snp.tab) > 0) {
        return(.snp.tab)
    }
    .snp.tab <- .snp.tab %>% mutate(ld.idx = .blk.idx)
    return(NULL)
}

.list.files <- function(chr) list.files(result.dir %&&% chr %&&% '/', pattern = '.zscore.gz', full.names = TRUE)

total.blks <- do.call(c, lapply(1:22, .list.files)) %>%
    sapply(gsub, pattern = '.zscore.gz', replacement = '')

if(!file.exists(trait.file)) {
    trait.tab <- do.call(rbind, lapply(total.blks, .take.trait.tab, lodds.cutoff = lodds.cutoff))
    log.msg('Read traits\n')
    write_tsv(trait.tab, path = trait.file)
    rm(trait.tab); gc();
}

if(!file.exists(var.file)) {
    var.tot.tab <- do.call(rbind, lapply(total.blks, .take.var.tab))
    log.msg('Read var\n')
    write_tsv(var.tot.tab, path = var.file)
    rm(var.tot.tab); gc();
}

if(!file.exists(snp.file)) {
    snp.tab <- do.call(rbind, lapply(total.blks, .take.snp.tab, lodds.cutoff = lodds.cutoff))
    log.msg('Read SNPs\n')
    write_tsv(snp.tab, path = snp.file)
}
