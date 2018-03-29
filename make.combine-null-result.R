#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
source('util.R')
library(dplyr)
library(readr)
library(tidyr)

.take.tab <- function(.blk, file.name = 'null') {

    .trait.file <- .blk %&&% '/' %&&% file.name %&&% '.trait-factor.gz'
    .var.file <- .blk %&&% '/' %&&% file.name %&&% '.var.gz'

    if(! all(sapply(c(.trait.file, .var.file), file.exists)) ) return(NULL)

    .trait.tab <- suppressMessages(read_tsv(.trait.file, col_names = TRUE, col_types = 'iiiccddd'))
    .var.tab <- suppressMessages(read_tsv(.var.file, col_names = TRUE, col_types = 'iiiccd'))

    ret <- suppressMessages(left_join(.var.tab, .trait.tab)) %>% as.data.frame()
    return(ret)
}

result.dir <- 'result/null/'
.list.files <- function(chr) list.files(result.dir %&&% chr %&&% '/', full.names = TRUE)
total.blks <- do.call(c, lapply(1:22, .list.files))
null.tab <- do.call(rbind, lapply(total.blks, .take.tab, file.name = 'null'))
write_tsv(null.tab, path = gzfile('result/ukbb-null.txt.gz'))
rm(null.tab); gc();

result.dir <- 'result/null-kron/'
.list.files <- function(chr) list.files(result.dir %&&% chr %&&% '/', full.names = TRUE)
total.blks <- do.call(c, lapply(1:22, .list.files))
null.tab <- do.call(rbind, lapply(total.blks, .take.tab, file.name = 'null-kron'))
write_tsv(null.tab, path = gzfile('result/ukbb-null-kron.txt.gz'))
rm(null.tab); gc();

