#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
source('util.R')
library(dplyr)
library(readr)

.safe.read <- function(...) {
    ret <- suppressMessages(read_tsv(...))
    if(nrow(ret) == 0) return(NULL)
    return(ret)
}

.read.conf <- function(hdr, chr, cutoff = log(0.25) - log(0.75)) {

    ## log.msg('Reading %s\n', hdr)

    lodds <- .safe.read(hdr %&&% '.conf-lodds.gz')
    conf.ind <- .safe.read(hdr %&&% '.conf-ind.gz')

    if(is.null(lodds) || is.null(conf.ind)) return(NULL)

    ret0 <- conf.ind %>% select(iid) %>% as.data.frame()
    ret1 <- conf.ind %>% select(-iid) %>% as.data.frame()

    .hdr <- chr %&&% '_' %&&% basename(hdr) %&&% '_'

    .cols <- .hdr %&&% names(ret1)
    names(ret1) <- .cols

    factors <- lodds %>% filter(conf.lodds > cutoff) %>%
        select(conf.factor) %>% unlist(use.names = FALSE)

    if(length(factors) < 1) return(ret0)
    .cols <- .hdr %&&% 'factor.' %&&% factors %>% sapply(as.name)
    ret1 <- ret1 %>% dplyr::select_(.dots = .cols)
    ret <- cbind(ret0, ret1 %>% signif(digits = 4))
    return(ret)
}

.join <- function(df1, df2) {
    if(is.null(df1) && is.null(df2)) return(NULL)
    if(is.null(df1)) return(df2)
    if(is.null(df2)) return(df1)
    left_join(df1, df2, by = 'iid')
}

.reduce.join <- function(...) Reduce(..., f = .join, init = NULL)

read.chr <- function(chr, K) {
    
    log.msg('Reading %d\n', chr)

    data.dir <- 'result/conf/' %&&% K %&&% '/' %&&% chr %&&% '/'
    
    file.hdr.list <- list.files(path = data.dir, pattern = '.conf-lodds.gz', full.names = TRUE) %>%
        gsub(pattern = '[.]conf-lodds.gz', replacement = '')
    
    conf.tab <- lapply(file.hdr.list, .read.conf, chr = chr) %>%
        .reduce.join()

    return(conf.tab %>% as.data.frame())
}

conf.chr.list <- lapply(1:22, read.chr, K = 50)

## leave one chromosome out
ret.loco.list <- lapply(1:22, function(chr) .reduce.join(conf.chr.list[-chr]))
loco.files <- 'result/conf/50/loco_' %&&% 1:22 %&&% '.txt.gz'

dir.create('result/conf/50/', recursive = TRUE)

temp <- sapply(1:22, function(chr) write_tsv(ret.loco.list[[chr]], loco.files[chr]))

