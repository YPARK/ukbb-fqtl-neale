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

.read.conf <- function(hdr, chr, cutoff = log(0.9) - log(0.1)) {

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

read.chr <- function(chr, ddir) {
    
    log.msg('Reading %d\n', chr)

    data.dir <- ddir %&&% '/' %&&% chr %&&% '/'
    
    file.hdr.list <- list.files(path = data.dir, pattern = '.conf-lodds.gz', full.names = TRUE) %>%
        gsub(pattern = '[.]conf-lodds.gz', replacement = '')
    
    conf.tab <- lapply(file.hdr.list, .read.conf, chr = chr) %>%
        .reduce.join()

    return(conf.tab %>% as.data.frame())
}

find.loco <- function(chr, data.list, p.val.cutoff) {

    .search.other <- function(other, self) {

        iid <- self[, 1]
        iid.other <- match(iid, other[, 1])
        ret0 <- other %c% 1

        stat <- calc.qtl.stat(other %r% iid.other %c% -1, self %c% -1) %>%
            filter(p.val < p.val.cutoff)

        if(nrow(stat) < 1) return(ret0)

        ret <- cbind(ret0, other %c% (1 + stat$x.col))
        log.msg('Found %d LOCO\n', ncol(ret) - 1)
        return(ret)
    }

    loco <- lapply(data.list[-chr], .search.other, self = data.list[[chr]]) %>%
        .reduce.join()

    log.msg('Constructed %d LOCO on chr%d\n', ncol(loco) - 1, chr)

    return(loco)
}


################################################################
conf.chr.list <- lapply(1:22, read.chr, ddir = 'result/conf/Q/30')
n.conf <- sum(sapply(conf.chr.list, ncol) - 1)
cutoff <- 0.05/n.conf

for(chr in 1:22) {
    loco.file <- 'result/conf/Q/30/loco_' %&&% chr %&&% '.txt.gz'
    if(!file.exists(loco.file)) {
        out <- find.loco(chr, data.list = conf.chr.list, p.val.cutoff = cutoff)
        write_tsv(out, loco.file)
    }
}

################################################################
conf.chr.list <- lapply(1:22, read.chr, ddir = 'result/conf/CC/30')
n.conf <- sum(sapply(conf.chr.list, ncol) - 1)
cutoff <- 0.05/n.conf

for(chr in 1:22) {
    loco.file <- 'result/conf/CC/30/loco_' %&&% chr %&&% '.txt.gz'
    if(!file.exists(loco.file)) {
        out <- find.loco(chr, data.list = conf.chr.list, p.val.cutoff = cutoff)
        write_tsv(out, loco.file)
    }
}

