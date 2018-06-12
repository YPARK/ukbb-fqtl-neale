#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 4) {
    q()
}

## create non-genetic confounder matrix for this LD block
ld.idx <- as.integer(argv[1])
ld.file <- argv[2]            # e.g., 'LD/fourier_ls-all.bed'
data.dir <- argv[3]           # e.g., 'result/20180606/factorization/20/'
out.file <- argv[4]

dir.create(dirname(out.file), recursive = TRUE, showWarnings = FALSE)
if(file.exists(out.file)) {
    q()
}
options(stringsAsFactors = FALSE)
source('util.R')

pval.cutoff <- 1e-4
lodds.cutoff <- 0

################################################################
source('util.fgwas.R')
library(dplyr)
library(readr)
ld.tab <- read.ld.info(ld.idx=NULL, ld.file)

################################################################
.read.conf <- function(ld.idx, data.dir, ...) {
    ## read factors on individuals
    ret <- suppressMessages(read.conf.ind(data.dir %&&% '/' %&&% ld.idx, ...))
    if(ncol(ret) < 2) { return(ret) }
    ## rename columns
    .names <- colnames(ret[, -1, drop = FALSE])
    colnames(ret) <- c('iid', ld.idx %&&% '_' %&&% .names)
    return(ret)
}

.join <- function(df1, df2) {
    if(is.null(df1) && is.null(df2)) return(NULL)
    if(is.null(df1)) return(df2)
    if(is.null(df2)) return(df1)
    left_join(df1, df2, by = 'iid')
}

.reduce.join <- function(...) Reduce(..., f = .join, init = NULL)

check.others <- function(other.idx, data.dir, Y1, lodds.cutoff, n.ctrl, pval.cutoff) {
    ret <- Y1[, 1, drop = FALSE]
    Y0 <- .read.conf(other.idx, data.dir, lodds.cutoff = lodds.cutoff)
    if(ncol(Y0) < 2) return(ret)
    iid.pos <- match(Y1$iid, Y0$iid)
    conf.cols <- find.cor.idx(Y1[, -1], Y0[iid.pos, -1], n.ctrl, pval.cutoff)
    if(length(conf.cols) < 1) return(ret)
    log.msg('Found %d correlated column(s) in %d\n', length(conf.cols), other.idx)
    ret <- Y0[iid.pos, c(1, 1 + conf.cols)]
    return(ret)
}

Y1 <- .read.conf(ld.idx, data.dir, lodds.cutoff = lodds.cutoff)

if(ncol(Y1) < 2) {
    if(ncol(Y1) > 1) {
        write_tsv(Y1 %>% select(iid), out.file)
    } else {
        write_tsv(data.frame(), out.file)
    }
    q()
}

others <- setdiff(ld.tab$ld.idx, ld.idx)

other.Y.list <- others %>% lapply(FUN = check.others,
                                  data.dir = data.dir,
                                  Y1 = Y1,
                                  lodds.cutoff = lodds.cutoff,
                                  n.ctrl = 3,
                                  pval.cutoff = pval.cutoff)

Y0 <- other.Y.list %>% .reduce.join() %>%
    (function(x) Y1 %>% select(iid) %>% left_join(x))

write_tsv(Y0, out.file)
