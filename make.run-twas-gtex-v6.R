#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) != 3) {
    q()
}

options(stringsAsFactors = FALSE)
source('util.R')
source('util.twas.R')
library(readr)
library(dplyr)
library(tidyr)
library(zqtl)

ld.idx <- as.integer(argv[1]) # e.g., ld.idx = 2
fgwas.hdr <- argv[2]          # e.g., fgwas.hdr = 'result/fgwas_nn/100/CC'
out.file <- argv[3]           # e.g., out.file = 'temp.txt.gz'

if(file.exists(out.file)) {
    log.msg('%s already exists\n', out.file)
    q()
}

dir.create(dirname(out.file), recursive = TRUE, showWarnings = FALSE)

################################################################
ld.file <- 'ldblocks/EUR/fourier_ls-all.bed'
ld.tab <- read_tsv(ld.file)[ld.idx, ]
chr.input <- ld.tab$chr %>% gsub(pattern = 'chr', replacement = '') %>% as.integer()
lb.input <- ld.tab$start %>% as.integer()
ub.input <- ld.tab$stop %>% as.integer()

fgwas.file <- fgwas.hdr %&&% '/' %&&% chr.input %&&% '/' %&&% ld.idx %&&% '.snp-factor.gz'
plink.hdr <- '1KG_EUR/chr' %&&% chr.input

temp.dir <- system('mkdir -p /broad/hptmp/ypp/ukbb-twas/' %&&%
                   out.file %&&% '; mktemp -d /broad/hptmp/ypp/ukbb-twas/' %&&%
                   out.file %&&% '/temp.XXXXXXXX',
                   intern = TRUE,
                   ignore.stderr = TRUE)

pip.cutoff <- 0.5
lodds.cutoff <- log(pip.cutoff) - log(1 - pip.cutoff)
fgwas.tab <- read_tsv(fgwas.file, col_types = 'iiiiiicidddicc')

n.sig.snps <- fgwas.tab %>% filter(lodds > lodds.cutoff) %>% nrow()
if(n.sig.snps < 1) {
    write_tsv(data.frame(), path = out.file)
    system('rm -r ' %&&% temp.dir)
    log.msg('No significant SNPs!\n')
    q()
}

dir.create(temp.dir, recursive = TRUE, showWarnings = FALSE)
plink <- subset.plink(plink.hdr, chr = chr.input,
                      plink.lb = lb.input, plink.ub = ub.input, temp.dir = temp.dir)

svd.out <- zqtl::take.ld.svd(plink$BED, options = list(eig.tol = 1e-2, do.stdize = TRUE))

x.bim <- plink$BIM %>%
    mutate(plink.pos = 1:n())

gtex.file <- 'FQTL-v6/chr' %&&% chr.input %&&% '/50/combined.txt.gz'
gtex.cols <- c('ensg', 'chr', 'tss', 'tis.idx', 'tis.name', 'tis.beta', 'tis.se', 'tis.lo',
               'snp.name', 'snp.beta', 'snp.se', 'snp.lo', 'factor', '.cutoff')

cis.dist <- 1e6
gtex.stat <- read_tsv(gtex.file, col_names = gtex.cols) %>%
    filter(chr == chr.input) %>%
        filter(tss > (lb.input - cis.dist), tss < (ub.input + cis.dist)) %>%
            filter(.cutoff == pip.cutoff)

if(nrow(gtex.stat) < 1) {
    write_tsv(data.frame(), path = out.file)
    system('rm -r ' %&&% temp.dir)
    log.msg('No GTEx SNPs!\n')
    q()
}

make.twas <- function(gg, gtex.valid.tab, gwas.z.tab, svd.out) {

    #############
    ## options ##
    #############

    n.cutoff <- 10
    n.perm <- 5e7
    n.blk <- 1024
    n.round <- ceiling(n.perm/n.blk)

    tss <- gtex.valid.tab[gg, 'tss'] %>% as.integer()
    info <- gtex.valid.tab %r% gg %>% select(ensg, factor) %>%
        rename(gtex.factor = factor)

    eqtl.tab <- (gtex.valid.tab %r% gg) %>%
        extract.gtex.v6() %>%
            mutate(eqtl.z = eqtl.beta / eqtl.se)

    matched <- gwas.z.tab %>% left_join(eqtl.tab) %>%
        mutate(eqtl.z = if_else((!is.na(eqtl.z)) & eqtl.a1 == plink.a1, -eqtl.z, eqtl.z))
    if(nrow(matched) < 1) return(NULL)

    eqtl.z.tab <- matched %>% filter(!is.na(eqtl.z))
    if(nrow(eqtl.z.tab) < 1) return(NULL)

    V.t <- svd.out$V.t %c% matched$plink.pos
    D <- svd.out$D

    V.t.obs <- svd.out$V.t %c% eqtl.z.tab$plink.pos
    gwas.z <- eqtl.z.tab$gwas.z
    eqtl.z <- eqtl.z.tab$eqtl.z

    gwas.z.tot <- matched$gwas.z
    obs.stat <- func.NWAS(eqtl.z, gwas.z, V.t.obs, D)
    blk.mat <- func.blk.ind(nrow(eqtl.z.tab), n.blk)

    z.perm.sum <- 0
    z.perm.sqsum <- 0

    ## adaptive permutation
    c.tot <- 0
    p.tot <- 0
    z.obs.abs <- abs(obs.stat[, 'z'])
    for(b in seq(1, n.round)){
        set.seed(b)

        stat.z <- func.NWAS.eqtl.perm(eqtl.z, gwas.z.tot, V.t, D, blk.mat) %>%
            select(z) %>% .unlist()

        c.tot <- c.tot + sum((abs(stat.z) + 1e-8) >= z.obs.abs)
        p.tot <- p.tot + length(stat.z)

        z.perm.sum <- z.perm.sum + sum(stat.z)
        z.perm.sqsum <- z.perm.sqsum + sum(stat.z^2)

        if(c.tot > n.cutoff){
            break;
        }
        cat('\n', c.tot, '/', p.tot, '... ')
    }

    z.perm.mean <- z.perm.sum / p.tot
    z.perm.var <- z.perm.sqsum / p.tot - z.perm.mean^2
    p.val <- (1 + c.tot)/(1 + p.tot)

    log.msg('Finished: %d, %s, %.2e\n', gg, info[1], p.val)

    ret <- bind_cols(info, obs.stat) %>%
        mutate(p.val,
               z.perm.mean,
               z.perm.sd = sqrt(z.perm.var))

    return(ret)
}

run.twas.factor <- function(ff, fgwas.tab, gwas.z.tab, plink, svd.out) {

    ## 1. generate factored summary GWAS
    fgwas.k <- fgwas.tab %>% filter(factor == ff) %>%
        select(snp.loc, a1, a2, theta, theta.se)

    effect.k <- x.bim %>% left_join(fgwas.k) %>% na.omit() %>%
        filter((a1 == plink.a2 & a2 == plink.a1) | (a1 == plink.a1 & a2 == plink.a2)) %>%
            mutate(theta = if_else(a1 == plink.a1, -theta, theta))

    ## 2. project onto all the SNPs
    .theta.z <- effect.k$theta
    .xx.k <- scale(plink$BED) %c% effect.k$plink.pos
    .xx.k[is.na(.xx.k)] <- 0
    .eta.k <- .xx.k %*% .theta.z

    .xx <- scale(plink$BED)
    .xx[is.na(.xx)] <- 0

    .zz <- t(.xx) %*% .eta.k / sqrt(nrow(.xx))
    .zz <- center(.zz)

    gwas.z.tab <-
        x.bim %>%
            mutate(gwas.z = as.numeric(.zz))

    valid.gtex.factors <- gtex.stat %>% extract.gtex.v6() %>%
        left_join(gwas.z.tab) %>% na.omit() %>%
            group_by(ensg, factor) %>%
                summarize(n.eqtl = n())

    gtex.valid <- valid.gtex.factors %>%
        left_join(gtex.stat)

    if(nrow(gtex.valid) < 1) return(NULL)

    gg.vec <- 1:nrow(gtex.valid)

    ret <- lapply(gg.vec, make.twas,
                  gtex.valid.tab = gtex.valid,
                  gwas.z.tab = gwas.z.tab,
                  svd.out = svd.out)

    ret <- ret %>% bind_rows() %>%
        mutate(gwas.factor = ff)
}

valid.gwas.factors <- fgwas.tab %>% filter(lodds > lodds.cutoff) %>%
    select(factor) %>% unique() %>% .unlist()

if(length(valid.gwas.factors) < 1) {
    write_tsv(data.frame(), path = out.file)
    system('rm -r ' %&&% temp.dir)
    log.msg('No GTEx factors!\n')
    q()
}

out.tab <- lapply(valid.gwas.factors,
                  run.twas.factor,
                  fgwas.tab = fgwas.tab,
                  gwas.z.tab = gwas.z.tab,
                  plink = plink,
                  svd.out = svd.out)

out.tab <- out.tab %>% bind_rows()

if(nrow(out.tab) > 0) {
    out.tab <- out.tab %>%
        mutate(chr = chr.input, LB = lb.input, UB = ub.input, ld.idx = ld.idx)
}

write_tsv(out.tab, path = out.file)
system('rm -r ' %&&% temp.dir)
log.msg('Successfully finished!\n')
