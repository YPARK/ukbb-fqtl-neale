#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) != 7) q()

in.dir <- argv[1]               # e.g., in.dir = 'data/1/1'
plink.hdr <- argv[2]            # e.g., plink.hdr = '1KG_EUR/chr1'
K.max <- as.integer(argv[3])    # e.g., K.max = 50
NON.NEG <- as.logical(argv[4])  # e.g., NON.NEG = FALSE
conf.file <- argv[5]            # e.g., conf.file = 'result/conf/50/loco_1.gz'
pheno.file <- argv[6]           # e.g., 'phenotypes/ukbb_pheno_Q.txt'
out.hdr <- argv[7]              # e.g., out.hdr = 'temp'

options(stringsAsFactors = FALSE)
source('util.R')
source('util.fgwas.R')
library(dplyr)
library(tidyr)
library(readr)
library(zqtl)

ld.idx <- in.dir %>% strsplit(split = '/') %>%
    (function(x) x[[1]]) %>%
        (function(x) x[length(x)]) %>% as.integer()

if(is.null(ld.idx)) {
    log.msg('Invalid in.dir : %s\n', in.dir)
    q()
}

dir.create(dirname(out.hdr), recursive = TRUE, showWarnings = FALSE)
snp.factor.file <- out.hdr %&&% '.snp-factor.gz'
trait.factor.file <- out.hdr %&&% '.trait-factor.gz'
zscore.file <- out.hdr %&&% '.zscore.gz'
var.file <- out.hdr %&&% '.var.gz'
conf.out.file <- out.hdr %&&% '.conf.gz'

.files <- c(snp.factor.file, trait.factor.file, zscore.file, var.file, conf.out.file)

if(all(sapply(.files, file.exists))) {
    log.msg('Files exists:\n%s\n', paste(.files, collapse = '\n'))
    q()
}

traits <- read_tsv(pheno.file) %>%
    select(Field.code) %>%
        unlist(use.names = FALSE) %>%
            unique() %>% sort()

stat.files <- in.dir %&&% '/' %&&% traits %&&% '.stat.gz'

if(!all(sapply(stat.files, file.exists))) {
    temp <- paste(stat.files[!sapply(stat.files, file.exists)], collapse = ', ')
    log.msg('Missing stat files : %s\n', temp)
    q()
}

LD.FILE <- 'ldblocks/EUR/fourier_ls-all.bed'
ld.tab <- read_tsv(LD.FILE) %>%
    rename(LB = start, UB = stop) %>%
    mutate(chr = gsub(chr, pattern = 'chr', replacement = ''))
ld.info <- ld.tab %r% ld.idx

################################################################
## read plink data
temp.dir <- system('mkdir -p /broad/hptmp/ypp/ukbb-fgwas/' %&&% out.hdr %&&%
                   '; mktemp -d /broad/hptmp/ypp/ukbb-fgwas/' %&&% out.hdr %&&%
                   '/temp.XXXXXXXX',
                   intern = TRUE,
                   ignore.stderr = TRUE)

dir.create(temp.dir, recursive = TRUE, showWarnings = FALSE)
plink <- subset.plink(plink.hdr, ld.info$chr, ld.info$LB, ld.info$UB, temp.dir)
system('[ -d ' %&&% temp.dir %&&% ' ] && rm -r ' %&&% temp.dir)

plink.fam <- plink$FAM
stat.tab <- read.stat.files(stat.files, plink)

if(nrow(stat.tab) < 1) {
    log.msg('Empty stat\n\n')
    write_tsv(data.frame(), path = gzfile(var.file))
    write_tsv(data.frame(), path = gzfile(conf.out.file))
    write_tsv(data.frame(), path = gzfile(snp.factor.file))
    write_tsv(data.frame(), path = gzfile(trait.factor.file))
    write_tsv(data.frame(), path = gzfile(zscore.file))
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

data <- construct.data.matrix(stat.tab, plink, traits)
gc()

if(is.null(data)) {
    log.msg('Too few number of rsids\n\n')
    write_tsv(data.frame(), path = gzfile(var.file))
    write_tsv(data.frame(), path = gzfile(conf.out.file))
    write_tsv(data.frame(), path = gzfile(snp.factor.file))
    write_tsv(data.frame(), path = gzfile(trait.factor.file))
    write_tsv(data.frame(), path = gzfile(zscore.file))
    q()
}

log.msg('Constructed data\n\n')

################################################################
## non-genetic confounder
if(!file.exists(conf.file)) {
    log.msg('Confounder file does not exist: %s\n', conf.file)
    q()
}

conf.tab <- read_tsv(conf.file)

C <- plink.fam %>% select(iid) %>%
    left_join(conf.tab) %>%
        select(-iid) %>%
            as.matrix()

stopifnot(all(!is.na(C)))

n <- nrow(data$X)
z.conf.0 <- t(data$X) %*% matrix(1, n, 1) / sqrt(n)

log.msg('Found %d covaraites\n\n', ncol(C))
Z.conf <- t(data$X) %*% C / sqrt(nrow(data$X))
Z.conf <- cbind(Z.conf, z.conf.0)

################################################################
K <- max(min(c(length(traits) - 1, ncol(data$X) - 1, K.max)), 1)

vb.opt <- list(pi.ub = -1/2, pi.lb = -3, tau = -5, do.hyper = TRUE,
               do.stdize = TRUE, eigen.tol = 1e-2, gammax = 1e4,
               svd.init = TRUE, do.rescale = TRUE,
               jitter = 0.1, vbiter = 5000, rate = 1e-2, decay = -1e-2,
               tol = 0, right.nn = NON.NEG, k = K)

z.out <- fit.zqtl(effect = data$beta, effect.se = data$se,
                  X = data$X, C.delta = Z.conf, factored = TRUE,
                  options = vb.opt)

log.msg('Finished fQTL estimation\n\n')

LD.info <- data$snps %>%
    summarize(chr = min(chr),
              snp.lb = min(snp.loc),
              snp.ub = max(snp.loc),
              ld.idx,
              LB = ld.info$LB,
              UB = ld.info$UB)

################################################################
var.tab <- estimate.variance(z.out, data$se)
var.tab <- data.frame(LD.info, var.tab)

right.tab <- melt.effect(z.out$param.right, traits, 1:K) %>%
    rename(theta.se = theta.var, trait = row, factor = col) %>%
        mutate(theta.se = signif(sqrt(theta.se), 2)) %>%
            as.data.frame()

right.tab <- data.frame(LD.info, right.tab)

left.tab <- melt.effect(z.out$param.left, data$zz$rs, 1:K) %>%
    rename(theta.se = theta.var, rs = row, factor = col) %>%
        mutate(theta.se = signif(sqrt(theta.se), 2)) %>%
            left_join(data$zz %>% select(rs, snp.loc, a1, a2)) %>%
                as.data.frame()

if(nrow(left.tab) > 0){
    left.tab <- data.frame(LD.info, left.tab)
} else {
    left.tab <- data.frame()
}

n.c <- ncol(Z.conf)

conf.assoc.tab <- melt.effect(z.out$conf.delta, 1:n.c, traits) %>%
    rename(theta.se = theta.var, conf = row, trait = col) %>%
        mutate(theta.se = signif(sqrt(theta.se), 2)) %>%
            as.data.frame()

write_tsv(conf.assoc.tab, path = gzfile(conf.out.file))
write_tsv(var.tab, path = gzfile(var.file))
write_tsv(left.tab, path = gzfile(snp.factor.file))
write_tsv(right.tab, path = gzfile(trait.factor.file))
write_tsv(data$zz, path = gzfile(zscore.file))

log.msg('Successfully finished fQTL\n\n')
