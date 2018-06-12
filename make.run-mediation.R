#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = FALSE)
source('util.R')
source('util.cammel.R')
library(readr)
library(dplyr)
library(tidyr)
library(zqtl)

ld.idx <- as.integer(argv[1])       # e.g., ld.idx = 2
fgwas.hdr <- argv[2]                # e.g., fgwas.hdr = 'result/fgwas_nn/100/CC'
gwas.stat.hdr <- argv[3]            # e.g., 'gwas_stat/pgc_mdd'
gammax.input <- as.numeric(argv[4]) # e.g., gammax.input = 1e4
eig.tol <- as.numeric(argv[5])      # e.g., eig.tol = 1e-2
out.file <- argv[6]                 # e.g., out.file = 'temp.txt.gz'

################################################################
ld.info <- read.ld.info(ld.idx)
chr.input <- ld.info$chr.input
lb.input <- ld.info$lb.input
ub.input <- ld.info$ub.input
gwas.name <- gwas.stat.hdr %>% basename()

################################################################
gwas.stat.file <- gwas.stat.hdr %&&% '_' %&&% ld.idx %&&% '.txt.gz'
gwas.tab <- read.gwas(gwas.stat.file)

fgwas.file <- fgwas.hdr %&&% '/' %&&% chr.input %&&% '/' %&&% ld.idx %&&% '.snp-factor.gz'
fgwas.tab <- read_tsv(fgwas.file, col_types = 'iiiiiicidddicc')

plink.hdr <- '1KG_EUR/chr' %&&% chr.input
temp.dir <- system('mkdir -p /broad/hptmp/ypp/ukbb-twas/' %&&%
                   out.file %&&% '; mktemp -d /broad/hptmp/ypp/ukbb-twas/' %&&%
                   out.file %&&% '/temp.XXXXXXXX',
                   intern = TRUE,
                   ignore.stderr = TRUE)
dir.create(temp.dir, recursive = TRUE, showWarnings = FALSE)
plink <- subset.plink(plink.hdr, chr = chr.input,
                      plink.lb = lb.input, plink.ub = ub.input, temp.dir = temp.dir)
system('[ ' %&&% temp.dir %&&% ' ] && rm -r ' %&&% temp.dir)

fgwas.med.tab <- fgwas.tab %>%
    rename(med.id = factor, qtl.a1 = a1, qtl.a2 = a2, qtl.beta = theta, qtl.se = theta.se) %>%
        select(chr, rs, snp.loc, med.id, qtl.a1, qtl.a2, qtl.beta, qtl.se, lodds) %>%
            mutate(qtl.z = qtl.beta / qtl.se)

.matched <- match.allele(gwas.tab, plink, fgwas.med.tab)
.data <- .matched %>% make.zqtl.data()

vb.opt <- list(pi.ub = -1/2, pi.lb = -2, tau = -5, do.hyper = TRUE, tol = 1e-8,
               gammax = gammax.input, nsingle = 100,
               vbiter = 3500, do.stdize = TRUE, eigen.tol = eig.tol,
               rate = 1e-2, decay = -1e-2, nsample = 11, print.interv = 500,
               weight = FALSE, do.rescale = TRUE,
               multivar.mediator = TRUE)

z.out <- .data %>%
    run.cammel(xx.gwas = plink$BED, xx = plink$BED, opt = vb.opt)

var.tab <- get.var.tab(z.out$var.decomp, .data$mediators) %>%
    select(med.id, var.mediated, var.direct.tot) %>%
        mutate(med.id = as.integer(med.id))

summary.tab <- .matched %>%
    na.omit() %>%
        group_by(med.id) %>%
            slice(which.max(abs(qtl.z))) %>%
                select(med.id, gwas.p, gwas.z) %>%
                    as.data.frame() %>%
                        mutate(med.id = as.integer(med.id))

out.tab <- melt.effect(z.out$param.mediated, .data$mediators, gwas.name) %>%
    rename(med.id = Var1, gwas = Var2) %>%
        mutate(med.id = as.integer(med.id)) %>%
            left_join(var.tab) %>%
                left_join(summary.tab)

out.tab <- out.tab %>%
    mutate(ld.idx = ld.idx, gwas.p.ld = min(.matched$gwas.p)) %>%
        rename(factor = med.id)

write_tsv(out.tab, out.file)
