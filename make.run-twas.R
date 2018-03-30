
options(stringsAsFactors = FALSE)
source('util.R')
library(readr)
library(dplyr)
library(zqtl)

ld.idx <- 2
eqtl.dir <- 'eqtl/mayo/'
fgwas.dir <- 'result/fgwas_nn/50_5/'
pip.cutoff <- 0.1
gwas.plink.hdr <- '1KG_EUR/chr'
eqtl.plink.hdr <- 'Mayo_GENO/mayo_imp_'
out.hdr <- 'temp'

ld.file <- 'ldblocks/EUR/fourier_ls-all.bed'

ld.tab <- read_tsv(ld.file)[ld.idx, ]

chr.input <- ld.tab$chr %>% gsub(pattern = 'chr', replacement = '') %>% as.integer()
lb.input <- ld.tab$start %>% as.integer()
ub.input <- ld.tab$stop %>% as.integer()

## run TWAS
eqtl.file <- eqtl.dir %&&% '/' %&&% ld.idx %&&% '_qtl.txt.gz'
fgwas.file <- fgwas.dir %&&% '/' %&&% chr.input %&&% '/' %&&% ld.idx %&&% '.snp-factor.gz'
lodds.cutoff <- log(pip.cutoff) - log(1 - pip.cutoff)

## Read fGWAS multivariate statistics
fgwas.tab <- read_tsv(fgwas.file, col_types = 'iiicidddicc')

n.sig.snps <- fgwas.tab %>% filter(lodds > lodds.cutoff) %>% nrow()

if(n.sig.snps < 1) {
    q()
}

################################################################
## Take PLINK on both sides
temp.dir <- system('mkdir -p /broad/hptmp/ypp/ukbb-twas/' %&&%
                   out.hdr %&&% '; mktemp -d /broad/hptmp/ypp/ukbb-twas/' %&&%
                   out.hdr %&&% '/temp.XXXXXXXX',
                   intern = TRUE,
                   ignore.stderr = TRUE)

dir.create(temp.dir, recursive = TRUE, showWarnings = FALSE)

plink.gwas <- subset.plink(gwas.plink.hdr %&&% chr.input, chr = chr.input,
                           plink.lb = lb.input, plink.ub = ub.input, temp.dir = temp.dir)

plink.eqtl <- subset.plink(eqtl.plink.hdr %&&% chr.input, chr = chr.input,
                           plink.lb = lb.input, plink.ub = ub.input, temp.dir = temp.dir)

plink.matched <- match.plink(plink.gwas, plink.eqtl)

plink.gwas <- plink.matched$gwas
plink.eqtl <- plink.matched$qtl
rm(plink.matched)
gc()

################################################################
## Read eQTL statistics and measure basic statistics
## chr rs snp.loc med.id qtl.a1 qtl.a2 qtl.beta qtl.z
## i   c  i       c      c      c      d        d
eqtl.tab <- read_tsv(eqtl.file, col_types = 'icicccdd')

## Find matched statistics only including fgwas factors with lodds > lodds.cutoff
inc.factors <- fgwas.tab %>% filter(lodds > lodds.cutoff) %>%
    select(factor) %>% unique() %>%
        unlist(use.names = FALSE)

matched.stat <- match.allele(fgwas.tab %>% filter(factor %in% inc.factors), plink.gwas, eqtl.tab)

.factor <- inc.factors[1]

.stat <- matched.stat %>% filter(factor == .factor)




## X.gwas <- plink$BED %c% matched.stat
## X.eqtl <- 


## check if eqtl data matches well with GWAS



## TODO: Reconcile two plink file sets




## expr




## Statistics = theta_ukbb * 

