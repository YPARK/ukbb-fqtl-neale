require(Matrix)
require(methods)

extract.gtex.v6 <- function(gtex) {
    eqtl.cols <- c('chr', 'snp.loc', 'eqtl.a1', 'eqtl.a2', 'remove')
    .split.bar <- function(s) strsplit(s, split = '[|]')
    ret <- gtex %>%
        select(ensg, factor, starts_with('snp')) %>%
            unnest(eqtl.name = .split.bar(snp.name),
                   eqtl.beta = .split.bar(snp.beta),
                   eqtl.se = .split.bar(snp.se)) %>%
                       select(ensg, factor, starts_with('eqtl')) %>%
                           separate('eqtl.name', eqtl.cols, sep = '_') %>%
                               select(-remove) %>%
                                   mutate(chr = as.integer(chr),
                                          snp.loc = as.integer(snp.loc),
                                          eqtl.beta = as.numeric(eqtl.beta),
                                          eqtl.se = as.numeric(eqtl.se))
    return(ret)
}

################################################################
## summary-based imputation
## 
## numerator :
##
##   z.GWAS (marginal)' * z.QTL (polygenic)
##
## denominator:
##
##   eta.QTL = D * V' * z.QTL
##   eta.GWAS = inv(D) * V' * z.GWAS
##
##   1. eta.QTL' * eta.QTL
##   2. eta.GWAS' * eta.GWAS   (These were too strong)
##   3. # SNPs                 (These were too strong)
##

func.NWAS <- function(.eqtl.z, .gwas.z, vt, D) {
    .num <- sum(.eqtl.z * .gwas.z)
    .eta.QTL <- sweep(sweep(vt, 2, .eqtl.z, `*`), 1, D, `*`)
    ## .eta.GWAS <- sweep(sweep(vt, 2, .gwas.z, `*`), 1, D, `/`)
    .denom <- sum(.eta.QTL^2)
    ## + sum(.eta.GWAS^2) + length(.eqtl.z)
    .sd <- pmax(sqrt(.denom), 1e-4)
    return(data.frame(z = .num/.sd, theta = .num, theta.se = .sd))
}

## Calculate SVD of reference genotype matrix X
func.LD.svd <- function(X, eig.tol = 1e-8, normalize = FALSE){
    ## use centered covariance matrix
    if(normalize){
        x.safe <- scale(X, scale = TRUE, center = TRUE)
    } else {
        x.safe <- scale(X, scale = FALSE, center = TRUE)
    }
    n.ind <- apply(!is.na(x.safe), 2, sum)
    x.safe <- sweep(x.safe, 2, sqrt(pmax(n.ind, 1)), `/`)
    x.safe[is.na(x.safe)] <- 0

    ## prevent underflow of Lapack dgesdd, (-5, 5) seems safe
    n.max <- 5 / (max(abs(x.safe)) + 1e-8)
    svd.out <- svd(x.safe * n.max)
    .valid <- which(svd.out$d^2 > eig.tol)

    ## scale back eigen values
    d <- svd.out$d[.valid] / n.max
    V.t <- t(svd.out$v[, .valid, drop = FALSE])

    log.msg('SVD finished: %d x %d V.t', dim(V.t)[1], dim(V.t)[2])

    return(list(V.t = V.t, d = d))
}

## Correct GWAS z-score bias
func.correct.bias <- function(z, V.t, D){
    V.t1 <- matrix(apply(V.t, 1, sum), ncol=1)
    num <- sweep(t(V.t), 2, D^2, `*`) %*% V.t1 * sum(z)
    denom <- sum((D * V.t1)^2)
    return(z - num/denom)
}

func.blk.ind <- function(n.qtl, n.blk) {

    ret <- sparseMatrix(i = as.integer(1:(n.qtl*n.blk)),
                        j = as.integer(as.vector(sapply(1:n.blk, function(x) rep(x, n.qtl)))),
                        x = 1)
    
    log.msg('Made %d x %d block index matrix\n', n.qtl*n.blk, n.blk)

    return(ret)
}

## faster QTL permutation
func.NWAS.eqtl.perm <- function(eqtl.z, gwas.z, V.t, D, blk.ind, DEBUG = FALSE) {

    n.blk <- dim(blk.ind)[2]
    n.tot.snp <- length(gwas.z)
    n.qtl <- length(eqtl.z)
    eqtl.idx <- sample(n.tot.snp, n.blk * n.qtl, replace = TRUE)
    
    stat.num <- sweep(matrix(gwas.z[eqtl.idx], nrow = n.qtl, byrow = FALSE), 1, eqtl.z, `*`)
    stat.num <- apply(stat.num, 2, sum)
    
    eta.EQTL.perm <- 
        sweep(sweep(V.t[, eqtl.idx, drop = FALSE], 1, D, `*`),
              2, rep(eqtl.z, n.blk), `*`)
    
    ## eta.GWAS.perm <- 
    ##     sweep(sweep(V.t[, eqtl.idx, drop = FALSE], 1, D, `/`),
    ##           2, gwas.z[eqtl.idx], `*`)
    
    stat.denom <- matrix(apply(eta.EQTL.perm^2, 2, sum), nrow = 1) %*% blk.ind
    ## apply(eta.GWAS.perm^2, 2, sum) %*% blk.ind +
    ##     n.qtl
    
    stat.num <- as.numeric(stat.num)
    stat.denom <- pmax(as.numeric(stat.denom), 1e-8)

    if(DEBUG){
        ## confirmed: slow permutation by permutation statistic
        .test <- matrix(eqtl.idx, nrow = n.qtl, byrow = FALSE)
        .func <- function(.idx) {
            func.NWAS(eqtl.z, gwas.z[.idx], V.t[, .idx, drop = FALSE], D)
        }
        test.temp <- do.call(rbind, apply(.test, 2, .func))
        cat(sum(abs(stat.num/sqrt(stat.denom) - test.temp[, 1])))
    }

    return(data.frame(z = stat.num/sqrt(stat.denom),
                      theta = stat.num,
                      theta.se = sqrt(stat.denom)))
}

################################################################

read.gene.info <- function() {

    gene.cols <- c('chr', 'lb', 'ub', 'strand', 'ensg', 'hgnc')
    gene.cols.types <- 'ciiccc_'
    gene.file <- 'data/coding.genes.txt.gz'

    gene.info <- read_tsv(gene.file, col_names = gene.cols, col_types = gene.cols.types) %>%
        mutate(tss = if_else(strand == '+', lb, ub)) %>%
            mutate(gene = 1:n())

    return(gene.info)
}
