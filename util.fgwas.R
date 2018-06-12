
################################################################
## read statistics
read.gwas <- function(stat.file) {
    require(dplyr)
    require(readr)
    stat.cols <- c('chr', 'rs', 'snp.loc', 'a1', 'a2', 'beta', 'se')
    stat.types <- 'iciccdd'

    trait <- basename(stat.file) %>%
        gsub(pattern='.stat.gz', replacement = '')

    ret <- read_tsv(stat.file, col_names = stat.cols, col_types = stat.types, skip = 1) %>%
        group_by(chr, rs, snp.loc, a1, a2) %>%
            summarize(beta = mean(beta), se = mean(se)) %>%
                mutate(trait = trait) %>%
                    as.data.frame()
    return(ret)
}

read.gwas.files <- function(stat.files, plink.obj) {

    plink.snps <- plink.obj$BIM %>% select(-rs) %>% mutate(plink.pos = 1:n())

    ## match effect directions
    stat.tab <- bind_rows(lapply(stat.files, read.gwas)) %>%
        left_join(plink.snps) %>% na.omit() %>%
            filter((a1 == plink.a2 & a2 == plink.a1) | (a1 == plink.a1 & a2 == plink.a2)) %>%
                mutate(beta = if_else(a1 == plink.a1, -beta, beta)) %>%
                    select(-a1, -a2, -missing)
}

read.conf.ind <- function(data.hdr, lodds.cutoff = 0) {
    conf.ind.file <- data.hdr %&&% '.conf-ind.gz'
    conf.loading.file <- data.hdr %&&% '.conf-trait.gz'
    conf.lodds.file <- data.hdr %&&% '.conf-lodds.gz'
    lodds <- read_tsv(conf.lodds.file)
    ret <- read_tsv(conf.ind.file)

    if(nrow(lodds) < 1) {
        ret <- ret %>% select(iid)
        return(ret)
    }

    lodds <- lodds %>%
        filter(conf.lodds > lodds.cutoff)

    if(nrow(lodds) < 1) {
        ret <- ret %>% select(iid)
        return(ret)
    }
    valid.factors <- lodds %>% select(conf.factor) %>% .unlist() %>%
        (function(x) 'factor.' %&&% x)
    ret <- ret %>% select_(.dots = c('iid', valid.factors))
    return(ret)
}

subset.plink <- function(plink.hdr, chr, plink.lb, plink.ub, temp.dir, maf = 0.05) {
    require(zqtl)
    require(dplyr)

    .error <- function(e) {
        log.msg('No QTL here!\n')
        return(NULL)
    }

    .subset <- function(plink.hdr, chr, plink.lb, plink.ub, temp.dir) {

        chr.num <- gsub(pattern = 'chr', replacement = '', chr) %>% as.integer()
        plink.cmd <- sprintf('./bin/plink --bfile %s --make-bed --geno 0.05 --maf %.2f --chr %d --from-bp %d --to-bp %d --out %s', plink.hdr, maf, chr.num, plink.lb, plink.ub, glue(temp.dir, '/plink'))
        system(plink.cmd)

        plink <- read.plink(glue(temp.dir, '/plink'))
        colnames(plink$BIM) <- c('chr', 'rs', 'missing', 'snp.loc', 'plink.a1', 'plink.a2')
        colnames(plink$FAM) <- c('fam', 'iid', 'father', 'mother', 'sex.code', '.pheno')
        plink$FAM <- plink$FAM %>% mutate(iid = sapply(iid, gsub, pattern = 'GTEX-', replacement = ''))

        if(any(is.logical(plink$BIM$plink.a1))) {
            plink$BIM$plink.a1 <- 'T'
        }

        if(any(is.logical(plink$BIM$plink.a2))) {
            plink$BIM$plink.a2 <- 'T'
        }
        return(plink)
    }

    plink <- tryCatch(.subset(plink.hdr, chr, plink.lb, plink.ub, temp.dir),
                      error = .error)
    return(plink)
}

construct.data.matrix <- function(stat.tab, plink.obj, traits) {

    require(dplyr)
    require(tidyr)

    ## construct beta and standard error matrix
    beta.tab <- stat.tab %>% dplyr::select(chr, rs, snp.loc, plink.pos, beta, trait) %>%
        spread(key = trait, value = beta) %>% arrange(snp.loc) %>% na.omit() %>%
            as.data.frame()

    se.tab <- stat.tab %>% dplyr::select(chr, rs, snp.loc, plink.pos, se, trait) %>%
        spread(key = trait, value = se) %>% arrange(snp.loc) %>% na.omit() %>%
            as.data.frame()

    zscore.tab <- stat.tab %>% dplyr::select(chr, rs, snp.loc, plink.pos, beta, se, trait) %>%
        mutate(z = signif(beta / se, 2)) %>% dplyr::select(-beta, -se) %>%
            spread(key = trait, value = z) %>% arrange(snp.loc) %>% na.omit() %>%
                as.data.frame()

    if(nrow(beta.tab) < 2) {
        log.msg('Too few number of rs\n\n')
        return(NULL)
    }

    X <- plink.obj$BED %c% beta.tab$plink.pos %>% scale()
    X[is.na(X)] <- 0

    plink.snps <- plink.obj$BIM %>% select(-rs) %>% mutate(plink.pos = 1:n()) %r% beta.tab$plink.pos

    tab2mat <- function(tab) {
        .loc <- match(traits, names(tab))
        ret <- tab %c% .loc  %>% as.matrix()
        return(ret)
    }

    beta.mat <- plink.snps %>% dplyr::select(chr, snp.loc) %>% left_join(beta.tab) %>% na.omit() %>%
        tab2mat()

    se.mat <- plink.snps %>% dplyr::select(chr, snp.loc) %>% left_join(se.tab) %>% na.omit() %>%
        tab2mat()

    zscore.tab <- plink.snps %>% select(chr, snp.loc) %>% left_join(zscore.tab) %>% na.omit()

    ret <- list(beta = beta.mat, se = se.mat, snps = plink.snps, X = X, zz = zscore.tab)
    return(ret)
}

estimate.variance <- function(z.out, se.mat) {

    require(dplyr)
    require(Matrix)
    require(methods)

    ## Variance of each factor
    theta.snp <- z.out$param.left$theta
    theta.trait <- z.out$param.right$theta

    take.var.k <- function(k) {

        theta.k <- (theta.snp %c% k) %*% t(theta.trait %c% k)
        theta.k <- theta.k * se.mat ## scale by SE
        eta.k <- sweep(z.out$Vt %*% theta.k, 1, sqrt(z.out$D2), `*`)
        var.k <- apply(eta.k^2, 2, sum)

        return(data.frame(trait = as.character(names(var.k)), factor = k, var = var.k))
    }

    theta.tot <- (theta.snp %*% t(theta.trait)) * se.mat ## scale by SE
    eta.tot <- sweep(z.out$Vt %*% theta.tot, 1, sqrt(z.out$D2), `*`)
    var.tot <- data.frame(trait = as.character(traits),
                          factor = 'total',
                          var = apply(eta.tot^2, 2, sum))
    var.tab <- bind_rows(lapply(1:K, take.var.k))
    rownames(var.tab) <- NULL
    var.tab <- rbind(var.tab, var.tot) %>%
        mutate(var = signif(var, 4))

    log.msg('Calculated Variance\n\n')
    return(var.tab)
}
