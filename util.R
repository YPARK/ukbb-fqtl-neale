options(stringsAsFactors = FALSE)

log.msg <- function(...) {
    cat('[', date() ,']', sprintf(...), file = stderr())
}

`%&&%` <- function(a, b) paste(a, b, sep = '')

glue <- function(...) paste(..., sep = '')

`%c%` <- function(a, b) a[, b, drop = FALSE]

`%r%` <- function(a, b) a[b, , drop = FALSE]

.NA <- function(nrow, ncol) {
    matrix(NA, nrow, ncol)
}

.rnorm <- function(nrow, ncol) {
    matrix(rnorm(nrow * ncol), nrow, ncol)
}

.sample <- function(nrow, ncol, pop = c(1,-1)) {
    matrix(sample(pop, nrow * ncol, TRUE), nrow, ncol)
}

melt.effect <- function(effect.obj, .rnames, .cnames) {
    library(dplyr)
    library(tidyr)
    .melt.effect <- function(mat.list, val.names, i) {
        mat <- signif(mat.list[[i]], 4)
        val <- val.names[[i]]
        colnames(mat) <- .cnames
        ret <- mat %>% as.data.frame() %>% mutate(row = .rnames) %>%
            gather_(key = 'col', value = 'value', gather_cols = .cnames) %>%
                rename_(.dots = setNames('value', val)) %>%
                    as.data.frame()
        return(ret)
    }

    melt.list <- lapply(seq_along(effect.obj), .melt.effect,
                        mat.list = effect.obj,
                        val.names = names(effect.obj))

    ret <- Reduce(function(...) left_join(..., by = c('row', 'col')), melt.list) %>%
        as.data.frame()
    return(ret)
}

################################################################
subset.plink <- function(plink.hdr, chr, plink.lb, plink.ub, temp.dir) {
    require(zqtl)
    require(dplyr)

    .error <- function(e) {
        print(e)
        log.msg('Failed to read plink!\n')
        return(NULL)
    }

    .subset <- function(plink.hdr, chr, plink.lb, plink.ub, temp.dir) {

        chr.num <- gsub(pattern = 'chr', replacement = '', chr) %>% as.integer()
        plink.cmd <- sprintf('./bin/plink --bfile %s --make-bed --geno 0.05 --maf 0.05 --chr %d --from-bp %d --to-bp %d --out %s', plink.hdr, chr.num, plink.lb, plink.ub, glue(temp.dir, '/plink'))
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

################################################################
## genrate theoretical null data

make.zqtl.null <- function(X, beta.mat, se.mat, eig.tol, is.indep = FALSE, stdize = TRUE) {
    require(zqtl)

    .rnorm <- function(nrow, ncol) matrix(rnorm(nrow * ncol), nrow, ncol)

    n.traits <- ncol(se.mat)

    ## Estimate LD matrix
    svd.out <- take.ld.svd(X, options = list(do.stdize = stdize, eigen.tol = eig.tol))
    
    if(is.indep || n.traits < 2) {
        ## R = V' D2 V
        ## (a) sample theta.tilde ~ N(0, D^2) <=> D * N(0, I)
        ## (b) sample beta.hat ~ se * V * theta.tilde
        d <- nrow(svd.out$D)
        theta.tilde <- sweep(.rnorm(d, n.traits), 1, svd.out$D, `*`)
        beta.null.mat <- (t(svd.out$V.t) %*% theta.tilde) * se.mat

    } else {        
        
        Vd <- sweep(t(svd.out$V.t), 2, svd.out$D, `*`) # R = Vd * Vd'
        d <- ncol(Vd)
        
        ## Estimate correlation between z-scores
        z.mat <- beta.mat / se.mat
        z.cov <- cor(z.mat, use = 'pairwise.complete.obs') # scale
        L <- chol(z.cov, pivot = FALSE) # C = L' * L

        z.null <- Vd %*% .rnorm(d, n.traits) %*% L
        beta.null.mat <- z.null * se.mat
    }

    for(k in 1:n.traits) {
        .val <- sort(beta.mat[, k])
        .o <- order(beta.null.mat[, k])
        beta.null.mat[.o, k] <- .val
    }

    return(beta.null.mat)
}

################################################################
fast.cov <- function(x, y) {
    n.obs <- crossprod(!is.na(x), !is.na(y))
    ret <- crossprod(replace(x, is.na(x), 0),
                     replace(y, is.na(y), 0)) / n.obs
    return(ret)
}

fast.z.cov <- function(x, y) {
    n.obs <- crossprod(!is.na(x), !is.na(y))
    ret <- crossprod(replace(x, is.na(x), 0),
                     replace(y, is.na(y), 0)) / sqrt(n.obs)
    return(ret)
}

## calculate univariate effect sizes and p-values
calc.qtl.stat <- function(xx, yy) {

    require(dplyr)
    require(tidyr)

    .xx <- scale(xx)
    .yy <- scale(yy)

    ## cross-product is much faster than covariance function
    n.obs <- crossprod(!is.na(.xx), !is.na(.yy))
    beta.mat <- crossprod(.xx %>% rm.na.zero(), .yy %>% rm.na.zero()) / n.obs

    log.msg('Computed cross-products')

    ## residual standard deviation
    resid.se.mat <- matrix(NA, ncol(.xx), ncol(.yy))

    for(k in 1:ncol(.yy)) {

        beta.k <- beta.mat[, k]
        yy.k <- .yy[, k]
        err.k <- sweep(sweep(.xx, 2, beta.k, `*`), 1, yy.k, `-`)
        se.k <- apply(err.k, 2, sd, na.rm = TRUE)

        log.msg('Residual on the column %d', k)
        resid.se.mat[, k] <- se.k + 1e-8
    }

    ## organize as consolidated table
    y.cols <- 1:ncol(yy)
    colnames(beta.mat) <- y.cols
    colnames(n.obs) <- y.cols
    colnames(resid.se.mat) <- y.cols

    beta.tab <- beta.mat %>%
        as.data.frame() %>%
            dplyr::mutate(x.col = 1:n()) %>%
                tidyr::gather(key = 'y.col', value = 'beta', y.cols)
    
    resid.se.tab <- resid.se.mat %>%
        as.data.frame() %>%
            dplyr::mutate(x.col = 1:n()) %>%
                tidyr::gather(key = 'y.col', value = 'resid.se', y.cols)
    
    nobs.tab <- n.obs %>%
        as.data.frame() %>%
            dplyr::mutate(x.col = 1:n()) %>%
                tidyr::gather(key = 'y.col', value = 'n', y.cols)
    
    out.tab <- beta.tab %>%
        left_join(nobs.tab) %>%
            left_join(resid.se.tab) %>%
                dplyr::mutate(se = resid.se/sqrt(n)) %>%
                    dplyr::mutate(p.val = zscore.pvalue(beta/se))
    
    out.tab <- out.tab %>%
        mutate(x.col = as.integer(x.col)) %>%
            mutate(y.col = as.integer(y.col))

    return(out.tab)
}

fast.cor <- function(x, y) {
    x.sd <- apply(x, 2, sd, na.rm = TRUE)
    y.sd <- apply(y, 2, sd, na.rm = TRUE)
    ret <- fast.cov(scale(x, scale = FALSE), scale(y, scale = FALSE))
    ret <- sweep(sweep(ret, 1, x.sd, `/`), 2, y.sd, `/`)
    return(ret)
}

################################################################
## Find most correlated (including zero values)
find.cor.idx <- function(Y1, Y0, n.ctrl, p.val.cutoff = 1) {

    colnames(Y1) <- 1:ncol(Y1)
    colnames(Y0) <- 1:ncol(Y0)

    require(dplyr)

    y01.stat <- get.marginal.qtl(Y0, Y1) %>%
        dplyr::rename(y0 = snp, y1 = gene) %>%
            dplyr::mutate(p.val = 2 * pnorm(abs(beta.z), lower.tail = FALSE)) %>%
                dplyr::filter(p.val < p.val.cutoff)

    ret <- y01.stat %>% dplyr::group_by(y1) %>%
        dplyr::top_n(n = -n.ctrl, wt = p.val)

    return(ret$y0)
}

################################################################
require(grid)
require(gridExtra)
require(gtable)
require(ggplot2)

match.widths.grob <- function(g.list) {

    max.width <- g.list[[1]]$widths[2:7]

    for(j in 2:length(g.list)) {
        max.width <- grid::unit.pmax(max.width, g.list[[j]]$widths[2:7])
    }

    for(j in 1:length(g.list)) {
        g.list[[j]]$widths[2:7] <- as.list(max.width)
    }
    return(g.list)
}

match.widths <- function(p.list) {
    g.list <- lapply(p.list, ggplotGrob)
    return(match.widths.grob(g.list))
}

grid.vcat <- function(p.list, ...) {
    g.list <- match.widths(p.list)
    ret <- grid.arrange(grobs = g.list, ncol = 1, newpage = FALSE, ...)
    return(ret)
}

match.heights.grob <- function(g.list, stretch = TRUE)  {
    max.height <- g.list[[1]]$heights[2:7]

    if(stretch) {
        for(j in 2:length(g.list)) {
            max.height <- grid::unit.pmax(max.height, g.list[[j]]$heights[2:7])
        }
    }

    for(j in 1:length(g.list)) {
        g.list[[j]]$heights[2:7] <- as.list(max.height)
    }

    return(g.list)
}

match.heights <- function(p.list, stretch = FALSE) {
    g.list <- lapply(p.list, ggplotGrob)
    return(match.heights(g.list, stretch))
}

grid.hcat <- function(p.list, ...) {
    g.list <- match.heights(p.list, stretch = TRUE)
    ret <- grid.arrange(grobs = g.list, nrow = 1, newpage = FALSE, ...)
    return(ret)
}


row.order <- function(mat) {
    require(cba)
    require(proxy)

    if(nrow(mat) < 3) {
        return(1:nrow(mat))
    }

    D <- proxy::dist(mat, method = function(a,b) 1 - cor(a,b, method = 'spearman'))
    D[!is.finite(D)] <- 0
    h.out <- hclust(D)
    o.out <- cba::order.optimal(D, h.out$merge)
    return(o.out$order)
}

gg.plot <- function(...) {
    ggplot(...) + theme_bw() + theme(plot.background = element_blank(),
                                     panel.background = element_blank(),
                                     strip.background = element_blank(),
                                     legend.background = element_blank())
}

order.pair <- function(pair.tab) {

    require(tidyr)
    require(dplyr)
    M <- pair.tab %>% tidyr::spread(key = col, value = weight, fill = 0)
    ro <- row.order(M %>% dplyr::select(-row) %>% as.matrix())
    co <- row.order(t(M %>% dplyr::select(-row) %>% as.matrix()))
    cc <- colnames(M)[-1]
    rr <- M[, 1]

    list(rows = rr[ro], cols = cc[co])
}
