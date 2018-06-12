################################################################
## match allele directions in statistics
## input
## 1. plink
## 2. gwas.tab
## 3. qtl.tab
match.allele <- function(gwas.tab, plink.obj, qtl.tab) {
    require(dplyr)

    x.bim <- plink.obj$BIM %>%
        dplyr::select(-rs, -missing) %>%
            dplyr::mutate(x.pos = 1:n())

    ## include all GWAS SNPs
    snp.tab <- x.bim %>%
        dplyr::left_join(gwas.tab, by = c('chr', 'snp.loc')) %>%
            na.omit()

    snp.tab <- snp.tab %>%
        dplyr::filter((plink.a1 == gwas.a1 & plink.a2 == gwas.a2) |
                          (plink.a1 == gwas.a2 & plink.a2 == gwas.a1))

    ## match with QTL tab (eQTL or hQTL ...)
    match.med <- function(g) {
        ret.g <- snp.tab %>%
            left_join(qtl.tab %>% dplyr::select(-rs) %>% dplyr::filter(med.id == g),
                      by = c('chr', 'snp.loc'))
        missing <- is.na(ret.g$med.id)
        ret.g[missing, 'med.id'] <- g
        ret.g[missing, 'qtl.beta'] <- 0
        ret.g[missing, 'qtl.se'] <- 1e-4
        ret.g[missing, 'qtl.z'] <- 0
        ret.g[missing, 'qtl.a1'] <- ret.g[missing, 'plink.a1']
        ret.g[missing, 'qtl.a2'] <- ret.g[missing, 'plink.a2']
        return(ret.g)
    }

    mediators <- qtl.tab$med.id %>% unique()

    ret <- lapply(mediators, match.med) %>% bind_rows()

    ret <- ret %>%
        dplyr::filter((plink.a1 == qtl.a1 & plink.a2 == qtl.a2) |
                          (plink.a1 == qtl.a2 & plink.a2 == qtl.a1))

    ret <- ret %>%
        dplyr::mutate(qtl.flip = if_else(plink.a1 == qtl.a1, 1.0, -1.0)) %>%
            dplyr::mutate(gwas.flip = if_else(plink.a1 == gwas.a1, 1.0, -1.0))

    ret <- ret %>%
        dplyr::mutate(gwas.z = gwas.flip * gwas.z, gwas.beta = gwas.flip * gwas.beta)

    ret <- ret %>%
        dplyr::mutate(qtl.z = qtl.flip * qtl.z, qtl.beta = qtl.flip * qtl.beta)

    ret <- ret %>%
        dplyr::select(-gwas.flip, -qtl.flip) %>%
            dplyr::rename(a1 = plink.a1, a2 = plink.a2)

    return(ret)
}

## match evething to plink.gwas
match.plink <- function(plink.gwas, plink.qtl) {

    if(is.null(plink.gwas)) return(NULL)
    if(is.null(plink.qtl)) return(NULL)

    ret.gwas <- plink.gwas
    ret.qtl <- plink.qtl

    gwas.bim <- plink.gwas$BIM %>%
        mutate(gwas.x.pos = 1:n()) %>%
            rename(gwas.plink.a1 = plink.a1,
                   gwas.plink.a2 = plink.a2) %>%
                       select(-missing)

    qtl.bim <- plink.qtl$BIM %>%
        mutate(qtl.x.pos = 1:n()) %>%
            rename(qtl.plink.a1 = plink.a1,
                   qtl.plink.a2 = plink.a2,
                   qtl.rs = rs) %>%
                       select(-missing)

    bim.matched <- gwas.bim %>%
        left_join(qtl.bim) %>%
            na.omit()

    if(nrow(bim.matched) < 1) return(NULL)

    bim.matched <- bim.matched %>%
        dplyr::filter(((gwas.plink.a1 == qtl.plink.a1) & (gwas.plink.a2 == qtl.plink.a2)) |
                          ((gwas.plink.a2 == qtl.plink.a1) & (gwas.plink.a1 == qtl.plink.a2))) %>%
                              arrange(chr, snp.loc)

    if(nrow(bim.matched) < 1) return(NULL)

    ret.gwas$BIM <- ret.gwas$BIM[bim.matched$gwas.x.pos, , drop = FALSE]
    ret.gwas$BED <- ret.gwas$BED[ , bim.matched$gwas.x.pos, drop = FALSE]

    ret.qtl$BIM <- ret.qtl$BIM[bim.matched$qtl.x.pos, , drop = FALSE]
    ret.qtl$BED <- ret.qtl$BED[ , bim.matched$qtl.x.pos, drop = FALSE]

    flip.tab <- ret.gwas$BIM %>% mutate(gwas.x.pos = 1:n()) %>%
        left_join(ret.qtl$BIM %>% mutate(qtl.x.pos = 1:n()),
                  by = c('chr', 'snp.loc'),
                  suffix = c('.gwas', '.qtl')) %>%                      
                      filter(plink.a1.gwas != plink.a1.qtl)

    ret.qtl$BIM[flip.tab$qtl.x.pos, ] <- ret.gwas$BIM[flip.tab$gwas.x.pos, ]

    flip.bed <- ret.qtl$BED[, flip.tab$qtl.x.pos]
    zero.idx <- flip.bed <= 0.5
    two.idx <- flip.bed >= 1.5
    flip.bed[two.idx] <- 0
    flip.bed[zero.idx] <- 2
    ret.qtl$BED[, flip.tab$qtl.x.pos] <- flip.bed

    return(list(gwas = ret.gwas, qtl = ret.qtl))
}

################################################################
## generate matrices from matched statistics
make.zqtl.data <- function(matched.stat, n.permuted = 0) {

    require(dplyr)
    require(tidyr)

    if(nrow(matched.stat) < 1) {
        return(NULL)
    }

    qtl.beta <- matched.stat %>% dplyr::select(med.id, x.pos, qtl.beta) %>%
        dplyr::group_by(med.id, x.pos) %>%
            dplyr::slice(which.max(abs(qtl.beta))) %>%
                tidyr::spread(key = med.id, value = qtl.beta)

    qtl.z <- matched.stat %>% dplyr::select(med.id, x.pos, qtl.z) %>%
        dplyr::group_by(med.id, x.pos) %>%
            dplyr::slice(which.max(abs(qtl.z))) %>%
                tidyr::spread(key = med.id, value = qtl.z)

    med.id <- colnames(qtl.beta)[-1]
    x.pos <- qtl.beta$x.pos

    .temp <- matched.stat %>% dplyr::select(x.pos, gwas.beta, gwas.se) %>%
        group_by(x.pos) %>% dplyr::slice(which.max(abs(gwas.beta / gwas.se)))

    gwas.beta <- qtl.beta %>%
        dplyr::select(x.pos) %>%
            left_join(.temp %>% dplyr::select(x.pos, gwas.beta))

    gwas.se <- qtl.beta %>%
        dplyr::select(x.pos) %>%
            left_join(.temp %>% dplyr::select(x.pos, gwas.se))

    .xx <- match(x.pos, qtl.beta$x.pos)
    .mm <- match(med.id, colnames(qtl.beta)) 
    qtl.beta <- as.matrix(qtl.beta[.xx, .mm])

    .xx <- match(x.pos, qtl.z$x.pos)
    .mm <- match(med.id, colnames(qtl.z)) 
    qtl.z <- as.matrix(qtl.z[.xx, .mm])

    qtl.se <- qtl.beta / qtl.z
    qtl.se[is.na(qtl.se)] <- mean(qtl.se, na.rm = TRUE)
    qtl.se <- pmax(qtl.se, 1e-8)

    ## augment permuted QTL effects
    if(ncol(qtl.beta) < n.permuted) {
        n.med <- ncol(qtl.beta)
        n.add <- n.permuted - n.med
        idx <- sample(n.med, n.add, replace = TRUE)
        n.snp <- nrow(qtl.beta)

        qtl.se.add <- qtl.se[, idx, drop = FALSE]
        qtl.beta.add <- qtl.beta[, idx, drop = FALSE]
        for(j in 1:n.add) {
            qtl.beta.add[, j] <- sample(qtl.beta.add[, j])
        }
        med.id <- c(med.id, paste('.perm', 1:n.add, sep='.'))
        qtl.beta <- cbind(qtl.beta, qtl.beta.add)
        qtl.se <- cbind(qtl.se, qtl.se.add)
    }

    .xx <- match(x.pos, gwas.beta$x.pos)
    gwas.beta <- as.matrix(gwas.beta[.xx, 'gwas.beta'])

    .xx <- match(x.pos, gwas.se$x.pos)
    gwas.se <- as.matrix(gwas.se[.xx, 'gwas.se'])

    gwas.se[is.na(gwas.se)] <- mean(gwas.se, na.rm = TRUE)

    ret <- list(x.pos = x.pos,
                mediators = med.id,
                qtl.beta = as.matrix(qtl.beta), qtl.se = as.matrix(qtl.se),
                gwas.beta = as.matrix(gwas.beta), gwas.se = as.matrix(gwas.se))
    return(ret)
}


################################################################
## Estimate variance
estimate.variance <- function(z.out, z.data, X, nn) {

    V.t <- z.out$Vt
    Y <- z.out$Y
    M <- z.out$M
    S.inv <- z.out$S.inv.y
    S <- 1/S.inv
    S.inv.m <- z.out$S.inv.m
    D2 <- z.out$D2
    D <- sqrt(z.out$D2)
    kk <- length(D)

    Vd <- sweep(t(V.t), 2, D, `*`)
    W <- sweep(t(V.t), 2, D, `/`)

    ## Estimate total variance using Shi et al.
    eta.gwas <- t(W) %*% z.data$gwas.beta
    gg <- sum(eta.gwas^2)
    var.shi <- (gg * nn - kk) / (nn - kk)

    ## Estimated true qtl effect
    ## theta.hat        ~ S R inv(S) (aa * bb)
    ## inv(S) theta.hat ~ R inv(S) (aa * bb)
    ## -- the mediated component
    aa <- sweep(W %*% sweep(M, 1, D, `/`), 1, S, `*`)
    bb <- z.out$param.mediated$theta
    ab <- aa %*% bb
    ## -- direct effect
    theta.dir <- z.out$param.direct$theta

    xx.std <- scale(X) %>% rm.na.zero()
    n.ref <- nrow(xx.std)

    ## residual variance
    r.hat <- z.out$resid$theta
    rr <- sweep(W %*% (t(W) %*% r.hat), 1, S, `*`)
    eta.r <- t(Vd) %*% rr
    var.resid <- sum(eta.r^2)

    ## Estimate variance using reference genotype matrix
    eta.tot <- xx.std %*% (theta.dir + ab)
    eta.dir <- xx.std %*% (theta.dir)
    eta.med.tot <- xx.std %*% (ab)

    var.ref.tot <- var(eta.tot, na.rm = TRUE) + var.resid
    var.ref.dir <- var(eta.dir, na.rm = TRUE)
    var.ref.med <- var(eta.med.tot, na.rm = TRUE)

    ## Estimate variance using covariance
    eta.tot <- t(Vd) %*% (theta.dir + ab)
    var.tot <- sum(eta.tot^2)
    eta.dir <- t(Vd) %*% (theta.dir)
    var.dir <- sum(eta.dir^2)
    eta.ab <- t(Vd) %*% ab
    var.med <- sum(eta.ab^2)

    ## each mediation effect
    var.each <- function(k) {
        ab.k <- (aa %c% k) %*% (bb %r% k)
        eta.ab.k <- t(Vd) %*% ab.k
        var.k <- sum(eta.ab.k^2)
        return(var.k)
    }

    n.med <- nrow(bb)
    var.med.vec <- sapply(1:n.med, var.each)

    return(list(shi = signif(var.shi, digits = 4),
                resid = signif(var.resid, 4),
                ref.tot = signif(var.ref.tot, 4),
                ref.dir = signif(var.ref.dir, 4),
                ref.med = signif(var.ref.med, 4),
                tot = signif(var.tot, 4),
                dir = signif(var.dir, 4),
                med.tot = signif(var.med, 4),
                med = signif(var.med.vec, 4),
                k = kk))
}

################################################################
melt.effect <- function(effect.obj, .rnames, .cnames) {
    require(reshape2)
    melt.list <- lapply(seq_along(effect.obj),
                        function(mat.list, val.names, i) {
                            mat <- signif(mat.list[[i]], 4)
                            val <- val.names[[i]]
                            rownames(mat) <- .rnames
                            colnames(mat) <- .cnames
                            melt(mat, value.name = val) },
                        mat.list = effect.obj, val.names = names(effect.obj))

    ret <- Reduce(function(...) left_join(..., by = c('Var1', 'Var2')), melt.list) %>%
        mutate(Var1 = as.character(Var1), Var2 = as.character(Var2))
    return(ret)
}

mat2tab <- function(mat, val.name) {
    require(dplyr)
    require(tidyr)
    ret <- as.matrix(mat) %>% as.data.frame()
    col.names <- 1:ncol(ret)
    colnames(ret) <- col.names
    ret <- ret %>%
        mutate(x.col = 1:n()) %>%
            gather_(key_col= 'y.col', value_col = val.name, col.names)
    return(ret)
}

effect2tab <- function(param) {
    require(dplyr)
    require(tidyr)
    .theta <- mat2tab(param$theta, val.name = 'theta')
    .theta.var <- mat2tab(param$theta.var, val.name = 'theta.var')
    .lodds <- mat2tab(param$lodds, val.name = 'lodds')
    ret <- .theta %>%
        left_join(.theta.var) %>%
            left_join(.lodds) %>%
                mutate(theta.se = sqrt(theta.var)) %>%
                    select(-theta.var) %>%
                        mutate(x.col = as.integer(x.col)) %>%
                            mutate(y.col = as.integer(y.col)) %>%
                                left_join(x.bim) %>%
                                    left_join(genes.Y1)
    ret <- ret %>%
        mutate(theta = signif(theta, 4),
               theta.se = signif(theta.se, 4),
               lodds = signif(lodds, 2))
}

simplify.ensg <- function(tab) {
    require(tidyr)
    tab %>% tidyr::separate(ensg, c('ensg', '.remove'), sep = '[.]') %>%
        select(-.remove)
}

run.cammel <- function(zqtl.data, xx.gwas, xx.med, opt) {

    if(is.null(zqtl.data)) return(NULL)

    xx.gwas.mat <- xx.gwas %c% zqtl.data$x.pos %>% as.matrix()
    xx.med.mat <- xx.med %c% zqtl.data$x.pos %>% as.matrix()

    z.out <- fit.med.zqtl(zqtl.data$gwas.beta, zqtl.data$gwas.se,
                          zqtl.data$qtl.beta, zqtl.data$qtl.se,
                          X.gwas = xx.gwas.mat,
                          X.med = xx.med.mat,                          
                          options = opt)
    return(z.out)
}

run.cammel.null <- function(zqtl.data, xx.gwas, xx.med, n.null, opt) {

    xx.gwas.mat <- xx.gwas %c% zqtl.data$x.pos %>% as.matrix()
    xx.med.mat <- xx.med %c% zqtl.data$x.pos %>% as.matrix()
    ret <- NULL
    for(r in 1:n.null) {

        gwas.null <- make.zqtl.null(xx.gwas.mat, zqtl.data$gwas.se, eig.tol = eig.tol)
        qtl.null <- make.zqtl.null(xx.med.mat, zqtl.data$qtl.se, eig.tol = eig.tol)

        z.out <- fit.med.zqtl(gwas.null, zqtl.data$gwas.se,
                              qtl.null, zqtl.data$qtl.se,
                              X.gwas = xx.gwas.mat, X.med = xx.med.mat,                              
                              options = opt)

        null.stat <- melt.effect(z.out$param.mediated, zqtl.data$mediators, r) %>%
            mutate(theta.var = sqrt(theta.var)) %>% rename(theta.se = theta.var) %>%
                mutate(log.var = as.numeric(signif(log(z.out$var.decomp$var.med.each), 2)))

        print(null.stat %>% filter(lodds > 0) %>% as.data.frame())
        ret <- bind_rows(ret, null.stat)
        log.msg('\nnull round = %d / %d\n', r, n.null)
    }
    ret <- ret %>% rename(med.id = Var1, null = Var2)
    return(ret)
}

get.var.tab <- function(var.decomp, mediators) {

    if(is.null(var.decomp)) return(NULL)

    ret <- data.frame(med.id = mediators,
                      var.mediated = signif(var.decomp$var.med.each, 2),
                      var.mediated.tot = signif(var.decomp$var.med.mean, 2),
                      var.mediated.tot.se = signif(sqrt(var.decomp$var.med.var), 2),
                      var.direct.tot = signif(var.decomp$var.direct.mean, 2),
                      var.direct.tot.se = signif(sqrt(var.decomp$var.direct.var), 2))
    return(ret)
}

## gene-level QTL and GWAS stat summary
get.summary.tab <- function(.gwas.tab, .qtl.tab, .plink.obj) {
    if(is.null(.gwas.tab)) return(NULL)

    if('rs' %in% colnames(.qtl.tab)){
        .qtl.tab <- .qtl.tab %>% select(-rs)
    }

    .temp <- .gwas.tab %>%
        match.allele(plink.obj = .plink.obj, qtl.tab = .qtl.tab)

    .temp.gwas <- .temp %>% group_by(med.id) %>%
        slice(which.min(gwas.p)) %>%
            select(med.id, rs, snp.loc, gwas.p, gwas.beta, gwas.se, qtl.z)

    .temp.qtl <- .temp %>% group_by(med.id) %>%
        slice(which.max(qtl.z)) %>%
            select(med.id, rs, snp.loc, gwas.p, gwas.beta, gwas.se, qtl.z)

    ret <- .temp.gwas %>%
        left_join(.temp.qtl, by = 'med.id', suffix = c('.by.gwas', '.by.qtl'))

    ret <- ret %>% dplyr::select_(.dots = sort(names(ret)))

    return(ret)
}

get.effect.tab <- function(z.out, z.data, gwas.tab, qtl.tab, data.name, plink.obj = NULL) {
    if(is.null(z.out)) return(NULL)
    z.effect <-
        melt.effect(z.out$param.mediated, z.data$mediators, data.name) %>%
            rename(med.id = Var1, gwas = Var2) %>%
                left_join(get.var.tab(z.out$var.decomp, z.data$mediators))

    if(!is.null(plink.obj)) {
        z.effect <- z.effect %>%
            left_join(get.summary.tab(gwas.tab, qtl.tab, plink.obj))
    }
    return(z.effect)
}

## read GWAS statistics
read.gwas <- function(gwas.file, se.reg = 1e-4) {
    require(dplyr)
    require(readr)

    ## expected columns
    expected.tab <- tibble(col.name = c('chr', 'rs', 'snp.loc', 'gwas.a1', 'gwas.a2', 'gwas.beta', 'gwas.se', 'gwas.p'), col.type = c('i', 'c', 'i', 'c', 'c', 'd', 'd', 'd'))

    ## check header
    hdr <- suppressMessages(read_tsv(gwas.file, n_max = 1, col_names = FALSE)) %>% unlist()
    if(length(hdr) < 1) return(NULL)
    hdr.tab <- tibble(col.name = hdr) %>%
        mutate(col.pos = 1:n())

    expected.hdr <- hdr.tab %>% left_join(expected.tab, by = 'col.name') %>%
        mutate(col.type = if_else(is.na(col.type), '_', col.type))

    cc <- expected.hdr %>% filter(col.type != '_') %>% select(col.name) %>% unlist()
    tt <- expected.hdr %>% select(col.type) %>% unlist() %>% paste(collapse = '')
        
    if(length(cc) < nrow(expected.tab)) return(NULL)

    ret <- read_tsv(gwas.file, col_names = cc, col_types = tt, skip = 1)
    ret <- ret %>% mutate(gwas.z = gwas.beta / (gwas.se + se.reg))

    return(ret)
}

################################################################
## read multiple eQTL effect sizes
read.multivar.eqtl <- function(eqtl.data, eqtl.data.files) {
    require(dplyr)
    require(readr)

    ## expected columns
    expected.tab <- tibble(col.name = c('chr', 'rs', 'snp.loc', 'med.id', 'qtl.a1', 'qtl.a2', 'qtl.beta', 'qtl.se'), col.type = c('i', 'c', 'i', 'c', 'c', 'c', 'd', 'd'))

    ## Read QTL statistics and measure basic statistics
    ## chr rs snp.loc med.id qtl.a1 qtl.a2 qtl.beta qtl.z
    read.eqtl <- function(ii) {

        ff <- eqtl.data.files[ii]

        ## a. read header
        hdr <- suppressMessages(read_tsv(ff, n_max = 1, col_names = FALSE)) %>% unlist()
        if(length(hdr) < 1) return(NULL)
        hdr.tab <- tibble(col.name = hdr) %>%
            mutate(col.pos = 1:n())

        expected.hdr <- hdr.tab %>% left_join(expected.tab, by = 'col.name') %>%
            mutate(col.type = if_else(is.na(col.type), '_', col.type))

        cc <- expected.hdr %>% filter(col.type != '_') %>% select(col.name) %>% unlist()
        tt <- expected.hdr %>% select(col.type) %>% unlist() %>% paste(collapse = '')
        
        if(length(cc) < nrow(expected.tab)) return(NULL)

        ret <- read_tsv(ff, col_names = cc, col_types = tt, skip = 1)

        if(nrow(ret) == 0) return(NULL)

        ret <- ret %>% mutate(data = eqtl.data[ii]) %>%
            mutate(rs = chr %&&% "_" %&&% snp.loc)

        ret <- ret %>%
            select(chr, rs, snp.loc, med.id, qtl.a1, qtl.a2, qtl.beta, qtl.se)

        ret <- ret %>%
            mutate(data = eqtl.data[ii])

        return(ret)
    }

    eqtl.tab <- 1:length(eqtl.data) %>%
        lapply(FUN = read.eqtl) %>%
            bind_rows()

    if(nrow(eqtl.tab) < 1) {
        return(data.frame())
    }

    eqtl.tab <- eqtl.tab %>%
        separate(col = med.id, into = c('med.id', 'factor'), sep = '@') %>%
            mutate(factor = if_else(is.na(factor), '0', factor)) %>%
                mutate(factor = as.integer(factor)) %>%
                    separate(col = med.id, into = c('med.id', 'remove'), sep = '[.]') %>%
                        select(-remove)

    eqtl.tab <- eqtl.tab %>%
        mutate(med.id = med.id %&&% '@' %&&% factor %&&% '@' %&&% data) %>%
            select(-factor, -data)

    eqtl.tab <- eqtl.tab %>%
        mutate(qtl.z = qtl.beta / qtl.se)

    return(eqtl.tab)
}

################################################################
## Separate z-scores by mediation
separate.zscore <- function(zqtl.out, x.bim) {
    require(dplyr)

    z.obs <- t(zqtl.out$Vt) %*% zqtl.out$Y %>% as.numeric()

    z.unmed.mat <- zqtl.out$resid.Z

    z.unmed <- apply(z.unmed.mat, 1, mean, na.rm = TRUE) %>%
        matrix(ncol=1) %*%
            zqtl.out$param.unmediated$theta %>%
                as.numeric()

    z.unmed.sd <- apply(z.unmed.mat, 1, sd, na.rm = TRUE) %>%
        matrix(ncol=1) %*%
            zqtl.out$param.unmediated$theta %>%
                as.numeric()

    z.med.hat <- t(zqtl.out$Vt) %*% zqtl.out$M %*% zqtl.out$param.mediated$theta %>%
        as.numeric()

    ret <- tibble(snp.loc = x.bim$snp.loc,
                  obs = z.obs,
                  unmed = z.unmed,
                  unmed.sd = z.unmed.sd,
                  med = z.med.hat)
}
