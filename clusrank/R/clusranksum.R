cluswilcox.test.ranksum.rgl <- function(x, cluster, group, stratum,
                                         alternative, exact, mu,
                                        DNAME = NULL, METHOD = NULL) {
    clus.grp <- lapply(split(cluster, group), unique)
    l.clus.grp <- length(unlist(clus.grp))
    l.clus <- length(unique(cluster))
    temp <- merge(table(cluster), cbind(cluster, group))
    temp <- temp[order(temp$Freq), ]
    temp <- split(temp$group, temp$Freq)
    check <- lapply(temp, function(x) length(unique(x)))
    ## check with each stratum in which clusters have the same cluster size 
    ## to see if they only are only assigned a single treatment.
    if(all(unlist(check) == 1))  warning("For each of the stratum of the data when data are
                      stratified by the cluster size,
                      there should be at least one stratum with both treatments ")
    arglist <- setNames(list(x, cluster, group, stratum, alternative,
                             mu, DNAME, METHOD, exact),
                        c("x", "cluster", "group", "stratum",
                          "alternative", "mu", "DNAME", "METHOD",
                          "exact"))
    if(l.clus != l.clus.grp) {
        result <- do.call("cluswilcox.test.ranksum.rgl.sub", c(arglist))
        return(result)
    } else {
        result <- do.call("cluswilcox.test.ranksum.rgl.clus", c(arglist))
        return(result)
    }
}

cluswilcox.test.ranksum.rgl.clus <- function(x, cluster, group,
                                             stratum, alternative, exact,
                                             mu, DNAME = NULL, METHOD = NULL) {
## The input data should be already in the increasing order with cluster.
    n.obs <- length(x)
    one <- rep(1, n.obs)
    x[which(group == 1)] <- x[which(group == 1)] - mu
    xrank <- rank(x)
    rksum <-stats::aggregate(xrank ~ cluster, FUN = sum)[, 2]
    csize <-stats::aggregate(one ~ cluster, FUN = sum)[, 2]
    grp <- stats::aggregate(group ~ cluster, FUN = mean)[, 2]
    strt <- stats::aggregate(stratum ~ cluster, FUN = mean)[, 2]
    clus <- unique(cluster)
    dat <- data.frame(clus, strt, grp, csize, rksum)
    bal <- (!(length(unique(table(cluster))) != 1L))

    csize.uniq <- unique(csize)
    strt.uniq <- unique(strt)
    l.csu <- length(csize.uniq)
    l.stu <- length(strt.uniq)

    dat.l <- split(dat, strt) ## Split the rank data by stratumum
    csize.split <- function(dat) {
        split(dat, dat$"csize")
    }
    dat.l <- lapply(dat.l, csize.split)
    if(is.null(exact)) {
        exact <- FALSE
    }
    if(exact == TRUE) {
        Wc <- sum(dat[dat$grp == 1, "rksum"])
        n.layer <- l.csu * l.stu
        mgv <- ngv <- rep(0, n.layer)
        ct <- 1
        rkx <- rkxc <- as.matrix(matrix(0, Wc + 1, n.layer))
        ## Small sample permutation test
        counter <- 1
        for( i in csize.uniq) {
            for( j in strt.uniq) {
                j.ch <- as.character(j)
                i.ch <- as.character(i)
                temp <- dat.l[[j.ch]][[i.ch]] ## data with csize i and stratum j
                temp.rksum <- temp[, "rksum"]
                temp.grp <- temp[, "grp"]
                mgv[counter] <- length(temp.grp[temp.grp == 1])
                ngv[counter] <- length(temp.grp)
                temp.x <- cumcrksum(Wc, mgv[counter], sort(temp.rksum))
                rkx[, counter] <- temp.x[, 1]
                rkxc[, counter] <- temp.x[, 2]
                counter <- counter + 1
             }
         }
        rkxi <- apply(rkxc, 1, function(x) all(x > 0))
        rkx <- as.matrix(rkx[rkxi, ])
        rkxc <- as.matrix(rkxc[rkxi, ])
        p.val.l <- pcrksum_str(Wc, rkx, rkxc, mgv, ngv, rep(nrow(rkx), n.layer))
        
        pval<- switch(alternative,
                      less = p.val.l,
                      greater = 1 - p.val.l,
                      two.sided = 2 * min(p.val.l, 1 - p.val.l))
        names(mu) <- "location"

        names(Wc) <- "rank sum statistic"

        result <- list(Rstat = Wc, p.value = pval,
                       null.value = mu, alternative = alternative,
                       data.name = DNAME, method = METHOD,
                       balance = bal, exact = exact)
        class(result) <- "ctest"
        return(result)



    } else {
        mgv <- ngv <- Ngv <- Rsumgv <- VRgv <- numeric(l.csu * l.stu)
    ## mgv: number of clusters under trt X with csize g in stratumum v
    ## ngv: number of clusters under trt Y with csize g in stratumum v
    ## Ngv: number of clusters with csize g in stratumum v
    ## Rsumgv: rank sum of clusters with size g in stratumum v
    ## VRgv: variance of rank of clusters with size g in stratumum v
    Csize <- rep(csize.uniq, l.stu)
    Str <- rep(strt.uniq, each = l.csu)
    summary.mat <- data.frame(mgv, ngv, Ngv, Rsumgv, VRgv, Csize, Str)
    for( i in csize.uniq) {
        for( j in strt.uniq) {
            j.ch <- as.character(j)
            i.ch <- as.character(i)
            temp <- dat.l[[j.ch]][[i.ch]] ## data with csize i and stratumum j
            temp.clus <- temp[, "clus"]
            temp.rksum <- temp[, "rksum"]
            temp.grp <- temp[, "grp"]
            summary.mat[which(Csize == i & Str == j), "mgv"] <- length(temp.grp[temp.grp == 1])
            summary.mat[which(Csize == i & Str == j), "Ngv"] <- length(temp.clus)
            summary.mat[which(Csize == i & Str == j), "Rsumgv"] <- sum(temp.rksum)
            summary.mat[which(Csize == i & Str == j), "VRgv"] <- var(temp.rksum) * (length(temp.clus) - 1)
        }
    }
        summary.mat[, "ngv"] <- summary.mat[, "Ngv"] - summary.mat[, "mgv"]
    ## Wc is the ranksum stat
        Wc <- sum(dat[dat$grp == 1, "rksum"])
        Ngv <- summary.mat[, "Ngv"]
        mgv <- summary.mat[, "mgv"]
        ngv <- summary.mat[, "ngv"]
        VRgv <- summary.mat[, "VRgv"]
        Rsumgv <- summary.mat[, "Rsumgv"]
        EWc <- sum(mgv * Rsumgv / Ngv)
        VarWc <- sum(mgv * ngv / (Ngv * (Ngv - 1)) * VRgv)
        Zc <- (Wc - EWc) / sqrt(VarWc)

        pval <- switch(alternative, less = pnorm(abs(Zc)),
                       greater = pnorm(abs(Zc), lower.tail = FALSE),
                       two.sided = 2 * min(pnorm(abs(Zc)),
                                           pnorm(abs(Zc), lower.tail = FALSE)))

        names(Wc) <- "rank sum statistic"
        names(EWc) <- "expected value of rank sum statistic"
        names(VarWc) <- "variance of rank sum statistic"
        names(Zc) <- "test statistic"
        names(mu) <- "difference in locations"
        result <- list(Rstat = Wc, ERstat = EWc,
                 VRstat = VarWc,
                 statistic = Zc, p.value = pval,
                 alternative = alternative,
                 null.value = mu,
                 data.name = DNAME,
                 method = METHOD,
                 balance = bal, exact = exact)
        class(result) <- "ctest"
        return(result)
    }

}



cluswilcox.test.ranksum.rgl.sub <- function(x, cluster, group, alternative,
                                       exact, mu, DNAME = NULL, METHOD = NULL, stratum) {
    ## The input data should be already arranged
    ## check balance of data.

    bal <- (length(table(table(cluster))) == 1)
    clus <- unique(cluster)
    cluster <- recoderFunc(cluster, clus, c(1 : length(clus)))
    n.clus <- length(clus)
    n.obs <- length(x)
    group <- recoderFunc(group, unique(group), c(1, 2))
    x[which(group == 1)] <- x[which(group == 1)] - mu

    temp <- order(cluster)
    x <- x[temp]
    cluster <- cluster[temp]
    group <- group[temp]

    
    if(bal == TRUE) {
        xrank <- rank(x)
        rksum <- stats::aggregate(xrank ~ cluster, FUN = sum)[, 2]
        csize <- (table(cluster))
        count.grp <- function(x) length(x[ x==1])
        ## q: obs under trt x in each cluster
        q <-stats::aggregate(group ~ cluster, FUN = count.grp)[, 2]
        
        ones <- rep(1, n.clus)
        nq.mat <- stats::aggregate(ones ~ q, FUN = sum)
        vrk <- stats::aggregate(xrank ~ cluster, FUN = var)[, 2]
        q.uniq <- nq.mat[, 1]
        nq <- nq.mat[, 2]
        varQ <- (sum(q.uniq^2 * nq) - (sum(q.uniq * nq)) ^ 2 / n.clus) / n.clus
        sb.2 <- sum((rksum - csize * (csize * n.clus + 1) / 2) ^ 2) / n.clus
        sw.2 <- sum(vrk) / n.clus
        Wc <- sum(as.numeric(group == 1) * xrank)
        EWc <- (n.clus * unique(csize) + 1) / 2 * sum(q.uniq * nq)
        G <- unique(csize)
        varb <- sum(nq * (G - q.uniq) * q.uniq * sw.2) / G
        VarWc <- n.clus * (n.clus / (n.clus - 1) * varQ * sb.2 / unique(csize) ^ 2 + varb / n.clus)
        Zc <- (Wc - EWc) / sqrt(VarWc)
        pval <- switch(alternative,
                       less = pnorm(abs(Zc)),
                       greater = pnorm(abs(Zc), lower.tail = FALSE),
                       two.sided = 2 * min(pnorm(abs(Zc)),
                                           pnorm(abs(Zc), lower.tail = FALSE)))
        names(Wc) <- "rank sum statistic"
        names(EWc) <- "expected value of rank sum statistic"
        names(VarWc) <- "variance of rank sum statistic"
        names(Zc) <- "test statistic"
        names(mu) <- "difference in locations"
        result <- list(Rstat = Wc, ERstat = EWc,
                 VRstat = VarWc,
                 statistic = Zc, p.value = pval,
                 alternative = alternative,
                 null.value = mu,
                 data.name = DNAME,
                 method = METHOD,
                 balance = bal)
        class(result) <- "ctest"
        return(result)

    } else {
        dat <- data.frame(x, cluster, group)
        ones <- rep(1, n.obs)
        csize <- stats::aggregate(ones ~ cluster, FUN = sum)
        n.csize <- length(table(csize[, "ones"]))
        colnames(csize)[2] <- "csize"
        dat <- merge(dat, csize, by = "cluster")
        dat.l <- split(dat, dat$csize)
        getWc <- function(dat) {
            xrank <- rank(dat[, "x"])
            sum(xrank * as.numeric(dat[, "group"] == 1))
        }
        Wc.l <- lapply(dat.l, getWc)
        Wc <- unlist(Wc.l)
        csize.grp <- as.numeric(names(dat.l))
        CountQplus <- function(dat) {
            sum(length(which(dat[, "group"] == 1 )))
        }
        Qplus <- lapply(dat.l, CountQplus)
        Qplus <- unlist(Qplus)

        CountN <- function(dat) {
            length(table(dat[, "cluster"]))
        }

        N.clus <- lapply(dat.l, CountN)
        N.clus <- unlist(N.clus)
        N.obs <- csize.grp * N.clus

        getRank <- function(dat) {
            rank(dat[, "x"])
        }
        rank.l <- lapply(dat.l, getRank)

        varR <- unlist(lapply(rank.l, var)) * (N.obs - 1) / N.obs
        varR <- unlist(varR)
        getGrp <- function(dat) {
            dat[, "group"]
        }
        grp.l <- lapply(dat.l, getGrp)


        ## theta: estimated prob of an obs under trt1 has a larger rank than an obs from trt 2
        ## deno: the denominator used in the estimation of intercluster correlation
        theta <- deno <- numeric(n.csize)
        trans.rk <- vector("list", n.csize)
        for( i in 1 : n.csize) {
            theta[i] <-( Wc[i] - Qplus[i] * (Qplus[i] + 1) / 2 ) / (Qplus[i] * (csize.grp[i] * N.clus[i] - Qplus[i]))
            deno[i] <- csize.grp[i] * N.clus[i] + 1
            trans.rk[[i]] <- qnorm(rank.l[[i]] / deno[i])
        }
        pscore <- unlist(trans.rk)
        aovout <- stats::aov(pscore ~ as.factor(cluster))
        aovout <- summary(aovout)

        errorms <- aovout[[1]][["Mean Sq"]][2]
        modelms <- aovout[[1]][["Mean Sq"]][1]
        sumn <- n.obs
        sumnsq <- sum(unlist(N.clus) * csize.grp ^ 2)

        g0 <- (sumn - sumnsq / sumn) / (n.clus - 1)
        sigsqa <- (modelms - errorms) / g0
        rhoh <- sigsqa / (sigsqa + errorms)
        if(rhoh < 0) {
            rhoh <- 0
        }
        rhoh <- rhoh * ( 1 + ( 1 - rhoh ^ 2) / (2 * (n.clus - 4)))
        rhor <- (6 / pi) * asin(rhoh / 2)
        varR <- unlist(varR)
        ssb <- csize.grp * varR * ( 1 + (csize.grp - 1) * rhor)
        ssw <- varR * ( 1 - rhor)

        getEQgQ <- function(dat) {
            tab.clus <- table(dat[, "cluster"])
            csize <- tab.clus[1]
            n.clus <- length(tab.clus)
            summary.q <- table(dat[, "cluster"] * as.numeric(dat[, "group"] == 1))
            sum.n <- names(summary.q)
            summary.q <- summary.q[sum.n[-which(sum.n == '0')]]
            sum.q <- table(summary.q)
            q.in <- as.numeric(names(sum.q))
            E <- sum(q.in * sum.q * (csize - q.in)) / n.clus
            return(E)
        }

        EQgQ <- lapply(dat.l, getEQgQ)
        EQgQ <- unlist(EQgQ)

        getVarQ <- function(dat) {
            tab.clus <- table(dat[, "cluster"])
            n.clus <- length(tab.clus)
            summary.q <- table(dat[, "cluster"] * as.numeric(dat[, "group"] == 1))
            sum.n <- names(summary.q)
            summary.q <- summary.q[sum.n[-which(sum.n == '0')]]
            sum.q <- table(summary.q)
            q.in <- as.numeric(names(sum.q))
            V <- (sum(q.in ^ 2 * sum.q) - sum(q.in * sum.q) ^ 2 / n.clus) / n.clus
            return(V)
        }


        VarQ <- lapply(dat.l, getVarQ)
        VarQ <- unlist(VarQ)
        
        VarTheta <- N.clus * (N.clus / (N.clus - 1) * as.numeric(N.clus > 1) *
                              VarQ * ssb / csize.grp ^ 2 + EQgQ * ssw / csize.grp) /
            (Qplus * (csize.grp * N.clus - Qplus)) ^ 2
        VarTheta[which(N.clus < 2)] <- 0
        wt <- 1 / VarTheta[which(VarTheta > 0)]

        theta <- theta[which(VarTheta > 0)]
        sumtheta <- sum(theta * wt)
        sumwt <- sum(wt)

        theta <- sumtheta / sumwt
        sdtheta <- sqrt(1 / sumwt)
        Zc <- (theta - 1/2) / sdtheta
        pval <- switch(alternative, less = pnorm(abs(Zc)),
                       greater = pnorm(abs(Zc), lower.tail = FALSE),
                       two.sided = 2 * min(pnorm(abs(Zc)),
                                           pnorm(abs(Zc), lower.tail = FALSE)))

        names(Zc) <- "Test statistic"
        
        names(mu) <- "difference in locations"
        result <- list( statistic = Zc, p.value = pval,
                       alternative = alternative, null.value = mu,
                       data.name = DNAME, method = METHOD,
                       balance = bal)
        class(result) <- "ctest"
        return(result)
    }
}

#' @importFrom MASS ginv

cluswilcox.test.ranksum.ds <- function(x, cluster, group,
                                       alternative,
                                       mu,
                                       DNAME, METHOD) {
    group.uniq <- length(unique(group))
    if(group.uniq == 1) {
        stop("invalid group variable, should contain at least 2 groups")
    }
    
    n.obs <- length(x)
    
    Fhat <- numeric(n.obs)
    order.c <- order(cluster)
    cluster <- cluster[order.c]
    x <- x[order.c]
    group <- group[order.c]
    cluster.uniq <- unique(cluster)
    M <- length(cluster.uniq)

    ni <- table(cluster)
    
    temp <- unique(sort(abs(x)))
    diff.min <- min(diff(temp)) / 10 ## A quntity 1 / 10 of the minimum difference between scores
    
    FHat <- ecdf(x)
    Fhat <- FHat(x) / 2 + FHat(x - diff.min) / 2
     
    
    F.prop <- Fprop(x, cluster, ni, M, n.obs)
    F.prop2 <- F.prop
    
    if(group.uniq == 2) {
#####calculate quantity 2 (using the pooled estimate of F)
        group <- recoderFunc(group, order(unique(group)), c(0, 1))
        x[which(group == 0)] <- x[which(group == 0)] - mu
        ni1 <- aggregate(group ~ cluster, FUN = sum)[, 2] # number of obs under trt 2 in each cluster
        if(all(ni1 / ni == 0.5)) {
            warning("The DS ranksum test is not reliable for colateral data where each cluster is equally split between the 2 treatments.")
        }
        ## Calculate S = E(W*|W, g)
        ##  S <- 0
        Ni <- rep.int(ni, times = ni)
        S <- sum(group / Ni * (1 + F.prop))
        ##       for( i in 1 : M) {
        ##         S <- S + sum(group[cluster == cluster.uniq[i]] /
        ##                      ni[i] * ( 1 + F.prop[cluster == cluster.uniq[i]]))
        ##    }
        ##   
        S <- S / (M + 1)
        
        ## Calculate E(S)
        ES <- sum(ni1 / ni) / 2
        
        ## Calculate var(S)
        
        W <- numeric(M)
        for(i in 1 : M) {
            Wi <- ((M - 1) * group[cluster == cluster.uniq[i]] -
                   sum(ni1[-i] / ni[-i])) *
                Fhat[cluster == cluster.uniq[i]]
            W[i] <- sum(Wi) / (ni[i] * (M+1))
        }
        a <- sum(ni1 / ni)
        EW <- M / (2 * (M + 1)) * (ni1 / ni - a / M)
        varS <- sum((W - EW)^2)
        Z <- (S - ES) / sqrt(varS)
        
        pval <- switch(alternative, less = pnorm(abs(Z)),
                       greater = pnorm(abs(Z), lower.tail = FALSE),
                       two.sided = 2 * min(pnorm(abs(Z)),
                                           pnorm(abs(Z), lower.tail = FALSE)))
        names(Z) <- "test statistic"
        names(mu) <- "difference in locations"

        result <- list(statistic = Z, p.value = pval,
                       alternative = alternative, null.value = mu,
                       data.name = DNAME, method = METHOD)

        class(result) <- "ctest"
        result
        
        
    } else {
#####calculate quantity 2 (using the pooled estimate of F)
        if(!is.null(mu) & mu != 0) {
            warning("comparison between m (m > 2) groups cannot set location shift parameter")
        }
        mu <- 0
        m <- length(unique(group))
        Sj <- numeric(m)
        group.uniq <- unique(group)
        for( j in 1 : m) {
            for( i in 1 : M) {
                gik.j <- ifelse(group == j, 1, 0)
                Sj[j] <- Sj[j] + sum(gik.j[cluster == cluster.uniq[i]] *
                                     (1 + F.prop[cluster == cluster.uniq[i]])) / ni[i]
            }
        }
        Sj <- Sj / (M + 1)
        
        nij <- matrix( 0, m, M)
        for( i in 1 : m) {
            nij[i, ] <- table(cluster[group == group.uniq[i]])
        }
        d <- apply(nij, 1, FUN = function(x) {x / ni})
        ESj <- apply(d, 2, sum) / 2
        
        What <- matrix(0, m, M)
        a <- t(d)
        for( i in 1 : M) {
            for( j in 1 : m) {
                gik.j <- ifelse(group[cluster == cluster.uniq[i]] == j, 1, 0)
                b <- 1 / (ni[i] * (M + 1))
                c <- gik.j * (M - 1)
                d <- sum(a[j, -i])
                What[j, i] <- b * sum((c - d) * Fhat[cluster == cluster.uniq[i]])
            }
        }
        EW <- matrix(0, m, M)
        for( j in 1 : m) {
            EW[j, ] <- (M / (2 * ( M + 1))) * (a[j, ] - sum(a[j, ]) / M)
        }

        dev.W <- What - EW
        V <- matrix(0, m, m)
        for( i in 1 : M) {
            V <- V + dev.W[, i] %*% t(dev.W[, i])
        }
        V <- V / M
        T <- t(Sj - ESj) %*% ginv(V) %*% (Sj - ESj) / M
        pval <- pchisq(T, df = (m - 1),lower.tail = F)
        names(T) <- "chi-square test statistic"
        ngrp <- length(group.uniq)
        df <- ngrp - 1
        names(ngrp) <- "number of groups: "
        names(df) <- "degree of freedom: "
        METHOD <- paste(METHOD, "using Chi-square test")
                                        #calculate the test statistic
        result <- list(statistic = T, p.value = pval, n.group = ngrp,
                       df = df,
                       data.name = DNAME, method = METHOD)
        class(result) <- "ctest"
        result
    }
}
