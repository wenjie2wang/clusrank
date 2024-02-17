##
## clusrank: Wilcoxon Rank Tests for Clustered Data
##
## Copyright (C) 2015-2024 Yujing Jiang, Mei-Ling Ting Lee, and Jun Yan
## Copyright (C) 2022-2024 Wenjie Wang
##
## This file is part of the R package clusrank.
##
## The R package clusrank is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package clusrank is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##


clusWilcox_test_ranksum_rgl <- function(x, cluster, group, stratum,
                                        alternative, exact, B, mu,
                                        DNAME = NULL, METHOD = NULL) {
    clus.grp <- lapply(split(cluster, group), unique)
    l.clus.grp <- length(unlist(clus.grp))
    l.clus <- length(unique(cluster))
    temp <- merge(table(cluster), cbind(cluster, group))
    temp <- temp[order(temp$Freq), ]
    temp <- split(temp$group, temp$Freq)
    check <- lapply(temp, function(x) length(unique(x)))
### check with each stratum in which clusters have the same cluster size
### to see if they only are only assigned a single treatment.
    if (all(unlist(check) == 1))
        warning("The two groups should contain clusters with the same size for at leaset one cluster size")
    n.obs <- length(x)
    one <- rep(1, n.obs)
    csize <- stats::aggregate(one ~ cluster, FUN = sum)[, 2]
    if (length(table(csize)) == length(unique(cluster))) {
        stop("For RGL method for rank-sum test, there should be at least a cluster size level that contains at least two clusters")
    }

    arglist <- setNames(list(x, cluster, group, stratum, alternative,
                             mu, DNAME, METHOD, exact, B),
                        c("x", "cluster", "group", "stratum",
                          "alternative", "mu", "DNAME", "METHOD",
                          "exact", "B"))
    if (l.clus != l.clus.grp) {
        result <- do.call("clusWilcox_test_ranksum_rgl_sub", c(arglist))
        return(result)
    } else {
        if (exact == FALSE | (exact == TRUE & B == 0)) {
            result <- do.call("clusWilcox_test_ranksum_rgl_clus", c(arglist))
            return(result)
        } else {
            result <- do.call("clusWilcox_test_ranksum_rgl_clus_exact", c(arglist))
            return(result)
        }

    }
}




clusWilcox_test_ranksum_rgl_clus_exact1 <- function(x, cluster, group,
                                                    stratum) {
    ## n.obs <- length(x)
    ## one <- rep(1, n.obs)
    ## x[which(group == 1)] <- x[which(group == 1)]
    xrank <- rank(x)
    rksum <-stats::aggregate(xrank ~ cluster, FUN = sum)[, 2]
    ## csize <-stats::aggregate(one ~ cluster, FUN = sum)[, 2]
    grp <- stats::aggregate(group ~ cluster, FUN = mean)[, 2]
    ## strt <- stats::aggregate(stratum ~ cluster, FUN = mean)[, 2]
    ## clus <- unique(cluster)
    ## n.clus <- length(clus)
    ## dat <- data.frame(clus, strt, grp, csize, rksum)
    ## bal <- (!(length(unique(table(cluster))) != 1L))

    ## csize.uniq <- unique(csize)
    ## strt.uniq <- unique(strt)
    ## l.csu <- length(csize.uniq)
    ## l.stu <- length(strt.uniq)

    ## dat.l <- split(dat, strt) ## Split the rank data by stratum
    ## csize.split <- function(dat) {
    ##     split(dat, dat$"csize")
    ## }
    ## dat.l <- lapply(dat.l, csize.split)
    ## sum(dat[dat$grp == 1, "rksum"])
    sum(rksum[grp == 1])
}


clusWilcox_test_ranksum_rgl_clus_exact <- function(x, cluster, group,
                                                   stratum, alternative,
                                                   exact, B, mu,
                                                   DNAME = NULL,
                                                   METHOD = NULL)
{
    METHOD <- paste0(METHOD, " (random permutation)")
    x <- x - mu
    n.obs <- length(x)
    one <- rep(1, n.obs)
    csize <-stats::aggregate(one ~ cluster, FUN = sum)
    x.csize <- merge(cbind(x, cluster), csize)
    csize <- csize[, 2]
    grp <- stats::aggregate(group ~ cluster, FUN = mean)[, 2]
    strt <- stats::aggregate(stratum ~ cluster, FUN = mean)[, 2]
    clus <- unique(cluster)
    n.clus <- length(clus)
    bal <- (!(length(unique(table(cluster))) != 1L))

    strt.uniq <- unique(strt)
    csize.uniq <- unique(csize)

    lsc <- length(strt.uniq) * length(csize.uniq)

    ind.l <- vector("list", lsc)
    ## Record the random group indicator, classified by cluster size and stratum

    ct <- 1
    for ( i in 1 : length(strt.uniq)) {
        for ( j in 1 : length(csize.uniq)) {
            temp <- (stratum == strt.uniq[i] & x.csize[, 3] == csize.uniq[j])
            if (all(temp == FALSE)) next
            ind.l[[ct]] <- temp
            ind.l[[ct]] <- cbind(x[temp], group[temp], cluster[temp],
                                 strt.uniq[i], csize.uniq[j])
            ct <- ct + 1

        }
    }

    ind.l <- ind.l[!unlist(lapply(ind.l, is.null))]

    samp.ind <- function(x) {
        grp <- x[, 2]
        cls <- x[, 3]
        temp.grp <- stats::aggregate(grp ~ cls, FUN = mean)[, 2]
        temp.grp <- sample(temp.grp, length(temp.grp))
        temp.grp <- rep(temp.grp, each = x[1, 5])
        x[, 2] <- temp.grp
        x
    }

    W <- clusWilcox_test_ranksum_rgl_clus_exact1(x, cluster, group, stratum)

    W.vec <- rep(NA, B)
    for ( i in 1 : B) {
        temp <- lapply(ind.l, samp.ind)
        temp1 <- NULL
        for ( j in 1 : length(temp)) {
            temp1 <- rbind(temp1, temp[[j]])
        }

        cluster.temp <- temp1[, 3]

        temp1 <- temp1[order(cluster.temp), ]

        x.temp <- temp1[, 1]
        group.temp <- temp1[, 2]
        cluster.temp <- temp1[, 3]

        str.temp <- temp1[, 4]




        W.vec[i] <- clusWilcox_test_ranksum_rgl_clus_exact1(x.temp,
                                                            cluster.temp,
                                                            group.temp,
                                                            str.temp)
    }
    pval <- perm_pvalue(W, W.vec, alternative)
    names(mu) <- "difference in locations"
    names(W) <- "W"
    result <- list(statistic = W, p.value = pval,
                   null.value = mu, alternative = alternative,
                   data.name = DNAME, method = METHOD,
                   balance = bal, exact = exact, B = B, nclus = n.clus,
                   nobs = n.obs)
    class(result) <- "ctest"
    return(result)
}



clusWilcox_test_ranksum_rgl_clus <- function(x, cluster, group,
                                             stratum, alternative,
                                             exact, B,
                                             mu, DNAME = NULL, METHOD = NULL) {
### The input data should be already in the increasing order with cluster.
    n.obs <- length(x)
    one <- rep(1, n.obs)
    x[which(group == 1)] <- x[which(group == 1)] - mu
    xrank <- rank(x)
    rksum <-stats::aggregate(xrank ~ cluster, FUN = sum)[, 2]
    csize <-stats::aggregate(one ~ cluster, FUN = sum)[, 2]
    grp <- stats::aggregate(group ~ cluster, FUN = mean)[, 2]
    strt <- stats::aggregate(stratum ~ cluster, FUN = mean)[, 2]
    clus <- unique(cluster)
    n.clus <- length(clus)
    dat <- data.frame(clus, strt, grp, csize, rksum)
    bal <- (!(length(unique(table(cluster))) != 1L))

    csize.uniq <- unique(csize)
    strt.uniq <- unique(strt)
    l.csu <- length(csize.uniq)
    l.stu <- length(strt.uniq)

    dat.l <- split(dat, strt) ## Split the rank data by stratum
    csize.split <- function(dat) {
        split(dat, dat$"csize")
    }
    dat.l <- lapply(dat.l, csize.split)

    if (exact == TRUE & B == 0) {
        METHOD <- paste0(METHOD, " (exact permutation)")
        if(length(table(cluster)) > 20)
            warning("Number of clusters exceeds 20 for RGL clustered rank exact test")
        W <- sum(dat[dat$grp == 1, "rksum"])
        n.layer <- l.csu * l.stu
        mgv <- ngv <- rep(0, n.layer)
        ct <- 1
        rkx <- rkxc <- as.matrix(matrix(0, W + 1, n.layer))
        counter <- 1
        for ( i in csize.uniq) {
            for ( j in strt.uniq) {
                j.ch <- as.character(j)
                i.ch <- as.character(i)
                temp <- dat.l[[j.ch]][[i.ch]] # data with csize i and stratum j
                temp.rksum <- temp[, "rksum"]
                temp.grp <- temp[, "grp"]
                mgv[counter] <- length(temp.grp[temp.grp == 1])
                ngv[counter] <- length(temp.grp)
                temp.x <- cumcrksum(W, mgv[counter], sort(temp.rksum), csize.uniq)
                rkx[, counter] <- temp.x[, 1]
                rkxc[, counter] <- temp.x[, 2]
                counter <- counter + 1
             }
         }
        rkxi <- apply(rkxc, 1, function(x) all(x > 0))
        rkx <- as.matrix(rkx[rkxi, ])
        rkxc <- as.matrix(rkxc[rkxi, ])
        p.val.l <- pcrksum_str(W, rkx, rkxc, mgv, ngv, rep(nrow(rkx), n.layer))

        pval<- switch(alternative,
                      less = p.val.l,
                      greater = 1 - p.val.l,
                      two.sided = 2 * min(p.val.l, 1 - p.val.l))
        names(mu) <- "location"

        names(W) <- "W"

        result <- list(statistic = W, p.value = pval,
                       null.value = mu, alternative = alternative,
                       data.name = DNAME, method = METHOD,
                       balance = bal, exact = exact, B = B, nclus = n.clus,
                       nobs = n.obs)
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
        for ( i in csize.uniq) {
            for ( j in strt.uniq) {
                j.ch <- as.character(j)
                i.ch <- as.character(i)
                temp <- dat.l[[j.ch]][[i.ch]] # data with csize i and stratumum j
                temp.clus <- temp[, "clus"]
                temp.rksum <- temp[, "rksum"]
                temp.grp <- temp[, "grp"]
                summary.mat[which(Csize == i & Str == j), "mgv"] <- length(temp.grp[temp.grp == 1])
                summary.mat[which(Csize == i & Str == j), "Ngv"] <- length(temp.clus)
                summary.mat[which(Csize == i & Str == j), "Rsumgv"] <- sum(temp.rksum)
                VRgv.temp <- 0
                if (length(temp.clus) > 1) {
                    VRgv.temp <- var(temp.rksum) * (length(temp.clus) - 1)
                }
                summary.mat[which(Csize == i & Str == j), "VRgv"] <- VRgv.temp
            }
        }
        summary.mat[, "ngv"] <- summary.mat[, "Ngv"] - summary.mat[, "mgv"]
        W <- sum(dat[dat$grp == 1, "rksum"]) # W is the ranksum stat
        Ngv <- summary.mat[, "Ngv"]
        mgv <- summary.mat[, "mgv"]
        ngv <- summary.mat[, "ngv"]
        VRgv <- summary.mat[, "VRgv"]
        Rsumgv <- summary.mat[, "Rsumgv"]
        EW <- sum((mgv * Rsumgv / Ngv)[Ngv > 0])
        VarW <- sum((mgv * ngv / (Ngv * (Ngv - 1)) * VRgv)[VRgv > 0])
        Z <- (W - EW) / sqrt(VarW)

        pval <- switch(alternative, less = pnorm(Z),
                       greater = pnorm(Z, lower.tail = FALSE),
                       two.sided = 2 * min(pnorm(abs(Z)),
                                           pnorm(abs(Z), lower.tail = FALSE)))

        names(W) <- "W"
        names(EW) <- "Expected value of W"
        names(VarW) <- "Variance of W"
        names(Z) <- "Z"
        names(mu) <- "difference in locations"
        result <- list(Rstat = W, ERstat = EW,
                       VRstat = VarW,
                       statistic = Z, p.value = pval,
                       alternative = alternative,
                       null.value = mu,
                       data.name = DNAME,
                       method = METHOD,
                       balance = bal, exact = exact,
                       B = B,
                       nobs = n.obs, nclus = n.clus)
        class(result) <- "ctest"
        return(result)
    }

}





clusWilcox_test_ranksum_rgl_sub_exact1 <- function(x, cluster, group,
                                                 stratum) {

    bal <- (length(table(table(cluster))) == 1) # check balance of data.
    clus <- unique(cluster)
    if (is.numeric(clus)) clus <- sort(clus)
    cluster <- recoderFunc(cluster, clus, c(1 : length(clus)))
    n.clus <- length(clus)
    n.obs <- length(x)
    group <- recoderFunc(group, unique(group), c(1, 2))
    temp <- order(cluster)
    x <- x[temp]
    cluster <- cluster[temp]
    group <- group[temp]

    if (bal == TRUE) {
        xrank <- rank(x)
        W <- sum(as.numeric(group == 1) * xrank)
        return(W)

    } else {
        dat <- data.frame(x, cluster, group)
        ones <- rep(1, n.obs)
        csize <- stats::aggregate(ones ~ cluster, FUN = sum)
        n.csize <- length(table(csize[, "ones"]))
        colnames(csize)[2] <- "csize"
        dat <- merge(dat, csize, by = "cluster")
        dat.l <- split(dat, dat$csize)
        getW <- function(dat) {
            xrank <- rank(dat[, "x"])
            sum(xrank * as.numeric(dat[, "group"] == 1))
        }
        W.l <- lapply(dat.l, getW)
        W <- unlist(W.l)
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
        for ( i in 1 : n.csize) {
            theta[i] <-( W[i] - Qplus[i] * (Qplus[i] + 1) / 2 ) / (Qplus[i] * (csize.grp[i] * N.clus[i] - Qplus[i]))
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
        Z <- (theta - 1/2) / sdtheta
        names(Z) <- "Z"
        return(Z)
    }
}


clusWilcox_test_ranksum_rgl_sub_exact <- function(x, cluster, group,
                                                  alternative, exact, B, mu,
                                                  DNAME = NULL,
                                                  METHOD = NULL,
                                                  stratum) {
    METHOD <- paste0(METHOD, " (Random Permutation)")
    bal <- (length(table(table(cluster))) == 1) # check balance of
                                        # data.
    x[which(group == 1)] <- x[which(group == 1)] - mu

    n.obs <- length(x)
    n.clus <- length(unique(cluster))
    one <- rep(1, n.obs)
    csize <- stats::aggregate(one ~ cluster, FUN = sum)
    colnames(csize)[2] <- "csize"
    group.1 <- ifelse(group == 1, 1, 0)

    q <- stats::aggregate(group.1 ~ cluster, FUN = sum)
    colnames(q)[2] <- "g1size"

    dat <- as.data.frame(merge(merge(cbind(x, cluster, group), csize), q))

    q.uniq <- unique(q[, 2])
    csize.uniq <- unique(csize[, 2])

    lqc <- length(q.uniq) * length(csize.uniq)
    ind.l <- vector("list", lqc)
    ## Record the random group indicator, classified by cluster size
    ## and number of cluster members from group 1

    ct <- 1
    for ( i in (q.uniq)) {
        for ( j in (csize.uniq)) {
            temp <- (dat["g1size"] == i & dat["csize"] == j)
            if (all(temp == FALSE)) next
            ind.l[[ct]] <- dat[temp, ]
            ct <- ct + 1
        }
    }
    ind.l <- ind.l[!unlist(lapply(ind.l, is.null))]

    samp.ind <- function(x) {
        q <- x["g1size"][[1]][1]
        csize <- x["csize"][[1]][1]
        q2 <- csize - q
        one <- rep.int(c(1, 2), times = c(q, q2))

        n <- length(unique(x["cluster"]))

        grp <- c(replicate(n, sample(one, csize)))
        x["group"] <- grp
        x
    }

    W <- clusWilcox_test_ranksum_rgl_sub_exact1(x, cluster,
                                                group, stratum)

    W.vec <- rep(NA, B)
    for ( i in 1 : B) {
        temp <- lapply(ind.l, samp.ind)
        temp1 <- NULL
        for ( j in 1 : length(temp)) {
            temp1 <- rbind(temp1, temp[[j]])
        }

        cluster.temp <- temp1[, 1]
        temp1 <- temp1[order(cluster.temp), ]
        x.temp <- temp1[, 2]
        group.temp <- temp1[, 3]
        cluster.temp <- temp1[, 1]

        W.vec[i] <- clusWilcox_test_ranksum_rgl_sub_exact1(x.temp,
                                                           cluster.temp,
                                                           group.temp, stratum)
    }
    pval <- perm_pvalue(W, W.vec, alternative)
    names(mu) <- "shift in location"
    if (bal)
        names(W) <- "W"
    else
        names(W) = "Z"
    result <- list(statistic = W, p.value = pval,
                   null.value = mu, alternative = alternative,
                   data.name = DNAME, method = METHOD,
                   balance = bal, exact = exact, B = B,  nclus = n.clus,
                   nobs = n.obs)
    class(result) <- "ctest"
    return(result)
}



clusWilcox_test_ranksum_rgl_sub <- function(x, cluster, group, alternative,
                                       exact, B, mu, DNAME = NULL, METHOD = NULL, stratum) {
### The input data should be already arranged
    if (exact == TRUE & B >= 1)
        return(clusWilcox_test_ranksum_rgl_sub_exact(x, cluster, group,
                                                    alternative, exact, B,
                                                    mu, DNAME, METHOD, stratum))
    if (exact == TRUE & B == 0)
        stop("Exact exactutation test is not available for RGL clustered rank-sum test when treatment is assigned at subunit level.")
    bal <- (length(table(table(cluster))) == 1) # check balance of data.
    clus <- unique(cluster)
    if (is.numeric(clus)) clus <- sort(clus)
    cluster <- recoderFunc(cluster, clus, c(1 : length(clus)))
    n.clus <- length(clus)
    n.obs <- length(x)
    group <- recoderFunc(group, unique(group), c(1, 2))
    x[which(group == 1)] <- x[which(group == 1)] - mu

    temp <- order(cluster)
    x <- x[temp]
    cluster <- cluster[temp]
    group <- group[temp]

    if (bal == TRUE) {
        xrank <- rank(x)
        rksum <- stats::aggregate(xrank ~ cluster, FUN = sum)[, 2]
        csize <- (table(cluster))
        count.grp <- function(x) length(x[ x==1])
        q <-stats::aggregate(group ~ cluster, FUN = count.grp)[, 2] # obs under trt x in each cluster

        ones <- rep(1, n.clus)
        nq.mat <- stats::aggregate(ones ~ q, FUN = sum)
        vrk <- stats::aggregate(xrank ~ cluster, FUN = var)[, 2]
        q.uniq <- nq.mat[, 1]
        nq <- nq.mat[, 2]
        varQ <- (sum(q.uniq^2 * nq) - (sum(q.uniq * nq)) ^ 2 / n.clus) / n.clus
        sb.2 <- sum((rksum - csize * (csize * n.clus + 1) / 2) ^ 2) / n.clus
        sw.2 <- sum(vrk) / n.clus
        W <- sum(as.numeric(group == 1) * xrank)
        EW <- (n.clus * unique(csize) + 1) / 2 * sum(q.uniq * nq)
        G <- unique(csize)
        varb <- sum(nq * (G - q.uniq) * q.uniq * sw.2) / G
        VarW <- n.clus * (n.clus / (n.clus - 1) * varQ * sb.2 / unique(csize) ^ 2 + varb / n.clus)
        Z <- (W - EW) / sqrt(VarW)
        pval <- switch(alternative,
                       less = pnorm(Z),
                       greater = pnorm(Z, lower.tail = FALSE),
                       two.sided = 2 * min(pnorm(abs(Z)),
                                           pnorm(abs(Z), lower.tail = FALSE)))
        names(W) <- "W"
        names(EW) <- "expected value of W"
        names(VarW) <- "variance of W"
        names(Z) <- "Z"
        names(mu) <- "difference in locations"
        result <- list(Rstat = W, ERstat = EW,
                 VRstat = VarW,
                 statistic = Z, p.value = pval,
                 alternative = alternative,
                 null.value = mu,
                 data.name = DNAME,
                 method = METHOD,
                 balance = bal, nobs = n.obs, nclus = n.clus)
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
        getW <- function(dat) {
            xrank <- rank(dat[, "x"])
            sum(xrank * as.numeric(dat[, "group"] == 1))
        }
        W.l <- lapply(dat.l, getW)
        W <- unlist(W.l)
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
        for ( i in 1 : n.csize) {
            theta[i] <-( W[i] - Qplus[i] * (Qplus[i] + 1) / 2 ) / (Qplus[i] * (csize.grp[i] * N.clus[i] - Qplus[i]))
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
        Z <- (theta - 1/2) / sdtheta
        pval <- switch(alternative, less = pnorm(Z),
                       greater = pnorm(Z, lower.tail = FALSE),
                       two.sided = 2 * min(pnorm(abs(Z)),
                                           pnorm(abs(Z), lower.tail = FALSE)))

        names(Z) <- "Z"

        names(mu) <- "difference in locations"
        result <- list(statistic = Z, p.value = pval,
                       alternative = alternative, null.value = mu,
                       data.name = DNAME, method = METHOD,
                       balance = bal, nobs = n.obs, nclus = n.clus)
        class(result) <- "ctest"
        return(result)
    }
}

##' @importFrom MASS ginv

clusWilcox_test_ranksum_ds_exact1 <- function(x, cluster, group, mu) {
    group.uniq <- length(unique(group))
    if (group.uniq == 1) {
        stop("invalid group variable, should contain at least 2 groups")
    }

    n.obs <- length(x)

    Fhat <- numeric(n.obs)
    order.c <- order(cluster)
    cluster <- cluster[order.c]
    x <- x[order.c]
    group <- group[order.c]
    cluster.uniq <- unique(cluster)
    n.clus <- M <- length(cluster.uniq)
    cluster <- recoderFunc(cluster, cluster.uniq, c(1 : M))



    ni <- table(cluster)

    temp <- unique(sort(abs(x)))
    diff.min <- min(diff(temp)) / 10 # A quntity 1 / 10 of the minimum difference between scores

    FHat <- ecdf(x)
    Fhat <- FHat(x) / 2 + FHat(x - diff.min) / 2


    F.prop <- Fprop(x, cluster, ni, M, n.obs)
    F.prop2 <- F.prop

    if (group.uniq == 2) {
###calculate quantity 2 (using the pooled estimate of F)
        group <- recoderFunc(group, sort(unique(group)), c(1, 0))
        x[which(group == 0)] <- x[which(group == 0)] - mu
        ni1 <- aggregate(group ~ cluster, FUN = sum)[, 2] # number of
                                        # obs under trt 2 in each
                                        # cluster
        if (all(ni1 / ni == 0.5)) {
            warning("The DS ranksum test is not reliable for colateral data where each cluster is equally split between the 2 treatments.")
        }
        ## Calculate S = E(W*|W, g)
        Ni <- rep.int(ni, times = ni)
        S <- sum(group / Ni * (1 + F.prop))
        S <- S / (M + 1)
        return(S)

    } else {
###calculate quantity 2 (using the pooled estimate of F)
        if (!is.null(mu) && mu != 0) {
            warning("comparison between m (m > 2) groups cannot set location shift parameter")
        }
        mu <- 0
        m <- length(unique(group))
        Sj <- numeric(m)
        group.uniq <- unique(group)
        for ( j in 1 : m) {
            for ( i in 1 : M) {
                gik.j <- ifelse(group == j, 1, 0)
                Sj[j] <- Sj[j] + sum(gik.j[cluster == i] *
                                     (1 + F.prop[cluster == i])) / ni[i]
            }
        }
        Sj <- Sj / (M + 1)

        nij <- matrix( 0, m, M)
        for ( i in 1 : m) {
            nij[i, ] <- table(factor(cluster[group == group.uniq[i]], levels = unique(cluster)))
        }
        d <- apply(nij, 1, FUN = function(x) {x / ni})
        ESj <- apply(d, 2, sum) / 2

        What <- matrix(0, m, M)
        a <- t(d)
        for ( i in 1 : M) {
            for ( j in 1 : m) {
                gik.j <- ifelse(group[cluster == i] == j, 1, 0)
                b <- 1 / (ni[i] * (M + 1))
                c <- gik.j * (M - 1)
                d <- sum(a[j, -i])
                What[j, i] <- b * sum((c - d) * Fhat[cluster == i])
            }
        }
        EW <- matrix(0, m, M)
        for ( j in 1 : m) {
            EW[j, ] <- (M / (2 * ( M + 1))) * (a[j, ] - sum(a[j, ]) / M)
        }

        dev.W <- What - EW
        V <- matrix(0, m, m)
        for ( i in 1 : M) {
            V <- V + dev.W[, i] %*% t(dev.W[, i])
        }
        V <- V / M
        T <- t(Sj - ESj) %*% ginv(V) %*% (Sj - ESj) / M
        return(T)
    }
}

clusWilcox_test_ranksum_ds_exact <- function(x, cluster, group,
                                            alternative,
                                            mu, exact, B,
                                            DNAME, METHOD) {
    x <- x - mu
    n.obs <- length(x)
    n.clus <- length(unique(cluster))
    W.vec <- rep(NA, B)
    W <- clusWilcox_test_ranksum_ds_exact1(x, cluster, group, 0)
    for ( i in 1 : B) {
        grp.temp <- sample(group, n.obs)
        W.vec[i] <- clusWilcox_test_ranksum_ds_exact1(x, cluster,
                                                      grp.temp, 0)
    }
    pval <- perm_pvalue(W, W.vec, alternative)
    METHOD <- paste0(METHOD, " (Random Permutation)")
    names(mu) <- "difference in locations"

    if (length(unique(group)) > 2) {
        names(W) <- "chi-square statistic"
        result <- list(statistic = W, p.value = pval,
                       alternative = alternative, null.value = mu,
                       data.name = DNAME, method = METHOD,
                       exact = exact, B = B,
                       nobs = n.obs, nclus = n.clus)
        class(result) <- "ctest"
    } else {
        names(W) <- "W"
        result <- list(statistic = W, p.value = pval,
                       alternative = alternative, null.value = mu,
                       data.name = DNAME, method = METHOD,
                       exact = exact, B = B,
                       nobs = n.obs, nclus = n.clus)
        class(result) <- "ctest"
    }
    result
}


##' @importFrom MASS ginv

clusWilcox_test_ranksum_ds <- function(x, cluster, group,
                                       alternative,
                                       mu, exact, B,
                                       DNAME, METHOD) {
    if (exact == TRUE)  {
        return(clusWilcox_test_ranksum_ds_exact(x, cluster, group,
                                               alternative,
                                               mu, exact, B, DNAME, METHOD))
    }
    uni_group <- unique(group)
    group.uniq <- length(uni_group)
    if (group.uniq == 1) {
        stop("invalid group variable, should contain at least 2 groups")
    }

    n.obs <- length(x)

    Fhat <- numeric(n.obs)
    order.c <- order(cluster)
    cluster <- cluster[order.c]
    x <- x[order.c]
    group <- group[order.c]
    cluster.uniq <- unique(cluster)
    n.clus <- M <- length(cluster.uniq)
    cluster <- recoderFunc(cluster, cluster.uniq, c(1 : M))



    ni <- table(cluster)

    temp <- unique(sort(abs(x)))
    diff.min <- min(diff(temp)) / 10 # A quntity 1 / 10 of the minimum difference between scores

    FHat <- ecdf(x)
    Fhat <- FHat(x) / 2 + FHat(x - diff.min) / 2


    F.prop <- Fprop(x, cluster, ni, M, n.obs)
    F.prop2 <- F.prop

    if (group.uniq == 2) {
###calculate quantity 2 (using the pooled estimate of F)
        group <- recoderFunc(group, sort(uni_group), c(1L, 0L))
        x[which(group == 0)] <- x[which(group == 0)] - mu
        ni1 <- tapply(group, cluster, sum) # number of obs under trt 2 in each cluster
        if (all(ni1 / ni == 0.5)) {
            warning("The DS ranksum test is not reliable for colateral data where each cluster is equally split between the 2 treatments.")
        }
        ## Calculate S = E(W*|W, g)
        Ni <- rep.int(ni, times = ni)
        S <- sum(group / Ni * (1 + F.prop))
        S <- S / (M + 1)

        ## Calculate E(S)
        ES <- sum(ni1 / ni) / 2

        ## Calculate var(S)

        W <- numeric(M)
        for (i in 1 : M) {
            Wi <- ((M - 1) * group[cluster == i] -
                   sum(ni1[-i] / ni[-i])) *
                Fhat[cluster == i]
            W[i] <- sum(Wi) / (ni[i] * (M+1))
        }
        a <- sum(ni1 / ni)
        EW <- M / (2 * (M + 1)) * (ni1 / ni - a / M)
        varS <- sum((W - EW)^2)
        Z <- (S - ES) / sqrt(varS)

        pval <- switch(alternative, less = pnorm(Z),
                       greater = pnorm(Z, lower.tail = FALSE),
                       two.sided = 2 * min(pnorm(abs(Z)),
                                           pnorm(abs(Z), lower.tail = FALSE)))
        names(Z) <- "Z"
        names(mu) <- "difference in locations"

        result <- list(statistic = Z, p.value = pval, S = S,
                       ES = ES, varS = varS,
                       alternative = alternative, null.value = mu,
                       data.name = DNAME, method = METHOD,
                       nobs = n.obs, nclus = n.clus)

        class(result) <- "ctest"
        return(result)


    } else {
###calculate quantity 2 (using the pooled estimate of F)
        if (!is.null(mu) && mu != 0) {
            warning("comparison between m (m > 2) groups cannot set location shift parameter")
        }
        mu <- 0
        m <- length(unique(group))
        Sj <- numeric(m)
        group.uniq <- unique(group)
        for ( j in 1 : m) {
            for ( i in 1 : M) {
                gik.j <- ifelse(group == j, 1, 0)
                Sj[j] <- Sj[j] + sum(gik.j[cluster == i] *
                                     (1 + F.prop[cluster == i])) / ni[i]
            }
        }
        Sj <- Sj / (M + 1)

        nij <- matrix( 0, m, M)
        for ( i in 1 : m) {
            nij[i, ] <- table(factor(cluster[group == group.uniq[i]], levels = unique(cluster)))
        }
        d <- apply(nij, 1, FUN = function(x) {x / ni})
        ESj <- apply(d, 2, sum) / 2

        What <- matrix(0, m, M)
        a <- t(d)
        for ( i in 1 : M) {
            for ( j in 1 : m) {
                gik.j <- ifelse(group[cluster == i] == j, 1, 0)
                b <- 1 / (ni[i] * (M + 1))
                c <- gik.j * (M - 1)
                d <- sum(a[j, -i])
                What[j, i] <- b * sum((c - d) * Fhat[cluster == i])
            }
        }
        EW <- matrix(0, m, M)
        for ( j in 1 : m) {
            EW[j, ] <- (M / (2 * ( M + 1))) * (a[j, ] - sum(a[j, ]) / M)
        }

        dev.W <- What - EW
        V <- matrix(0, m, m)
        for ( i in 1 : M) {
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
        result <- list(statistic = T, p.value = pval, ngroup = ngrp,
                       df = df, data.name = DNAME, method = METHOD,
                       nobs = n.obs, nclus = n.clus)
        class(result) <- "ctest"
        result
    }
}



clusWilcox_test_ranksum_dd <- function(x, cluster, group,
                                       alternative,
                                       mu, exact, B,
                                       DNAME, METHOD) {
    group <- recoderFunc(group, sort(unique(group)), c(1, 0))
     if (exact == TRUE)  {
        stop("No exact test is available for the DD ranksum test.")
    }
    group.uniq <- length(unique(group))
    if (group.uniq == 1) {
        stop("invalid group variable, should contain 2 groups")
    }

     cgrp0 <- aggregate(group == 0, list(cluster), sum)
     cgrp1 <- aggregate(group == 1, list(cluster), sum)
     cid <- cgrp0[, 1]
     cid <- cid[which(cgrp0[, 2] > 0)]  # Take out clusters with at least one member from group 0
     cid10 <- cid[which((cgrp0[, 2] > 0) & cgrp1[, 2] == 0)] # Further take out clusters with no member from group 1

     ## Take out clusters which have obs from group 1
     data <- cbind(cluster, x, group)
     data <- data[cluster %in% cid, ]
     m <- length(unique(data[, 1]))
     n.obs <- nrow(data)

     rn <- function(dv) {
         cx <- dv[1]
         x <- dv[2]
         ds1 <- data[data[, 3] == 0, ]
         vs1 <- (ds1[, 2] < x) + (ds1[, 2] <= x)
         sl1 <- aggregate(vs1, list(ds1[, 1]), mean)[, 2]
         ds2 <- data[data[, 3] == 1, ]

         if (length(cid10) > 0) {
             ds2 <- rbind(ds2, cbind(cid10, 0, 2))
         }

         vs2 <- (ds2[, 2] < x) + (ds2[, 2] <= x)
         sl2 <- aggregate(vs2, list(ds2[, 1]), mean)[, 2]

         id <- cx %in% cid10
         fg <- (id == FALSE) * (sl1 + sl2) / 2 + (id == TRUE) * (sl1)
         fg[cx] <- 0
         return(fg)
     }

      rst <- function(il) {
        ly <- sum(mat[-which(d0[, 1] == il), -il])
        return(ly)
      }

     rst.add <- function(il) {
         ly <- sum(idadd[-which(d0[, 1] == il)])
         return(ly)
     }

     d0 <- data[data[, 3] == 0, ]
     cd0 <- (d0[, 1])
     nv <- as.vector(table(cd0)[match(cd0, names(table(cd0)))])
     mat <- t((apply(cbind(d0[, 1:2]), 1, rn)))/ (nv * 2)
     idmul <- ((!(unique(cd0) %in% cid10)) + 2 * cgrp1[, 2] * (unique(cd0) %in% cid10)) / 2
     mat <- t(t(mat) * idmul)
     idadd <- (!(cd0 %in% cid10)) / (2 * nv) + (cd0 %in% cid10) / nv

     v1 <- sum(mat) + sum(idadd)
     vd <- apply(cbind(seq(1, m)), 1, rst) + apply(cbind(seq(1, m)), 1, rst.add)
     S <- v1
     ES <- 0.25 * (m + 1) * (m + length(unique(cid10)))
     h <- 1
     test <- (m/m^h) * v1 - ((m - 1)/(m - 1)^h) * vd
     v.test <- var(test)
     v_hat <- (((m^h)^2)/(m - 1)) * v.test
     varS <- ifelse(v_hat == 0, 1e-08, v_hat)
     Z <- (S - ES)/sqrt(varS)

     pval <- switch(alternative, less = pnorm(Z),
                   greater = pnorm(Z, lower.tail = FALSE),
                   two.sided = 2 * min(pnorm(abs(Z)),
                                       pnorm(abs(Z), lower.tail = FALSE)))
    names(Z) <- "Z"
    names(mu) <- "difference in locations"

    result <- list(statistic = Z, p.value = pval, S = S,
                   ES = ES, varS = varS,
                   alternative = alternative, null.value = mu,
                   data.name = DNAME, method = METHOD,
                   nobs = n.obs, nclus = m)

    class(result) <- "ctest"
    result
}
