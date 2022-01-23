################################################################################
##
## clusrank: Wilcoxon Rank Tests for Clustered Data
## Copyright (C) 2015-2022  Yujing Jiang, Mei-Ling Ting Lee, and Jun Yan
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
################################################################################
clusWilcox.test.signedrank.rgl.exact <- function(x, cluster,
                                                alternative,
                                                mu, B, DNAME, METHOD) {
    METHOD <- paste0(METHOD, " (random permutation)")
    x <- x - mu
    data <- data.frame(x, cluster)
    data <- data[x != 0, ]
    cluster.size <- table(data$cluster)
    m <- length(cluster.size)
    n <- nrow(data)
    if (length(table(cluster.size)) != 1) {
        balance <- FALSE
    } else {
        balance <- TRUE
    }

    xrank <- rank(abs(data$x))
    data <- cbind(data, xrank)
    signrank <- ifelse(data$x > 0, 1, -1) * data$xrank
    data <- cbind(data, signrank)
    colnames(data)[4] <- "signrank"

    srksum <-  stats::aggregate(signrank ~ cluster, FUN = sum)[, 2]

    T <- sum(data$signrank)

    ind <- replicate(B, sample(c(1, -1), m, TRUE))

    T.ecdf <- ecdf(colSums(ind * srksum))

    pval <- switch(alternative,
                   less = T.ecdf(T),
                   greater = 1 - T.ecdf(T),
                   two.sided = 2 * min(T.ecdf(T), 1 - T.ecdf(T)))

    names(T) <- "T"
    names(n) <- "total number of observations"
    names(m) <- "total number of clusters"
    names(mu) <- "shift in location"
    result <- list(statistic = T,
                   p.value = pval, nobs = n, nclus = m, null.value = mu,
                   alternative = alternative,
                   data.name = DNAME, method = METHOD)
    class(result) <- "ctest"
    return(result)
}





clusWilcox.test.signedrank.rgl <- function(x, cluster, alternative,
                                           mu, exact, B, DNAME, METHOD) {
### Ties are dropped
    if (exact == TRUE && B >= 1)
        return(clusWilcox.test.signedrank.rgl.exact(x, cluster, alternative,
                                                   mu, B, DNAME, METHOD))
    x <- x - mu
    data <- data.frame(x, cluster)
    data <- data[x != 0, ]
    cluster.size <- table(data$cluster)
    m <- length(cluster.size)
    n <- nrow(data)
    if (length(table(cluster.size)) != 1) {
        balance <- FALSE
    } else {
        balance <- TRUE
    }

    xrank <- rank(abs(data$x))
    data <- cbind(data, xrank)
    signrank <- ifelse(data$x > 0, 1, -1) * data$xrank
    data <- cbind(data, signrank)
    colnames(data)[4] <- "signrank"

    if (exact == TRUE) {
        METHOD <- paste0(METHOD, " (exact)")
        if (length(table(cluster)) > 40)
            warning("Number of clusters exceeds 40",
                    " for RGL clustered signed-rank test,",
                    " the exact signed rank test may not work due to overflow.")
        T <- sum(data$signrank)
        srksum <-  stats::aggregate(signrank ~ cluster, FUN = sum)[, 2]
        T.pos <- sum(srksum[srksum > 0])
        p.val.l <- psrkg(T.pos, sort(abs(srksum)))
        pval <- switch(alternative,
                       less = p.val.l,
                       greater = 1 - p.val.l,
                       two.sided = 2 * min(p.val.l, 1 - p.val.l))

        names(T) <- "T"
        names(n) <- "total number of observations"
        names(m) <- "total number of clusters"
        names(mu) <- "shift in location"
        result <- list(statistic = T,
                       p.value = pval, nobs = n, nclus = m, null.value = mu,
                       alternative = alternative,
                       data.name = DNAME, method = METHOD)
        class(result) <- "ctest"
        return(result)

    } else {
        if (balance == TRUE){
            T <- sum(data$signrank)
            sumrank <- c(by(data$signrank, data$cluster, sum))
            sumsq <- sum(sumrank ^ 2)
            VarT <- sumsq
            Z <- T / sqrt(VarT)
            pval <- switch(alternative,
                           less = pnorm(Z),
                           greater = pnorm(Z, lower.tail = FALSE),
                           two.sided = 2 * min(pnorm(abs(Z)),
                                               pnorm(abs(Z), lower.tail = FALSE)))
            ADJUST <- FALSE
            names(T) <- "T"
            names(Z) <- "Z"
            names(VarT) <- "variance of T"
            names(n) <- "total number of observations"
            names(m) <- "total number of clusters"
            names(VarT) <- paste("variance of ", names(T))
            names(mu) <- "shift in location"
            result <- list(Rstat = T, VRstat = VarT, statistic = Z,
                           p.value = pval, nobs = n, nclus = m, null.value = mu,
                           alternative = alternative,
                           data.name = DNAME, method = METHOD,
                           adjusted = ADJUST)
            class(result) <- "ctest"
            return(result)
        } else {
            sumclusterrank <- c(by(data$signrank, data$cluster, sum))
            sumsq <- sum(sumclusterrank ^ 2)
            meansumrank <- sumclusterrank / cluster.size
            sumsqi <- sum(cluster.size ^ 2)
            ## calculate intraclass correlation between signed ranks within the same cluster
            data$cluster.f <- as.factor(data$cluster)
            mod <- lm(signrank ~ cluster.f, data, y = TRUE)
            errordf <- mod$df.residual
            errorss <- sum(mod$residuals ^ 2)
            modeldf <- n - errordf - 1
            modelss <- sum((mod$y - mean(mod$y)) ^ 2) - errorss
            sumi <- n
            m0 <- (sumi - (sumsqi / sumi)) / (m - 1)
            totalss <- errorss + modelss
            totaldf <- errordf + modeldf
            wthms <- errorss / errordf
            betms <- modelss / modeldf
            vars <- totalss / totaldf
            s2b <- (betms - wthms) / m0
            s2w <- wthms
            rosglm <- s2b / (s2b + s2w)
            if (rosglm < 0) {
                rosglm = 0
            }
            ros <- rosglm
            roscor <- ros * (1 + (1 - ros ^ 2) / (m - 2.5))
            if (roscor > 1) {
                roscor = 1
            }
            wi <- cluster.size / (vars * (1 + (cluster.size - 1) * roscor))
            T <- sum(meansumrank * wi)
            sqweightsum <- sum(wi ^ 2 * meansumrank ^ 2)
            VarT <- sqweightsum
            Z <-  T / (sqrt(VarT))
            pval <- switch(alternative,
                           less = pnorm(Z),
                           greater = pnorm(Z, lower.tail = FALSE),
                           two.sided = 2 * min(pnorm(abs(Z)),
                                               pnorm(abs(Z), lower.tail = FALSE)))

            names(T) <- "T"
            names(Z) <- "Z"

            ADJUST <- TRUE
            names(mu) <- "shift in location"
            names(n) <- "total number of observations"
            names(m) <- "total number of clusters"
            names(VarT) <- paste("variance of ", names(T))
            result <- list(Rstat = T, VRstat = VarT, statistic = Z,
                           p.value = pval, nobs = n, nclus = m,
                           alternative = alternative,
                           null.value = mu,
                           data.name = DNAME, method = METHOD,
                           adjusted = ADJUST)
            class(result) <- "ctest"
            return(result)

        }
    }
}

clusWilcox.test.signedrank.ds.exact.1 <- function(x, cluster) {
    order.c <- order(cluster)
    x <- x[order.c]
    cluster <- cluster[order.c]
    csize <- as.vector(table(cluster))
    m <- length(csize)
    cid <- as.numeric(names(csize))
    n <- sum(csize)
    plus <- as.numeric(x > 0)
    minus <- as.numeric(x < 0)
    niplus <- aggregate(plus ~ cluster, FUN = sum)[, 2]
    niminus <- aggregate(minus ~ cluster, FUN = sum)[, 2]
    ni <- table(cluster)


    csize.cum <- cumsum(csize)
    csize.cum <- c(0, csize.cum)
    cluster <- recoderFunc(cluster, unique(cluster), c(1 : m))


    Ftot.vec  <- Ftot_vec(abs(x), cluster, csize, n, m)
    Fi.vec <- Fi_vec(abs(x), cluster, csize, n, m)
    Fcom.vec <- Fcom_vec(abs(x), cluster, csize, n, m)

    T <- sum((niplus - niminus) / ni) +
        sum((sign(x) * (Ftot.vec - Fi.vec)) / rep.int(ni, times = ni))
    T

}


clusWilcox.test.signedrank.ds.exact <- function(x, cluster, alternative,
                                               B,
                                               mu, DNAME, METHOD) {
    METHOD <- paste0(METHOD, " (random permutation)")
    T.vec <- rep(NA, B)
    x <- x - mu
    n.obs <- length(x)
    n.clus <- length(unique(cluster))
    for ( i in 1 : B) {
        sgn.samp <- sample(c(-1, 1), n.obs, TRUE)
        x.samp <- abs(x) * sgn.samp
        T.vec[i] <- clusWilcox.test.signedrank.ds.exact.1(x.samp, cluster)
    }
    T <- clusWilcox.test.signedrank.ds.exact.1(x, cluster)
    t.ecdf <- ecdf(T.vec)

    pval <- switch(alternative, less = t.ecdf(T),
                   greater = 1 - t.ecdf(T),
                   two.sided = 2 * min(t.ecdf(T), 1 - t.ecdf(T)))

    names(n.obs) <- "total number of observations"
    names(n.clus) <- "total number of clusters"
    names(T) <- "T"
    names(mu) <- "shift in location"
    result <- list(statistic = T,
                   p.value = pval, nobs = n.obs, nclus = n.clus,
                   alternative = alternative,
                   null.value = mu,
                   data.name = DNAME, method = METHOD)
    class(result) <- "ctest"
    return(result)


}


clusWilcox.test.signedrank.ds <- function(x, cluster, alternative, exact, B,
                                          mu, DNAME, METHOD) {


    if (exact == TRUE & B >= 1){
        return(clusWilcox.test.signedrank.ds.exact(x, cluster, alternative,
                                                  B, mu, DNAME, METHOD))
    }

    if (exact == TRUE & B == 0) {
        warning("Exact test is not provided for DS method for signed rank test, large-sample test will be carried out")
    }

    Xij <- x - mu
    ni <- as.vector(table(cluster))


    g <- length(ni)
    n <- sum(ni)
    cni <- cumsum(ni)
    cni <- c(0, cni)
    Fi <- function(x, i) {
        Xi <- Xij[(cni[i] + 1):(cni[i + 1])]
        (sum(abs(Xi) <= x) + sum(abs(Xi) < x))/(2 * ni[i])
        }
    Ftot <- function(x) {
        st <- 0
        for (i in 1:g) st <- st + Fi(x, i)
        return(st)
    }
    Fcom <- function(x) {
        st <- 0
        for (i in 1:g) st <- st + Fi(x, i) * ni[i]
        return(st/n)
    }
    TS <- VTS <- 0
    for (i in 1:g) {
        Xi <- Xij[(cni[i] + 1):(cni[i + 1])]
        first <- (sum(Xi > 0) - sum(Xi < 0))/length(Xi)
        second <- 0
        third <- 0
        for (x in Xi) {
            second <- second + sign(x) * (Ftot(abs(x)) -
                                          Fi(abs(x), i))
                third <- third + sign(x) * Fcom(abs(x))
        }
        TS <- TS + first + second/length(Xi)
        VTS <- VTS + (first + (g - 1) * third/length(Xi))^2
    }

    Z <- TS / sqrt(VTS)
    pval <- switch(alternative, less = pnorm(Z),
                   greater = pnorm(Z, lower.tail = FALSE),
                   two.sided = 2 * min(pnorm(abs(Z)),
                                       pnorm(abs(Z), lower.tail = FALSE)))

    names(n) <- "total number of observations"
    names(g) <- "total number of clusters"
    names(Z) <- "Z"
    names(mu) <- "shift in location"
    result <- list(statistic = Z,
                   srstat = T,
                   vsrstat = VTS,
                   p.value = pval, nobs = n, nclus = g,
                   alternative = alternative,
                   null.value = mu,
                   data.name = DNAME, method = METHOD)
    class(result) <- "ctest"
    return(result)
}
