################################################################################
##
##   R package clusrank by Mei-Ling Ting Lee, Jun Yan, and Yujing Jiang
##   Copyright (C) 2016
##
##   This file is part of the R package clusrank.
##
##   The R package clusrank is free software: you can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package clusrank is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package clusrank. If not, see <http://www.gnu.org/licenses/>.
##
################################################################################
cluswilcox.test.signedrank.rgl <- function(x, cluster, alternative,
                                           mu, exact, DNAME, METHOD) {
  ## Drop the ties
  x <- x - mu
  data <- data.frame(x, cluster)
  data <- data[x != 0, ]
  cluster.size <- table(data$cluster)
  m <- length(cluster.size)
  n <- nrow(data)
  if (length(table(cluster.size)) != 1) {
    balance = FALSE
  } else {
    balance = TRUE
  }

  xrank <- rank(abs(data$x))
  data <- cbind(data, xrank)
  signrank <- ifelse(data$x > 0, 1, -1) * data$xrank
  data <- cbind(data, signrank)
  colnames(data)[4] <- "signrank"
    if(is.null(exact)) {
        exact <- FALSE
    }

    if(exact == TRUE) {
        if(balance == TRUE ) {            
            Tc <- sum(data$signrank)
            srksum <-  stats::aggregate(signrank ~ cluster, FUN = sum)[, 2]
            Tc.pos <- sum(srksum[srksum > 0])
            p.val.l <- psrkg(Tc.pos, sort(abs(srksum)))
            pval <- switch(alternative,
                          less = p.val.l,
                          greater = 1 - p.val.l,
                          two.sided = 2 * min(p.val.l, 1 - p.val.l))
        } else {
            sumclusterrank <- c(by(data$signrank, data$cluster, sum))
            sumsq <- sum(sumclusterrank ^ 2)
            meansumrank <- sumclusterrank / cluster.size
            sumsqi <- sum(cluster.size ^ 2)
            
                                        # calculate intraclass correlation between signed ranks within the same cluster
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
            
            Tc <- sum(meansumrank * wi)
            msr.pos <- msr.w <-  meansumrank * wi
            msr.pos <- msr.pos[msr.pos > 0]
            Tc.pos <- sum(msr.pos)

            msr.dif <- min(lag(sort(abs(msr.w))))
            msr.pos.int <- as.integer(sort(abs(msr.w) / msr.dif))
            Tc.pos.int <- as.integer(Tc.pos / msr.dif)
            p.val.l <- psrkg(Tc.pos.int, msr.pos.int)
            pval <- switch(alternative,
                          less = p.val.l,
                          greater = 1 - p.val.l,
                          two.sided = 2 * min(p.val.l, 1 - p.val.l))
        }
        names(Tc) <- "rank statistic"
        names(n) <- "total number of observations"
        names(m) <- "total number of clusters"
        names(mu) <- "shift in location"
        result <- list(Rstat = Tc,
                       p.value = pval, n = n, cn = m, null.value = mu,
                       alternative = alternative,
                       data.name = DNAME, method = METHOD)
        class(result) <- "ctest"
        return(result)
        
    } else {
        if(balance == TRUE){
            T_c <- sum(data$signrank)
            sumrank <- c(by(data$signrank, data$cluster, sum))
            sumsq <- sum(sumrank ^ 2)
            Var_t <- sumsq
            W_c <- T_c / sqrt(Var_t)
            P_val <- 2 * (1 - pnorm(abs(W_c)))       
            ADJUST <- FALSE
            names(T_c) <- "rank statistic"
            names(W_c) <- "test statistic"
            names(Var_t) <- "variance of rank statistic"       
            names(n) <- "total number of observations"
            names(m) <- "total number of clusters"
            names(Var_t) <- paste("variance of ", names(T_c))
            names(mu) <- "shift in location"
            result <- list(Rstat = T_c, VRstat = Var_t, statistic = W_c,
                           p.value = P_val, n = n, cn = m, null.value = mu,
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
                                        # calculate intraclass correlation between signed ranks within the same cluster
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
            T_c <- sum(meansumrank * wi)
            sqweightsum <- sum(wi ^ 2 * meansumrank ^ 2)
            Var_t <- sqweightsum
            W_c <-  T_c / (sqrt(Var_t))
            P_val <- switch(alternative,
                            less = pnorm(abs(W_c)),
                            greater = pnorm(abs(W_c), lower.tail = FALSE),
                            two.sided = 2 * min(pnorm(abs(W_c)),
                                                pnorm(abs(W_c), lower.tail = FALSE)))
            
            names(T_c) <- "adjusted rank statistic"
            names(W_c) <- "test statistic"
            
            ADJUST <- TRUE
            names(mu) <- "shift in location"
            
            names(n) <- "total number of observations"
            names(m) <- "total number of clusters"
            names(Var_t) <- paste("variance of ", names(T_c))
            result <- list(Rstat = T_c, VRstat = Var_t, statistic = W_c,
                           p.value = P_val, n = n, cn = m,
                           alternative = alternative,
                           null.value = mu,
                           data.name = DNAME, method = METHOD,
                               adjusted = ADJUST)
            class(result) <- "ctest"
            return(result)
            
        }
    }
}

cluswilcox.test.signedrank.ds <- function(x, cluster, alternative,
                                          mu, DNAME, METHOD) {
    x <- x - mu
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
    temp <- sign(x) * Fcom.vec

    temp <- aggregate(temp ~ cluster, FUN = sum)[, 2]
    VTS <- sum(((niplus - niminus) / ni + (m - 1) / ni * temp) ^ 2 )
    Z <- T / sqrt(VTS)
    P_val <- switch(alternative, less = pnorm(abs(Z)),
                    greater = pnorm(abs(Z), lower.tail = FALSE),
                    two.sided = 2 * min(pnorm(abs(Z)),
                                        pnorm(abs(Z), lower.tail = FALSE)))
    
    names(n) <- "total number of observations"
    names(m) <- "total number of clusters"
    names(Z) <- "test Statistic"
    names(mu) <- "shift in location"
    result <- list(statistic = Z,
                   p.value = P_val, n = n, cn = m,
                   alternative = alternative,
                   null.value = mu,
                 data.name = DNAME, method = METHOD)
    class(result) <- "ctest"
    return(result)
}
