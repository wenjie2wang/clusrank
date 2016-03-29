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

cluswilcox.test.ranksum.rgl <- function(x, cluster, group, strats, 
           alternative, mu, DNAME = NULL, METHOD = NULL) {
  ## Check if any cluster contains subjects from both groups
  ## Make use of the number of group, use methematical way to check
  ## Assuming the index of group is in sequence
  n.clus <- length(unique(cluster))
  n.group <- length(unique(group))
  
  if(n.group > 2) {
    stop("RGL's ranksum test only allows comparison between 2 groups")
  }
  
  n.check <- length(unique(cluster * n.group + group))
  if(n.clus != n.check) {
    ## If there is at least a cluster containing subjects from both groups
    ## Use the rgl's subunit method
    cluswilcox.test.ranksum.rgl.sub(x, cluster, group,
                                alternative, mu, DNAME = NULL, METHOD = NULL)
  } else {
    ## Otherwise, use the rgl's method which allows stratification
    cluswilcox.test.ranksum.rgl.str(x, cluster, group,
                                    alternative, mu, DNAME = NULL, METHOD = NULL)
  }
}

cluswilcox.test.ranksum.rgl.str <- function(x, cluster, group, strats, 
                                        alternative, mu, DNAME = NULL, METHOD = NULL) {

  data <- as.data.frame(cbind(x, cluster, group, strats))
  ## group is assumed to take value 1 and 2; could be made
  ## more general
  group.uniq <- unique(data$group) ## Record possible groups
  
  data <- data[with(data, order(cluster)), ]
  ## Reorder the observations by the order of id.
  xrank <- rank(data$x, na.last = NA)
  ## Compute the rank of x score
  
  data <- cbind(data, xrank)
  ## Add z score to the dataframe crd1
  
  ### calculate ranksum within each cluster within each strats
  
  cluster.size <- table(data$cluster)
  ## cluster.size is the cluster size
  sumrank <- c(by(data$xrank, data$cluster, sum))
  ## Compute rank sum within each cluster
  
  strats <- c(by(data$strats, data$cluster, mean))
  ## Compute strats  of each cluster
  
  group <- c(by(data$group, data$cluster, mean))
  ## Compute group  of each cluster
  
  
  #count number of subunits within cluster size group within each strats
  
  n.csize.strats <- as.data.frame(table(cluster.size, strats))
  ## Count the objects for each group and each strats
  n.csize.xy <- as.data.frame(table(cluster.size, strats, group))
  n.csize.x <- subset(n.csize.xy, group == 1)
  ## Take out the count for group 1
  n.csize.y <- subset(n.csize.xy, group == 2)
  
  ## Ngv is the number of clusters which have the same
  ## subunit, controlled by stratum.
  colnames(n.csize.strats)[3] <- "Ngv"
  
  ## mgv is number of clusters which have the same
  ## subunit in group x, ontrolled by stratum.
  colnames(n.csize.x)[4] <- "mgv"
  
  psumrnk <- cbind(cluster.size, strats, group, sumrank)
  
  str.uniq <- sort(unique(strats))
  n.str.uniq <- length(str.uniq)
  cluser.size.uniq <- unique(cluster.size)
  n.cluser.size.uniq <- length(cluser.size.uniq)
  
  if(n.cluser.size.uniq > 1L) {
    balance <- FALSE
  } else {
    balance <- TRUE
  }
  
  psumrnk_int <- matrix(0,  n.str.uniq * cluser.size.uniq , 3)
  
  k <- 1
  for(i in 1:n.str.uniq ){
    for(j in 1:n.cluser.size.uniq){
      ind <- strats == str.uniq[i] & cluster.size == cluser.size.uniq[j]
      foo <- sum(psumrnk[ind, "sumrank"])
      if (foo) {
        psumrnk_int[k,] <- c(str.uniq[i], cluser.size.uniq[j], foo)
        k <- k + 1
      }
    }
  }
  colnames(psumrnk_int) <- c("strats", "cluster.size", "psumrank")
  psumrnk_int <- as.data.frame(psumrnk_int[psumrnk_int[, "strats"] != 0, , drop = FALSE])
  
  WC <- sum(psumrnk[psumrnk[,"group"] == 1,"sumrank"])
  
  expwc <- merge(merge(n.csize.x, psumrnk_int), n.csize.strats)
  ExpWc <- sum(expwc[ , "mgv"] * expwc[ , "psumrank"] / expwc[, "Ngv"])
  
  varwc <- merge(merge(merge(psumrnk, n.csize.xy), n.csize.strats), psumrnk_int)
  varwc_int <- (varwc[ ,"sumrank"] - varwc[ ,"psumrank"] / varwc[, "Ngv"])^2
  VarWc <- cbind(varwc, varwc_int)
  varwc_pre_final <- VarWc[,"Freq"] * (VarWc[, "Ngv"]- VarWc[, "Freq"]) / (VarWc[,"Ngv"] * (VarWc[,"Ngv"] -1)) * VarWc[,"varwc_int"]
  varwc_final <- sum(varwc_pre_final[VarWc[,"Ngv"] > 1])
  ## Drop the value where there is only 1 cluster in that strats,
  ## and that cluster contains only 1
  
  zc <- (WC - ExpWc)/sqrt(varwc_final)
  pval <- switch(alternative,
                 less = pnorm(abs(zc)),
                 greater = pnorm(abs(zc), lower.tail = FALSE),
                 two.sided = 2 * min(pnorm(abs(zc)),
                                     pnorm(abs(zc), lower.tail = FALSE)))
  
  names(WC) <- "Rank sum statistic"
  names(ExpWc) <- "Expected value of rank sum statistic"
  names(varwc_final) <- "Variance of rank sum statistic"
  names(zc) <- "Test statistic"
  names(mu) <- "difference in locations"
  result <- list(rstatistic = WC, erstatistic = ExpWc,
                 vrstatistic = varwc_final,
                 statistic = zc, p.value = pval,
                 alternative = alternative,
                 null.value = mu,
                 data.name = DNAME,
                 method = METHOD,
                 balance = balance)
  class(result) <- "ctest"
  result
}


cluswilcox.test.ranksum.rgl.sub <- function(x, cluster, group, 
                                            alternative, mu, DNAME = NULL, METHOD = NULL) {
  
  n.obs <- length(x)
  n.cluster <- length(unique(cluster))
  ##recode group into 0 and 1 to follow the paper
  group <- recoderFunc(group, unique(group), c(0, 1))
  data <- as.data.frame(cbind(x, cluster, group))

  group.uniq <- c(0, 1) ## Record possible groups
  
  data <- data[with(data, order(cluster)), ]
  ## Reorder the observations by the order of id.


  ### calculate ranksum within each cluster within each strats
  
  cluster.size <- table(data$cluster)
  gmax <- max(as.numeric(names(cluster.size)))
  if(unique(cluster.size) > 1) {
    balance <- FALSE
  } else {
    balance <- TRUE
  }
  
  if(balance == TRUE) {
    
    xrank <- rank(data$x, na.last = NA)
    ## Compute the rank of x score
    
    data <- cbind(data, xrank)
    ## cluster.size is the cluster size
    sumrank <- c(by(data$xrank, data$cluster, sum))
    ## Compute rank sum within each cluster
    
    meanrank <- c(by(data$xrank, data$cluster, mean))
    ## Compute rank mean within each cluster
    
    cluster.gr1 <- c(by(data$group, data$cluster, function(x) length(x[x == 0])))
    ## Count number of subject belongs to group 1 in each cluster
      
    gmax <- g
    n.cl <- numeric(g)
    n.cl.gr1 <- numeric(g)
    ## Record number of clusters with each cluster size
    ## Record number of clusters with all possible
    ## number of members belongs to group1
    
    n.cluster.gr1 <- table(cluster.gr1)
    n.cl.gr1[as.numeric(names(n.cluster.gr1))] <- n.cluster.gr1
    
    n.cl[gmax] <- n.cluster
    ## Since this is a balance design, all clusters have the same size.
    
    wc <- sum(data$xrank * group)
    
    sum.Ngr1 <- sum(which(n.cl.gr1 != 0) * n.cl.gr1[n.cl.gr1 != 0])
    ## Count the total number of subjects in group 0
    
    sum2.Ngr1 <- sum(which(n.cl.gr1 != 0) ^ 2 * n.cl.gr1[n.cl.gr1 != 0])
    
    EWc <- (gmax * n.obs + 1) / 2 * sum.Ngr1
    VarQ <- (sum2.Ngr1 - sum.Ngr1 ^ 2 / n.obs) / n.obs
    sumssb <- sum((xrank - gmax * (gmax * n.obs + 1) / 2) ^ 2)
    ssb <- sumssb / n.obs
    
    varwc <- n.obs * (n.obs / (n.obs - 1) * VarQ * ssb / g.max ^ 2)
    stwc <- sqrt(varwc)
    zc <- (wc - EWc) / stwc
    pval <- switch(alternative,
                   less = pnorm(abs(zc)),
                   greater = pnorm(abs(zc), lower.tail = FALSE),
                   two.sided = 2 * min(pnorm(abs(zc)),
                                       pnorm(abs(zc), lower.tail = FALSE)))
    
    names(wc) <- "Subunit-specific Clustered Wilcoxon RankSum Statistic"
    names(EWc) <- "Expected value of rank sum statistic"
    names(stwc) <- "Standard Deviation of rank sum statistic"
    names(zc) <- "Test statistic"
    names(mu) <- "difference in locations"
    result <- list(rstatistic = WC, erstatistic = ExpWc,
                   vrstatistic = varwc_final,
                   statistic = zc, p.value = pval,
                   alternative = alternative,
                   null.value = mu,
                   data.name = DNAME,
                   method = METHOD,
                   balance = balance)
    class(result) <- "ctest"
    result
  } else {
    ## Unbalanced design
    cluster.s <- c(by(rep(1, n.obs), data$cluster, length))
    ## Count cluster size
    cluster.gr1 <- c(by(data$group, data$cluster, function(x) length(x[x == 0])))
    ## Count number of subject belongs to group 1 in each cluster
    
    n.qplus <- n.cluster <- numeric(max(as.numeric(names(cluster.size))))
    ## names(cluster.size) records the all the possibilities of cluster size
    n.cluster[as.numeric(names(cluster.size))] <- cluster.size
    ## Rcord the number of clusters of each size
    n.qplus0 <- c(by(cluster.gr1, cluster.s, sum))
    n.qplus[as.numeric(names(n.qplus0))] <- n.qplus0
    ## Record the total number of exposed subjects within each cluster size
    
    gNg <- which(n.cluster != 0) * n.cluster[n.cluster != 0]
    ## Number of obs with a specific cluster size
    deno <- gNg + 1
    deno[deno == 1] <- 0
    ## Denominator for some later quantity
    
    clus.s <- cbind(as.numeric(names(cluster.size)), cluster.size)
    colnames(clus.s) <- c("cluster", "cluster.size")
    data <- merge(data, clus.s, by = "cluster")
    
    ## Order data by cluster.size
    data <- data[with(data, order(data$cluster.size)), ]
    ## Assign rank of observations within each cluster size
    ## xrank is a list which contains rank for each cluster size.
    xrank <- c(by(data$x, data$cluster.size, rank))
    meanrank <- lapply(xrank, mean)
    meanrank <- unlist(meanrank)
    ## Prepare meanrnak to merge with data
    meanrank <- cbind(as.numeric(names(meanrank)), meanrank)
    colnames(meanrank) <- c("cluster.size", "meanrank")
    data <- merge(data, meanrank, by = "cluster.size")
    
    xrank <- unlist(xrank)
    data <- cbind(data, xrank)
    
    sumvarr <- varr <- numeric(gmax)
    for( i in 1 : gmax) {
      sumvarr[i] <- sum(((data$xrank - data$meanrank) ^ 2)
                        [(data$cluster.size == i)])
    }
    sumvarr[!is.numeric(sumvarr)] <- 0
    varr[gNg != 0] <- (sumvarr / gNg)[gNg != 0]
    
    
    wc <- numeric(gmax)
    ## Record the wc stat for each cluster size
    for( i in 1 : gmax) {
      wc[i] <- sum((data$xrank * data$group)[data$cluster.size == i])
    }
    wc[!is.numeric(wc)] <- 0
    
    deno0 <- rep(deno, times = gNg)
    pscore <- qnorm(data$xrank / deno0)
    data <- cbind(data, pscore)
    aovout <- aov(pscore ~ as.factor(data$group))
    aovout <- summary(aovout)
    
    errorms <- aovout[[1]][["Mean Sq"]][2]
    modelms <- aovout[[1]][["Mean Sq"]][1]
    
    ## Compute theta
    uc <- theta <- numeric(gmax)
    uc <- wc - n.qplus * (n.qplus + 1) / 2
    theta <- uc / (n.qplus * (gNg - n.qplus))
    
    sumn <- n.obs
    sumnsq <- sum(n.cluster * c(1 : gmax) ^ 2)
    
    g0 <- (sumn - sumnsq / sumn) / (n.obs - 1)
    sigsqa <- (modelms - errorms) / g0
    rhoh <- sigsqa / (sigsqa + errorms)
    if(rhoh < 0) {
      rhoh <- 0
    }
    rhoh <- rhoh * ( 1 + ( 1 - rhoh ^ 2) / (2 * (n.obs - 4)))
    rhor <- (6 / pi) * asin(rhoh / 2)
    
    ssb <- c(1 : gmax) * varr * ( 1 + (c(1 : gmax) - 1) * rhor)
    ssw <= varr * (1 - rhor)
    
    sume0 <- (cluster.s - cluster.gr1) * cluster.gr1
    sume <- e <-  numeric(gmax)
    for( i in 1 : gmax) {
      sume[i] <- sum(sume0[cluster.s == i])
      e[i] <- sume / n.cluster[i]
    }
    sume[!is.numeric(sume)] <- e[!is.numeric(e)] <- 0
    
    ## Compute variance of observations from group 1
    sumq <- sumqsq <- numeric(gmax)
    for( i in 1 : gmax) {
      sumq[i] <- sum(cluster.gr1[cluster.s == i])
      if(n.cluster[i] != 0) {
        sumqsq[i] <- sumq[i] ^ 2 / n.cluster[i]
      }
    }
    
    sumqq <- varq <- numeric(gmax)
    for( i in 1 : gmax) {
      sumqq[i] <- sum((cluster.gr1[cluster.s == i]) ^ 2)
      if(n.cluster[i] != 0) {
        varq[i] <- (sumqq[i] - sumqsq[i]) / n.cluster[i]
      }
    }
    
    vartheta <- numeric(gmax)
    for( i in 1 : gmax) {
      if(n.cluster[i] >= 2) {
        vartheta[i] <- n.cluster[i] * 
          ((n.cluster[i] / (n.cluster[i] - 1)) * 
             varq[i] * ssb[i] / i ^ 2 + e[i] * ssw[i] / i) /
          (n.qplus[i] * (gNg[i] - n.qplus[i])) ^ 2
      } else if( n.cluster[i] == 1) {
        vartheta[i] <- n.cluster[i] * (e[i] * ssw[i] / i) /
          (n.qplus[i] *(gNg[i] - n.qplus[i])) ^ 2
      }
    }
    
    weight <- numeric(gmax)
    weight[vartheta != 0] <- 1 / vartheta[vartheta != 0]
    
    sumwtheta <- sum((weight * theta)[n.cluster != 0])
    sumwt <- sum(weight[n.cluster != 0])
    
    theta <- sumwtheta / sumwt
    sdtheta <- sqrt(1 / sumwt)
    zstat <- (theta - 1 / 2) / sdtheta
    pval <- switch(alternative,
                   less = pnorm(abs(zstat)),
                   greater = pnorm(abs(zstat), lower.tail = FALSE),
                   two.sided = 2 * min(pnorm(abs(zstat)),
                                       pnorm(abs(zstat), lower.tail = FALSE)))    
    names(theta) <- "Estimated Probability That a Random 
    Exposed Subunit Will Have a Higher Score Than a Random Unexposed Subunit"
    names(sdtheta) <- "Standard Deviation of Theta"
    names(zstat) <- "Z statistic for Subunit-specific Clustered 
    Wilcoxon RankSum Statistic"
    names(pval) <- "P-value for SUbunit-specific Clustered
    Wilcoxon RankSum Z Statistic"
    
 } 
  

}


