
cluswilcox.test.ranksum <- function(x, cluster, group, strats, method,
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




