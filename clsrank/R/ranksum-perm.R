cluswilcox.permutation <-
  function(formula, data = parent.frame(), subset = NULL,
           id = NULL, stratum = NULL, contrasts = NULL,
           na.action = na.omit, n.rep = 1000, ...) {
    # Incoporating clustering effects for the WilcoxonRank Sum Test 
    # for stratified balanced or unbalanced designs for small 
    # sample size, where permutation test is used to simulate 
    # the distribution of test statistics. In addition,
    # one can control for confounding by forming strata which are 
    # defined based on cross-classfication of one or more categorical 
    # covariates which are defined in terms of a single categorical
    # variable denoted by strata.
    # 
    #
    # Requirement:
    #   An ASCII data file with 4 variables per record. The data file
    #   does not have to be sorted in ID order.
    #   The 4 variables need to be space selimited and arranged in the 
    #   following order:
    #   1. id
    #   2. score
    #   3. group (X, Y) indicators : need to be 1 and 2.
    #   4. stratum
    #
    # Args:
    #   1. name of data set.
    #   2. number of stratum.
    #
    # Returns:
    #   1. Clustered Wilcoxon RankSum Statistic
    #   2. Expected Clustered Wilcoxon RankSum Statistic
    #   3. Variance of clustered Wilcoxon RankSum Statistic
    #   4. Z statistic for Clustered Wilcoxon RankSum Statistic
    #   5. P-value for CLustered Wilcoxon RankSum Z Statistic
    
    ## preparation
    ## data <- na.omit(data) ## data[complete.cases(data),]
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "id", "stratum", "na.action"), 
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    z <- model.response(mf, "numeric")
    x <- model.matrix(mt, mf, contrasts)
    group <- x[,2]
    id <- model.extract(mf, "id")
    if (is.null(id)) id <- 1:length(z) ## non-clustered data
    stratum <- model.extract(mf, "stratum")
    if (is.null(stratum)) stratum <- rep(1, length(z))
    
    crd <- data.frame(z, group, id, stratum)
    
    ## group is assumed to take value 1 and 2; could be made
    ## more general
    gr.uniq <- unique(crd$group) ## Record possible groups
    
    crd1 <- crd[with(crd, order(id)), ]
    ## Reorder the observations by the order of id.
    zrank <- rank(crd1$z, na.last = NA)
    ## Compute the rank of z score
    
    crd2 <- cbind(subset(crd1, z != "NA"), zrank)
    ## Add z score to the dataframe crd1
    
    #calculate ranksum within each cluster within each stratum
    
    g <- table(crd2$id)
    ## g is the cluster size
    sumrank <- c(by(crd2$zrank, crd2$id, sum))
    ## Compute rank sum within each cluster
    
    stratum <- c(by(crd2$stratum, crd2$id, mean))
    ## Compute stratum mean of each cluster
    
    group <- c(by(crd2$group, crd2$id, mean))
    ## Compute group mean of each cluster
    
    
    #count number of subunits within cluster size group within each stratum
    
    
    ng.stratum <- as.data.frame(table(g, stratum))
    ## Count the objects for each group and each stratum 
    ng.xy <- as.data.frame(table(g, stratum, group))
    ng.x <- subset(ng.xy, group == 1)
    ## Take out the count for group 1
    ng.y <- subset(ng.xy, group == 2)
    
    colnames(ng.stratum)[3] <- "Ngv"
    colnames(ng.x)[4] <- "mgv"
    
    psumrnk <- cbind(g, stratum, group, sumrank)
    
    str_unique <- sort(unique(stratum))
    str_unique_n <- length(str_unique)
    g_unique <- unique(g)
    g_unique_n <- length(g_unique)
    psumrnk_int <- matrix(0,str_unique_n*g_unique_n, 3)
    
    k <- 1
    for(i in 1:str_unique_n){
      for(j in 1:g_unique_n){
        foo <- sum(subset(psumrnk, stratum == str_unique[i] & g == g_unique[j])[,"sumrank"])
        if(foo){
          psumrnk_int[k,] <- c(str_unique[i], g_unique[j], foo)
          k <- k + 1
        }
      }
    }
    colnames(psumrnk_int) <- c("stratum", "g", "psumrank")
    psumrnk_int <- psumrnk_int[psumrnk_int[,"stratum"]!=0,]
    
    WC <- sum(psumrnk[psumrnk[,"group"]==1,"sumrank"])
  
    ## Matrix to sample permuation from
    sample.base <- as.data.frame(cbind(stratum, g, group, sumrank))
    W.samp <- rep(0, n.rep)
    for(stm in unique(stratum)) {
      for(gz in unique(g)) {
        index <- which(sample.base$stratum == stm & sample.base$g == gz)
        x.size <- ng.x[which(ng.x$stratum == stm & ng.x$g == gz), ]$mgv
        index.samp <- replicate(n.rep, sample(index, x.size))
        x.samp <- sample.base[index.samp, "sumrank"]
        W.samp <- W.samp + colSums(matrix(x.samp, nrow = x.size))
        
      }
    } 
    pval <- 2 * min(ecdf(W.samp)(WC), 1 - ecdf(W.samp)(WC), 0.5)
    result <- list(WC = WC, p.value = pval)
    result  
  }
