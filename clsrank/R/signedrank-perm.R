#'The Wilcoxon Signed Rank Test for Paired Comparisons of Clustered Data
#'
#'This is the signed rank test for clustered data where there are 
#'changed score for each subject per cluster. The hypothesis to tset
#'is that the distribution of change in score 
#'is symmetric about zero. Each subunit (instead of the most basic unit)
#'it the unit of change.
#'
#'@param z  The name of score for each observation 
#'in the data
#'@param cluster  The name of cluster of cluster for each 
#'observation in the data. If not specified, each score has 
#'its own cluster, which reduces to non-clustereddata.
#'@param data  An optional data frame, list or environment (or object
#'coercible by \code{as.data.frame} to a data frame) containing the
#'variables in the model.  If not found in \code{data}, the 
#'variables are taken from \code{environment(formula)}, typically
#'the environment from which this function is called.
#'
#'@return  a list with the following components
#'\item{T_c}{Clustered Wilcoxon signed ranksum statistic}
#'\item{Var_t}{Variance of clustered Wilcoxon signed ranksum statistic}
#'\item{W_c}{Standardized clustered Wilcoxon signed ranksum statistic}
#'\item{p.value}{P-value for clustered Wilcoxon ranksum test statistic W_c}
#'\item{n}{Total number of observations}
#'\item{m}{Total number of clusters}
#'
#'@note  The function will determine if the dataset is balanced or 
#'unbalanced, and treat them accordingly
#' Zero values can be included in a dataset but are not used in 
#' the analysis. 
#'@examples
#'data(data)
#'clusignrank(z, cluster, data = data)
#'data(data.unb)
#'clusignrank(z, cluster, data = data.unb)
#'
#'@references
#'Bernard Rosner, Robert J. Glynn, Mei-Ling Ting Lee(2006) 
#'\emph{The Wilcoxon Signed Rank Test for Paired Comparisons of
#' Clustered Data.} Biometrics, \bold{62}, 185-192.

cluswilcox.test.signedrank.permutation <- 
  function(z, cluster,  alternative, n.rep = 1000,
           DNAME = NULL, METHOD = NULL){
    
    #Calculate number of observations per cluster
    
    
    data <- data.frame(z, cluster)
    cluster.size <- table(data$cluster)
    m <- length(cluster.size)
    n <- nrow(data)
    if (length(table(cluster.size)) != 1) {
      balance = FALSE
    } else {
      balance = TRUE
    }
    zrank <- rank(abs(data$z))
    data <- cbind(data, zrank)
    signrank <- ifelse(data$z > 0, 1, -1) * data$zrank
    data <- cbind(data, signrank)
    colnames(data)[4] <- "signrank"
    
    if(balance == TRUE){
      T_c <- sum(data$signrank)
      sumrank <- c(by(data$signrank, data$cluster, sum))
      delta <- replicate(n.rep, sample(c(-1, 1), m, TRUE))
      T_c.sim <- colSums(delta * sumrank)
      ecdf.tc <- ecdf(T_c.sim)
      
      
      P_val <- switch(alternative,
                      less = ecdf.tc(abs(T_c)), 
                      greater = 1 - ecdf.tc(abs(T_c)), 
                      two.sided = 2 * min(ecdf.tc(abs(T_c)),     
                                          1 - ecdf.tc(abs(T_c)), 
                                          0.5))
      
      
      names(T_c) <- "rank statistic"
      names(n) <- "total number of observations"
      names(m) <- "total number of clusters"
      result <- list(rstatistic = T_c,  
                     p.value = P_val, 
                     n = n,  cn = m, permutation = TRUE,
                     method = METHOD, data.name = DNAME)
      class(result) <- "ctest"
      return(result)
    } else {
      if (balance == FALSE) {
        sumclusterrank <- c(by(data$signrank, data$cluster, sum))
        sumsq <- sum(sumclusterrank ^ 2)
        meansumrank <- sumclusterrank / cluster.size
        sumsqi <- sum(cluster.size ^ 2)
        
        # calculate intraclass correlation between signed ranks within the same cluster
        data$cluster.f <- as.factor(data$cluster)
        mod <- lm(signrank ~ cluster.f, data, y = TRUE)
        errordf <- mod$df.resclusterual
        errorss <- sum(mod$resclusteruals ^ 2)
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
        wi <- cluster.size / (vars * (1 + (g - 1) * roscor))
        T_c <- sum(meansumrank * wi)
        
        delta <- replicate(n.rep, sample(c(-1, 1), m, TRUE))
        T_c.sim <- colSums(delta * c(meansumrank) * c(wi))
        
        ecdf.tc <- ecdf(T_c.sim)
        
        
        P_val <- switch(alternative,
                        less = ecdf.tc(abs(T_c)), 
                        greater = 1 - ecdf.tc(abs(T_c)), 
                        two.sided = 2 * min(ecdf.tc(abs(T_c)),     
                                            1 - ecdf.tc(abs(T_c)), 
                                            0.5))
        

        names(T_c) <- "rank statistic"
        names(n) <- "total number of observations"
        names(m) <- "total number of clusters"
        result <- list(rstatistic = T_c,  
                       p.value = P_val, 
                       n = n,  cn = m, permutation = TRUE,
                       method = METHOD, data.name = DNAME)
        class(result) <- "ctest"
        return(result)
      }
    }
    
    
    
  }