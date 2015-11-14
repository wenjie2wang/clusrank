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
#'@param id  The name of id of cluster for each 
#'observation in the data. If not specified, each score has 
#'its own id, which reduces to non-clustereddata.
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
#'data(crsd)
#'clusignrank(z, id, data = crsd)
#'data(crsd.unb)
#'clusignrank(z, id, data = crsd.unb)
#'
#'@references
#'Bernard Rosner, Robert J. Glynn, Mei-Ling Ting Lee(2006) 
#'\emph{The Wilcoxon Signed Rank Test for Paired Comparisons of
#' Clustered Data.} Biometrics, \bold{62}, 185-192.

clusignrank.s <- 
  function(z, id = NULL, data = parent.frame(), n.sim = 1000){
    
    #Calculate number of observations per cluster
    pars <- as.list(match.call()[-1])
    data <- data[, as.character(pars$z) !=0]
    z <- data[,as.character(pars$z)]
    id <- pars$id
    if (is.null(id)) {
      id <- 1:length(z)
    } else {
      id <- data[,as.character(id)]
    }
    crsd <- data.frame(z, id)
    g <- table(crsd$id)
    m <- length(g)
    n <- nrow(crsd)
    
    zrank <- rank(abs(crsd$z))
    crsd1 <- cbind(crsd, zrank)
    signrank <- ifelse(crsd1$z > 0, 1, -1) * crsd1$zrank
    crsd1 <- cbind(crsd1, signrank)
    colnames(crsd1)[4] <- "signrank"  
      T_c <- sum(crsd1$signrank)
      sumrank <- c(by(crsd1$signrank, crsd1$id, sum))
      delta <- replicate(n.sim, sample(c(-1, 1), m, TRUE))
      T_c.sim <- colSums(delta * sumrank)
      P_val <- 2 * min(1 - ecdf(T_c.sim)(T_c), ecdf(T_c.sim)(T_c), 0.5)
      result <- list(T_C = T_c,  
                     p.value = P_val, n = n, m = m)
      return(result)
    
    
    
  }