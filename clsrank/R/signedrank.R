################################################################################
##
##   R package clusrank by Mei-Ling Ting Lee, Jun Yan, and Yujing Jiang
##   Copyright (C) 2015
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
##   along with the R package reda. If not, see <http://www.gnu.org/licenses/>.
##
################################################################################

#' The Wilcoxon Signed Rank Test for Paired Comparisons of Clustered Data
#' 
#' This is the signed rank test for clustered data where there are 
#' changed score for each subject per cluster. The hypothesis to tset
#' is that the distribution of change in score 
#' is symmetric about zero. Each subunit (instead of the most basic unit)
#' it the unit of change.
#' 
#' @param z  a vector, contains the difference between 
#' score of observations, or the score of objects before treatment.
#' @param y  a vector, the score of objects after treatment if is 
#' not NULL.
#' @param cluster the cluster cluster for each 
#' observation. If not specified, each observation will 
#' be assigned a distinct cluster, i.e., no cluster in the data.
#' @param data  an optional data frame
#' 
#' @return  a list with the following components
#' \item{rstatistic}{Clustered Wilcoxon signed rank statistic.}
#' \item{vrstatistic}{Variance of clustered Wilcoxon signed rank statistic.}
#' \item{statistic}{Standardized clustered Wilcoxon signed rank statistic.}
#' \item{p.value}{P-value for clustered Wilcoxon rank test statistic W_c.}
#' \item{n}{Total number of observations.}
#' \item{m}{Total number of clusters.}
#' \item{data.name}{a character string giving the names of the data.}
#' \item{method}{the type of test applied.}
#' \item{adjusted}{indicator of whether adjusted signed rank statistic is used.}
#' @note This function is able to deal with data with 
#' clusterentitical or variable cluster size. When the data
#' is unbalanced, adjusted signed rank statistic is used.
#' Ties are dropped in the test. 
#' @examples
#' data(data)
#' clusignrank(z, cluster, data = data)
#' data(data.unb)
#' clusignrank(z, cluster, data = data.unb)
#' @author Yujing Jiang 
#' @references
#' Bernard Rosner, Robert J. Glynn, Mei-Ling Ting Lee(2006) 
#' \emph{The Wilcoxon Signed Rank Test for Paired Comparisons of
#'  Clustered Data.} Biometrics, \bold{62}, 185-192.

cluswilcox.test.signedrank <- 
  function(z, cluster,
           alternative = c("two.sided", "less", "greater"),
           DNAME = NULL, METHOD = NULL) {
    #Calculate number of observations per cluster
    
    METHOD <- "Wilcoxon signed rank test for clutered data"

    
    ## Drop the ties

    data <- data.frame(z, cluster)
    data <- data[z != 0, ]
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
      names(Var_t) <- paste("Variance of ", names(T_c))
      result <- list(rstatistic = T_c, vrstatistic = Var_t, statistic = W_c, 
                     p.value = P_val, n = n, cn = m, 
                     data.name = DNAME, method = METHOD, 
                     adjusted = ADJUST)
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

        names(T_c) <- "T_cs"
        names(W_c) <- "W_cs"
        ADJUST <- TRUE
        
        names(n) <- "total number of observations"
        names(m) <- "total number of clusters"
        names(Var_t) <- paste("Variance of ", names(T_c))
        result <- list(rstatistic = T_c, vrstatistic = Var_t, statistic = W_c, 
                       p.value = P_val, n = n, cn = m, 
                       data.name = DNAME, method = METHOD, 
                       adjusted = ADJUST)
        class(result) <- "ctest"
        return(result)
      }
    }
    
    
  }

