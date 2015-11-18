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
#' @param id the cluster id for each 
#' observation. If not specified, each observation will 
#' be assigned a distinct id, i.e., no cluster in the data.
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
#' identitical or variable cluster size. When the data
#' is unbalanced, adjusted signed rank statistic is used.
#' Ties are dropped in the test. 
#' @examples
#' data(crsd)
#' clusignrank(z, id, data = crsd)
#' data(crsd.unb)
#' clusignrank(z, id, data = crsd.unb)
#' @author Yujing Jiang 
#' @references
#' Bernard Rosner, Robert J. Glynn, Mei-Ling Ting Lee(2006) 
#' \emph{The Wilcoxon Signed Rank Test for Paired Comparisons of
#'  Clustered Data.} Biometrics, \bold{62}, 185-192.

clusignrank <- 
  function(z, y = NULL, id = NULL, data = parent.frame()) {
    #Calculate number of observations per cluster
    
    METHOD <- "Wilcoxon signed rank test for clutered data"
    pars <- as.list(match.call()[-1])
    
    if(!is.null(pars$data)) {
      z <- data[, as.character(pars$z)]
      DNAME <- paste(pars$z, "from",
                     pars$data)
      
      if(!is.null(y)) {
        y <- data[, as.character(pars$y)]
        if(!is.numeric(y))
          stop("'y' vector must be numeric")
        if(length(z) != length(y))
          stop("lengths of 'z' and 'y' are not the same")
        DNAME <- paste(pars$z, "and",
                       pars$y, "from",
                       pars$data)
        
      } else {
        y <- 0
        DNAME <- paste(pars$z, "from",
                       pars$data)
      }
      z <- z - y
      if(!is.null(id)) {
        id <- data[, as.character(pars$id)]
        data <- data.frame(z, id)
      } else {
        id <- 1 : length(z)
        data <- data.frame(z, id)
      }
    } else {
      if(!is.null(y)) {
        if(!is.numeric(y))
          stop("'y' vector must be numeric")
        if(length(z) != length(y))
          stop("lengths of 'z' and 'y' are not the same")
        z <- y - z
        DNAME <- paste(pars$z, "and",
                       pars$y)
      }
      
      if(is.null(y)) {
        DNAME <- paste(pars$z)
      }
      if(!is.null(id)) {
        data <- data.frame(z, id)
      } else {
        id <- 1 : length(z)
        data <- data.frame(z, id)
      }
    }
    
    
    
    if(!is.numeric(z))
      stop("'z' vector must be numeric")
    
    if(length(z) < 1L) 
      stop("Not enough observations")
    
    
    
    ## Drop the ties
    data <- data[z != 0, ]
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
    if (length(table(g)) != 1) {
      balance = FALSE
    } else {
      balance = TRUE
    }
    
    zrank <- rank(abs(crsd$z))
    crsd1 <- cbind(crsd, zrank)
    signrank <- ifelse(crsd1$z > 0, 1, -1) * crsd1$zrank
    crsd1 <- cbind(crsd1, signrank)
    colnames(crsd1)[4] <- "signrank"  
    if(balance == TRUE){
      T_c <- sum(crsd1$signrank)
      sumrank <- c(by(crsd1$signrank, crsd1$id, sum))
      sumsq <- sum(sumrank ^ 2)
      Var_t <- sumsq
      W_c <- T_c / sqrt(Var_t)
      P_val <- 2 * (1 - pnorm(abs(W_c)))
      
      ADJUST <- FALSE
      names(T_c) <- "T_c"
      names(W_c) <- "W_c"
      names(n) <- "total number of observation"
      names(m) <- "total number of cluster"
      names(Var_t) <- paste("Variance of ", names(T_c))
      result <- list(rstatistic = T_c, vrstatistic = Var_t, statistic = W_c, 
                     p.value = P_val, n = n, m = m, 
                     data.name = DNAME, method = METHOD, 
                     adjusted = ADJUST)
      class(result) <- "ctest"
      return(result)
    } else {
      if (balance == FALSE) {
        sumidrank <- c(by(crsd1$signrank, crsd1$id, sum))
        sumsq <- sum(sumidrank ^ 2)
        meansumrank <- sumidrank / g
        sumsqi <- sum(g ^ 2)
        
        # calculate intraclass correlation between signed ranks within the same cluster
        crsd1$id.f <- as.factor(crsd1$id)
        mod <- lm(signrank ~ id.f, crsd1, y = TRUE)
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
        wi <- g / (vars * (1 + (g - 1) * roscor))
        T_c <- sum(meansumrank * wi)
        sqweightsum <- sum(wi ^ 2 * meansumrank ^ 2)
        Var_t <- sqweightsum
        W_c <-  T_c / (sqrt(Var_t))
        P_val <- 2 * (1 - pnorm(abs(W_c)))
        
        names(T_c) <- "T_cs"
        names(W_c) <- "W_cs"
        ADJUST <- TRUE
        
        names(n) <- "total number of observation"
        names(m) <- "total number of cluster"
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
#' @param id the cluster id for each 
#' observation. If not specified, each observation will 
#' be assigned a distinct id, i.e., no cluster in the data.
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
#' identitical or variable cluster size. When the data
#' is unbalanced, adjusted signed rank statistic is used.
#' Ties are dropped in the test. 
#' @examples
#' data(crsd)
#' clusignrank(z, id, data = crsd)
#' data(crsd.unb)
#' clusignrank(z, id, data = crsd.unb)
#' @author Yujing Jiang 
#' @references
#' Bernard Rosner, Robert J. Glynn, Mei-Ling Ting Lee(2006) 
#' \emph{The Wilcoxon Signed Rank Test for Paired Comparisons of
#'  Clustered Data.} Biometrics, \bold{62}, 185-192.

clusignrank <- 
  function(z, y = NULL, id = NULL, data = parent.frame()){
    #Calculate number of observations per cluster
    
    METHOD <- "Wilcoxon signed rank test for clutered data"
    pars <- as.list(match.call()[-1])
    
    if(!is.null(pars$data)) {
      z <- data[, as.character(pars$z)]
      DNAME <- paste(pars$z, "from",
                     pars$data)
      
      if(!is.null(y)) {
        y <- data[, as.character(pars$y)]
        if(!is.numeric(y))
          stop("'y' vector must be numeric")
        if(length(z) != length(y))
          stop("lengths of 'z' and 'y' are not the same")
        DNAME <- paste(pars$z, "and",
                       pars$y, "from",
                       pars$data)
        
      } else {
        y <- 0
        DNAME <- paste(pars$z, "from",
                       pars$data)
      }
      z <- z - y
      if(!is.null(id)) {
        id <- data[, as.character(pars$id)]
        data <- data.frame(z, id)
      } else {
        id <- 1 : length(z)
        data <- data.frame(z, id)
      }
    } else {
      if(!is.null(y)) {
        if(!is.numeric(y))
          stop("'y' vector must be numeric")
        if(length(z) != length(y))
          stop("lengths of 'z' and 'y' are not the same")
        z <- y - z
        DNAME <- paste(pars$z, "and",
                       pars$y)
      }
      
      if(is.null(y)) {
        DNAME <- paste(pars$z)
      }
      if(!is.null(id)) {
        data <- data.frame(z, id)
      } else {
        id <- 1 : length(z)
        data <- data.frame(z, id)
      }
    }
    
    
    
    if(!is.numeric(z))
      stop("'z' vector must be numeric")
    
    if(length(z) < 1L) 
      stop("Not enough observations")
    
    
    
    ## Drop the ties
    data <- data[z != 0, ]
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
    if (length(table(g)) != 1) {
      balance = FALSE
    } else {
      balance = TRUE
    }
    
    zrank <- rank(abs(crsd$z))
    crsd1 <- cbind(crsd, zrank)
    signrank <- ifelse(crsd1$z > 0, 1, -1) * crsd1$zrank
    crsd1 <- cbind(crsd1, signrank)
    colnames(crsd1)[4] <- "signrank"  
    if(balance == TRUE){
      T_c <- sum(crsd1$signrank)
      sumrank <- c(by(crsd1$signrank, crsd1$id, sum))
      sumsq <- sum(sumrank ^ 2)
      Var_t <- sumsq
      W_c <- T_c / sqrt(Var_t)
      P_val <- 2 * (1 - pnorm(abs(W_c)))
      
      ADJUST <- FALSE
      names(T_c) <- "T_c"
      names(W_c) <- "W_c"
      names(n) <- "total number of observation"
      names(m) <- "total number of cluster"
      names(Var_t) <- paste("Variance of ", names(T_c))
      result <- list(rstatistic = T_c, vrstatistic = Var_t, statistic = W_c, 
                     p.value = P_val, n = n, m = m, 
                     data.name = DNAME, method = METHOD, 
                     adjusted = ADJUST)
      class(result) <- "ctest"
      return(result)
    } else {
      if (balance == FALSE) {
        sumidrank <- c(by(crsd1$signrank, crsd1$id, sum))
        sumsq <- sum(sumidrank ^ 2)
        meansumrank <- sumidrank / g
        sumsqi <- sum(g ^ 2)
        
        # calculate intraclass correlation between signed ranks within the same cluster
        crsd1$id.f <- as.factor(crsd1$id)
        mod <- lm(signrank ~ id.f, crsd1, y = TRUE)
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
        wi <- g / (vars * (1 + (g - 1) * roscor))
        T_c <- sum(meansumrank * wi)
        sqweightsum <- sum(wi ^ 2 * meansumrank ^ 2)
        Var_t <- sqweightsum
        W_c <-  T_c / (sqrt(Var_t))
        P_val <- 2 * (1 - pnorm(abs(W_c)))
        
        names(T_c) <- "T_cs"
        names(W_c) <- "W_cs"
        ADJUST <- TRUE
        
        names(n) <- "total number of observation"
        names(m) <- "total number of cluster"
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
