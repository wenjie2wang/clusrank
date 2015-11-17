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
##   along with the R package clusrank. If not, see <http://www.gnu.org/licenses/>.
##
################################################################################


#'Wilcoxon Rank Sum Test for Clustered Data
#'This is the sum rank test to compare the means of scores from 
#'two groups for clustered data. The cluster size can be either
#'identitical or variable. Effect of stratification on the test 
#'is also adjusted for if in presence.
#'@param formula   an object of class \code{"formula"} in the
#'form of score ~ group,
#'@param data a optional data frame
#'@param subset an optional vector specifying a 
#'subset of observations to be used.
#'@param id a numeric vector indicating the cluster ids of the scores. 
#'If not specified, each score has its own id, i.e., there is no 
#'cluster in the data.
#'@param stratum a numeric vector indicating the stratum ids of 
#'the scores. The default is that the data has no stratum.
#'@param constrasts  the same as that in \code{lm}. Currently unused.
#'@param na.action a function which indicates what should happen 
#'when the data contains NAs. The  default action is to omit them.
#'@param ... additional arguments, currently ignored.
#'@return  a list with the following components
#'\item{Wc}{Clustered Wilcoxon ranksum statistic.}
#'\item{ExpWc}{Expected value clustered Wilcoxon ranksum statistic.}
#'\item{VarWc}{Variance of clustered Wilcoxon ranksum statistic.}
#'\item{zc}{z statistic for clustered Wilcoxon ranksum statistic.}
#'\item{p.value}{P-value for clustered Wilcoxon ranksum z statistic}
#'
#'@examples
#'data(crd)
#'cluswilcox(z ~ group, data = crd, id = id)
#'data(crd.str)
#'cluswilcox(z ~ group, data = crd.str, id = id, stratum = stratum)
#'
#'@author Yujing Jiang \email{yujing.jiang@uconn.edu}
#'@references
#'Bernard Rosner, Robert J. Glynn, Mei-Ling Ting Lee(2003)
#' \emph{Incorporation of Clustering Effects for the Wilcoxon Rank 
#' Sum Test: A Large-Sample Approach.} Biometrics, \bold{59}, 1089-1098.



cluswilcox <-
  function(formula, data = parent.frame(), subset = NULL,
           id = NULL, stratum = NULL, contrasts = NULL,
           na.action = na.omit, ...) {
  # Incoporating clustering effects for the WilcoxonRank Sum Test 
  # for stratified balanced or unbalanced designs. In addition,
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
    
    METHOD <- "Wilcoxon rank sum test for clutered data"
    
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  
  DNAME <- paste("from", mf$data)
  m <- match(c("formula", "data", "subset", "id", "stratum", "na.action"), 
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  
  DNAME <- paste(names(mf)[1], DNAME, ", with group id:", names(mf)[2])
#   if(!is.null(stratum)) {
#     DNAME <- paste(DNAME, ", with stratum id:", as.character(cl$stratum))
#   }
  mt <- attr(mf, "terms")
  z <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf, contrasts)
  group <- x[,2]
  id <- model.extract(mf, "id")
  if (is.null(id)) id <- 1:length(z) ## non-clustered data
  stratum <- model.extract(mf, "stratum")
  if (is.null(stratum)) stratum <- rep(1, length(z))
  OK <- complete.cases(z)
  crd <- data.frame(z, group, id, stratum)[OK,]
  if (nrow(crd) < 1L) 
    stop("There is not enough observations")
  if (!is.numeric(z))
    stop("The score has to be numeric")
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
  psumrnk_int <- (psumrnk_int[psumrnk_int[,"stratum"]!=0,])
  if(nrow(psumrnk_int) == 1) {
    psumrnk_int <- t(as.data.frame(psumrnk_int))
  } else {
    psumrnk_int <- (as.data.frame(psumrnk_int))
  }
  
  WC <- sum(psumrnk[psumrnk[,"group"]==1,"sumrank"])
  
  expwc <- merge(merge(ng.x, psumrnk_int), ng.stratum)
  ExpWc <- sum(expwc[ , "mgv"] * expwc[ , "psumrank"] / expwc[, "Ngv"])
  
  varwc <- merge(merge(merge(psumrnk, ng.xy), ng.stratum), psumrnk_int)
  varwc_int <- (varwc[ ,"sumrank"] - varwc[ ,"psumrank"] / varwc[, "Ngv"])^2
  VarWc <- cbind(varwc, varwc_int)  
  varwc_pre_final <- VarWc[,"Freq"] * (VarWc[, "Ngv"]-VarWc[, "Freq"]) / (VarWc[,"Ngv"] * (VarWc[,"Ngv"] -1)) * VarWc[,"varwc_int"]
  varwc_final <- sum(varwc_pre_final[VarWc[,"Ngv"] > 1])
  ## Drop the value where there is only 1 cluster in that stratum,
  ## and that cluster contains only 1
  
  zc <- (WC - ExpWc)/sqrt(varwc_final)
  pval <- 2*(1-pnorm(abs(zc)))
  names(WC) <- "Wc"
  names(ExpWc) <- "ExpWc" 
  names(varwc_final) <- "VarWc"
  names(zc) <- "Zc"
  result <- list(rstatistic = WC, erstatistic = ExpWc, 
                 vrstatistic = varwc_final,
                 statistic = zc, p.value = pval,
                 data.name = DNAME, method = METHOD, 
                 balance = balance)
  class(result) <- "ctest"
  result  
}
