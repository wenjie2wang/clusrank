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
##    with the R package clusrank. If not, see <http://www.gnu.org/licenses/>.
##
################################################################################


#'Wilcoxon Rank Sum Test for Clustered Data
#'
#'This is the rank sum test for clustered data.
#' The cluster size can be either
#'identitical or variable. Effect of stratification on the test
#'is also adjusted for if in presence.
#'@param x a numeric vector of data values.
#'@param cluster a numeric vector indicating the cluster ids of the scores.
#'If not specified, each score has its own id, i.e., there is no
#'cluster in the data.
#'@param group a numeric vector indicating group id.
#'@param strats a numeric vector indicatg the strats ids of
#'the scores. The default is that the data has no stratum.
#'@param alternative a character string specifying the
#' alternative hypothesis, must be one of "two.sided" (default),
#'  "greater" or "less". You can specify just the initial letter.
#'@param mu null value of the hypothesis
#'@param DNAME a character string, inheritated from
#' \code{cluswilcox.test.model}, for result output.
#' @param METHOD a character string, inheritated from
#' \code{cluswilcox.test.model}, for result output.

#'@return  a list with the following components
#'\item{rstatistic}{Clustered Wilcoxon ranksum statistic.}
#'\item{erstatistics}{Expected value clustered Wilcoxon ranksum statistic.}
#'\item{vrstatistics}{Variance of clustered Wilcoxon ranksum statistic.}
#'\item{statistics}{the value of the test statistic.}
#'\item{p.value}{the p-value for the test}
#'\item{data.name}{a character string giving the names of the data.}
#'\item{method}{the name of the method}
#'\item{balance}{a logical, indicating if the data is balanced.}
#'@importFrom  stats pnorm
#'@examples
#'data(crd)
#'cluswilcox.test(z ~ group(group) + cluster(id), data = crd)
#'data(crdStr)
#'cluswilcox.test(z ~ group(group) + cluster(id) + stratum(stratum), data = crdStr)
#'
#'@author Yujing Jiang
#'@references
#'Bernard Rosner, Robert J. Glynn, Mei-Ling Ting Lee(2003)
#' \emph{Incorporation of Clustering Effects for the Wilcoxon Rank
#' Sum Test: A Large-Sample Approach.} Biometrics, \bold{59}, 1089-1098.



cluswilcox.test.ranksum <-
  function(x, cluster, group, strats, method,
           alternative, mu, DNAME = NULL, METHOD = NULL) {
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
    #   4. strats
    #
    # Args:
    #   1. name of data set.
    #   2. number of strats.
    #
    # Returns:
    #   1. Clustered Wilcoxon RankSum Statistic
    #   2. Expected Clustered Wilcoxon RankSum Statistic
    #   3. Variance of clustered Wilcoxon RankSum Statistic
    #   4. Z statistic for Clustered Wilcoxon RankSum Statistic
    #   5. P-value for CLustered Wilcoxon RankSum Z Statistic




  METHOD <- "Wilcoxon rank sum test for clutered data"
  if(method == "rgl") {
    cluswilcox.test.ranksum.rgl(x, cluster, group, strats,
             alternative, mu, DNAME, METHOD)
  } 
  if(method == "ds") {
    cluswilcox.test.ranksum.ds(x, cluster, group, 
                               alternative, mu, DNAME, METHOD)
  }
 

  
}
