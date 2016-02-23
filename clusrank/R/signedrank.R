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
##   along with the R package reda. If not, see <http://www.gnu.org/licenses/>.
##
################################################################################

#' The Wilcoxon Signed Rank Test for Clustered Data
#'
#' Performs one-sample Wilcoxon test on vectors of data using
#' large sample.
#'
#' @param x  numeric vector of data values. Non-finite (e.g.,
#' infinite or missing) values will be omitted.
#' @param cluster numeric or charater vector, the id of clusters.
#'  If not specified, each observation will
#' be assigned a distinct cluster, i.e., no cluster in the data.
#' @param alternative a character string specifying the
#' alternative hypothesis, must be one of "two.sided" (default),
#'  "greater" or "less". You can specify just the initial letter.
#'@param mu null value of the hypothesis
#' @param DNAME a character string, inheritated from
#' \code{cluswilcox.test.numeric}, for result output.
#' @param METHOD a character string, inheritated from
#' \code{cluswilcox.test.numeric}, for result output.
#' @param ...  Further arguments to be passed to or from methods.
#' @return  a list with class "\code{ctest}" containing
#' the following components:
#' \item{rstatistic}{the value of the signed rank statistic
#'  with a name describing it.}
#' \item{vrstatistic}{Variance of \code{rstatistic}.}
#' \item{statistic}{the value of the test statistic.}
#' \item{p.value}{the p-value for the test.}
#' \item{n}{Total number of observations.}
#' \item{cn}{Total number of clusters.}
#' \item{data.name}{a character string giving the names of the data.}
#' \item{method}{the type of test applied.}
#' \item{adjusted}{indicator of whether adjusted signed rank statistic is used.}
#'  @note This function is able to deal with data with
#' clusterentitical or variable cluster size. When the data
#' is unbalanced, adjusted signed rank statistic is used.
#' Ties are dropped in the test.
#' @importFrom  stats pnorm lm
#' @examples
#' data(crsd)
#' cluswilcox.test(z, cluster = id, data = crsd)
#' data(crsdUnb)
#' cluswilcox.test(z, cluster = id, data = crsdUnb)
#' @author Yujing Jiang
#'
#' @references
#' Bernard Rosner, Robert J. Glynn, Mei-Ling Ting Lee(2006)
#' \emph{The Wilcoxon Signed Rank Test for Paired Comparisons of
#'  Clustered Data.} Biometrics, \bold{62}, 185-192.
#'  

cluswilcox.test.signedrank <-
  function(x, cluster, method,
           alternative = c("two.sided", "less", "greater"),
           mu,
           DNAME = NULL, METHOD = NULL, ...) {
    names(mu) <- "location shift"
    testfunc <- switch(method, rgl = cluswilcox.test.signedrank.rgl,
                       ds = cluswilcox.test.signedrank.ds)
    testfunc(x, cluster, althernative, mu, DNAME, METHOD)
  }

