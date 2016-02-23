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
cluswilcox.test.signedrank.ds <- function(x, cluster, alternative,
                                          mu, DNAME, METHOD) {
  METHOD <- paste(METHOD, "using Datta-Satten method")
  Xij <- x
  cluster.size <- as.vector(table(cluster))
  m <- length(cluster.size)
  n <- sum(cluster.size)
  cum.csize <- cumsum(cluster.size)
  cum.csize <- c(0,cum.csize)
  ## H_i-hat function in the paper.
  Fi <- function(x,i) { 
    Xi <- Xij[(cum.csize[i] + 1):(cum.csize[i + 1])]
    (sum(abs(Xi) <= x) + sum(abs(Xi)<x)) / (2 * cluster.size[i])
  }
  ## Function to compute D_i-hat function in the paper.
  Ftot <- function(x) {
    st <- 0
    for (i in 1:m) st <- st + Fi(x,i)
    return(st)
  }
  ## H-hat function in the paper.
  Fcom <- function(x) {
    st <- 0
    for (i in 1:m) st <- st + Fi(x,i) * cluster.size[i]
    return(st/n)
  }
  
  # SIGNED RANK TEST STATISTIC
  TS <- VTS <- 0
  for (i in 1:m) {
    Xi <- Xij[(cum.csize[i] + 1):(cum.csize[i + 1])]
    first <- (sum(Xi > 0) - sum(Xi < 0)) / length(Xi)
    second <- 0
    third <- 0
    for (x in Xi) { second <- second + sign(x) * (Ftot(abs(x)) - Fi(abs(x),i));
    third <- third + sign(x) * Fcom(abs(x))}
    
    TS <- TS + first + second / length(Xi)
    VTS <- VTS + (first + (m - 1) * third / length(Xi)) ^ 2
  }
  
  Z <- TS / sqrt(VTS)
  P_val <- switch(alternative,
                  less = pnorm(abs(Z)),
                  greater = pnorm(abs(Z), lower.tail = FALSE),
                  two.sided = 2 * min(pnorm(abs(Z)),
                                      pnorm(abs(Z), lower.tail = FALSE)))
  
  names(n) <- "total number of observations"
  names(m) <- "total number of clusters"
  names(Z) <- "Test Statistic"
  result <- list(statistic = Z,
                 p.value = P_val, n = n, cn = m,
                 alternative = alternative,
                 null.value = mu,
                 data.name = DNAME, method = METHOD)
  class(result) <- "ctest"
  return(result)
}