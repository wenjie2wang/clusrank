##
## clusrank: Wilcoxon Rank Tests for Clustered Data
## Copyright (C) 2016-2021  Yujing Jiang, Mei-Ling Ting Lee, and Jun Yan
##
## This file is part of the R package clusrank.
##
## The R package clusrank is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package clusrank is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##


##' Pairwise Wilcoxon Rank Sum and Signed Rank Tests for Clustered Data
##'
##' Performs pairwise comparisons between group levels with corrections for
##' multiple testing.
##'
##' @inheritParams clusWilcox.test
##'
##' @param x A numeric vector of data values.
##' @param p.adjust.method Method for adjusting p values.  \code{?p.adjust} for
##'     details.
##' @param ... Other arguments passed to \code{clusWilcox.test}.
##'
##' @return A \code{pairwise.htest} object.
##'
##' @examples
##' library(clusrank)
##' data(crd)
##'
##' ## for demonstration purpose, create random groups
##' set.seed(123)
##' g <- sample(seq_len(3), size = nrow(crd), replace = TRUE)
##' with(crd, pairwise.clusWilcox.test(z, group = g, cluster = id,
##'                                    method = "ds", p.adjust.method = "fdr"))
##' @importFrom stats p.adjust.methods p.adjust pairwise.table
##' @export
pairwise.clusWilcox.test <- function(x,
                                     group,
                                     cluster,
                                     p.adjust.method = p.adjust.methods,
                                     ...)
{
    ## base::deparse1 since R 4.0.0
    deparse_ <- function(expr, collapse = " ", width.cutoff = 500L, ...) {
        paste(deparse(expr, width.cutoff, ...), collapse = collapse)
    }
    p.adjust.method <- match.arg(p.adjust.method)
    DNAME <- paste(deparse_(substitute(x)), "and", deparse_(substitute(group)))
    group <- factor(group)
    METHOD <- NULL
    compare.levels <- function(i, j) {
        is_i <- as.integer(group) == i
        is_j <- as.integer(group) == j
        xi <- x[is_i]
        xj <- x[is_j]
        cluster_i <- cluster[is_i]
        cluster_j <- cluster[is_j]
        if (is.null(METHOD)) {
            wt <- clusWilcox.test.default(
                c(xi, xj),
                group = c(rep(i, length(xi)),
                          rep(j, length(xj))),
                cluster = c(cluster_i, cluster_j),
                ...
            )
            METHOD <<- wt$method
            wt$p.value
        }
        else {
            clusWilcox.test.default(
                c(xi, xj),
                group = c(rep(i, length(xi)),
                          rep(j, length(xj))),
                cluster = c(cluster_i, cluster_j),
                ...
            )$p.value
        }
    }
    PVAL <- pairwise.table(compare.levels, levels(group), p.adjust.method)
    ans <- list(method = METHOD,
                data.name = DNAME,
                p.value = PVAL,
                p.adjust.method = p.adjust.method)
    class(ans) <- "pairwise.htest"
    ans
}
