################################################################################
##
## clusrank: Wilcoxon Rank Tests for Clustered Data
## Copyright (C) 2015-2021  Yujing Jiang, Mei-Ling Ting Lee, and Jun Yan
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
################################################################################


#' @keywords internal
#' @export
print.ctest <- function (x, digits = getOption("digits"), prefix = "\t", ...)
{
    cat("\n")
    cat(strwrap(x$method, prefix = prefix), sep = "\n")
    cat("\n")
    cat("data:  ", x$data.name, "\n", sep = "")
    out <- character()
    cat("number of observations: ", x$nobs, ";  ", sep = "")
    cat("number of clusters: ", x$nclus, "\n", sep = "")
    if (!is.null(x$ngroup))
        cat("number of groups: ", x$ngroup, "\n", sep = "")
    if (!is.null(x$statistic))
        out <- c(out, paste(names(x$statistic), "=", format(signif(x$statistic,
            max(1L, digits - 2L)))))
    if (!is.null(x$parameter))
        out <- c(out, paste(names(x$parameter), "=", format(signif(x$parameter,
            max(1L, digits - 2L)))))
    if (!is.null(x$p.value)) {
        fp <- format.pval(x$p.value, digits = max(1L, digits -
            3L))
        out <- c(out, paste("p-value", if (substr(fp, 1L, 1L) ==
            "<") fp else paste("=", fp)))
    }
    cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
    if (!is.null(x$alternative)) {
        cat("alternative hypothesis: ")
        if (!is.null(x$null.value)) {
            if (length(x$null.value) == 1L) {
                alt.char <- switch(x$alternative, two.sided = "not equal to",
                  less = "less than", greater = "greater than")
                cat("true ", names(x$null.value), " is ", alt.char,
                  " ", x$null.value, "\n", sep = "")
            }
            else {
                cat(x$alternative, "\nnull values:\n", sep = "")
                print(x$null.value, digits = digits, ...)
            }
        }
        else cat(x$alternative, "\n", sep = "")
    }
    if (!is.null(x$conf.int)) {
        cat(format(100 * attr(x$conf.int, "conf.level")), " percent confidence interval:\n",
            " ", paste(format(c(x$conf.int[1L], x$conf.int[2L])),
                collapse = " "), "\n", sep = "")
    }
    if (!is.null(x$estimate)) {
        cat("sample estimates:\n")
        print(x$estimate, digits = digits, ...)
    }
    cat("\n")
    invisible(x)
}
