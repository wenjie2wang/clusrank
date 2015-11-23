#' Identify clusters
#'
#' This is a special function used in the context of formula
#' used for Wilcoxon sum rank test for clustered data.
#' It identifies the cluster id of observations, and is used
#' on the right hand side of a formula.
#'
#' @param x A numeric variable of cluster id.
#'
#' @details THe function's only action is semantic, to mark
#' a variable as the cluster indicator. If not supplied,
#' will assume no cluster in the data.
#' @return x
#' @seealso cluswilcox.test.formula
#'
#' @examples
#' data(crd)
#' cluswilcox.test(z ~ cluster(id) + group(group), data = crd)
#' @export
cluster <- function(x) {x}

#' Identify stratums.
#'
#' This is a special function used in the context of formula
#' used for Wilcoxon sum rank test for clustered data.
#' It identifies the stratum id of observations, and is used
#' on the right hand side of a formula.
#'
#' @param x A numeric variable of stratum id.
#'
#' @details THe function's only action is semantic, to mark
#' a variable as the stratum indicator. If not supplied,
#' will assume no stratification in the data.
#' @seealso cluswilcox.test.formula
#'
#' @examples
#' data(crdStr)
#' cluswilcox.test(z ~ cluster(id) + group(group) + stratum(stratum), data = crdStr)
#' @export
stratum <- function(x) {x}

#' Identify treatment groups.
#'
#' This is a special function used in the context of formula
#' used for Wilcoxon sum rank test for clustered data.
#' It identifies the treatment group id of observations, and is used
#' on the right hand side of a formula.
#'
#' @param x A numeric variable of cluster id.
#'
#' @details THe function's only action is semantic, to mark
#' a variable as the group indicator. Must be supplied.
#' @seealso cluswilcox.test.formula
#'
#' @examples
#' data(crd)
#' cluswilcox.test(z ~ cluster(id) + group(group), data = crd)
#' @export
group <- function(x) {x}
