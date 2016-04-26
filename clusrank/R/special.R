#' Identify treatment groups.
#'
#' This is a special function used in the context of formula
#' used for Wilcoxon sum rank test for clustered data.
#' It identifies the treatment group id of observations, and is used
#' on the right hand side of a formula.
#'
#' @param x A numeric variable of cluster id.
#'
#' @details The only action of this function is semantic, to mark
#' a variable as the group indicator. Must be supplied.
#' @seealso \code{\link{cluswilcox.test.formula}}
#'
#' @examples
#' data(crd)
#' cluswilcox.test(z ~ cluster(id) + group(group), data = crd)
#' @keywords internal
#' @export
group <- function(x) {x}


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
#' @seealso \code{\link{cluswilcox.test.formula}}
#'
#' @examples
#' data(crd)
#' cluswilcox.test(z ~ cluster(id) + group(group), data = crd)
#' @keywords internal
#' @export
cluster <- function(x) {x}

#' Identify strata.
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
#' @keywords internal
#' @export
stratum <- function(x) {x}


untangle.specials <- function (tt, special, order = 1) 
{
  spc <- attr(tt, "specials")[[special]]
  if (length(spc) == 0) 
    return(list(vars = character(0), terms = numeric(0)))
  facs <- attr(tt, "factors")
  fname <- dimnames(facs)
  ff <- apply(facs[spc, , drop = FALSE], 2, sum)
  list(vars = (fname[[1]])[spc], terms = seq(ff)[ff & match(attr(tt, 
                                                                 "order"), order, nomatch = 0)])
}
