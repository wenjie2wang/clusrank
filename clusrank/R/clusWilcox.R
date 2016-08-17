##' Wilcoxon Rank Sum and Signed Rank Test for Clustered Data
##'
##' Performs one-sample and two-sample Wilcoxon test for clutered data
##' on vectors of data.
##'
##' @param x A numeric vector of data values or a formula. Non-finite (e.g.,
##'     infinite or missing) values will be omitted
##' @param y An optional numeric vector of data values: as with
##'     \code{x} non-finite values will be omitted
##' @param cluster A optional numeric vector of cluster id. Does not
##'     apply to the formula interface
##' @param group A optional numeric vector of treatment id. Does not
##'     apply to the formula interface
##' @param stratum A opptional numeric vector of stratum id. Only
##'     applies for the \code{"rgl"} method. Does not apply to the
##'     formula interface.
##' @param method A character string specifying the method of
##'     clustered wilcoxon rank test used, must be one of \code{"rgl"}
##'     or \code{"ds"}
##' @param paired A logical indicating whether you want a paired test
##' @param exact A logical indicating whether an exact p-value should
##'     be computed. Not recommended currently.
##' @param formula A formula of the form \code{lhs ~ rhs} where the
##'     \code{lhs} is the observed values and
##'     the \code{rhs} of the form group + \code{cluster}(x1) +
##'     \code{stratum}(x2) for clustered rank sum test, where
##'     \code{x1} and \code{x2} are cluster id and stratum id in the
##'     data.  For clustered signed rank test, only \code{cluster}(x1)
##'     should appear on the right hand side.
##' @param data An optional matrix or dataframe of data used in the
##'     formula.
##' @param subset An optional vector specifying a subset of
##'     observations to be used.
##' @param na.action A function which indicates what should happen
##'     when the data contain NAs. Defaults to getOption("na.action").
##' @param alternative A character string specifying the alternative
##'     hypothesis, must be one of \code{"two sided"} (default),
##'     \code{"greater"} or \code{"less"}. You can specify just the
##'     initial letter
##' @param mu A number specifying an optional parameter used to form
##'     the null hypothesis. See 'Details'
##' @param ... Further arguments to be passed to or from methods
##' @details The formula interface is to both clustered signed rank
##'     test and clustered rank sum test.
##'
##' Given the cluster id, if both \code{x} and \code{y} are provided
##' or only \code{x} is provided and \code{paired} is \code{TRUE}, a
##' clustered Wilcoxon signed rank test of the null hypothesis that
##' the distribution of \code{x - y} or of \code{x} is symmetric about
##' \code{mu} is performed.
##'
##' Otherwise, if only \code{x} is given and \code{paired} is
##' \code{FALSE}, a Wilcoxon rank sum test is carried out. In this
##' case, the \code{group} variable is required. If the \code{method}
##' is \code{"rgl"} (default), the null hypothesis is that the
##' distributions of values from the two groups differ by a location
##' shift of \code{mu} and the alternative is that they differ by some
##' other location shift. If the \code{method} is \code{"ds"}, when
##' the \code{group} has 2 levels, the null and hypothesis are the
##' same as for \code{"rgl"} test; when \code{group} has more than 2
##' levels, the null hypothesis is that the locations
##' are the same for data in all groups and the alternative is that
##' they are not all the same.
##'
##' 
##' If \code{cluster} is not provided, the default is that there is no
##' clutering in the data. Both \code{"rgl"} and \code{"ds"} method
##' support balanced and unbalanced data (cluster size is identical or
##' varied) and individual level and cluster level treatment assignment.
##'
##' If \code{method} is \code{"rgl"}, then a strafication variable,
##' \code{stratum}, is allowed for the clustered Wilcoxon rank sum
##' test.
##'
##' The exact test is still under development and is only available
##' for ranksum test when treatment is assigned at cluster level and
##' signed rank test for \code{rgl} method and can be applied when
##' number of cluster is small (under 50). It will be slow if number
##' of cluster is too large.
##'
##' There is also a formula interface for both tests. For details look
##' at the examples.
##' 
##' @return A list with class \code{"ctest"}.
##'
##' \item{Rstat}{The value of the rank statistic with a name discribing it}
##' \item{ERstat}{The expectation of the rank statistic}
##' \item{VRstat}{The variance of the rank statistic}
##' \item{statistic}{The value of the test statistic}
##' \item{p.value}{The p-value for the test}
##' \item{null.value}{The location parameter \code{mu}}
##' \item{method}{The type of test applied}
##' \item{data.name}{A character string giving the names of the data}
##' \item{balance}{A logical indicating whether the data are balanced}
##' \item{n.group}{Number of treatment, is returned when there are more than 2 treatment groups}
##' \item{df}{Degrees of freedom of chi-square distribution, is returned when there are more than 2 treatment groups}

#' @examples
#' ## Clustered signed rank test using RGL method.
#' data(crsd)
#' clusWilcox.test(z, cluster = id, data = crsd, paired = TRUE)
#' ## or
#' clusWilcox.test(z ~ cluster(id), data = crsd, paired = TRUE)
#' \dontrun{clusWilcox.test(z, cluster = id, data = crsd)
#' ## Default is rank sum test. The group variable is required.}
#' ## Clustered rank sum test using RGL method.
#' data(crd)
#' clusWilcox.test(z ~ group + cluster(id), data = crd)
#' ## or
#' clusWilcox.test(z, cluster = id, group = group, data = crd)
#' @author Yujing Jiang
#' @references
#' Bernard Rosner, Robert J. Glynn, Mei-Ling T. Lee (2006)
#' \emph{The Wilcoxon Signed Rank Test for Paired Comparisons of
#'  Clustered Data}. Biometrics, \bold{62}, 185-192.
#'
#' Bernard Rosner, Robert J. Glynn, Mei-Ling T. Lee (2003)
#' \emph{Incorporation of Clustering Effects for the Wilcoxon Rank
#' Sum Test: A Large-Sample Approach}. Biometrics, \bold{59}, 1089-1098.
#' 
#' Bernard Rosner, Robert J. Glynn, Mei-Ling T. Lee (2006)
#' \emph{Extension of the Rank Sum Test for Clustered Data: 
#' Two-Group Comparisons with Group}. Biometrics, \bold{62}, 1251-1259.
#' 
#' Somnath Datta, Glen A. Satten (2005) \emph{Rank-Sum Tests for Clustered Data}.
#' Journal of the American Statistical Association, \bold{100}, 908-915.
#' 
#' Somath Datta, Glen A. Satten (2008) \emph{A Signed-Rank test for Clustered Data}.
#' Biometric, \bold{64}, 501-507.
#'
#' @note Exact tests are not recommended in the current version of package.
#' @importFrom  stats complete.cases na.omit terms complete.cases model.extract aggregate
#' @importFrom stats lm ecdf pnorm qnorm var  pchisq setNames lag
#' @importFrom MASS ginv
#' @importFrom Rcpp evalCpp
#' @useDynLib clusrank 
#' @export

clusWilcox.test <- function(x, ...) {
    pars <- as.list(match.call()[-1])
    if (!is.null(pars$data)) {
        data.temp <- eval(pars$data, parent.frame())
    } else {
        data.temp <- NULL
    }
    if (!is.null(data.temp) && length(pars$x) == 1) {
        if (is.data.frame(data.temp) &&
            any(as.character(pars$x) %in% names(data.temp))) {
            x <- data.temp[, as.character(pars$x)]
        } else if(is.matrix(data.temp) &&
                  any(as.character(pars$x) %in% colnames(data.temp))) {
            x <- data.temp[, as.character(pars$x)]
        }
    }
    UseMethod("clusWilcox.test", x)
}


#' @method clusWilcox.test formula
#' @describeIn clusWilcox.test \code{S3} method for class 'formula'
#' @export


clusWilcox.test.formula <- function(formula, data = parent.frame(), subset = NULL,
                                    na.action = na.omit,
                                    alternative = c("two.sided", "less", "greater"),
                                    mu = 0, paired = FALSE, exact = FALSE,
                                    method = c("rgl", "ds"),  ...) {
    if (missing(formula) || (length(formula) != 3L)) {
        stop("'formula' missing or incorrect")
    }
    m <- match.call(expand.dots = FALSE)
    

    if (is.matrix(eval(m$data, parent.frame()))) {
        m$data <- as.data.frame(data)
    }
    special <- c("stratum", "cluster")
    m[[1L]] <- quote(stats::model.frame)
    m$... <- NULL
    m$formula <- if(missing(data)) terms(formula, special)
                 else terms(formula, special, data = data)

    m <- if (missing(data)) m <- m[1 : 2]
         else m <- m[1 : 3]
    mf <- eval(m, parent.frame())
    Terms <- terms(mf) ## Delete 

    x.name <- rownames(attr(m$formula, "factors"))[1]
    DNAME <- paste0(x.name, ";")
    response <- attr(terms(mf), "response")
    x <- mf[[response]]
    n.obs <- length(x)
    
    group <- extractTerm("group", mf, n.obs, paired)
    DNAME <- paste0(DNAME, group[["name"]])
    group <- group[["var"]]

    cluster <- extractTerm("cluster", mf, n.obs, paired)
    DNAME <- paste0(DNAME, cluster[["name"]])
    cluster <- cluster[["var"]]

    stratum <- extractTerm("stratum", mf, n.obs, paired)
    DNAME <- paste0(DNAME, stratum[["name"]])
    stratum <- stratum[["var"]]
    
    
    ## please use tab key to keep the indentation correct
    if (!missing(data)) {
        DNAME <- paste(DNAME, "(from", paste0(m$data, ")"))
    }
    
    y <- do.call("clusWilcox.test.default",
                 c( list(x = x, cluster = cluster,
                         group = group, stratum = stratum,
                         DNAME = DNAME, paired = paired,
                         alternative = alternative,
                         method = method, mu = mu, exact = exact),
                   list(...)))
    return(y)
}


#' @method clusWilcox.test default
#' @describeIn clusWilcox.test Default \code{S3} method.
#' @export

clusWilcox.test.default <- function(x, y = NULL, cluster = NULL,
                                    group = NULL, stratum = NULL, data = NULL,
                                    alternative = c("two.sided", "less", "greater"),
                                    mu = 0, paired = FALSE, exact = FALSE,
                                    method = c("rgl", "ds"), ...) {
    alternative <- match.arg(alternative)
    method <- match.arg(method)
    pars <- as.list(match.call()[-1])
    
    if (!missing(mu) && ((length(mu) > 1L) || !is.finite(mu)))
        stop("'mu' must be a single number")
    DNAME <- list(...)$"DNAME"
    if (is.null(DNAME))  {
        if (!is.null(pars$data)) {
            x <- extractVar("x", pars, data)
            DNAME <- extractName("x", pars)

            y <- extractVar("y", pars, data)
            DNAME <- paste0(DNAME, extractName("y", pars))

            group <- extractVar("group", pars, data)
            DNAME <- paste0(DNAME, extractName("group", pars))

            cluster <- extractVar("cluster", pars, data)
            DNAME <- paste0(DNAME, extractName("cluster", pars))

            stratum <- extractVar("stratum", pars, data)
            DNAME <- paste0(DNAME, extractName("stratum", pars))

            DNAME <- paste0(DNAME, " (from ", paste0(pars$data, ")"))
        } else {
            DNAME <- (pars$x)
            DNAME <- paste0(DNAME, extractName("y", pars))
            DNAME <- paste0(DNAME, extractName("group", pars))
            DNAME <- paste0(DNAME, extractName("cluster", pars))
            DNAME <- paste0(DNAME, extractName("stratum", pars))    
        }
    } 


    if (!is.null(y)) {
        if(!is.numeric(y)) {
            stop("'y' must be numeric")
        }
        if(length(y) != length(x)) {
            stop("'y' must have the same length as 'x' for the clustered signed rank test")
        }
        paired <- TRUE
        x <- x - y
    }

    if (is.null(cluster)) {
        stop("'cluster' is required")
    }

    if (is.null(group) && paired == FALSE) {
        stop("'group' is required for the clustered rank sum test")
    }
    
    if (!is.null(group) && paired == TRUE) {
        warning("'group' will be ignored for the clustered signed rank test")
    }

    if (is.null(stratum)) {
        stratum <- rep(1, length(x))
    }
    if (is.null(DNAME)) {
        if(is.null(y)) {
            DNAME  <-  deparse(substitute(x))
        } else {
            DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
        }
    }
    OK <- complete.cases(x, cluster, group, stratum) & is.finite(x)
    x <- x[OK]
    cluster <- cluster[OK]
    group <- group[OK]
    stratum <- stratum[OK]
    
    if (length(x) < 1L) {
        stop("not enough (finite) 'x' observation")
    }
    
    if (paired == TRUE) {
        if (length(table(stratum)) > 1L) {
            warning("'stratum' will be ignored for the clustered signed rank test")
        }
        x <- x - mu
        METHOD <- "Clustered Wilcoxon signed rank test"
        if (method == "rgl") {
            METHOD <- paste(METHOD, "using Rosner-Glynn-Lee method", sep = " ")
            arglist <- setNames(list(x, cluster, alternative, mu, METHOD,
                                     DNAME, exact),
                                c("x", "cluster", "alternative",
                                  "mu",
                                  "METHOD", "DNAME",  "exact"))
            result <- do.call("clusWilcox.test.signedrank.rgl", c(arglist))
            return(result)
        }
        
        if (method == "ds") {
            METHOD <- paste(METHOD, "using Datta-Satten method", sep = " ")
            arglist <- setNames(list(x, cluster, alternative, mu, METHOD, DNAME),
                                c("x", "cluster", "alternative",
                                  "mu", "METHOD", "DNAME"))
            result <-  do.call("clusWilcox.test.signedrank.ds",
                               c(arglist))
            return(result)
        } else {
            stop("Method should be one of 'rgl' and 'ds'")
        }
        
    } else {
        METHOD <- "Clustered Wilcoxon rank sum test"
        if ((method) == "rgl") {
            METHOD <- paste( METHOD, "using Rosner-Glynn-Lee method", sep = " ")
            arglist <- setNames(list(x, cluster, group, stratum,
                                     alternative, mu, DNAME,
                                     METHOD, exact),
                                c("x", "cluster", "group", "stratum",
                                  "alternative", "mu", "DNAME", "METHOD",
                                  "exact"))
            result <- do.call("clusWilcox.test.ranksum.rgl", c(arglist))
            return(result)
        }
        
        if ((method) == "ds") {
            METHOD <- paste(METHOD, "using Datta-Satten method", sep = " ")
            ## FIXME: The length of the list and the name do not match!
            arglist <- setNames(list(x, cluster, group, mu, alternative,
                                     DNAME, METHOD),
                                c("x", "cluster", "group", "mu",
                                  "alternative", "DNAME", "METHOD"))
            if (exact == TRUE) {
                warning(" No exact test is provided for 'ds' method")
            }
            
            if (length(table(stratum)) > 1L) {
                warning("'stratum' will be ignored for the clustered rank sum test, 'ds' method")
            }
            result <- do.call("clusWilcox.test.ranksum.ds", c(arglist))
            return(result)
        }
        
        else {
            stop("Method should be one of 'rgl' and 'ds'")
        }
    }
}
