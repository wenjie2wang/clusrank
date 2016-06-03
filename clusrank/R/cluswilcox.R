##' Wilcoxon Rank Sum and Signed Rank Test for Clustered Data
##'
##' Performs one-sample and two-sample Wilcoxon test for clutered data
##' on vectors of data.
##'
##' @param x A numeric vector of data values. Non-finite (e.g.,
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
##'     \code{lhs} is a numeric variable giving the data values and
##'     the \code{rhs} of the form with special term \code{cluster(x1)
##'     + group(x2) + stratum(x3)}, where \code{x1, x2, x3} are the
##'     corresponding variables.
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
##' @param DNAME An optional character string for data name when
##'     printing the result
##' @param ... Further arguments to be passed to or from methods
##' @details The formula interface is only applicable for the
##'     {m}-sample rank sum tests \eqn{m \ge 2}. If the data are saved
##'     in a data frame where the observation, cluster id, group id
##'     and stratum id are saved as \code{z}, \code{id}, \code{grp}
##'     and \code{strat} respectively, then the formula should be
##'     written as \code{z ~ cluster(id) + group(grp) +
##'     stratum(strat)}. The \code{group} variable is required.
##'
##' If both \code{x} and \code{y} are given or only \code{x} is given
##' and \code{paired} is \code{TRUE}, a clustered Wilcoxon signed rank
##' test of the null that the distribution of \code{x - y} or of
##' \code{x} is symmetric about \code{mu} is performed.
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
##' for ranksum test and signed rank test when the \code{method} is
##' \code{"rgl"}. Currently the exact test is not recommended.
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
#' cluswilcox.test(z, cluster = id, data = crsd, paired = TRUE)
#' \dontrun{cluswilcox.test(z, cluster = id, data = crsd)
#' ## Default is rank sum test. The group variable is required.}
#' ## Clustered rank sum test using RGL method.
#' data(crd)
#' cluswilcox.test(z ~ cluster(id) + group(group), data = crd)
#' ## or
#' cluswilcox.test(z, cluster = id, group = group, data = crd)
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
#'@note Exact tests are not recommended in the current version of package.
#' @importFrom  stats complete.cases na.omit terms complete.cases model.extract aggregate
#' @importFrom stats lm ecdf pnorm qnorm var  pchisq setNames lag
#' @importFrom MASS ginv
#' @importFrom Rcpp evalCpp
#' @useDynLib clusrank 
#' @export

cluswilcox.test <- function(x, ...) {
    pars <- as.list(match.call()[-1])
    if(!is.null(pars$data)) {
        data.temp <- eval(pars$data, parent.frame())
    }
    if(!is.null(data.temp) & length(pars$x) == 1) {
        if(is.data.frame(data.temp) & any(as.character(pars$x)
            %in% names(data.temp))) {
            x <- data.temp[, as.character(pars$x)]
        } else if(is.matrix(data.temp) & any(as.character(pars$x)
            %in% colnames(data.temp))) {
            x <- data.temp[, as.character(pars$x)]
        }
    }
    UseMethod("cluswilcox.test", x)
}


#' @method cluswilcox.test formula
#' @describeIn cluswilcox.test \code{S3} method for class 'formula'
#' @export


cluswilcox.test.formula <- function(formula, data = NULL, subset = NULL, na.action = na.omit, ...)
{
    if(missing(formula) ||
       (length(formula) != 3L)) {
        stop("'formula' missing or incorrect")
    }
    m <- match.call(expand.dots = FALSE)
    if(!missing(data)) {
        DNAME <- paste("from", m$data)
    } else {
        DNAME <- NULL
    }

    if(is.matrix(eval(m$data, parent.frame()))) {
        m$data <- as.data.frame(data)
    }
    special <- c("stratum", "cluster", "group")
    m[[1L]] <- quote(stats::model.frame)
    m$... <- NULL
    m$formula <- if(missing(data)) terms(formula, special)
                 else terms(formula, special, data = data)

    mf <- eval(m, parent.frame())
    Terms <- terms(mf)

    x.name <- rownames(attr(m$formula, "factors"))[1]
    DNAME <- paste0(paste(x.name, "from", m$data), ",")
    response <- attr(attr(mf, "terms"), "response")
    x <- mf[[response]]
    n.obs <- length(x)

    group <- attr(Terms, "specials")$group
    if(length(group)) {
        gtemp <- untangle.specials(Terms, "group", 1)
        group.name <- gsub("[\\(\\)]", "",
                           regmatches(gtemp$vars,
                                      gregexpr("\\(.*?\\)", gtemp$vars))[[1]])
        DNAME <- paste0(DNAME, " group: ", group.name, ",")
        

        if(length(gtemp$vars) == 1) {
            group.keep <- mf[[gtemp$vars]]
        } else {
            stop("more than one variable are set as the group id")
        }
        group.uniq <- unique(group.keep)
        group.uniq.l <- length(group.uniq)
        
        if(!is.character(group.uniq) && !is.numeric(group.uniq)) {
            stop("group id has to be numeric or character")
        }
        
        group <- recoderFunc(group.keep, group.uniq, c(1 : group.uniq.l))
    }
    
    cluster <- attr(attr(mf, "terms"), "specials")$cluster
    if(length(cluster)) {
    ctemp <- untangle.specials(Terms, "cluster", 1)
    cluster.name <- gsub("[\\(\\)]", "",
                         regmatches(ctemp$vars,
                                    gregexpr("\\(.*?\\)", ctemp$vars))[[1]])
    DNAME <- paste0(DNAME, " cluster: ", cluster.name, ",")


    if(length(ctemp$vars) == 1) {
        cluster.keep <- mf[[ctemp$vars]]
    } else {
        stop("more than one variable are set as the cluster id")
    }
    cluster.uniq <- unique(cluster.keep)
    cluster.uniq.l <- length(cluster.uniq)
    
    if(is.character(cluster.uniq)) {
        cluster <- recoderFunc(cluster.keep, cluster.uniq, c(1 : cluster.uniq.l))
    } else {
        if(!is.numeric(cluster.uniq)) {
            stop("cluster id should be numeric or character")
        }
        cluster <- cluster.keep
    }
    } else {
        cluster <- c(1 : n.obs)
    }


    
    stratum <- attr(attr(mf, "terms"), "specials")$stratum
    if(!is.null(stratum)) {
        stratum <- mf[[stratum]]
    } else {
        stratum <- rep(1, n.obs)
    }
    
    y <- do.call("cluswilcox.test.default",
                 c( list(x = x, cluster = cluster,
                         group = group, stratum = stratum,
                         DNAME = DNAME),
                   list(...)))
    return(y)
}


#' @method cluswilcox.test default
#' @describeIn cluswilcox.test Default \code{S3} method.
#' @export

cluswilcox.test.default <- function(x, y = NULL, cluster = NULL,
            group = NULL, stratum = NULL, data = parent.frame(),
            alternative = c("two.sided", "less", "greater"),
            mu = 0, paired = FALSE, exact = FALSE,
            method = c("rgl", "ds"), DNAME = NULL, ...) {
    alternative <- match.arg(alternative)
    method <- match.arg(method)
      pars <- as.list(match.call()[-1])

    if (!missing(mu) && ((length(mu) > 1L) || !is.finite(mu)))
        stop("'mu' must be a single number")
    if(is.null(DNAME))  {
 if(!is.null(pars$data)) {
            x <- data[, as.character(pars$x)]
            DNAME <- (pars$x)
              if(!is.null(pars$y)) {
                  y <- data[, as.character(pars$y)]
                  DNAME <- paste(DNAME, "and", pars$y)
                } else {
                    y <- NULL
                  }
            if(!is.null(pars$cluster)) {
                cluster <- data[, as.character(pars$cluster)]
                DNAME <- paste0(DNAME, ", cluster: ", pars$cluster)
              } else {
                  cluster <- NULL
              }
            if(!is.null(pars$group)) {
              group <- data[, as.character(pars$group)]
              DNAME <- paste0(DNAME, ", group: ", pars$group)
            } else {
              group <- NULL
            }
            if(!is.null(pars$stratum)) {
              stratum <- data[, as.character(pars$stratum)]
              DNAME <- paste0(DNAME, ", stratum: ", pars$stratum)
            } else {
              stratum <- NULL
            }
            DNAME <- paste0(DNAME, " from ", pars$data)
          } else {
              DNAME <- (pars$x)
              if(!is.null(y)) {
                  DNAME <- paste(DNAME, "and", (pars$y))
                }

               if(!is.null(cluster)) {
                    DNAME <- paste0(DNAME, ", cluster id: ", pars$cluster)
               }
              if(!is.null(pars$group)) {
                DNAME <- paste0(DNAME, ", group: ", pars$group)
              }
              if(!is.null(pars$stratum)) {
                DNAME <- paste0(DNAME, ", stratum: ", pars$stratum)
              }
            }
    }
    
      if(!is.null(pars$data)) {
            x <- data[, as.character(pars$x)]
              if(!is.null(pars$y)) {
                  y <- data[, as.character(pars$y)]

                } else {
                    y <- NULL
                  }
            if(!is.null(pars$cluster)) {
                cluster <- data[, as.character(pars$cluster)]

              } else {
                  cluster <- NULL
              }
            if(!is.null(pars$group)) {
              group <- data[, as.character(pars$group)]
            } else {
              group <- NULL
            }
            if(!is.null(pars$stratum)) {
              stratum <- data[, as.character(pars$stratum)]

            } else {
              stratum <- NULL
            }

          }
    if(!is.numeric(x)){
        stop("'x' must be numeric")
    }



    if(!is.null(y)) {
        if(!is.numeric(y)) {
            stop("'y' must be numeric")
        }
        if(length(y) != length(x)) {
            stop("'y' must have the same length as 'x' for the clustered signed rank test")
        }
        paired <- TRUE
        x <- x - y
    }

    if(is.null(cluster)) {
        stop("'cluster' is required")
    }

    if(is.null(group) & paired == FALSE) {
        stop("'group' is required for the clustered rank sum test")
    }

    if(!is.null(group) & paired == TRUE) {
        warning("'group' will be ignored for the clustered signed rank test")
    }

    if(is.null(stratum)) {
        stratum <- rep(1, length(x))
    }
    if(is.null(DNAME)) {
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

    if(length(x) < 1L) {
        stop("not enough (finite) 'x' observation")
    }

    if(paired == TRUE) {
        if(length(table(stratum)) > 1L) {
            warning("'stratum' will be ignored for the clustered signed rank test")
        }
        x <- x - mu
        METHOD <- "Clustered Wilcoxon signed rank test"
       if(toupper(method) == "RGL") {
           METHOD <- paste(METHOD, "using RGL method", sep = " ")
           arglist <- setNames(list(x, cluster, alternative, mu, METHOD,
                                    DNAME, exact),
                            c("x", "cluster", "alternative",
                              "mu",
                              "METHOD", "DNAME",  "exact"))
            result <- do.call("cluswilcox.test.signedrank.rgl", c(arglist))
            return(result)
        }

        if(toupper(method) == "DS") {
            METHOD <- paste(METHOD, "using DS method", sep = " ")
             arglist <- setNames(list(x, cluster, alternative, mu, METHOD, DNAME),
                            c("x", "cluster", "alternative",
                              "mu",
                              "METHOD", "DNAME"))
           result <-  do.call("cluswilcox.test.signedrank.ds",
                              c(arglist))
           return(result)
        } else {
            stop("Method should be one of 'rgl' and 'ds'")
        }
        
    } else {
        METHOD <- "Clustered Wilcoxon rank sum test"
        if(toupper(method) == "RGL") {
            METHOD <- paste( METHOD, "using Rosner-Glynn-Lee method", sep = " ")
             arglist <- setNames(list(x, cluster, group, stratum, alternative,
                                 mu, DNAME, METHOD, exact),
                            c("x", "cluster", "group", "stratum",
                              "alternative", "mu", "DNAME", "METHOD",
                              "exact"))
            result <- do.call("cluswilcox.test.ranksum.rgl", c(arglist))
            return(result)
        }

        if(toupper(method) == "DS") {
            METHOD <- paste(METHOD, "using Datta-Satten method", sep = " ")
            ## FIXME: The length of the list and the name do not match!
            arglist <- setNames(list(x, cluster, group, mu, alternative,
                                     DNAME, METHOD),
                                c("x", "cluster", "group", "mu",
                                  "alternative", "DNAME", "METHOD"))
            if(exact == TRUE) {
                warning(" No exact test is provided for 'ds' method")
            }

            if(length(table(stratum)) > 1L) {
                warning("'stratum' will be ignored for the clustered rank sum test, 'ds' method")
            }
            result <- do.call("cluswilcox.test.ranksum.ds", c(arglist))
            return(result)
        }
        
        else {
            stop("Method should be one of 'rgl' and 'ds'")
        }
    }
}





