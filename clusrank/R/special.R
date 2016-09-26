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
#' @seealso \code{\link{clusWilcox.test.formula}}
#'
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
#' @seealso clusWilcox.test.formula
#'
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
    list(vars = (fname[[1]])[spc],
         terms = seq(ff)[ff & match(attr(tt, "order"), order, nomatch = 0)])
}


extractTerm <- function(term, mf, nobs, paired) {
    if (term == "group") {
        term.lab <- attr(terms(mf), "term.labels")
        term.mf <- term.lab[!grepl("[\\(\\)]", term.lab)]
    } else {
        term.mf <- attr(attr(mf, "terms"), "specials")[[term]]
    }
    
    if (length(term.mf) == 0) {
        if (term == "cluster") var <- c(1 : nobs)
        if (term == "stratum") var <- rep(1, nobs)
        if (term == "group") {
            if (!paired) {
                stop("group variable is missing")
            } else {
                var <- NULL
            }
        }
        name <- NULL
    } else {
        if (term == "group") {
            temp <- name <- term.mf
        } else {
            temp <- untangle.specials(terms(mf), term, 1)
            name <- gsub("[\\(\\)]", "",
                         regmatches(temp$vars,
                                    gregexpr("\\(.*?\\)", temp$vars))[[1]])
            temp <- temp$vars
        }
        name <- paste0(" ", term, ": ", name, ";")

        if (length(temp) == 1) {
            keep <- mf[[temp]]
            if (is.null(keep)) {
                stop(paste(term, "is missing from the data"))
            }
        } else {
            stop(paste("more than one variable are set as the",
                       term, "id"))
        }
        uniq <- unique(keep)
        uniq.l <- length(uniq)

        if ((term == "group") & (uniq.l == 1)) {
            stop("group must contain at least two levels")
        }

        var <- keep
        
        if (is.numeric(uniq) | is.character(uniq)) {
            var <- keep
            if (is.character(uniq)) var <- recoderFunc(keep, uniq, c(1 : uniq.l)) 
        } else {
            stop(paste(term, "id should be numeric or character"))
        }
    }

    return(list(name = name, var = var))
}


extractVar <- function(var, pars, data) {
    if (!is.null(pars[[var]])) {
        return(data[, as.character(pars[[var]])])
    } 
}

extractName <- function(var, pars) {
    if (!is.null(pars[[var]])) {
        if (var == "y") {
            return(paste0(" and ", pars[[var]], ";"))
        } else {
            if (var == "x") return(paste0(pars[[var]], ";"))
            return(paste0(" ", var, ": ", pars[[var]], ";"))
        } 
    } else {
        return(NULL)
    }
}

