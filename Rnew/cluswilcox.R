cluswilcox <- function(x, ...) {
    UseMethod("wilcox.test")
}

cluswilcox.formula <- function(formula, data, subset, na.action, ...)
{
    if(missing(formula) ||
      (length(formula) != 3L)) {
        stop("'formula' missing or incorrect")
    }
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame()))) {
        m$data <- as.data.frame(data)
    }
    special <- c("stratum", "cluster", "group")
    m[[1L]] <- quote(stats::model.frame)
    m$... <- NULL
    m$formula <- terms(formula, special, data = data)
    mf <- eval(m, parent.frame())

    x.name <- rownames(attr(temp$formula, "factors"))[1]
    DNAME <- paste0(paste(x.name, "from", m$data), ",")
    response <- attr(attr(mf, "terms"), "response")
    group <- attr(attr(mf, "terms"), "specials")$group
    g <- factor(mf[[group]])
    cluster <- attr(attr(mf, "terms"), "specials")$cluster
    if(is.null(cluster)) {
        stop("cluster id is missing")
    }
    clus <- mf[[cluster]]

    stratum <- attr(attr(mf, "terms"), "specials")$stratum
    if(!is.null(stratum)) {
        strat <- mf[[stratum]]
    }
    
    if(nlevels(g) == 2L) {
        DATA <- setNames(split(mf[[response]], g), c("x", "y"))
        CLUS <- setNames(split(clus, g), c("xclus", "yclus"))
        if(!is.null(stratum)) {
            STRAT <- setNames(split(clus, g), c("xstrat", "ystrat"))
        } else {
            STRAT <- NULL
        }
        y <- do.call("cluswilcox.test", c(DATA, CLUS, STRAT, list(...)))
        y$data.name <- DNAME
        return(y)
    } else if (nlevels(g) > 2L) {
        method <- match.arg(method)
        if(toupper(method) != "DS") {
            stop("RGL method cannot handle comparison among m groups, try DS method")
        } else {
            DATA <- list(x = mf[[response]])
            CLUS <- list(xclus = clus)
            if(!is.null(stratum)) {
                stop("DS method cannot handel stratified data")
            }
            y <- do.call("cluswilcox.test.ranksum.ds.mgroup", c(DATA, CLUS, list(...)))
        }
    } else {
        stop("grouping factor must have 2 or more levels")
    }
}


cluswilcox.test.default <- function(x, y = NULL, xclus = NULL, yclus = NULL,
                                    xstrat = NULL, ystrat = NULL,
                                    alternative = c("two.sided", "less", "greater"),
                                    mu = 0, paired = FALSE, exact = NULL,
                                    method = c("rgl", "ds"), ...) {
    if(!is.numeric(x)){
        stop("'x' must be numeric")
    }
    x.len <- length(x)
    if(!is.null(xclus)) {
        if(length(xclus) != x.len) {
            stop("'xclus' and 'x' must have the same length")
        }
    } else {
        xclus <- c(1:x.len)
    }

    ## Check stratum
    if(is.null(y) | paired == TRUE | toupper(method) == DS) {
        if(!is.null(xtrat) | !is.null(ystrat)) {
            stop(" Only accept stratified data for RGL rank sum test")
        }
    }
    
    if(!is.null(y)) {
        if(!is.numeric(y)){
            stop("'y' must be numeric")
        }
            if(!is.null(yclus)) {
                if(length(yclus) != length(y)) {
                    stop("'yclus' and 'y' must have the same length")
                } else {
                    yclus <- xclus
                }
            }
            if(paired){
                if(x.len != length(y))
                    stop("'x' and 'y' must have the same length")
                OK <- complete.cases(x, y)
                x <- x[OK] - y[OK]
                xclus <- xclus[OK]
                y <- NULL
            } else {
                x.OK <- is.finite(x)
                x <- x[x.OK]
                xclus <- xclus[x.OK]
                y.OK <- is.finite(y)
                y <- y[y.OK]
                yclus <- yclus[y.OK]
                if(!is.null(xstrat))
                    xstrat <- xstrat[x.OK]
                if(!is.null(ystrat))
                    ystrat <- ystrat[y.OK]
            }
    } else {
        DNAME <- deparse(substitute(x))
        if(paired)
            stop("'y' is missing for paired test")
        x.OK <- is.finite(x)
        x <- x[x.OK]
        xclus <- xclus[x.OK]       
    }

    if(length(x) < 1L)
        stop("not enough (finite) 'x' observation")
    if(is.null(y)) {
        METHOD <- "Cluster Wilcoxon signed rank test"
        arglist <- setNames(list(x, xclus, alternative, mu, DNAME, METHOD, exact),
                            c("x", "cluster", "alternative",
                              "mu", "DNAME", 
                              "METHOD", "exact"))
        if(toupper(method) == "RGL") {
            result <- do.call("cluswilcox.test.signedrank.rgl", c(arglist))
            return(result)
        }
        if(toupper(method) == "DS") {
           result <-  do.call("cluswilcox.test.signedrank.ds",
                              c(arglist))
           return(result)
        }
        else {
            stop("Method should be one of 'rgl' and 'ds'")
        }
    } else {
        METHOD <- "Cluster Wilcoxon rank sum test"
        x <- c(x, y)
        cluster <- c(xclus, yclus)
        group <- c(rep(1, length(x)), rep(2, length(y)))
        strats <- c(xstrat, ystrat)
        arglist <- setNames(list(x, cluster, group, strats, alternative,
                                 mu, DNAME, METHOD, exact),
                            c("x", "cluster", "group", "strats",
                              "alternative", "mu", "DNAME", "METHOD",
                              "exact"))
        if(toupper(method) == "RGL") {
            result <- do.call("cluswilcox.test.ranksum.rgl", c(arglist))
            return(result)
        }
        if(toupper(method) == "DS") {
            result <- do.call("cluswilcox.test.ranksum.ds", c(arglist))
            return(result)
        }
        else {
            stop("Method should be one of 'rgl' and 'ds'")
        }
    }
}
        
        
    
