cluswilcox <- function(x, ...) {
    UseMethod("wilcox.test")
}

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
    m$formula <- if(missing(data))
    terms(formula, special)
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


    


cluswilcox.test.default <- function(x, y = NULL, cluster = NULL,
                                    group = NULL, stratum = NULL,
                                    alternative = c("two.sided", "less", "greater"),
                                    mu = 0, paired = FALSE, exact = NULL,
                                    method = c("rgl", "ds"), DNAME = NULL) {
    alternative <- match.arg(alternative)
    method <- match.arg(method)
      pars <- as.list(match.call()[-1])

    if (!missing(mu) && ((length(mu) > 1L) || !is.finite(mu)))
    stop("'mu' must be a single number")
    if(!is.numeric(x)){
        stop("'x' must be numeric")
    }
    if(is.null(DNAME)) {
         DNAME <- (pars$x)

    if(!is.null(pars$y)) {
      DNAME <- paste(DNAME, "and", pars$y)
    }
    if(!is.null(pars$cluster)) {
      DNAME <- paste0(DNAME, ", cluster: ", pars$cluster)
    }
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
             arglist <- setNames(list(x, cluster, alternative, METHOD, DNAME),
                            c("x", "cluster", "alternative",
                              "mu", 
                              "METHOD", "DNAME"))
           result <-  do.call("cluswilcox.test.signedrank.ds",
                              c(arglist))
           return(result)
        }
        
        else {
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
            METHOD <- paste( METHOD, "using Datta-Satten method", sep = " ")
             arglist <- setNames(list(x, cluster, group, alternative,
                                 DNAME, METHOD, exact),
                            c("x", "cluster", "group", "stratum",
                              "alternative", "DNAME", "METHOD",
                              "exact"))
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
        
        
    
