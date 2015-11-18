cluswilcox.test <- function(x, ...) {
  UseMethod("cluswilcox.test")
}

cluswilcox.test.default <- function(x, y = NULL, group = NULL,
                                    id = NULL, stratum = NULL, 
                                    data = parent.frame(),
                                    group.x = NULL,
                                    alternative = c("two.sided", "less", "greater"),
                                    mu = 0, paired = FALSE, exact = NULL, 
                                    correct = FALSE, conf.int = FALSE,
                                    conf.level = 0.95, ...) {
  pars <- as.list(match.call()[-1])
  
  if(!is.null(pars$data)) {
    x <- data[, as.character(pars$x)]
    if(!is.null(pars$y)) {
      y <- data[, as.character(pars$y)]
    } else {
      y <- NULL
    }
  }
  
  if(!is.null(group)) {
    if (!isTRUE(all(group == floor(group))) ||
        !(is.character(group))) 
      stop("'group' must only contain integer values or
           characters")
    if(length(unique(group)) > 2)
      stop("'group' can only have no more than two level.")
  }
  
  
  alternative <- match.arg(alternative)
  if(!missing(mu) && ((length(mu) > 1L) || !is.finite(mu))) {
    stop("'mu' must be a single number")
  }
  if(conf.int) {
    if (!((length(conf.level) == 1L) && is.finite(conf.level) && 
          (conf.level > 0) && (conf.level < 1))) 
      stop("'conf.level' must be a single number between 0 and 1")
  }
  
  if (!is.numeric(x)) 
    stop("'x' must be numeric")
  if (!is.null(y)) {
    if (!is.numeric(y)) 
      stop("'y' must be numeric")
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    if (paired) {
      if (length(x) != length(y)) 
        stop("'x' and 'y' must have the same length")
      OK <- complete.cases(x, y)
      x <- x[OK] - y[OK]
      y <- NULL
    }
    else {
      x <- x[is.finite(x)]
      y <- y[is.finite(y)]
    }
  } else {
    if(!is.null(group.x)) {
      group[group == group.x] <- 1
      group[group != group.x] <- 2
    }
    ## If paired is TRUE, assume x vector is the difference vector
    DNAME <- deparse(substitute(x))
    x <- x[is.finite(x)]
  }
  if (length(x) < 1L) 
    stop("not enough (finite) 'x' observations")
  CORRECTION <- 0
  
  
  
  ## Main body of tests.
  
  if(!is.null(pars$group)) {
    
  }
  
  if(!is.null(y)) {
    
  }
  
  
}