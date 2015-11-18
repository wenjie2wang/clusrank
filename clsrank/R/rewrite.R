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
  
  ## If data name existed, take out the x (and y) observations,
  ## group id, cluster id, stratum id, otherwise, no need to 
  ## take values from a data frame. 

  if(!is.null(pars$data)) {
    x <- data[, as.character(pars$x)]
    DNAME <- (pars$x)
    
    if(!is.null(pars$y)) {
      y <- data[, as.character(pars$y)]
      DNAME <- paste(DNAME, "and", pars$y)
    } else {
      y <- NULL
    }
    if(!is.null(pars$group)) {
      group <- data[, as.character(pars$group)]
      DNAME <- paste0(DNAME, ", group id: ", pars$group)
    } else {
      group <- NULL
    }
    if(!is.null(pars$id)) {
      id <- data[, as.character(pars$id)]
      DNAME <- paste0(DNAME, ", cluster id: ", pars$id)
    } else {
      id <- NULL
    }
    if(!is.null(pars$stratum)) {
      stratum <- data[, as.character(pars$stratum)]
      DNAME <- paste0(DNAME, ", stratum id: ", pars$stratum)
    } else {
      stratum <- NULL
    }
    DNAME <- paste0(DNAME, "from", pars$data)
  } else {
    DNAME <- deparse(substitute(x))
    if(!is.null(y)) {
      DNAME <- paste(DNAME, "and", deparse(substitute(y)))
    } 
    if(!is.null(group)) {
      DNAME <- paste0(DNAME, ", group id: ", deparse(substitute(group)))
    } 
    if(!is.null(id)) {
      DNAME <- paste0(DNAME, ", cluster id: ", deparse(substitute(id)))
    } 
    if(!is.null(stratum)) {
      DNAME <- paste0(DNAME, ", stratum id: ", deparse(substitute(stratum)))
    } 
  }
  
  ## If group variable provided, then x comes in 
  ## under two treatments, i.e., it is not that
  ## x is from treatment 1 and y from treatment 2.
  ## Or x consistis of data from pre-treatment and
  ## post-treatment. Need "paired" option to decide.
  ## Will give warning if both "group" and "y"
  ## are available.
  
  
  if(!is.null(group)) {
    if (!isTRUE(all(group == floor(group))) ||
        !(is.character(group))) 
      stop("'group' must only contain integers or
           characters")
    
    group.uniq <- unique(group)
    group.uniq.freq <- table(group)
    
    if(length(group.uniq) > 2)
      stop("'group' can only have no more than two levels.")
    
    if(length(group) != length(x))
      stop("'group' must have the same length as 'x'")
    
    if(!is.null(y) || length(group.uniq) != 1) {
      warning("Both 'y' and 'group' are provided, only need one,
              'y' is ignored.")
    }
    if(paired == TRUE && is.null(y) 
       && group.uniq.freq[1] != group.uniq.freq[2]) {
      stop("For paired test, unequal sample sizes
           from pre-treatment and post-treatment")
    }
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
    
    ## If the group.x is not provided, use the smaller
    ## id as the x group, either integer or character.
    
    if(is.null(group.x)) {
      group.x <- min(unique(group))
    }
    group[group == group.x] <- 1
    group[group != group.x] <- 2
    
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