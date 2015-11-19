cluswilcox.test.data <- function(x, y, group,
                            id, stratum, 
                            data,
                            group.x,
                            alternative,
                            mu, paired, exact, 
                            correct, conf.int,
                            conf.level, ...) {
  ## Process the input arguments before feeding them to
  ## main test functions. Assign a class (better to 
  ## be S4) to the processed arguments for them to be
  ## sent to the corresponding functions.
  ##
  ## Inputs:
  ##   The same as cluswilcox.test.default.
  ##   x: numeric vector of data values. Non-finite 
  ##     (e.g., infinite or missing) values will be omitted.
  ##     If to carry out wilcoxon sum rank test, x is a vector
  ##     contains observations from both treatments, the
  ##     treatment id is the 'group' variable; If to 
  ##     carry out signed rank test, x is either a vector contains
  ##     the difference between pre and post treatments; or a 
  ##     vector contains 
  ##     observations from both pre and post treatments, and the
  ##     treatment id is the 'group' variable; or contains 
  ##     observations from pre treatment while 'y' contains
  ##     observations from post treatment. 
  ##
  ##   y: an optional numeric vector of data values: 
  ##     as with x non-finite values will be omitted.
  ##     Contains observations from post treatment, is not 
  ##     NULL only for signed rank test.
  ##
  ##   group: a vector of either integers or characters.
  ##     Cannot take more than 2 values. Treatment id, 
  ##     Id for group X is
  ##     chosen by 'group.x' if provided, otherwise, will 
  ##     choose the minimum of two distinct values as the 
  ##     id for group X, for both integers and characters. 
  ##
  ##   id:  an integer vector. Cluster id.If not provided, 
  ##     assume there is no cluster.
  ##
  ##   data: an optional matrix or data frame 
  ##     (or similar: see model.frame) containing the variables.
  ##     By default the variables are taken from environment(formula).
  ##  
  ##   stratum: an integer vector. Stratum id. If not provided, 
  ##     assume there is no stratification.
  ##
  ##   group.x: an integer or a character, chosen as the
  ##     treatment X for rank sum test or pre treatment
  ##     for signed rank test.
  ##     
  ##   alternative: a character string specifying the 
  ##     alternative hypothesis, must be one of 
  ##    "two.sided" (default), "greater" or "less". 
  ##     You can specify just the initial letter.
  ##
  ##   mu: a number specifying an optional parameter 
  ##       used to form the null hypothesis.
  ##
  ##   paired: a logical indicating whether you want a paired test.
  ##
  ##   exact:	a logical indicating whether an exact 
  ##     p-value should be computed.
  ##
  ##   conf.int: a logical indicating whether a confidence 
  ##     interval should be computed.
  ##   
  ##   conf.level: confidence level of the interval.
  ##
   
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
    
    ## If the group.x is not provided, use the smaller
    ## id as the x group, either integer or character.
    
    if(is.null(group.x)) {
      group.x <- min(unique(group))
    }
    group[group == group.x] <- 1
    group[group != group.x] <- 2
    
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
  
  
  ## Initialize id, group, stratum if not provided.
  
  l.x <- length(x)
  
  if( is.null(id)) {
    id <- c(1 : l.x)
  }
  
  if( is.null(group)) {
    if ( !paired) {
      stop("'group' variable missing for two sample rank sum test")
    }
    group <- rep(1, l.x)
  }
  
  if( is.null(stratum)) {
    stratum <- rep(1, l.x)
  }
  
  
  
  
  if (!is.numeric(x)) 
    stop("'x' must be numeric")
  if (!is.null(y)) {
    if (!is.numeric(y)) 
      stop("'y' must be numeric")
    if (paired) {
      if (length(x) != length(y)) 
        stop("'x' and 'y' must have the same length")
      OK <- complete.cases(x, y, id, )
      x <- x[OK] - y[OK]
      id <- id[OK]
      y <- NULL
    }
    else {
      ## If not for paired test, "y" vector will be ignored. 
      warning("'y' vector is ignored for two sample rank sum test")
      OK <- complete.cases(x, id)
      x <- x[OK]
      id <- id[OK]
      x <- x[is.finite(x)]
    }
  } else {
    
    
    
    ## If paired is TRUE, assume x vector is the difference vector
    DNAME <- deparse(substitute(x))
    x <- x[is.finite(x)]
  }
  if (length(x) < 1L) 
    stop("not enough (finite) 'x' observations")
}