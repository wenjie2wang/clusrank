cluswilcox.test.default <- function(x, y, group,
                                    id, stratum, 
                                    data,
                                    group.x,
                                    alternative,
                                    mu, paired, exact, 
                                    ...) {
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
  
  
  ## Check and initialize id if not given,
  ## transform it to numeric if given as characters.
  
  l.x <- length(x)
  
  if( is.null(id)) {
    id <- c(1 : l.x)
  } else {
    if(!is.numeric(id)) {
      if(!is.character(id)) {
        stop("'id' has to be numeric or characters")
      }
      if(length(id) != l.x) {
        stop("'id' and 'x' must have the same lengths")
      }
      uniq.id <- unique(id)
      l.uniq.id <- length(uniq.id)
      id <- as.numeric(recoderFunc(id, uniq.id, c(1 : l.uniq.id)))
    }
  }
  
  
  
  ## Check x.
  if ( !is.numeric(x)) 
    stop("'x' must be numeric")
  
  ## If y is present, carry out paired test.
  
  if( !is.null(y)) {
    paired <- TRUE
  }
  
  ## Check data for paired test, paired test
  ## do not deal with stratified data.
  
  if( paired == TRUE) {
    
    if(!is.null(stratum)) {
      warning("Signed rank test does not handle stratified data,
              stratum id will be ignored.")
    }
    
    if( !is.null(y)) {
      
      if (!is.numeric(y)) 
        stop("'y' must be numeric")
      
      l.y <- length(y)
      
      if( l.y != l.x) {
        stop("'x' and 'y' must have the same 
             lengths for signed rank test.")
      }
      
      OK <- complete.cases(x, y, id)
      x <- x[OK] - y[OK]
      id <- id[OK]
      finite.x <- is.finite(x)
      x <- x[finite.x]
      id <- id[finite.x]
      
      if(length(x) < 1L) {
                  stop("not enough (finite) 'x' observations")
      }
      
      } else {
        
        ## If only x is given, it is the difference score.
        
        OK <- complete.cases(x, id)
        x <- x[OK]
        id <- id[OK]
        finite.x <- is.finite(x)
        x <- x[finite.x]
        id <- id[finite.x]
        if(length(x) < 1L) {
          stop("not enough (finite) 'x' observations")
        }
      }
    
    if(!is.null(group)) {
      warning("'group' variable is ignored for signed rank test.")
    }
    
    } else {
      ## Rank sum test.
      if( is.null(group)) {
        stop( "'group id' required for rank sum test.")
      }
      
      l.group <- length(group) 
      if( l.group != l.x) {
        stop( "'x' and 'group' must have the same lengths for rank
              sum test.")
      }
      
      
      if (!isTRUE(all(group == floor(group))) ||
          !(is.character(group))) 
        stop("'group' must only contain integers or
             characters")
      
      group.uniq <- unique(group)
      group.uniq.freq <- table(group)
      
      if(length(group.uniq) > 2)
        stop("'group' can only have no more than two levels.")
      
      ## If the group.x is not provided, use the smaller
      ## id as the x group, either integer or character.
      
      if(is.null(group.x)) {
        group.x <- min(unique(group))
      }
      group[group == group.x] <- 1
      group[group != group.x] <- 2
      
      if( is.null(stratum)) {
        stratum <- rep(1, l.x)
      }
   
      }
  
  
  
  alternative <- match.arg(alternative)
  if(!missing(mu) && ((length(mu) > 1L) || !is.finite(mu))) {
    stop("'mu' must be a single number")
  }  
  
  if(paired == TRUE) {
    result <- list(x = x, id = id, alternative = alternative, 
                   mu = mu, exact = exact)
    class(result) <- "signedrank"
    return(result)
  } else {
    result <- list(x = x, id = id, group = group, stratum = stratum,
                   alternative = alternative, 
                   mu = mu, exact = exact)
    class(result) <- "ranksum"
    return(result)
  }

  
  
  

  
  
  
  
}