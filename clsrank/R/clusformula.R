cluswilcox.test.formula <- function(formula, data, subset, na.action, 
                                    group.x = NULL, permutation = FALSE,
                                    ...) {
  ## This function is only used for rank sum test. 
  ## Mainly to process data.
  ## Inputs:
  ##  formula: a formula of the form "lhs ~ rhs",
  ##    where the lhs is a numeric variable giving the 
  ##    data values, the observed score and the rhs 
  ##    is of the form with special terms: 
  ##    cluster(x1) + group(x2) + stratum(x3), where
  ##    x1, x2, x3 are the corresponding variables.
  ##  data: an optional matrix or dataframe of data 
  ##    used in the formula.
  ##  subset: an optional vector specifying
  ##   a subset of observations to be used.
  ##  na.action: a function which indicates what should happen when 
  ##    the data contain NAs. Defaults to getOption("na.action").
  ##
  ##
  ##
  
  METHOD <- "Wilcoxon rank sum test for clutered data"
  
  Call <- match.call()
  
  DNAME <- paste("from", Call$data)
  
  
  indx <- match(c("formula", "data", "subset", "na.action"),
                names(Call), nomatch = 0)
  if(indx[1] == 0)
    stop(" A formula argument is required for clustered rank sum test")
  temp <- Call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")
  special <- c("strata", "cluster", "group")
  temp$formula <- if(missing(data))
    terms(formula, special)
  else terms(formula, special, data = data)
  
  cluster <- function(x) {x}
  
  stratum <- function(x) {x}
  
  group <- function(x) {x}
  
  
  mf <- eval(temp, parent.frame())
  
  if(nrow(mf) == 0) 
    stop("No (non-missing) observations")
  
  Terms <- terms(mf)
  extraArgs <- list(...)
  
  x <- model.extract(mf, "response")
  if(is.vector(x)) {
    data.n <- length(x)
  } else {
    data.n <- nrow(x)
  }
  
  
  strats <- attr(Terms, "specials")$strata
  
  if(length(strats)) {
    DNAME <- paste("stratum id:", strats, DNAME)
    stemp <- untangle.specials(Terms, "strata", 1)
    if(length(stemp$vars) == 1) {
      strata.keep <- mf[[stemp$vars]]
    } else {
      stop("more than one variable are set as the stratum id")
    }
    strata.uniq <- unique(strata.keep)
    strata.uniq.l <- length(strata.uniq)
    
    if(is.character(strata.uniq)) {
      strata <- recoderFunc(strata.keep, strata.uniq, c(1 : strata.uniq.l))
    } else {
      if(!is.numeric(strata.uniq)) {
        stop("stratum id should be numeric or character")
      }
      strata <- strata.keep
    }
  } else {
    strats <- rep(1, data.n)
  }
  
  cluster <- attr(Terms, "specials")$cluster
  if(length(cluster)) {
    DNAME <- paste("cluster id:", cluster, DNAME)
    ctemp <- untangle.specials(Terms, "cluster", 1)
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
    cluster <- c(1 : data.n)
  }
  
  group <- attr(Terms, "specials")$group
  if(length(group)) {
    DNAME <- paste("group id:", group, DNAME)
    gtemp <- untangle.specials(Terms, "group", 1)
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
    
    if(group.uniq.l > 2L) {
      stop("can only handle 2 groups in rank sum test")
    }
    if(!is.null(group.x)) {
      if(! group.x %in% group.uniq) {
        stop("label of group x is not valid")
      }
    } else {
      group.x <- min(group.uniq)
    }
    group.y <- group.uniq[group.uniq != group.x]
    group <- recoderFunc(group.keep, c(group.x, group.y), c(1, 2))
  } else {
    stop("no group id for rank sum test")
  }
  
  OK <- complete.cases(x)
  finite.x <- is.finite(x)
  x <- x[OK && finite.x]
  group <- group[OK && finite.x]
  stratum <- group[OK && finite.x]
  cluster <- cluster[OK && finite.x]
  
  if(permutation == FALSE) {
    return(cluswilcox.test.ranksum(x,  cluster, group, strats))
  } else {
    return(cluswilcox.test.ranksum.permutation(x,  cluster, group, strats, ...))
  }
  
  

} 
