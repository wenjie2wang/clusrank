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
  if(is.null(pars$data)) {
    data <- NULL
  }
  
  
  ##  Process the data use cluswilcox.test.data 
  
  
  CORRECTION <- 0
  
  
  
  ## Main body of tests.
  
  if(!is.null(pars$group)) {
    
  }
  
  if(!is.null(y)) {
    
  }
  
  
}