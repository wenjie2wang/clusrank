cluswilcox.test.signedrank.rgl <- function(x, cluster, alternative,
                                           mu, DNAME, METHOD) {
  METHOD <- paste(METHOD, "using Rosner-Glynn-Lee method")
  ## Drop the ties
  data <- data.frame(x, cluster)
  data <- data[x != 0, ]
  cluster.size <- table(data$cluster)
  m <- length(cluster.size)
  n <- nrow(data)
  if (length(table(cluster.size)) != 1) {
    balance = FALSE
  } else {
    balance = TRUE
  }
  
  xrank <- rank(abs(data$x))
  data <- cbind(data, xrank)
  signrank <- ifelse(data$x > 0, 1, -1) * data$xrank
  data <- cbind(data, signrank)
  colnames(data)[4] <- "signrank"
  if(balance == TRUE){
    T_c <- sum(data$signrank)
    sumrank <- c(by(data$signrank, data$cluster, sum))
    sumsq <- sum(sumrank ^ 2)
    Var_t <- sumsq
    W_c <- T_c / sqrt(Var_t)
    P_val <- 2 * (1 - pnorm(abs(W_c)))
    
    ADJUST <- FALSE
    names(T_c) <- "rank statistic"
    names(W_c) <- "test statistic"
    names(Var_t) <- "variance of rank statistic"
    
    names(n) <- "total number of observations"
    names(m) <- "total number of clusters"
    names(Var_t) <- paste("Variance of ", names(T_c))
    result <- list(rstatistic = T_c, vrstatistic = Var_t, statistic = W_c,
                   p.value = P_val, n = n, cn = m, null.value = mu,
                   alternative = alternative,
                   data.name = DNAME, method = METHOD,
                   adjusted = ADJUST)
    class(result) <- "ctest"
    return(result)
  } else {
    if (balance == FALSE) {
      sumclusterrank <- c(by(data$signrank, data$cluster, sum))
      sumsq <- sum(sumclusterrank ^ 2)
      meansumrank <- sumclusterrank / cluster.size
      sumsqi <- sum(cluster.size ^ 2)
      
      # calculate intraclass correlation between signed ranks within the same cluster
      data$cluster.f <- as.factor(data$cluster)
      mod <- lm(signrank ~ cluster.f, data, y = TRUE)
      errordf <- mod$df.residual
      errorss <- sum(mod$residuals ^ 2)
      modeldf <- n - errordf - 1
      modelss <- sum((mod$y - mean(mod$y)) ^ 2) - errorss
      sumi <- n
      
      m0 <- (sumi - (sumsqi / sumi)) / (m - 1)
      totalss <- errorss + modelss
      totaldf <- errordf + modeldf
      wthms <- errorss / errordf
      betms <- modelss / modeldf
      vars <- totalss / totaldf
      s2b <- (betms - wthms) / m0
      s2w <- wthms
      rosglm <- s2b / (s2b + s2w)
      if (rosglm < 0) {
        rosglm = 0
      }
      ros <- rosglm
      roscor <- ros * (1 + (1 - ros ^ 2) / (m - 2.5))
      if (roscor > 1) {
        roscor = 1
      }
      wi <- cluster.size / (vars * (1 + (cluster.size - 1) * roscor))
      T_c <- sum(meansumrank * wi)
      sqweightsum <- sum(wi ^ 2 * meansumrank ^ 2)
      Var_t <- sqweightsum
      W_c <-  T_c / (sqrt(Var_t))
      P_val <- switch(alternative,
                      less = pnorm(abs(W_c)),
                      greater = pnorm(abs(W_c), lower.tail = FALSE),
                      two.sided = 2 * min(pnorm(abs(W_c)),
                                          pnorm(abs(W_c), lower.tail = FALSE)))
      
      names(T_c) <- "adjusted rank statistic"
      names(W_c) <- "test statistic"
      
      ADJUST <- TRUE
      
      names(n) <- "total number of observations"
      names(m) <- "total number of clusters"
      names(Var_t) <- paste("Variance of ", names(T_c))
      result <- list(rstatistic = T_c, vrstatistic = Var_t, statistic = W_c,
                     p.value = P_val, n = n, cn = m,
                     alternative = alternative,
                     null.value = mu,
                     data.name = DNAME, method = METHOD,
                     adjusted = ADJUST)
      class(result) <- "ctest"
      return(result)
    }
}