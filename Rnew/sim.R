library(mvtnorm)

#' @title Generating clustered data for Wilcoxon test
#' @param nclus the number of clusters
#'
#' @return a data.frame with columns: score, group, id
datgen <- function(nclus, maxclsize, delta = 0.,
                   rho = c(0.9, -.1), rate = 1,
                   grplevel = c("cluster", "individual")) {
    nn <- nclus * maxclsize
    Sigma1 <- diag(1 - rho[1], nrow = maxclsize) + rho[1]
    x <- c(t(rmvnorm(nclus, sigma = Sigma1)))
    Sigma2 <- diag(1 - rho[2], nrow = maxclsize) + rho[2]
    y <- c(t(rmvnorm(nclus, sigma = Sigma2)))
    group <- rep(c(1, 2), each = nn)
    if (grplevel == "individual")
        group <- sample(group, nn, FALSE)
    id <- rep(1:(2 * nclus), each = maxclsize)
    score <- exp(c(x, y)) + delta * group
    dat <- data.frame(score = score, group = group, id = id)
    keep <- sort(sample(1:(2 * nn), size = rate * (2 * nn), FALSE))
    if (rate == 1) dat else dat[keep,]
}


simpower <- function(nrep, level, delta, rho, rate, grplevel, nclus, maxclsize) {
    do1rep <- function() {
        dat <- datgen(nclus, maxclsize, delta, rho, rate, grplevel);
        p.rgl <- cluswilcox.test(score, cluster = id, group = grp, data = dat, method = "rgl")$p.value
        p.ds  <- cluswilcox.test(score, cluster = id, group = grp, data = dat, method = "ds" )$p.value
        c(rgl = p.rgl, ds = p.ds)
    }
    sim <- t(replicate(nrep, do1rep()))
    apply(sim, 2, function(x) mean(x < level))
}
