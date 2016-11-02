

ex  <- function(dim, rho) {
    diag(1 - rho, dim) + matrix(rho, dim, dim)
}
ar1 <- function(dim, rho) {
    rho ^ outer(1:dim, 1:dim, function(x, y) abs(x - y))
}



library(mvtnorm)
datgen.sum <- function(nclus, maxclsize, delta = 0., rho = c(0.1, 0.1),
                       corr = ex, misrate = 0., clusgrp = TRUE) {
    nn <- nclus * maxclsize
    Sigma1 <- corr(maxclsize, rho[1])
    Sigma2 <- corr(maxclsize, rho[2])
    y1 <- c(t(rmvnorm(nclus, sigma = Sigma1)))
    y2 <- c(t(rmvnorm(nclus, sigma = Sigma2)))
    group <- rep(c(0, 1), each = nn)
    if (!clusgrp) group  <- sample(group, nn, FALSE)
    cid <- rep(1:(2 * nclus), each = maxclsize)
    x <- exp(c(y1, y2)) + delta * group
    dat <- data.frame(x = x, grp = group, cid = cid)
    drop <-  sort(sample(1:(2 * nn), size = misrate * (2 * nn), FALSE))
    if (misrate == 0.) dat else dat[-drop, ]
}


datgen.sgn <- function(nclus, maxclsize, delta = 0., rho = 0.1,
                       corr = ex, misrate = 0.) {
    nn <- nclus * maxclsize
    Sigma <- corr(maxclsize, rho)
    z <- delta + c(t(rmvnorm(nclus, sigma = Sigma)))
    x <- sign(z) * exp(abs(z))
    cid <- rep(1:nclus, each = maxclsize)
    dat <- data.frame(x = x, cid = cid)
    drop <- sort(sample(1:nn, size = (1 - misrate) * nn, FALSE))
    if (misrate == 0.) dat else dat[-drop,]
}



library(clusrank)
simpower <- function(nrep, level, paired, nclus, maxclsize,
                     delta, rho, corr, misrate, ...) {
    do1rep <- function() {
        datgen <- if (paired) datgen.sgn else datgen.sum
        formula <- if (paired) x ~ cluster(cid)
                   else x ~ cluster(cid) + grp
        dat <- datgen(nclus, maxclsize, delta, rho, corr, misrate, ...)
        p.rgl <- clusWilcox.test(formula, paired = paired,
                                 data = dat, method = "rgl")$p.value
        p.ds  <- clusWilcox.test(formula, paired = paired,
                                 data = dat, method = "ds" )$p.value
        c(rgl = p.rgl, ds = p.ds)
    }
    sim <- t(replicate(nrep, do1rep()))
    apply(sim, 2, function(x) mean(x < level))
}






rs.b.e.11.2 <- marix(0, 3, 6)
rs.b.e.11.2[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 2, 0, c(0.1, 0.1), ex, 0)
rs.b.e.11.2[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 2, 0.2, c(0.1, 0.1), ex, 0)
rs.b.e.11.2[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 2, 0.5, c(0.1, 0.1), ex, 0)
rs.b.e.11.2[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 2, 0, c(0.1, 0.1), ex, 0)
rs.b.e.11.2[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 2, 0.2, c(0.1, 0.1), ex, 0)
rs.b.e.11.2[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 2, 0.5, c(0.1, 0.1), ex, 0)
rs.b.e.11.2[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 2, 0, c(0.1, 0.1), ex, 0)
rs.b.e.11.2[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 2, 0.2, c(0.1, 0.1), ex, 0)
rs.b.e.11.2[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 2, 0.5, c(0.1, 0.1), ex, 0)





rs.b.e.55.2 <- marix(0, 3, 6)
rs.b.e.55.2[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 2, 0, c(0.5, 0.5), ex, 0)
rs.b.e.55.2[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 2, 0.2, c(0.5, 0.5), ex, 0)
rs.b.e.55.2[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 2, 0.5, c(0.5, 0.5), ex, 0)
rs.b.e.55.2[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 2, 0, c(0.5, 0.5), ex, 0)
rs.b.e.55.2[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 2, 0.2, c(0.5, 0.5), ex, 0)
rs.b.e.55.2[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 2, 0.5, c(0.5, 0.5), ex, 0)
rs.b.e.55.2[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 2, 0, c(0.5, 0.5), ex, 0)
rs.b.e.55.2[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 2, 0.2, c(0.5, 0.5), ex, 0)
rs.b.e.55.2[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 2, 0.5, c(0.5, 0.5), ex, 0)







rs.b.e.19.2 <- marix(0, 3, 6)
rs.b.e.19.2[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 2, 0, c(-0.1, 0.9), ex, 0)
rs.b.e.19.2[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 2, 0.2, c(-0.1, 0.9), ex, 0)
rs.b.e.19.2[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 2, 0.5, c(-0.1, 0.9), ex, 0)
rs.b.e.19.2[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 2, 0, c(-0.1, 0.9), ex, 0)
rs.b.e.19.2[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 2, 0.2, c(-0.1, 0.9), ex, 0)
rs.b.e.19.2[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 2, 0.5, c(-0.1, 0.9), ex, 0)
rs.b.e.19.2[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 2, 0, c(-0.1, 0.9), ex, 0)
rs.b.e.19.2[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 2, 0.2, c(-0.1, 0.9), ex, 0)
rs.b.e.19.2[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 2, 0.5, c(-0.1, 0.9), ex, 0)







rs.b.e.11.5 <- marix(0, 3, 6)
rs.b.e.11.5[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 5, 0, c(0.1, 0.1), ex, 0)
rs.b.e.11.5[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 5, 0.2, c(0.1, 0.1), ex, 0)
rs.b.e.11.5[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 5, 0.5, c(0.1, 0.1), ex, 0)
rs.b.e.11.5[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 5, 0, c(0.1, 0.1), ex, 0)
rs.b.e.11.5[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 5, 0.2, c(0.1, 0.1), ex, 0)
rs.b.e.11.5[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 5, 0.5, c(0.1, 0.1), ex, 0)
rs.b.e.11.5[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 5, 0, c(0.1, 0.1), ex, 0)
rs.b.e.11.5[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 5, 0.2, c(0.1, 0.1), ex, 0)
rs.b.e.11.5[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 5, 0.5, c(0.1, 0.1), ex, 0)





rs.b.e.55.5 <- marix(0, 3, 6)
rs.b.e.55.5[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 5, 0, c(0.5, 0.5), ex, 0)
rs.b.e.55.5[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 5, 0.2, c(0.5, 0.5), ex, 0)
rs.b.e.55.5[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 5, 0.5, c(0.5, 0.5), ex, 0)
rs.b.e.55.5[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 5, 0, c(0.5, 0.5), ex, 0)
rs.b.e.55.5[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 5, 0.2, c(0.5, 0.5), ex, 0)
rs.b.e.55.5[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 5, 0.5, c(0.5, 0.5), ex, 0)
rs.b.e.55.5[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 5, 0, c(0.5, 0.5), ex, 0)
rs.b.e.55.5[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 5, 0.2, c(0.5, 0.5), ex, 0)
rs.b.e.55.5[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 5, 0.5, c(0.5, 0.5), ex, 0)











rs.b.e.19.5 <- marix(0, 3, 6)
rs.b.e.19.5[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 5, 0, c(-0.1, 0.9), ex, 0)
rs.b.e.19.5[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 5, 0.2, c(-0.1, 0.9), ex, 0)
rs.b.e.19.5[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 5, 0.5, c(-0.1, 0.9), ex, 0)
rs.b.e.19.5[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 5, 0, c(-0.1, 0.9), ex, 0)
rs.b.e.19.5[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 5, 0.2, c(-0.1, 0.9), ex, 0)
rs.b.e.19.5[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 5, 0.5, c(-0.1, 0.9), ex, 0)
rs.b.e.19.5[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 5, 0, c(-0.1, 0.9), ex, 0)
rs.b.e.19.5[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 5, 0.2, c(-0.1, 0.9), ex, 0)
rs.b.e.19.5[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 5, 0.5, c(-0.1, 0.9), ex, 0)








rs.ub.e.11.10 <- marix(0, 3, 6)
rs.ub.e.11.10[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 10, 0, c(0.1, 0.1), ex, 0.5)
rs.ub.e.11.10[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 10, 0.2, c(0.1, 0.1), ex, 0.5)
rs.ub.e.11.10[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 10, 0.5, c(0.1, 0.1), ex, 0.5)
rs.ub.e.11.10[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 10, 0, c(0.1, 0.1), ex, 0.5)
rs.ub.e.11.10[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 10, 0.2, c(0.1, 0.1), ex, 0.5)
rs.ub.e.11.10[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 10, 0.5, c(0.1, 0.1), ex, 0.5)
rs.ub.e.11.10[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 10, 0, c(0.1, 0.1), ex, 0.5)
rs.ub.e.11.10[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 10, 0.2, c(0.1, 0.1), ex, 0.5)
rs.ub.e.11.10[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 10, 0.5, c(0.1, 0.1), ex, 0.5)










rs.ub.e.55.10 <- marix(0, 3, 6)
rs.ub.e.55.10[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 10, 0, c(0.5, 0.5), ex, 0.5)
rs.ub.e.55.10[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 10, 0.2, c(0.5, 0.5), ex, 0.5)
rs.ub.e.55.10[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 10, 0.5, c(0.5, 0.5), ex, 0.5)
rs.ub.e.55.10[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 10, 0, c(0.5, 0.5), ex, 0.5)
rs.ub.e.55.10[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 10, 0.2, c(0.5, 0.5), ex, 0.5)
rs.ub.e.55.10[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 10, 0.5, c(0.5, 0.5), ex, 0.5)
rs.ub.e.55.10[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 10, 0, c(0.5, 0.5), ex, 0.5)
rs.ub.e.55.10[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 10, 0.2, c(0.5, 0.5), ex, 0.5)
rs.ub.e.55.10[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 10, 0.5, c(0.5, 0.5), ex, 0.5)







rs.ub.e.19.10 <- marix(0, 3, 6)
rs.ub.e.19.10[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 10, 0, c(-0.1, 0.9), ex, 0.5)
rs.ub.e.19.10[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 10, 0.2, c(-0.1, 0.9), ex, 0.5)
rs.ub.e.19.10[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 10, 0.5, c(-0.1, 0.9), ex, 0.5)
rs.ub.e.19.10[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 10, 0, c(-0.1, 0.9), ex, 0.5)
rs.ub.e.19.10[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 10, 0.2, c(-0.1, 0.9), ex, 0.5)
rs.ub.e.19.10[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 10, 0.5, c(-0.1, 0.9), ex, 0.5)
rs.ub.e.19.10[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 10, 0, c(-0.1, 0.9), ex, 0.5)
rs.ub.e.19.10[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 10, 0.2, c(-0.1, 0.9), ex, 0.5)
rs.ub.e.19.10[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 10, 0.5, c(-0.1, 0.9), ex, 0.5)





rs.b.e.11.2.sub <- marix(0, 3, 6)
rs.b.e.11.2.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 2, 0, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.2.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 2, 0.2, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.2.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 2, 0.5, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.2.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 2, 0, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.2.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 2, 0.2, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.2.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 2, 0.5, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.2.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 2, 0, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.2.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 2, 0.2, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.2.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 2, 0.5, c(0.1, 0.1), ex, 0, FALSE)





rs.b.e.55.2.sub <- marix(0, 3, 6)
rs.b.e.55.2.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 2, 0, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.2.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 2, 0.2, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.2.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 2, 0.5, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.2.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 2, 0, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.2.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 2, 0.2, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.2.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 2, 0.5, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.2.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 2, 0, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.2.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 2, 0.2, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.2.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 2, 0.5, c(0.5, 0.5), ex, 0, FALSE)







rs.b.e.19.2.sub <- marix(0, 3, 6)
rs.b.e.19.2.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 2, 0, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.2.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 2, 0.2, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.2.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 2, 0.5, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.2.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 2, 0, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.2.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 2, 0.2, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.2.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 2, 0.5, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.2.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 2, 0, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.2.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 2, 0.2, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.2.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 2, 0.5, c(-0.1, 0.9), ex, 0, FALSE)







rs.b.e.11.5.sub <- marix(0, 3, 6)
rs.b.e.11.5.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 5, 0, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.5.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 5, 0.2, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.5.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 5, 0.5, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.5.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 5, 0, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.5.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 5, 0.2, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.5.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 5, 0.5, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.5.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 5, 0, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.5.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 5, 0.2, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.5.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 5, 0.5, c(0.1, 0.1), ex, 0, FALSE)





rs.b.e.55.5.sub <- marix(0, 3, 6)
rs.b.e.55.5.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 5, 0, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.5.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 5, 0.2, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.5.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 5, 0.5, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.5.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 5, 0, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.5.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 5, 0.2, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.5.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 5, 0.5, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.5.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 5, 0, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.5.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 5, 0.2, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.5.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 5, 0.5, c(0.5, 0.5), ex, 0, FALSE)











rs.b.e.19.5.sub <- marix(0, 3, 6)
rs.b.e.19.5.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 5, 0, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.5.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 5, 0.2, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.5.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 5, 0.5, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.5.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 5, 0, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.5.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 5, 0.2, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.5.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 5, 0.5, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.5.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 5, 0, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.5.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 5, 0.2, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.5.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 5, 0.5, c(-0.1, 0.9), ex, 0, FALSE)








rs.ub.e.11.10.sub <- marix(0, 3, 6)
rs.ub.e.11.10.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 10, 0, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.e.11.10.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 10, 0.2, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.e.11.10.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 10, 0.5, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.e.11.10.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 10, 0, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.e.11.10.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 10, 0.2, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.e.11.10.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 10, 0.5, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.e.11.10.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 10, 0, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.e.11.10.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 10, 0.2, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.e.11.10.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 10, 0.5, c(0.1, 0.1), ex, 0.5, FALSE)










rs.ub.e.55.10.sub <- marix(0, 3, 6)
rs.ub.e.55.10.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 10, 0, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.e.55.10.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 10, 0.2, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.e.55.10.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 10, 0.5, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.e.55.10.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 10, 0, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.e.55.10.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 10, 0.2, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.e.55.10.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 10, 0.5, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.e.55.10.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 10, 0, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.e.55.10.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 10, 0.2, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.e.55.10.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 10, 0.5, c(0.5, 0.5), ex, 0.5, FALSE)







rs.ub.e.19.10.sub <- marix(0, 3, 6)
rs.ub.e.19.10.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 10, 0, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.e.19.10.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 10, 0.2, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.e.19.10.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 10, 0.5, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.e.19.10.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 10, 0, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.e.19.10.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 10, 0.2, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.e.19.10.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 10, 0.5, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.e.19.10.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 10, 0, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.e.19.10.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 10, 0.2, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.e.19.10.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 10, 0.5, c(-0.1, 0.9), ex, 0.5, FALSE)









%%%%%%%%%%%%%% Below are simulation for signrank ex










sr.b.e.1.2 <- marix(0, 3, 6)
sr.b.e.1.2[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 2, 0, 0.1, ex, 0)
sr.b.e.1.2[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 2, 0.2, 0.1, ex, 0)
sr.b.e.1.2[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 2, 0.5, 0.1, ex, 0)
sr.b.e.1.2[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 2, 0, 0.1, ex, 0)
sr.b.e.1.2[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 2, 0.2, 0.1, ex, 0)
sr.b.e.1.2[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 2, 0.5, 0.1, ex, 0)
sr.b.e.1.2[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 2, 0, 0.1, ex, 0)
sr.b.e.1.2[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 2, 0.2, 0.1, ex, 0)
sr.b.e.1.2[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 2, 0.5, 0.1, ex, 0)





sr.b.e.5.2 <- marix(0, 3, 6)
sr.b.e.5.2[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 2, 0, 0.5, ex, 0)
sr.b.e.5.2[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 2, 0.2, 0.5, ex, 0)
sr.b.e.5.2[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 2, 0.5, 0.5, ex, 0)
sr.b.e.5.2[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 2, 0, 0.5, ex, 0)
sr.b.e.5.2[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 2, 0.2, 0.5, ex, 0)
sr.b.e.5.2[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 2, 0.5, 0.5, ex, 0)
sr.b.e.5.2[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 2, 0, 0.5, ex, 0)
sr.b.e.5.2[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 2, 0.2, 0.5, ex, 0)
sr.b.e.5.2[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 2, 0.5, 0.5, ex, 0)


sr.b.e.9.2 <- marix(0, 3, 6)
sr.b.e.9.2[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 2, 0, 0.9, ex, 0)
sr.b.e.9.2[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 2, 0.2, 0.9, ex, 0)
sr.b.e.9.2[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 2, 0.5, 0.9, ex, 0)
sr.b.e.9.2[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 2, 0, 0.9, ex, 0)
sr.b.e.9.2[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 2, 0.2, 0.9, ex, 0)
sr.b.e.9.2[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 2, 0.5, 0.9, ex, 0)
sr.b.e.9.2[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 2, 0, 0.9, ex, 0)
sr.b.e.9.2[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 2, 0.2, 0.9, ex, 0)
sr.b.e.9.2[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 2, 0.5, 0.9, ex, 0)








sr.b.e.1.10 <- marix(0, 3, 6)
sr.b.e.1.10[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 10, 0, 0.1, ex, 0)
sr.b.e.1.10[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 10, 0.2, 0.1, ex, 0)
sr.b.e.1.10[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 10, 0.5, 0.1, ex, 0)
sr.b.e.1.10[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 10, 0, 0.1, ex, 0)
sr.b.e.1.10[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 10, 0.2, 0.1, ex, 0)
sr.b.e.1.10[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 10, 0.5, 0.1, ex, 0)
sr.b.e.1.10[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 10, 0, 0.1, ex, 0)
sr.b.e.1.10[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 10, 0.2, 0.1, ex, 0)
sr.b.e.1.10[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 10, 0.5, 0.1, ex, 0)





sr.b.e.5.10 <- marix(0, 3, 6)
sr.b.e.5.10[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 10, 0, 0.5, ex, 0)
sr.b.e.5.10[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 10, 0.2, 0.5, ex, 0)
sr.b.e.5.10[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 10, 0.5, 0.5, ex, 0)
sr.b.e.5.10[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 10, 0, 0.5, ex, 0)
sr.b.e.5.10[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 10, 0.2, 0.5, ex, 0)
sr.b.e.5.10[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 10, 0.5, 0.5, ex, 0)
sr.b.e.5.10[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 10, 0, 0.5, ex, 0)
sr.b.e.5.10[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 10, 0.2, 0.5, ex, 0)
sr.b.e.5.10[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 10, 0.5, 0.5, ex, 0)


sr.b.e.9.10 <- marix(0, 3, 6)
sr.b.e.9.10[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 10, 0, 0.9, ex, 0)
sr.b.e.9.10[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 10, 0.2, 0.9, ex, 0)
sr.b.e.9.10[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 10, 0.5, 0.9, ex, 0)
sr.b.e.9.10[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 10, 0, 0.9, ex, 0)
sr.b.e.9.10[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 10, 0.2, 0.9, ex, 0)
sr.b.e.9.10[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 10, 0.5, 0.9, ex, 0)
sr.b.e.9.10[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 10, 0, 0.9, ex, 0)
sr.b.e.9.10[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 10, 0.2, 0.9, ex, 0)
sr.b.e.9.10[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 10, 0.5, 0.9, ex, 0)






sr.ub.e.1.5 <- marix(0, 3, 6)
sr.ub.e.1.5[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 5, 0, 0.1, ex, 0.5)
sr.ub.e.1.5[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 5, 0.2, 0.1, ex, 0.5)
sr.ub.e.1.5[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 5, 0.5, 0.1, ex, 0.5)
sr.ub.e.1.5[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 5, 0, 0.1, ex, 0.5)
sr.ub.e.1.5[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 5, 0.2, 0.1, ex, 0.5)
sr.ub.e.1.5[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 5, 0.5, 0.1, ex, 0.5)
sr.ub.e.1.5[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 5, 0, 0.1, ex, 0.5)
sr.ub.e.1.5[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 5, 0.2, 0.1, ex, 0.5)
sr.ub.e.1.5[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 5, 0.5, 0.1, ex, 0.5)





sr.ub.e.5.5 <- marix(0, 3, 6)
sr.ub.e.5.5[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 5, 0, 0.5, ex, 0.5)
sr.ub.e.5.5[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 5, 0.2, 0.5, ex, 0.5)
sr.ub.e.5.5[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 5, 0.5, 0.5, ex, 0.5)
sr.ub.e.5.5[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 5, 0, 0.5, ex, 0.5)
sr.ub.e.5.5[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 5, 0.2, 0.5, ex, 0.5)
sr.ub.e.5.5[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 5, 0.5, 0.5, ex, 0.5)
sr.ub.e.5.5[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 5, 0, 0.5, ex, 0.5)
sr.ub.e.5.5[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 5, 0.2, 0.5, ex, 0.5)
sr.ub.e.5.5[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 5, 0.5, 0.5, ex, 0.5)


sr.ub.e.9.5 <- marix(0, 3, 6)
sr.ub.e.9.5[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 5, 0, 0.9, ex, 0.5)
sr.ub.e.9.5[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 5, 0.2, 0.9, ex, 0.5)
sr.ub.e.9.5[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 5, 0.5, 0.9, ex, 0.5)
sr.ub.e.9.5[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 5, 0, 0.9, ex, 0.5)
sr.ub.e.9.5[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 5, 0.2, 0.9, ex, 0.5)
sr.ub.e.9.5[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 5, 0.5, 0.9, ex, 0.5)
sr.ub.e.9.5[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 5, 0, 0.9, ex, 0.5)
sr.ub.e.9.5[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 5, 0.2, 0.9, ex, 0.5)
sr.ub.e.9.5[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 5, 0.5, 0.9, ex, 0.5)










sr.ub.e.1.10 <- marix(0, 3, 6)
sr.ub.e.1.10[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 10, 0, 0.1, ex, 0.5)
sr.ub.e.1.10[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 10, 0.2, 0.1, ex, 0.5)
sr.ub.e.1.10[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 10, 0.5, 0.1, ex, 0.5)
sr.ub.e.1.10[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 10, 0, 0.1, ex, 0.5)
sr.ub.e.1.10[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 10, 0.2, 0.1, ex, 0.5)
sr.ub.e.1.10[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 10, 0.5, 0.1, ex, 0.5)
sr.ub.e.1.10[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 10, 0, 0.1, ex, 0.5)
sr.ub.e.1.10[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 10, 0.2, 0.1, ex, 0.5)
sr.ub.e.1.10[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 10, 0.5, 0.1, ex, 0.5)





sr.ub.e.5.10 <- marix(0, 3, 6)
sr.ub.e.5.10[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 10, 0, 0.5, ex, 0.5)
sr.ub.e.5.10[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 10, 0.2, 0.5, ex, 0.5)
sr.ub.e.5.10[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 10, 0.5, 0.5, ex, 0.5)
sr.ub.e.5.10[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 10, 0, 0.5, ex, 0.5)
sr.ub.e.5.10[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 10, 0.2, 0.5, ex, 0.5)
sr.ub.e.5.10[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 10, 0.5, 0.5, ex, 0.5)
sr.ub.e.5.10[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 10, 0, 0.5, ex, 0.5)
sr.ub.e.5.10[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 10, 0.2, 0.5, ex, 0.5)
sr.ub.e.5.10[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 10, 0.5, 0.5, ex, 0.5)


sr.ub.e.9.10 <- marix(0, 3, 6)
sr.ub.e.9.10[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 10, 0, 0.9, ex, 0.5)
sr.ub.e.9.10[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 10, 0.2, 0.9, ex, 0.5)
sr.ub.e.9.10[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 10, 0.5, 0.9, ex, 0.5)
sr.ub.e.9.10[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 10, 0, 0.9, ex, 0.5)
sr.ub.e.9.10[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 10, 0.2, 0.9, ex, 0.5)
sr.ub.e.9.10[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 10, 0.5, 0.9, ex, 0.5)
sr.ub.e.9.10[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 10, 0, 0.9, ex, 0.5)
sr.ub.e.9.10[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 10, 0.2, 0.9, ex, 0.5)
sr.ub.e.9.10[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 10, 0.5, 0.9, ex, 0.5)













rs.b.a.11.2 <- marix(0, 3, 6)
rs.b.a.11.2[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 2, 0, c(0.1, 0.1), ar1, 0)
rs.b.a.11.2[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 2, 0.2, c(0.1, 0.1), ar1, 0)
rs.b.a.11.2[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 2, 0.5, c(0.1, 0.1), ar1, 0)
rs.b.a.11.2[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 2, 0, c(0.1, 0.1), ar1, 0)
rs.b.a.11.2[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 2, 0.2, c(0.1, 0.1), ar1, 0)
rs.b.a.11.2[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 2, 0.5, c(0.1, 0.1), ar1, 0)
rs.b.a.11.2[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 2, 0, c(0.1, 0.1), ar1, 0)
rs.b.a.11.2[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 2, 0.2, c(0.1, 0.1), ar1, 0)
rs.b.a.11.2[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 2, 0.5, c(0.1, 0.1), ar1, 0)





rs.b.a.55.2 <- marix(0, 3, 6)
rs.b.a.55.2[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 2, 0, c(0.5, 0.5), ar1, 0)
rs.b.a.55.2[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 2, 0.2, c(0.5, 0.5), ar1, 0)
rs.b.a.55.2[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 2, 0.5, c(0.5, 0.5), ar1, 0)
rs.b.a.55.2[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 2, 0, c(0.5, 0.5), ar1, 0)
rs.b.a.55.2[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 2, 0.2, c(0.5, 0.5), ar1, 0)
rs.b.a.55.2[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 2, 0.5, c(0.5, 0.5), ar1, 0)
rs.b.a.55.2[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 2, 0, c(0.5, 0.5), ar1, 0)
rs.b.a.55.2[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 2, 0.2, c(0.5, 0.5), ar1, 0)
rs.b.a.55.2[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 2, 0.5, c(0.5, 0.5), ar1, 0)







rs.b.a.19.2 <- marix(0, 3, 6)
rs.b.a.19.2[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 2, 0, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.2[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 2, 0.2, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.2[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 2, 0.5, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.2[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 2, 0, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.2[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 2, 0.2, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.2[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 2, 0.5, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.2[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 2, 0, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.2[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 2, 0.2, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.2[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 2, 0.5, c(-0.1, 0.9), ar1, 0)







rs.b.a.11.5 <- marix(0, 3, 6)
rs.b.a.11.5[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 5, 0, c(0.1, 0.1), ar1, 0)
rs.b.a.11.5[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 5, 0.2, c(0.1, 0.1), ar1, 0)
rs.b.a.11.5[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 5, 0.5, c(0.1, 0.1), ar1, 0)
rs.b.a.11.5[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 5, 0, c(0.1, 0.1), ar1, 0)
rs.b.a.11.5[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 5, 0.2, c(0.1, 0.1), ar1, 0)
rs.b.a.11.5[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 5, 0.5, c(0.1, 0.1), ar1, 0)
rs.b.a.11.5[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 5, 0, c(0.1, 0.1), ar1, 0)
rs.b.a.11.5[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 5, 0.2, c(0.1, 0.1), ar1, 0)
rs.b.a.11.5[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 5, 0.5, c(0.1, 0.1), ar1, 0)





rs.b.a.55.5 <- marix(0, 3, 6)
rs.b.a.55.5[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 5, 0, c(0.5, 0.5), ar1, 0)
rs.b.a.55.5[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 5, 0.2, c(0.5, 0.5), ar1, 0)
rs.b.a.55.5[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 5, 0.5, c(0.5, 0.5), ar1, 0)
rs.b.a.55.5[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 5, 0, c(0.5, 0.5), ar1, 0)
rs.b.a.55.5[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 5, 0.2, c(0.5, 0.5), ar1, 0)
rs.b.a.55.5[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 5, 0.5, c(0.5, 0.5), ar1, 0)
rs.b.a.55.5[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 5, 0, c(0.5, 0.5), ar1, 0)
rs.b.a.55.5[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 5, 0.2, c(0.5, 0.5), ar1, 0)
rs.b.a.55.5[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 5, 0.5, c(0.5, 0.5), ar1, 0)











rs.b.a.19.5 <- marix(0, 3, 6)
rs.b.a.19.5[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 5, 0, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.5[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 5, 0.2, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.5[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 5, 0.5, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.5[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 5, 0, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.5[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 5, 0.2, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.5[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 5, 0.5, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.5[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 5, 0, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.5[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 5, 0.2, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.5[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 5, 0.5, c(-0.1, 0.9), ar1, 0)








rs.ub.a.11.10 <- marix(0, 3, 6)
rs.ub.a.11.10[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 10, 0, c(0.1, 0.1), ar1, 0.5)
rs.ub.a.11.10[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 10, 0.2, c(0.1, 0.1), ar1, 0.5)
rs.ub.a.11.10[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 10, 0.5, c(0.1, 0.1), ar1, 0.5)
rs.ub.a.11.10[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 10, 0, c(0.1, 0.1), ar1, 0.5)
rs.ub.a.11.10[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 10, 0.2, c(0.1, 0.1), ar1, 0.5)
rs.ub.a.11.10[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 10, 0.5, c(0.1, 0.1), ar1, 0.5)
rs.ub.a.11.10[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 10, 0, c(0.1, 0.1), ar1, 0.5)
rs.ub.a.11.10[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 10, 0.2, c(0.1, 0.1), ar1, 0.5)
rs.ub.a.11.10[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 10, 0.5, c(0.1, 0.1), ar1, 0.5)










rs.ub.a.55.10 <- marix(0, 3, 6)
rs.ub.a.55.10[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 10, 0, c(0.5, 0.5), ar1, 0.5)
rs.ub.a.55.10[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 10, 0.2, c(0.5, 0.5), ar1, 0.5)
rs.ub.a.55.10[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 10, 0.5, c(0.5, 0.5), ar1, 0.5)
rs.ub.a.55.10[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 10, 0, c(0.5, 0.5), ar1, 0.5)
rs.ub.a.55.10[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 10, 0.2, c(0.5, 0.5), ar1, 0.5)
rs.ub.a.55.10[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 10, 0.5, c(0.5, 0.5), ar1, 0.5)
rs.ub.a.55.10[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 10, 0, c(0.5, 0.5), ar1, 0.5)
rs.ub.a.55.10[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 10, 0.2, c(0.5, 0.5), ar1, 0.5)
rs.ub.a.55.10[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 10, 0.5, c(0.5, 0.5), ar1, 0.5)







rs.ub.a.19.10 <- marix(0, 3, 6)
rs.ub.a.19.10[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 10, 0, c(-0.1, 0.9), ar1, 0.5)
rs.ub.a.19.10[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 10, 0.2, c(-0.1, 0.9), ar1, 0.5)
rs.ub.a.19.10[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 10, 0.5, c(-0.1, 0.9), ar1, 0.5)
rs.ub.a.19.10[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 10, 0, c(-0.1, 0.9), ar1, 0.5)
rs.ub.a.19.10[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 10, 0.2, c(-0.1, 0.9), ar1, 0.5)
rs.ub.a.19.10[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 10, 0.5, c(-0.1, 0.9), ar1, 0.5)
rs.ub.a.19.10[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 10, 0, c(-0.1, 0.9), ar1, 0.5)
rs.ub.a.19.10[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 10, 0.2, c(-0.1, 0.9), ar1, 0.5)
rs.ub.a.19.10[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 10, 0.5, c(-0.1, 0.9), ar1, 0.5)






%%%%%%%%%%%%%%%%%%%%%%%%%% Above are simulation for ranksum ar1, cluster level










%%%%%%%%%%%%%%%%%%%%%%%%% Below are simulation for ranksum ar1, subunit level





rs.b.a.11.2 <- marix(0, 3, 6)
rs.b.a.11.2[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 2, 0, c(0.1, 0.1), ar1, 0, FALSE)
rs.b.a.11.2[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 2, 0.2, c(0.1, 0.1), ar1, 0, FALSE)
rs.b.a.11.2[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 2, 0.5, c(0.1, 0.1), ar1, 0, FALSE)
rs.b.a.11.2[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 2, 0, c(0.1, 0.1), ar1, 0, FALSE)
rs.b.a.11.2[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 2, 0.2, c(0.1, 0.1), ar1, 0, FALSE)
rs.b.a.11.2[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 2, 0.5, c(0.1, 0.1), ar1, 0, FALSE)
rs.b.a.11.2[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 2, 0, c(0.1, 0.1), ar1, 0, FALSE)
rs.b.a.11.2[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 2, 0.2, c(0.1, 0.1), ar1, 0, FALSE)
rs.b.a.11.2[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 2, 0.5, c(0.1, 0.1), ar1, 0, FALSE)





rs.b.a.55.2 <- marix(0, 3, 6)
rs.b.a.55.2[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 2, 0, c(0.5, 0.5), ar1, 0, FALSE)
rs.b.a.55.2[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 2, 0.2, c(0.5, 0.5), ar1, 0, FALSE)
rs.b.a.55.2[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 2, 0.5, c(0.5, 0.5), ar1, 0, FALSE)
rs.b.a.55.2[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 2, 0, c(0.5, 0.5), ar1, 0, FALSE)
rs.b.a.55.2[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 2, 0.2, c(0.5, 0.5), ar1, 0, FALSE)
rs.b.a.55.2[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 2, 0.5, c(0.5, 0.5), ar1, 0, FALSE)
rs.b.a.55.2[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 2, 0, c(0.5, 0.5), ar1, 0, FALSE)
rs.b.a.55.2[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 2, 0.2, c(0.5, 0.5), ar1, 0, FALSE)
rs.b.a.55.2[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 2, 0.5, c(0.5, 0.5), ar1, 0, FALSE)







rs.b.a.19.2 <- marix(0, 3, 6)
rs.b.a.19.2[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 2, 0, c(-0.1, 0.9), ar1, 0, FALSE)
rs.b.a.19.2[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 2, 0.2, c(-0.1, 0.9), ar1, 0, FALSE)
rs.b.a.19.2[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 2, 0.5, c(-0.1, 0.9), ar1, 0, FALSE)
rs.b.a.19.2[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 2, 0, c(-0.1, 0.9), ar1, 0, FALSE)
rs.b.a.19.2[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 2, 0.2, c(-0.1, 0.9), ar1, 0, FALSE)
rs.b.a.19.2[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 2, 0.5, c(-0.1, 0.9), ar1, 0, FALSE)
rs.b.a.19.2[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 2, 0, c(-0.1, 0.9), ar1, 0, FALSE)
rs.b.a.19.2[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 2, 0.2, c(-0.1, 0.9), ar1, 0, FALSE)
rs.b.a.19.2[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 2, 0.5, c(-0.1, 0.9), ar1, 0, FALSE)







rs.b.a.11.5 <- marix(0, 3, 6)
rs.b.a.11.5[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 5, 0, c(0.1, 0.1), ar1, 0, FALSE)
rs.b.a.11.5[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 5, 0.2, c(0.1, 0.1), ar1, 0, FALSE)
rs.b.a.11.5[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 5, 0.5, c(0.1, 0.1), ar1, 0, FALSE)
rs.b.a.11.5[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 5, 0, c(0.1, 0.1), ar1, 0, FALSE)
rs.b.a.11.5[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 5, 0.2, c(0.1, 0.1), ar1, 0, FALSE)
rs.b.a.11.5[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 5, 0.5, c(0.1, 0.1), ar1, 0, FALSE)
rs.b.a.11.5[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 5, 0, c(0.1, 0.1), ar1, 0, FALSE)
rs.b.a.11.5[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 5, 0.2, c(0.1, 0.1), ar1, 0, FALSE)
rs.b.a.11.5[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 5, 0.5, c(0.1, 0.1), ar1, 0, FALSE)





rs.b.a.55.5 <- marix(0, 3, 6)
rs.b.a.55.5[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 5, 0, c(0.5, 0.5), ar1, 0, FALSE)
rs.b.a.55.5[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 5, 0.2, c(0.5, 0.5), ar1, 0, FALSE)
rs.b.a.55.5[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 5, 0.5, c(0.5, 0.5), ar1, 0, FALSE)
rs.b.a.55.5[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 5, 0, c(0.5, 0.5), ar1, 0, FALSE)
rs.b.a.55.5[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 5, 0.2, c(0.5, 0.5), ar1, 0, FALSE)
rs.b.a.55.5[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 5, 0.5, c(0.5, 0.5), ar1, 0, FALSE)
rs.b.a.55.5[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 5, 0, c(0.5, 0.5), ar1, 0, FALSE)
rs.b.a.55.5[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 5, 0.2, c(0.5, 0.5), ar1, 0, FALSE)
rs.b.a.55.5[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 5, 0.5, c(0.5, 0.5), ar1, 0, FALSE)











rs.b.a.19.5 <- marix(0, 3, 6)
rs.b.a.19.5[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 5, 0, c(-0.1, 0.9), ar1, 0, FALSE)
rs.b.a.19.5[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 5, 0.2, c(-0.1, 0.9), ar1, 0, FALSE)
rs.b.a.19.5[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 5, 0.5, c(-0.1, 0.9), ar1, 0, FALSE)
rs.b.a.19.5[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 5, 0, c(-0.1, 0.9), ar1, 0, FALSE)
rs.b.a.19.5[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 5, 0.2, c(-0.1, 0.9), ar1, 0, FALSE)
rs.b.a.19.5[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 5, 0.5, c(-0.1, 0.9), ar1, 0, FALSE)
rs.b.a.19.5[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 5, 0, c(-0.1, 0.9), ar1, 0, FALSE)
rs.b.a.19.5[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 5, 0.2, c(-0.1, 0.9), ar1, 0, FALSE)
rs.b.a.19.5[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 5, 0.5, c(-0.1, 0.9), ar1, 0, FALSE)








rs.ub.a.11.10 <- marix(0, 3, 6)
rs.ub.a.11.10[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 10, 0, c(0.1, 0.1), ar1, 0.5, FALSE)
rs.ub.a.11.10[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 10, 0.2, c(0.1, 0.1), ar1, 0.5, FALSE)
rs.ub.a.11.10[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 10, 0.5, c(0.1, 0.1), ar1, 0.5, FALSE)
rs.ub.a.11.10[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 10, 0, c(0.1, 0.1), ar1, 0.5, FALSE)
rs.ub.a.11.10[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 10, 0.2, c(0.1, 0.1), ar1, 0.5, FALSE)
rs.ub.a.11.10[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 10, 0.5, c(0.1, 0.1), ar1, 0.5, FALSE)
rs.ub.a.11.10[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 10, 0, c(0.1, 0.1), ar1, 0.5, FALSE)
rs.ub.a.11.10[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 10, 0.2, c(0.1, 0.1), ar1, 0.5, FALSE)
rs.ub.a.11.10[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 10, 0.5, c(0.1, 0.1), ar1, 0.5, FALSE)










rs.ub.a.55.10 <- marix(0, 3, 6)
rs.ub.a.55.10[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 10, 0, c(0.5, 0.5), ar1, 0.5, FALSE)
rs.ub.a.55.10[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 10, 0.2, c(0.5, 0.5), ar1, 0.5, FALSE)
rs.ub.a.55.10[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 10, 0.5, c(0.5, 0.5), ar1, 0.5, FALSE)
rs.ub.a.55.10[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 10, 0, c(0.5, 0.5), ar1, 0.5, FALSE)
rs.ub.a.55.10[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 10, 0.2, c(0.5, 0.5), ar1, 0.5, FALSE)
rs.ub.a.55.10[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 10, 0.5, c(0.5, 0.5), ar1, 0.5, FALSE)
rs.ub.a.55.10[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 10, 0, c(0.5, 0.5), ar1, 0.5, FALSE)
rs.ub.a.55.10[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 10, 0.2, c(0.5, 0.5), ar1, 0.5, FALSE)
rs.ub.a.55.10[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 10, 0.5, c(0.5, 0.5), ar1, 0.5, FALSE)







rs.ub.a.19.10 <- marix(0, 3, 6)
rs.ub.a.19.10[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 10, 0, c(-0.1, 0.9), ar1, 0.5, FALSE)
rs.ub.a.19.10[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 10, 0.2, c(-0.1, 0.9), ar1, 0.5, FALSE)
rs.ub.a.19.10[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 10, 0.5, c(-0.1, 0.9), ar1, 0.5, FALSE)
rs.ub.a.19.10[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 10, 0, c(-0.1, 0.9), ar1, 0.5, FALSE)
rs.ub.a.19.10[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 10, 0.2, c(-0.1, 0.9), ar1, 0.5, FALSE)
rs.ub.a.19.10[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 10, 0.5, c(-0.1, 0.9), ar1, 0.5, FALSE)
rs.ub.a.19.10[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 10, 0, c(-0.1, 0.9), ar1, 0.5, FALSE)
rs.ub.a.19.10[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 10, 0.2, c(-0.1, 0.9), ar1, 0.5, FALSE)
rs.ub.a.19.10[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 10, 0.5, c(-0.1, 0.9), ar1, 0.5, FALSE)







sr.b.a.1.2 <- marix(0, 3, 6)
sr.b.a.1.2[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 2, 0, 0.1, ar1, 0)
sr.b.a.1.2[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 2, 0.2, 0.1, ar1, 0)
sr.b.a.1.2[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 2, 0.5, 0.1, ar1, 0)
sr.b.a.1.2[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 2, 0, 0.1, ar1, 0)
sr.b.a.1.2[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 2, 0.2, 0.1, ar1, 0)
sr.b.a.1.2[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 2, 0.5, 0.1, ar1, 0)
sr.b.a.1.2[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 2, 0, 0.1, ar1, 0)
sr.b.a.1.2[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 2, 0.2, 0.1, ar1, 0)
sr.b.a.1.2[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 2, 0.5, 0.1, ar1, 0)





sr.b.a.5.2 <- marix(0, 3, 6)
sr.b.a.5.2[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 2, 0, 0.5, ar1, 0)
sr.b.a.5.2[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 2, 0.2, 0.5, ar1, 0)
sr.b.a.5.2[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 2, 0.5, 0.5, ar1, 0)
sr.b.a.5.2[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 2, 0, 0.5, ar1, 0)
sr.b.a.5.2[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 2, 0.2, 0.5, ar1, 0)
sr.b.a.5.2[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 2, 0.5, 0.5, ar1, 0)
sr.b.a.5.2[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 2, 0, 0.5, ar1, 0)
sr.b.a.5.2[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 2, 0.2, 0.5, ar1, 0)
sr.b.a.5.2[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 2, 0.5, 0.5, ar1, 0)


sr.b.a.9.2 <- marix(0, 3, 6)
sr.b.a.9.2[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 2, 0, 0.9, ar1, 0)
sr.b.a.9.2[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 2, 0.2, 0.9, ar1, 0)
sr.b.a.9.2[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 2, 0.5, 0.9, ar1, 0)
sr.b.a.9.2[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 2, 0, 0.9, ar1, 0)
sr.b.a.9.2[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 2, 0.2, 0.9, ar1, 0)
sr.b.a.9.2[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 2, 0.5, 0.9, ar1, 0)
sr.b.a.9.2[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 2, 0, 0.9, ar1, 0)
sr.b.a.9.2[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 2, 0.2, 0.9, ar1, 0)
sr.b.a.9.2[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 2, 0.5, 0.9, ar1, 0)








sr.b.a.1.10 <- marix(0, 3, 6)
sr.b.a.1.10[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 10, 0, 0.1, ar1, 0)
sr.b.a.1.10[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 10, 0.2, 0.1, ar1, 0)
sr.b.a.1.10[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 10, 0.5, 0.1, ar1, 0)
sr.b.a.1.10[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 10, 0, 0.1, ar1, 0)
sr.b.a.1.10[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 10, 0.2, 0.1, ar1, 0)
sr.b.a.1.10[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 10, 0.5, 0.1, ar1, 0)
sr.b.a.1.10[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 10, 0, 0.1, ar1, 0)
sr.b.a.1.10[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 10, 0.2, 0.1, ar1, 0)
sr.b.a.1.10[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 10, 0.5, 0.1, ar1, 0)





sr.b.a.5.10 <- marix(0, 3, 6)
sr.b.a.5.10[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 10, 0, 0.5, ar1, 0)
sr.b.a.5.10[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 10, 0.2, 0.5, ar1, 0)
sr.b.a.5.10[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 10, 0.5, 0.5, ar1, 0)
sr.b.a.5.10[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 10, 0, 0.5, ar1, 0)
sr.b.a.5.10[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 10, 0.2, 0.5, ar1, 0)
sr.b.a.5.10[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 10, 0.5, 0.5, ar1, 0)
sr.b.a.5.10[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 10, 0, 0.5, ar1, 0)
sr.b.a.5.10[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 10, 0.2, 0.5, ar1, 0)
sr.b.a.5.10[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 10, 0.5, 0.5, ar1, 0)


sr.b.a.9.10 <- marix(0, 3, 6)
sr.b.a.9.10[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 10, 0, 0.9, ar1, 0)
sr.b.a.9.10[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 10, 0.2, 0.9, ar1, 0)
sr.b.a.9.10[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 10, 0.5, 0.9, ar1, 0)
sr.b.a.9.10[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 10, 0, 0.9, ar1, 0)
sr.b.a.9.10[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 10, 0.2, 0.9, ar1, 0)
sr.b.a.9.10[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 10, 0.5, 0.9, ar1, 0)
sr.b.a.9.10[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 10, 0, 0.9, ar1, 0)
sr.b.a.9.10[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 10, 0.2, 0.9, ar1, 0)
sr.b.a.9.10[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 10, 0.5, 0.9, ar1, 0)






sr.ub.a.1.5 <- marix(0, 3, 6)
sr.ub.a.1.5[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 5, 0, 0.1, ar1, 0.5)
sr.ub.a.1.5[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 5, 0.2, 0.1, ar1, 0.5)
sr.ub.a.1.5[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 5, 0.5, 0.1, ar1, 0.5)
sr.ub.a.1.5[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 5, 0, 0.1, ar1, 0.5)
sr.ub.a.1.5[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 5, 0.2, 0.1, ar1, 0.5)
sr.ub.a.1.5[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 5, 0.5, 0.1, ar1, 0.5)
sr.ub.a.1.5[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 5, 0, 0.1, ar1, 0.5)
sr.ub.a.1.5[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 5, 0.2, 0.1, ar1, 0.5)
sr.ub.a.1.5[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 5, 0.5, 0.1, ar1, 0.5)





sr.ub.a.5.5 <- marix(0, 3, 6)
sr.ub.a.5.5[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 5, 0, 0.5, ar1, 0.5)
sr.ub.a.5.5[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 5, 0.2, 0.5, ar1, 0.5)
sr.ub.a.5.5[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 5, 0.5, 0.5, ar1, 0.5)
sr.ub.a.5.5[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 5, 0, 0.5, ar1, 0.5)
sr.ub.a.5.5[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 5, 0.2, 0.5, ar1, 0.5)
sr.ub.a.5.5[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 5, 0.5, 0.5, ar1, 0.5)
sr.ub.a.5.5[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 5, 0, 0.5, ar1, 0.5)
sr.ub.a.5.5[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 5, 0.2, 0.5, ar1, 0.5)
sr.ub.a.5.5[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 5, 0.5, 0.5, ar1, 0.5)


sr.ub.a.9.5 <- marix(0, 3, 6)
sr.ub.a.9.5[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 5, 0, 0.9, ar1, 0.5)
sr.ub.a.9.5[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 5, 0.2, 0.9, ar1, 0.5)
sr.ub.a.9.5[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 5, 0.5, 0.9, ar1, 0.5)
sr.ub.a.9.5[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 5, 0, 0.9, ar1, 0.5)
sr.ub.a.9.5[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 5, 0.2, 0.9, ar1, 0.5)
sr.ub.a.9.5[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 5, 0.5, 0.9, ar1, 0.5)
sr.ub.a.9.5[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 5, 0, 0.9, ar1, 0.5)
sr.ub.a.9.5[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 5, 0.2, 0.9, ar1, 0.5)
sr.ub.a.9.5[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 5, 0.5, 0.9, ar1, 0.5)










sr.ub.a.1.10 <- marix(0, 3, 6)
sr.ub.a.1.10[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 10, 0, 0.1, ar1, 0.5)
sr.ub.a.1.10[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 10, 0.2, 0.1, ar1, 0.5)
sr.ub.a.1.10[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 10, 0.5, 0.1, ar1, 0.5)
sr.ub.a.1.10[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 10, 0, 0.1, ar1, 0.5)
sr.ub.a.1.10[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 10, 0.2, 0.1, ar1, 0.5)
sr.ub.a.1.10[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 10, 0.5, 0.1, ar1, 0.5)
sr.ub.a.1.10[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 10, 0, 0.1, ar1, 0.5)
sr.ub.a.1.10[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 10, 0.2, 0.1, ar1, 0.5)
sr.ub.a.1.10[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 10, 0.5, 0.1, ar1, 0.5)





sr.ub.a.5.10 <- marix(0, 3, 6)
sr.ub.a.5.10[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 10, 0, 0.5, ar1, 0.5)
sr.ub.a.5.10[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 10, 0.2, 0.5, ar1, 0.5)
sr.ub.a.5.10[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 10, 0.5, 0.5, ar1, 0.5)
sr.ub.a.5.10[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 10, 0, 0.5, ar1, 0.5)
sr.ub.a.5.10[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 10, 0.2, 0.5, ar1, 0.5)
sr.ub.a.5.10[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 10, 0.5, 0.5, ar1, 0.5)
sr.ub.a.5.10[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 10, 0, 0.5, ar1, 0.5)
sr.ub.a.5.10[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 10, 0.2, 0.5, ar1, 0.5)
sr.ub.a.5.10[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 10, 0.5, 0.5, ar1, 0.5)


sr.ub.a.9.10 <- marix(0, 3, 6)
sr.ub.a.9.10[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 10, 0, 0.9, ar1, 0.5)
sr.ub.a.9.10[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 10, 0.2, 0.9, ar1, 0.5)
sr.ub.a.9.10[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 10, 0.5, 0.9, ar1, 0.5)
sr.ub.a.9.10[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 10, 0, 0.9, ar1, 0.5)
sr.ub.a.9.10[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 10, 0.2, 0.9, ar1, 0.5)
sr.ub.a.9.10[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 10, 0.5, 0.9, ar1, 0.5)
sr.ub.a.9.10[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 10, 0, 0.9, ar1, 0.5)
sr.ub.a.9.10[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 10, 0.2, 0.9, ar1, 0.5)
sr.ub.a.9.10[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 10, 0.5, 0.9, ar1, 0.5)





sum.tab.rs.e <- rbind(rs.b.e.11.2, rs.b.e.55.2, rs.b.e.19.2,
                    rs.b.e.11.5, rs.b.e.55.5, rs.b.e.19.5,
                    rs.ub.e.11.10, rs.ub.e.55.10, rs.ub.e.19.10,
                    rs.b.e.11.2.sub, rs.b.e.55.2.sub, rs.b.e.19.2.sub,
                    rs.b.e.11.5.sub, rs.b.e.55.5.sub, rs.b.e.19.5.sub,
                    rs.ub.e.11.10.sub, rs.ub.e.55.10.sub, rs.ub.e.19.10.sub)
sum.tab.rs.e <- sum.tab.rs.e * 100





sum.tab.rs.a <- rbind(rs.b.a.11.2, rs.b.a.55.2, rs.b.a.19.2,
                    rs.b.a.11.5, rs.b.a.55.5, rs.b.a.19.5,
                    rs.ub.a.11.10, rs.ub.a.55.10, rs.ub.a.19.10,
                    rs.b.a.11.2.sub, rs.b.a.55.2.sub, rs.b.a.19.2.sub,
                    rs.b.a.11.5.sub, rs.b.a.55.5.sub, rs.b.a.19.5.sub,
                    rs.ub.a.11.10.sub, rs.ub.a.55.10.sub, rs.ub.a.19.10.sub)
sum.tab.rs.a <- sum.tab.rs.a * 100






sum.tab.sr.e <- rbind(sr.b.e.1.2, sr.b.e.5.2, sr.b.e.9.2,
                      sr.b.e.1.10, sr.b.e.5.10, sr.b.e.9.10,
                      sr.ub.e.1.5, sr.ub.e.5.5, sr.ub.e.9.5,
                      sr.ub.e.1.10, sr.ub.e.5.10, sr.ub.e.9.10) * 100


sum.tab.sr.a <- rbind(sr.b.a.1.2, sr.b.a.5.2, sr.b.a.9.2,
                      sr.b.a.1.10, sr.b.a.5.10, sr.b.a.9.10,
                      sr.ub.a.1.5, sr.ub.a.5.5, sr.ub.a.9.5,
                      sr.ub.a.1.10, sr.ub.a.5.10, sr.ub.a.9.10) * 100




