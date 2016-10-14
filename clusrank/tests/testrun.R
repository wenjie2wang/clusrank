

library(clusrank)
data(crsd)
clusWilcox.test(z, cluster = id, data = crsd, paired = TRUE, method = "rgl")
clusWilcox.test(z, cluster = id, data = crsd, paired = TRUE, method = "rgl", alternative = "greater", mu = 1)
clusWilcox.test(z ~ cluster(id), data = crsd, paired = TRUE, method = "ds")
clusWilcox.test(z ~ cluster(id), data = crsd, paired = TRUE, method = "rgl", exact = TRUE)
data(crsdUnb)
clusWilcox.test(z, cluster = id, data = crsdUnb, paired = TRUE, method = "rgl")
clusWilcox.test(z, cluster = id, data = crsdUnb, paired = TRUE, method = "ds")
data(crd)
clusWilcox.test(z ~ group + cluster(id), data = crd)
crd1 <- crd[c(1:20, 141:160), ]
clusWilcox.test(z ~ group + cluster(id), data = crd1, method = "rgl", exact = TRUE)
data(crdStr)
clusWilcox.test(z ~ group + cluster(id) + stratum(stratum), data = crdStr)
clusWilcox.test(z ~ group + cluster(id) + stratum(stratum), data = crdStr, method = "ds")



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

set.seed(1234)
dat.cl <- datgen.sum(10, 3, 0, c(.9, .9), ex, 0, TRUE)
clusWilcox.test(x ~ grp + cluster(cid), dat.cl, method = "rgl", exact = TRUE)
clusWilcox.test(x, group =  grp, cluster = cid, data = dat.cl, method = "ds")

dat.cl$strat <- rep(rep(1:2, each = 15), 2)
clusWilcox.test(x ~ grp + cluster(cid) + stratum(strat), dat = dat.cl,
                method = "rgl")
dat.cl$grp <- rep(1:4, each = 15)
clusWilcox.test(x ~ grp + cluster(cid), dat = dat.cl, method = "ds")
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

dat.sgn <- datgen.sgn(10, 3, cor = ex, rho = 0.5)
clusWilcox.test(x ~ cluster(cid),  dat.sgn, paired = TRUE, method = "rgl")

clusWilcox.test(x ~ cluster(cid),  dat.sgn, paired = TRUE, method = "ds")
data(amd)
clusWilcox.test(CARMS ~ Variant + cluster(ID), data = amd,
                subset = CARMS %in% c(1, 2, 3, 4), method = "rgl")
clusWilcox.test(CARMS ~ Variant + cluster(ID), data = amd,
                subset = CARMS %in% c(1, 2, 3, 4), method = "ds")
clusWilcox.test(CARMS ~ Variant + cluster(ID), data = amd, method = "rgl",
                subset = CARMS %in% c(1, 2, 3, 5))
clusWilcox.test(CARMS ~ Variant + cluster(ID), data = amd, method = "ds",
                subset = CARMS %in% c(1, 2, 3, 5))
