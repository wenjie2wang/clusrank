### R code from vignette source 'clusrank.Rnw'

###################################################
### code chunk number 1: knitr
###################################################
library(knitr)
opts_chunk$set(cache=TRUE, autodep=TRUE)


###################################################
### code chunk number 2: args-cc
###################################################
library(clusrank)
args(clusWilcox.test)


###################################################
### code chunk number 3: args-cd
###################################################
args(getS3method("clusWilcox.test", "default"))


###################################################
### code chunk number 4: formula
###################################################
args(getS3method("clusWilcox.test", "formula"))


###################################################
### code chunk number 5: corstr
###################################################
ex  <- function(dim, rho) {
    diag(1 - rho, dim) + matrix(rho, dim, dim)
}
ar1 <- function(dim, rho) {
    rho ^ outer(1:dim, 1:dim, function(x, y) abs(x - y))
}


###################################################
### code chunk number 6: datgen
###################################################
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


###################################################
### code chunk number 7: seed
###################################################
set.seed(1234)


###################################################
### code chunk number 8: dat-clus
###################################################
dat.cl <- datgen.sum(10, 3, 0, c(.9, .9), ex, 0, TRUE)


###################################################
### code chunk number 9: dat.cl
###################################################
cbind(head(dat.cl, 6), "head / tail" = "     ", tail(dat.cl, 6))


###################################################
### code chunk number 10: clus-rs-rgl-asymp
###################################################
clusWilcox.test(x ~ grp + cluster(cid), dat.cl, method = "rgl")


###################################################
### code chunk number 11: clus-rs-rgl-exact
###################################################
clusWilcox.test(x ~ grp + cluster(cid), dat.cl, method = "rgl",
                exact = TRUE, B = 0)


###################################################
### code chunk number 12: clus-rs-ds
###################################################
clusWilcox.test(x, group =  grp, cluster = cid, data = dat.cl, method = "ds")


###################################################
### code chunk number 13: clus-rs-rgl-str
###################################################
dat.cl$strat <- rep(rep(1:2, each = 15), 2)
cbind(head(dat.cl, 6), "head / tail" = "     ", tail(dat.cl, 6))
clusWilcox.test(x ~ grp + cluster(cid) + stratum(strat), dat = dat.cl,
                method = "rgl")


###################################################
### code chunk number 14: clus-rs-ds-mgrp
###################################################
dat.cl$grp <- rep(1:4, each = 15)
cbind(head(dat.cl, 6), "head / tail" = "     ", tail(dat.cl, 6))
clusWilcox.test(x ~ grp + cluster(cid), dat = dat.cl, method = "ds")


###################################################
### code chunk number 15: datagen-sr
###################################################
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


###################################################
### code chunk number 16: dat-sr
###################################################
dat.sgn <- datgen.sgn(10, 3, cor = ex, rho = 0.5)
cbind(head(dat.cl, 6), "head / tail" = "     ", tail(dat.cl, 6))


###################################################
### code chunk number 17: sr-rgl
###################################################
clusWilcox.test(x ~ cluster(cid),  dat.sgn, paired = TRUE, method = "rgl")


###################################################
### code chunk number 18: sr-ds
###################################################
clusWilcox.test(x ~ cluster(cid),  dat.sgn, paired = TRUE, method = "ds")


###################################################
### code chunk number 19: amd1234-both
###################################################
data(amd)
clusWilcox.test(CARMS ~ Variant + cluster(ID), data = amd,
                subset = CARMS %in% c(1, 2, 3, 4), method = "rgl")
clusWilcox.test(CARMS ~ Variant + cluster(ID), data = amd,
                subset = CARMS %in% c(1, 2, 3, 4), method = "ds")


###################################################
### code chunk number 20: amd1234-rgl
###################################################
clusWilcox.test(CARMS ~ Variant + cluster(ID) + stratum(Agesex), data = amd,
                subset = CARMS %in% c(1, 2, 3, 4))


###################################################
### code chunk number 21: amd12355-both
###################################################
clusWilcox.test(CARMS ~ Variant + cluster(ID), data = amd, method = "rgl",
                subset = CARMS %in% c(1, 2, 3, 5))
clusWilcox.test(CARMS ~ Variant + cluster(ID), data = amd, method = "ds",
                subset = CARMS %in% c(1, 2, 3, 5))


###################################################
### code chunk number 22: amd12355-rgl
###################################################
clusWilcox.test(CARMS ~ Variant + cluster(ID) + stratum(Agesex), data = amd,
                subset = CARMS %in% c(1, 2, 3, 5), method = "rgl")


###################################################
### code chunk number 23: simpower
###################################################
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


###################################################
### code chunk number 24: sim-clus-rs
###################################################
simpower(1000, 0.05, FALSE, 20, 3, 0.0, c(0.5, 0.5), ex, 0., clusgrp = TRUE)


###################################################
### code chunk number 25: sim-sr
###################################################
simpower(1000, 0.05, TRUE,  20, 3, 0.0, 0.5, ar1, 0.)


