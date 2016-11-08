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
                perm = TRUE, B = 1)


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


rs.b.e.11.2 <- matrix(0, 3, 6)
rs.b.e.11.2[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 2, 0, c(0.1, 0.1), ex, 0)
rs.b.e.11.2[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 2, 0.2, c(0.1, 0.1), ex, 0)
rs.b.e.11.2[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 2, 0.5, c(0.1, 0.1), ex, 0)
rs.b.e.11.2[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 2, 0, c(0.1, 0.1), ex, 0)
rs.b.e.11.2[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 2, 0.2, c(0.1, 0.1), ex, 0)
rs.b.e.11.2[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 2, 0.5, c(0.1, 0.1), ex, 0)
rs.b.e.11.2[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 2, 0, c(0.1, 0.1), ex, 0)
rs.b.e.11.2[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 2, 0.2, c(0.1, 0.1), ex, 0)
rs.b.e.11.2[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 2, 0.5, c(0.1, 0.1), ex, 0)





rs.b.e.55.2 <- matrix(0, 3, 6)
rs.b.e.55.2[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 2, 0, c(0.5, 0.5), ex, 0)
rs.b.e.55.2[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 2, 0.2, c(0.5, 0.5), ex, 0)
rs.b.e.55.2[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 2, 0.5, c(0.5, 0.5), ex, 0)
rs.b.e.55.2[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 2, 0, c(0.5, 0.5), ex, 0)
rs.b.e.55.2[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 2, 0.2, c(0.5, 0.5), ex, 0)
rs.b.e.55.2[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 2, 0.5, c(0.5, 0.5), ex, 0)
rs.b.e.55.2[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 2, 0, c(0.5, 0.5), ex, 0)
rs.b.e.55.2[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 2, 0.2, c(0.5, 0.5), ex, 0)
rs.b.e.55.2[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 2, 0.5, c(0.5, 0.5), ex, 0)







rs.b.e.19.2 <- matrix(0, 3, 6)
rs.b.e.19.2[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 2, 0, c(-0.1, 0.9), ex, 0)
rs.b.e.19.2[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 2, 0.2, c(-0.1, 0.9), ex, 0)
rs.b.e.19.2[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 2, 0.5, c(-0.1, 0.9), ex, 0)
rs.b.e.19.2[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 2, 0, c(-0.1, 0.9), ex, 0)
rs.b.e.19.2[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 2, 0.2, c(-0.1, 0.9), ex, 0)
rs.b.e.19.2[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 2, 0.5, c(-0.1, 0.9), ex, 0)
rs.b.e.19.2[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 2, 0, c(-0.1, 0.9), ex, 0)
rs.b.e.19.2[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 2, 0.2, c(-0.1, 0.9), ex, 0)
rs.b.e.19.2[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 2, 0.5, c(-0.1, 0.9), ex, 0)







rs.b.e.11.5 <- matrix(0, 3, 6)
rs.b.e.11.5[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 5, 0, c(0.1, 0.1), ex, 0)
rs.b.e.11.5[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 5, 0.2, c(0.1, 0.1), ex, 0)
rs.b.e.11.5[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 5, 0.5, c(0.1, 0.1), ex, 0)
rs.b.e.11.5[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 5, 0, c(0.1, 0.1), ex, 0)
rs.b.e.11.5[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 5, 0.2, c(0.1, 0.1), ex, 0)
rs.b.e.11.5[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 5, 0.5, c(0.1, 0.1), ex, 0)
rs.b.e.11.5[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 5, 0, c(0.1, 0.1), ex, 0)
rs.b.e.11.5[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 5, 0.2, c(0.1, 0.1), ex, 0)
rs.b.e.11.5[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 5, 0.5, c(0.1, 0.1), ex, 0)





rs.b.e.55.5 <- matrix(0, 3, 6)
rs.b.e.55.5[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 5, 0, c(0.5, 0.5), ex, 0)
rs.b.e.55.5[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 5, 0.2, c(0.5, 0.5), ex, 0)
rs.b.e.55.5[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 5, 0.5, c(0.5, 0.5), ex, 0)
rs.b.e.55.5[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 5, 0, c(0.5, 0.5), ex, 0)
rs.b.e.55.5[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 5, 0.2, c(0.5, 0.5), ex, 0)
rs.b.e.55.5[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 5, 0.5, c(0.5, 0.5), ex, 0)
rs.b.e.55.5[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 5, 0, c(0.5, 0.5), ex, 0)
rs.b.e.55.5[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 5, 0.2, c(0.5, 0.5), ex, 0)
rs.b.e.55.5[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 5, 0.5, c(0.5, 0.5), ex, 0)











rs.b.e.19.5 <- matrix(0, 3, 6)
rs.b.e.19.5[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 5, 0, c(-0.1, 0.9), ex, 0)
rs.b.e.19.5[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 5, 0.2, c(-0.1, 0.9), ex, 0)
rs.b.e.19.5[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 5, 0.5, c(-0.1, 0.9), ex, 0)
rs.b.e.19.5[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 5, 0, c(-0.1, 0.9), ex, 0)
rs.b.e.19.5[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 5, 0.2, c(-0.1, 0.9), ex, 0)
rs.b.e.19.5[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 5, 0.5, c(-0.1, 0.9), ex, 0)
rs.b.e.19.5[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 5, 0, c(-0.1, 0.9), ex, 0)
rs.b.e.19.5[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 5, 0.2, c(-0.1, 0.9), ex, 0)
rs.b.e.19.5[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 5, 0.5, c(-0.1, 0.9), ex, 0)








rs.ub.e.11.10 <- matrix(0, 3, 6)
rs.ub.e.11.10[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 10, 0, c(0.1, 0.1), ex, 0.5)
rs.ub.e.11.10[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 10, 0.2, c(0.1, 0.1), ex, 0.5)
rs.ub.e.11.10[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 10, 0.5, c(0.1, 0.1), ex, 0.5)
rs.ub.e.11.10[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 10, 0, c(0.1, 0.1), ex, 0.5)
rs.ub.e.11.10[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 10, 0.2, c(0.1, 0.1), ex, 0.5)
rs.ub.e.11.10[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 10, 0.5, c(0.1, 0.1), ex, 0.5)
rs.ub.e.11.10[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 10, 0, c(0.1, 0.1), ex, 0.5)
rs.ub.e.11.10[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 10, 0.2, c(0.1, 0.1), ex, 0.5)
rs.ub.e.11.10[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 10, 0.5, c(0.1, 0.1), ex, 0.5)










rs.ub.e.55.10 <- matrix(0, 3, 6)
rs.ub.e.55.10[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 10, 0, c(0.5, 0.5), ex, 0.5)
rs.ub.e.55.10[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 10, 0.2, c(0.5, 0.5), ex, 0.5)
rs.ub.e.55.10[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 10, 0.5, c(0.5, 0.5), ex, 0.5)
rs.ub.e.55.10[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 10, 0, c(0.5, 0.5), ex, 0.5)
rs.ub.e.55.10[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 10, 0.2, c(0.5, 0.5), ex, 0.5)
rs.ub.e.55.10[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 10, 0.5, c(0.5, 0.5), ex, 0.5)
rs.ub.e.55.10[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 10, 0, c(0.5, 0.5), ex, 0.5)
rs.ub.e.55.10[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 10, 0.2, c(0.5, 0.5), ex, 0.5)
rs.ub.e.55.10[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 10, 0.5, c(0.5, 0.5), ex, 0.5)







rs.ub.e.19.10 <- matrix(0, 3, 6)
rs.ub.e.19.10[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 10, 0, c(-0.1, 0.9), ex, 0.5)
rs.ub.e.19.10[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 10, 0.2, c(-0.1, 0.9), ex, 0.5)
rs.ub.e.19.10[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 10, 0.5, c(-0.1, 0.9), ex, 0.5)
rs.ub.e.19.10[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 10, 0, c(-0.1, 0.9), ex, 0.5)
rs.ub.e.19.10[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 10, 0.2, c(-0.1, 0.9), ex, 0.5)
rs.ub.e.19.10[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 10, 0.5, c(-0.1, 0.9), ex, 0.5)
rs.ub.e.19.10[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 10, 0, c(-0.1, 0.9), ex, 0.5)
rs.ub.e.19.10[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 10, 0.2, c(-0.1, 0.9), ex, 0.5)
rs.ub.e.19.10[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 10, 0.5, c(-0.1, 0.9), ex, 0.5)





rs.b.e.11.2.sub <- matrix(0, 3, 6)
rs.b.e.11.2.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 2, 0, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.2.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 2, 0.2, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.2.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 2, 0.5, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.2.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 2, 0, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.2.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 2, 0.2, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.2.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 2, 0.5, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.2.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 2, 0, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.2.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 2, 0.2, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.2.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 2, 0.5, c(0.1, 0.1), ex, 0, FALSE)





rs.b.e.55.2.sub <- matrix(0, 3, 6)
rs.b.e.55.2.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 2, 0, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.2.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 2, 0.2, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.2.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 2, 0.5, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.2.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 2, 0, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.2.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 2, 0.2, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.2.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 2, 0.5, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.2.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 2, 0, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.2.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 2, 0.2, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.2.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 2, 0.5, c(0.5, 0.5), ex, 0, FALSE)







rs.b.e.19.2.sub <- matrix(0, 3, 6)
rs.b.e.19.2.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 2, 0, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.2.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 2, 0.2, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.2.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 2, 0.5, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.2.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 2, 0, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.2.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 2, 0.2, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.2.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 2, 0.5, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.2.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 2, 0, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.2.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 2, 0.2, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.2.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 2, 0.5, c(-0.1, 0.9), ex, 0, FALSE)







rs.b.e.11.5.sub <- matrix(0, 3, 6)
rs.b.e.11.5.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 5, 0, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.5.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 5, 0.2, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.5.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 5, 0.5, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.5.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 5, 0, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.5.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 5, 0.2, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.5.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 5, 0.5, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.5.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 5, 0, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.5.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 5, 0.2, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.5.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 5, 0.5, c(0.1, 0.1), ex, 0, FALSE)





rs.b.e.55.5.sub <- matrix(0, 3, 6)
rs.b.e.55.5.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 5, 0, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.5.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 5, 0.2, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.5.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 5, 0.5, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.5.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 5, 0, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.5.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 5, 0.2, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.5.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 5, 0.5, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.5.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 5, 0, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.5.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 5, 0.2, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.5.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 5, 0.5, c(0.5, 0.5), ex, 0, FALSE)











rs.b.e.19.5.sub <- matrix(0, 3, 6)
rs.b.e.19.5.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 5, 0, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.5.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 5, 0.2, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.5.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 5, 0.5, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.5.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 5, 0, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.5.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 5, 0.2, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.5.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 5, 0.5, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.5.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 5, 0, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.5.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 5, 0.2, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.5.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 5, 0.5, c(-0.1, 0.9), ex, 0, FALSE)








rs.ub.e.11.10.sub <- matrix(0, 3, 6)
rs.ub.e.11.10.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 10, 0, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.e.11.10.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 10, 0.2, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.e.11.10.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 10, 0.5, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.e.11.10.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 10, 0, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.e.11.10.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 10, 0.2, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.e.11.10.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 10, 0.5, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.e.11.10.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 10, 0, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.e.11.10.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 10, 0.2, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.e.11.10.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 10, 0.5, c(0.1, 0.1), ex, 0.5, FALSE)










rs.ub.e.55.10.sub <- matrix(0, 3, 6)
rs.ub.e.55.10.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 10, 0, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.e.55.10.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 10, 0.2, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.e.55.10.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 10, 0.5, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.e.55.10.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 10, 0, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.e.55.10.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 10, 0.2, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.e.55.10.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 10, 0.5, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.e.55.10.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 10, 0, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.e.55.10.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 10, 0.2, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.e.55.10.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 10, 0.5, c(0.5, 0.5), ex, 0.5, FALSE)







rs.ub.e.19.10.sub <- matrix(0, 3, 6)
rs.ub.e.19.10.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 10, 0, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.e.19.10.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 10, 0.2, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.e.19.10.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 10, 0.5, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.e.19.10.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 10, 0, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.e.19.10.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 10, 0.2, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.e.19.10.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 10, 0.5, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.e.19.10.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 10, 0, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.e.19.10.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 10, 0.2, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.e.19.10.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 10, 0.5, c(-0.1, 0.9), ex, 0.5, FALSE)















rs.b.a.11.2.sub <- matrix(0, 3, 6)
rs.b.a.11.2.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 2, 0, c(0.1, 0.1), ex, 0, FALSE)
rs.b.a.11.2.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 2, 0.2, c(0.1, 0.1), ex, 0, FALSE)
rs.b.a.11.2.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 2, 0.5, c(0.1, 0.1), ex, 0, FALSE)
rs.b.a.11.2.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 2, 0, c(0.1, 0.1), ex, 0, FALSE)
rs.b.a.11.2.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 2, 0.2, c(0.1, 0.1), ex, 0, FALSE)
rs.b.a.11.2.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 2, 0.5, c(0.1, 0.1), ex, 0, FALSE)
rs.b.a.11.2.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 2, 0, c(0.1, 0.1), ex, 0, FALSE)
rs.b.a.11.2.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 2, 0.2, c(0.1, 0.1), ex, 0, FALSE)
rs.b.a.11.2.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 2, 0.5, c(0.1, 0.1), ex, 0, FALSE)





rs.b.a.55.2.sub <- matrix(0, 3, 6)
rs.b.a.55.2.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 2, 0, c(0.5, 0.5), ex, 0, FALSE)
rs.b.a.55.2.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 2, 0.2, c(0.5, 0.5), ex, 0, FALSE)
rs.b.a.55.2.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 2, 0.5, c(0.5, 0.5), ex, 0, FALSE)
rs.b.a.55.2.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 2, 0, c(0.5, 0.5), ex, 0, FALSE)
rs.b.a.55.2.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 2, 0.2, c(0.5, 0.5), ex, 0, FALSE)
rs.b.a.55.2.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 2, 0.5, c(0.5, 0.5), ex, 0, FALSE)
rs.b.a.55.2.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 2, 0, c(0.5, 0.5), ex, 0, FALSE)
rs.b.a.55.2.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 2, 0.2, c(0.5, 0.5), ex, 0, FALSE)
rs.b.a.55.2.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 2, 0.5, c(0.5, 0.5), ex, 0, FALSE)







rs.b.a.19.2.sub <- matrix(0, 3, 6)
rs.b.a.19.2.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 2, 0, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.a.19.2.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 2, 0.2, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.a.19.2.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 2, 0.5, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.a.19.2.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 2, 0, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.a.19.2.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 2, 0.2, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.a.19.2.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 2, 0.5, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.a.19.2.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 2, 0, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.a.19.2.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 2, 0.2, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.a.19.2.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 2, 0.5, c(-0.1, 0.9), ex, 0, FALSE)







rs.b.a.11.5.sub <- matrix(0, 3, 6)
rs.b.a.11.5.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 5, 0, c(0.1, 0.1), ex, 0, FALSE)
rs.b.a.11.5.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 5, 0.2, c(0.1, 0.1), ex, 0, FALSE)
rs.b.a.11.5.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 5, 0.5, c(0.1, 0.1), ex, 0, FALSE)
rs.b.a.11.5.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 5, 0, c(0.1, 0.1), ex, 0, FALSE)
rs.b.a.11.5.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 5, 0.2, c(0.1, 0.1), ex, 0, FALSE)
rs.b.a.11.5.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 5, 0.5, c(0.1, 0.1), ex, 0, FALSE)
rs.b.a.11.5.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 5, 0, c(0.1, 0.1), ex, 0, FALSE)
rs.b.a.11.5.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 5, 0.2, c(0.1, 0.1), ex, 0, FALSE)
rs.b.a.11.5.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 5, 0.5, c(0.1, 0.1), ex, 0, FALSE)





rs.b.a.55.5.sub <- matrix(0, 3, 6)
rs.b.a.55.5.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 5, 0, c(0.5, 0.5), ex, 0, FALSE)
rs.b.a.55.5.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 5, 0.2, c(0.5, 0.5), ex, 0, FALSE)
rs.b.a.55.5.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 5, 0.5, c(0.5, 0.5), ex, 0, FALSE)
rs.b.a.55.5.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 5, 0, c(0.5, 0.5), ex, 0, FALSE)
rs.b.a.55.5.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 5, 0.2, c(0.5, 0.5), ex, 0, FALSE)
rs.b.a.55.5.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 5, 0.5, c(0.5, 0.5), ex, 0, FALSE)
rs.b.a.55.5.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 5, 0, c(0.5, 0.5), ex, 0, FALSE)
rs.b.a.55.5.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 5, 0.2, c(0.5, 0.5), ex, 0, FALSE)
rs.b.a.55.5.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 5, 0.5, c(0.5, 0.5), ex, 0, FALSE)











rs.b.a.19.5.sub <- matrix(0, 3, 6)
rs.b.a.19.5.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 5, 0, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.a.19.5.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 5, 0.2, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.a.19.5.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 5, 0.5, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.a.19.5.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 5, 0, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.a.19.5.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 5, 0.2, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.a.19.5.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 5, 0.5, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.a.19.5.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 5, 0, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.a.19.5.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 5, 0.2, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.a.19.5.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 5, 0.5, c(-0.1, 0.9), ex, 0, FALSE)








rs.ub.a.11.10.sub <- matrix(0, 3, 6)
rs.ub.a.11.10.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 10, 0, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.a.11.10.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 10, 0.2, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.a.11.10.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 10, 0.5, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.a.11.10.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 10, 0, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.a.11.10.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 10, 0.2, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.a.11.10.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 10, 0.5, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.a.11.10.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 10, 0, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.a.11.10.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 10, 0.2, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.a.11.10.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 10, 0.5, c(0.1, 0.1), ex, 0.5, FALSE)










rs.ub.a.55.10.sub <- matrix(0, 3, 6)
rs.ub.a.55.10.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 10, 0, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.a.55.10.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 10, 0.2, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.a.55.10.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 10, 0.5, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.a.55.10.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 10, 0, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.a.55.10.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 10, 0.2, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.a.55.10.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 10, 0.5, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.a.55.10.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 10, 0, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.a.55.10.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 10, 0.2, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.a.55.10.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 10, 0.5, c(0.5, 0.5), ex, 0.5, FALSE)







rs.ub.a.19.10.sub <- matrix(0, 3, 6)
rs.ub.a.19.10.sub[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 10, 0, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.a.19.10.sub[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 10, 0.2, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.a.19.10.sub[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 10, 0.5, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.a.19.10.sub[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 10, 0, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.a.19.10.sub[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 10, 0.2, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.a.19.10.sub[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 10, 0.5, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.a.19.10.sub[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 10, 0, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.a.19.10.sub[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 10, 0.2, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.a.19.10.sub[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 10, 0.5, c(-0.1, 0.9), ex, 0.5, FALSE)







sr.b.e.1.2 <- matrix(0, 3, 6)
sr.b.e.1.2[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 2, 0, 0.1, ex, 0)
sr.b.e.1.2[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 2, 0.2, 0.1, ex, 0)
sr.b.e.1.2[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 2, 0.5, 0.1, ex, 0)
sr.b.e.1.2[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 2, 0, 0.1, ex, 0)
sr.b.e.1.2[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 2, 0.2, 0.1, ex, 0)
sr.b.e.1.2[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 2, 0.5, 0.1, ex, 0)
sr.b.e.1.2[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 2, 0, 0.1, ex, 0)
sr.b.e.1.2[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 2, 0.2, 0.1, ex, 0)
sr.b.e.1.2[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 2, 0.5, 0.1, ex, 0)





sr.b.e.5.2 <- matrix(0, 3, 6)
sr.b.e.5.2[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 2, 0, 0.5, ex, 0)
sr.b.e.5.2[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 2, 0.2, 0.5, ex, 0)
sr.b.e.5.2[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 2, 0.5, 0.5, ex, 0)
sr.b.e.5.2[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 2, 0, 0.5, ex, 0)
sr.b.e.5.2[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 2, 0.2, 0.5, ex, 0)
sr.b.e.5.2[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 2, 0.5, 0.5, ex, 0)
sr.b.e.5.2[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 2, 0, 0.5, ex, 0)
sr.b.e.5.2[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 2, 0.2, 0.5, ex, 0)
sr.b.e.5.2[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 2, 0.5, 0.5, ex, 0)


sr.b.e.9.2 <- matrix(0, 3, 6)
sr.b.e.9.2[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 2, 0, 0.9, ex, 0)
sr.b.e.9.2[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 2, 0.2, 0.9, ex, 0)
sr.b.e.9.2[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 2, 0.5, 0.9, ex, 0)
sr.b.e.9.2[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 2, 0, 0.9, ex, 0)
sr.b.e.9.2[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 2, 0.2, 0.9, ex, 0)
sr.b.e.9.2[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 2, 0.5, 0.9, ex, 0)
sr.b.e.9.2[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 2, 0, 0.9, ex, 0)
sr.b.e.9.2[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 2, 0.2, 0.9, ex, 0)
sr.b.e.9.2[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 2, 0.5, 0.9, ex, 0)








sr.b.e.1.10 <- matrix(0, 3, 6)
sr.b.e.1.10[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 10, 0, 0.1, ex, 0)
sr.b.e.1.10[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 10, 0.2, 0.1, ex, 0)
sr.b.e.1.10[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 10, 0.5, 0.1, ex, 0)
sr.b.e.1.10[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 10, 0, 0.1, ex, 0)
sr.b.e.1.10[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 10, 0.2, 0.1, ex, 0)
sr.b.e.1.10[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 10, 0.5, 0.1, ex, 0)
sr.b.e.1.10[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 10, 0, 0.1, ex, 0)
sr.b.e.1.10[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 10, 0.2, 0.1, ex, 0)
sr.b.e.1.10[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 10, 0.5, 0.1, ex, 0)





sr.b.e.5.10 <- matrix(0, 3, 6)
sr.b.e.5.10[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 10, 0, 0.5, ex, 0)
sr.b.e.5.10[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 10, 0.2, 0.5, ex, 0)
sr.b.e.5.10[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 10, 0.5, 0.5, ex, 0)
sr.b.e.5.10[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 10, 0, 0.5, ex, 0)
sr.b.e.5.10[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 10, 0.2, 0.5, ex, 0)
sr.b.e.5.10[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 10, 0.5, 0.5, ex, 0)
sr.b.e.5.10[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 10, 0, 0.5, ex, 0)
sr.b.e.5.10[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 10, 0.2, 0.5, ex, 0)
sr.b.e.5.10[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 10, 0.5, 0.5, ex, 0)


sr.b.e.9.10 <- matrix(0, 3, 6)
sr.b.e.9.10[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 10, 0, 0.9, ex, 0)
sr.b.e.9.10[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 10, 0.2, 0.9, ex, 0)
sr.b.e.9.10[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 10, 0.5, 0.9, ex, 0)
sr.b.e.9.10[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 10, 0, 0.9, ex, 0)
sr.b.e.9.10[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 10, 0.2, 0.9, ex, 0)
sr.b.e.9.10[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 10, 0.5, 0.9, ex, 0)
sr.b.e.9.10[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 10, 0, 0.9, ex, 0)
sr.b.e.9.10[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 10, 0.2, 0.9, ex, 0)
sr.b.e.9.10[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 10, 0.5, 0.9, ex, 0)






sr.ub.e.1.5 <- matrix(0, 3, 6)
sr.ub.e.1.5[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 5, 0, 0.1, ex, 0.5)
sr.ub.e.1.5[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 5, 0.2, 0.1, ex, 0.5)
sr.ub.e.1.5[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 5, 0.5, 0.1, ex, 0.5)
sr.ub.e.1.5[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 5, 0, 0.1, ex, 0.5)
sr.ub.e.1.5[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 5, 0.2, 0.1, ex, 0.5)
sr.ub.e.1.5[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 5, 0.5, 0.1, ex, 0.5)
sr.ub.e.1.5[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 5, 0, 0.1, ex, 0.5)
sr.ub.e.1.5[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 5, 0.2, 0.1, ex, 0.5)
sr.ub.e.1.5[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 5, 0.5, 0.1, ex, 0.5)





sr.ub.e.5.5 <- matrix(0, 3, 6)
sr.ub.e.5.5[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 5, 0, 0.5, ex, 0.5)
sr.ub.e.5.5[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 5, 0.2, 0.5, ex, 0.5)
sr.ub.e.5.5[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 5, 0.5, 0.5, ex, 0.5)
sr.ub.e.5.5[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 5, 0, 0.5, ex, 0.5)
sr.ub.e.5.5[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 5, 0.2, 0.5, ex, 0.5)
sr.ub.e.5.5[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 5, 0.5, 0.5, ex, 0.5)
sr.ub.e.5.5[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 5, 0, 0.5, ex, 0.5)
sr.ub.e.5.5[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 5, 0.2, 0.5, ex, 0.5)
sr.ub.e.5.5[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 5, 0.5, 0.5, ex, 0.5)


sr.ub.e.9.5 <- matrix(0, 3, 6)
sr.ub.e.9.5[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 5, 0, 0.9, ex, 0.5)
sr.ub.e.9.5[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 5, 0.2, 0.9, ex, 0.5)
sr.ub.e.9.5[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 5, 0.5, 0.9, ex, 0.5)
sr.ub.e.9.5[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 5, 0, 0.9, ex, 0.5)
sr.ub.e.9.5[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 5, 0.2, 0.9, ex, 0.5)
sr.ub.e.9.5[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 5, 0.5, 0.9, ex, 0.5)
sr.ub.e.9.5[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 5, 0, 0.9, ex, 0.5)
sr.ub.e.9.5[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 5, 0.2, 0.9, ex, 0.5)
sr.ub.e.9.5[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 5, 0.5, 0.9, ex, 0.5)










sr.ub.e.1.10 <- matrix(0, 3, 6)
sr.ub.e.1.10[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 10, 0, 0.1, ex, 0.5)
sr.ub.e.1.10[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 10, 0.2, 0.1, ex, 0.5)
sr.ub.e.1.10[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 10, 0.5, 0.1, ex, 0.5)
sr.ub.e.1.10[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 10, 0, 0.1, ex, 0.5)
sr.ub.e.1.10[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 10, 0.2, 0.1, ex, 0.5)
sr.ub.e.1.10[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 10, 0.5, 0.1, ex, 0.5)
sr.ub.e.1.10[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 10, 0, 0.1, ex, 0.5)
sr.ub.e.1.10[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 10, 0.2, 0.1, ex, 0.5)
sr.ub.e.1.10[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 10, 0.5, 0.1, ex, 0.5)





sr.ub.e.5.10 <- matrix(0, 3, 6)
sr.ub.e.5.10[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 10, 0, 0.5, ex, 0.5)
sr.ub.e.5.10[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 10, 0.2, 0.5, ex, 0.5)
sr.ub.e.5.10[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 10, 0.5, 0.5, ex, 0.5)
sr.ub.e.5.10[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 10, 0, 0.5, ex, 0.5)
sr.ub.e.5.10[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 10, 0.2, 0.5, ex, 0.5)
sr.ub.e.5.10[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 10, 0.5, 0.5, ex, 0.5)
sr.ub.e.5.10[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 10, 0, 0.5, ex, 0.5)
sr.ub.e.5.10[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 10, 0.2, 0.5, ex, 0.5)
sr.ub.e.5.10[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 10, 0.5, 0.5, ex, 0.5)


sr.ub.e.9.10 <- matrix(0, 3, 6)
sr.ub.e.9.10[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 10, 0, 0.9, ex, 0.5)
sr.ub.e.9.10[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 10, 0.2, 0.9, ex, 0.5)
sr.ub.e.9.10[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 10, 0.5, 0.9, ex, 0.5)
sr.ub.e.9.10[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 10, 0, 0.9, ex, 0.5)
sr.ub.e.9.10[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 10, 0.2, 0.9, ex, 0.5)
sr.ub.e.9.10[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 10, 0.5, 0.9, ex, 0.5)
sr.ub.e.9.10[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 10, 0, 0.9, ex, 0.5)
sr.ub.e.9.10[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 10, 0.2, 0.9, ex, 0.5)
sr.ub.e.9.10[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 10, 0.5, 0.9, ex, 0.5)













rs.b.a.11.2 <- matrix(0, 3, 6)
rs.b.a.11.2[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 2, 0, c(0.1, 0.1), ar1, 0)
rs.b.a.11.2[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 2, 0.2, c(0.1, 0.1), ar1, 0)
rs.b.a.11.2[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 2, 0.5, c(0.1, 0.1), ar1, 0)
rs.b.a.11.2[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 2, 0, c(0.1, 0.1), ar1, 0)
rs.b.a.11.2[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 2, 0.2, c(0.1, 0.1), ar1, 0)
rs.b.a.11.2[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 2, 0.5, c(0.1, 0.1), ar1, 0)
rs.b.a.11.2[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 2, 0, c(0.1, 0.1), ar1, 0)
rs.b.a.11.2[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 2, 0.2, c(0.1, 0.1), ar1, 0)
rs.b.a.11.2[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 2, 0.5, c(0.1, 0.1), ar1, 0)





rs.b.a.55.2 <- matrix(0, 3, 6)
rs.b.a.55.2[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 2, 0, c(0.5, 0.5), ar1, 0)
rs.b.a.55.2[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 2, 0.2, c(0.5, 0.5), ar1, 0)
rs.b.a.55.2[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 2, 0.5, c(0.5, 0.5), ar1, 0)
rs.b.a.55.2[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 2, 0, c(0.5, 0.5), ar1, 0)
rs.b.a.55.2[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 2, 0.2, c(0.5, 0.5), ar1, 0)
rs.b.a.55.2[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 2, 0.5, c(0.5, 0.5), ar1, 0)
rs.b.a.55.2[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 2, 0, c(0.5, 0.5), ar1, 0)
rs.b.a.55.2[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 2, 0.2, c(0.5, 0.5), ar1, 0)
rs.b.a.55.2[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 2, 0.5, c(0.5, 0.5), ar1, 0)







rs.b.a.19.2 <- matrix(0, 3, 6)
rs.b.a.19.2[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 2, 0, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.2[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 2, 0.2, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.2[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 2, 0.5, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.2[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 2, 0, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.2[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 2, 0.2, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.2[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 2, 0.5, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.2[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 2, 0, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.2[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 2, 0.2, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.2[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 2, 0.5, c(-0.1, 0.9), ar1, 0)







rs.b.a.11.5 <- matrix(0, 3, 6)
rs.b.a.11.5[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 5, 0, c(0.1, 0.1), ar1, 0)
rs.b.a.11.5[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 5, 0.2, c(0.1, 0.1), ar1, 0)
rs.b.a.11.5[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 5, 0.5, c(0.1, 0.1), ar1, 0)
rs.b.a.11.5[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 5, 0, c(0.1, 0.1), ar1, 0)
rs.b.a.11.5[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 5, 0.2, c(0.1, 0.1), ar1, 0)
rs.b.a.11.5[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 5, 0.5, c(0.1, 0.1), ar1, 0)
rs.b.a.11.5[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 5, 0, c(0.1, 0.1), ar1, 0)
rs.b.a.11.5[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 5, 0.2, c(0.1, 0.1), ar1, 0)
rs.b.a.11.5[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 5, 0.5, c(0.1, 0.1), ar1, 0)





rs.b.a.55.5 <- matrix(0, 3, 6)
rs.b.a.55.5[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 5, 0, c(0.5, 0.5), ar1, 0)
rs.b.a.55.5[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 5, 0.2, c(0.5, 0.5), ar1, 0)
rs.b.a.55.5[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 5, 0.5, c(0.5, 0.5), ar1, 0)
rs.b.a.55.5[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 5, 0, c(0.5, 0.5), ar1, 0)
rs.b.a.55.5[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 5, 0.2, c(0.5, 0.5), ar1, 0)
rs.b.a.55.5[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 5, 0.5, c(0.5, 0.5), ar1, 0)
rs.b.a.55.5[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 5, 0, c(0.5, 0.5), ar1, 0)
rs.b.a.55.5[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 5, 0.2, c(0.5, 0.5), ar1, 0)
rs.b.a.55.5[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 5, 0.5, c(0.5, 0.5), ar1, 0)











rs.b.a.19.5 <- matrix(0, 3, 6)
rs.b.a.19.5[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 5, 0, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.5[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 5, 0.2, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.5[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 5, 0.5, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.5[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 5, 0, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.5[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 5, 0.2, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.5[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 5, 0.5, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.5[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 5, 0, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.5[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 5, 0.2, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.5[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 5, 0.5, c(-0.1, 0.9), ar1, 0)








rs.ub.a.11.10 <- matrix(0, 3, 6)
rs.ub.a.11.10[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 10, 0, c(0.1, 0.1), ar1, 0.5)
rs.ub.a.11.10[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 10, 0.2, c(0.1, 0.1), ar1, 0.5)
rs.ub.a.11.10[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 10, 0.5, c(0.1, 0.1), ar1, 0.5)
rs.ub.a.11.10[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 10, 0, c(0.1, 0.1), ar1, 0.5)
rs.ub.a.11.10[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 10, 0.2, c(0.1, 0.1), ar1, 0.5)
rs.ub.a.11.10[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 10, 0.5, c(0.1, 0.1), ar1, 0.5)
rs.ub.a.11.10[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 10, 0, c(0.1, 0.1), ar1, 0.5)
rs.ub.a.11.10[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 10, 0.2, c(0.1, 0.1), ar1, 0.5)
rs.ub.a.11.10[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 10, 0.5, c(0.1, 0.1), ar1, 0.5)










rs.ub.a.55.10 <- matrix(0, 3, 6)
rs.ub.a.55.10[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 10, 0, c(0.5, 0.5), ar1, 0.5)
rs.ub.a.55.10[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 10, 0.2, c(0.5, 0.5), ar1, 0.5)
rs.ub.a.55.10[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 10, 0.5, c(0.5, 0.5), ar1, 0.5)
rs.ub.a.55.10[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 10, 0, c(0.5, 0.5), ar1, 0.5)
rs.ub.a.55.10[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 10, 0.2, c(0.5, 0.5), ar1, 0.5)
rs.ub.a.55.10[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 10, 0.5, c(0.5, 0.5), ar1, 0.5)
rs.ub.a.55.10[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 10, 0, c(0.5, 0.5), ar1, 0.5)
rs.ub.a.55.10[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 10, 0.2, c(0.5, 0.5), ar1, 0.5)
rs.ub.a.55.10[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 10, 0.5, c(0.5, 0.5), ar1, 0.5)







rs.ub.a.19.10 <- matrix(0, 3, 6)
rs.ub.a.19.10[1, 1:2] <- simpower(4000, 0.05, FALSE, 20, 10, 0, c(-0.1, 0.9), ar1, 0.5)
rs.ub.a.19.10[1, 3:4] <- simpower(4000, 0.05, FALSE, 20, 10, 0.2, c(-0.1, 0.9), ar1, 0.5)
rs.ub.a.19.10[1, 5:6] <- simpower(4000, 0.05, FALSE, 20, 10, 0.5, c(-0.1, 0.9), ar1, 0.5)
rs.ub.a.19.10[2, 1:2] <- simpower(4000, 0.05, FALSE, 50, 10, 0, c(-0.1, 0.9), ar1, 0.5)
rs.ub.a.19.10[2, 3:4] <- simpower(4000, 0.05, FALSE, 50, 10, 0.2, c(-0.1, 0.9), ar1, 0.5)
rs.ub.a.19.10[2, 5:6] <- simpower(4000, 0.05, FALSE, 50, 10, 0.5, c(-0.1, 0.9), ar1, 0.5)
rs.ub.a.19.10[3, 1:2] <- simpower(4000, 0.05, FALSE, 100, 10, 0, c(-0.1, 0.9), ar1, 0.5)
rs.ub.a.19.10[3, 3:4] <- simpower(4000, 0.05, FALSE, 100, 10, 0.2, c(-0.1, 0.9), ar1, 0.5)
rs.ub.a.19.10[3, 5:6] <- simpower(4000, 0.05, FALSE, 100, 10, 0.5, c(-0.1, 0.9), ar1, 0.5)





























sr.b.a.1.2 <- matrix(0, 3, 6)
sr.b.a.1.2[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 2, 0, 0.1, ar1, 0)
sr.b.a.1.2[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 2, 0.2, 0.1, ar1, 0)
sr.b.a.1.2[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 2, 0.5, 0.1, ar1, 0)
sr.b.a.1.2[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 2, 0, 0.1, ar1, 0)
sr.b.a.1.2[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 2, 0.2, 0.1, ar1, 0)
sr.b.a.1.2[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 2, 0.5, 0.1, ar1, 0)
sr.b.a.1.2[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 2, 0, 0.1, ar1, 0)
sr.b.a.1.2[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 2, 0.2, 0.1, ar1, 0)
sr.b.a.1.2[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 2, 0.5, 0.1, ar1, 0)





sr.b.a.5.2 <- matrix(0, 3, 6)
sr.b.a.5.2[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 2, 0, 0.5, ar1, 0)
sr.b.a.5.2[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 2, 0.2, 0.5, ar1, 0)
sr.b.a.5.2[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 2, 0.5, 0.5, ar1, 0)
sr.b.a.5.2[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 2, 0, 0.5, ar1, 0)
sr.b.a.5.2[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 2, 0.2, 0.5, ar1, 0)
sr.b.a.5.2[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 2, 0.5, 0.5, ar1, 0)
sr.b.a.5.2[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 2, 0, 0.5, ar1, 0)
sr.b.a.5.2[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 2, 0.2, 0.5, ar1, 0)
sr.b.a.5.2[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 2, 0.5, 0.5, ar1, 0)


sr.b.a.9.2 <- matrix(0, 3, 6)
sr.b.a.9.2[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 2, 0, 0.9, ar1, 0)
sr.b.a.9.2[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 2, 0.2, 0.9, ar1, 0)
sr.b.a.9.2[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 2, 0.5, 0.9, ar1, 0)
sr.b.a.9.2[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 2, 0, 0.9, ar1, 0)
sr.b.a.9.2[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 2, 0.2, 0.9, ar1, 0)
sr.b.a.9.2[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 2, 0.5, 0.9, ar1, 0)
sr.b.a.9.2[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 2, 0, 0.9, ar1, 0)
sr.b.a.9.2[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 2, 0.2, 0.9, ar1, 0)
sr.b.a.9.2[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 2, 0.5, 0.9, ar1, 0)








sr.b.a.1.10 <- matrix(0, 3, 6)
sr.b.a.1.10[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 10, 0, 0.1, ar1, 0)
sr.b.a.1.10[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 10, 0.2, 0.1, ar1, 0)
sr.b.a.1.10[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 10, 0.5, 0.1, ar1, 0)
sr.b.a.1.10[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 10, 0, 0.1, ar1, 0)
sr.b.a.1.10[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 10, 0.2, 0.1, ar1, 0)
sr.b.a.1.10[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 10, 0.5, 0.1, ar1, 0)
sr.b.a.1.10[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 10, 0, 0.1, ar1, 0)
sr.b.a.1.10[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 10, 0.2, 0.1, ar1, 0)
sr.b.a.1.10[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 10, 0.5, 0.1, ar1, 0)





sr.b.a.5.10 <- matrix(0, 3, 6)
sr.b.a.5.10[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 10, 0, 0.5, ar1, 0)
sr.b.a.5.10[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 10, 0.2, 0.5, ar1, 0)
sr.b.a.5.10[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 10, 0.5, 0.5, ar1, 0)
sr.b.a.5.10[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 10, 0, 0.5, ar1, 0)
sr.b.a.5.10[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 10, 0.2, 0.5, ar1, 0)
sr.b.a.5.10[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 10, 0.5, 0.5, ar1, 0)
sr.b.a.5.10[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 10, 0, 0.5, ar1, 0)
sr.b.a.5.10[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 10, 0.2, 0.5, ar1, 0)
sr.b.a.5.10[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 10, 0.5, 0.5, ar1, 0)


sr.b.a.9.10 <- matrix(0, 3, 6)
sr.b.a.9.10[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 10, 0, 0.9, ar1, 0)
sr.b.a.9.10[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 10, 0.2, 0.9, ar1, 0)
sr.b.a.9.10[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 10, 0.5, 0.9, ar1, 0)
sr.b.a.9.10[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 10, 0, 0.9, ar1, 0)
sr.b.a.9.10[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 10, 0.2, 0.9, ar1, 0)
sr.b.a.9.10[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 10, 0.5, 0.9, ar1, 0)
sr.b.a.9.10[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 10, 0, 0.9, ar1, 0)
sr.b.a.9.10[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 10, 0.2, 0.9, ar1, 0)
sr.b.a.9.10[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 10, 0.5, 0.9, ar1, 0)






sr.ub.a.1.5 <- matrix(0, 3, 6)
sr.ub.a.1.5[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 5, 0, 0.1, ar1, 0.5)
sr.ub.a.1.5[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 5, 0.2, 0.1, ar1, 0.5)
sr.ub.a.1.5[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 5, 0.5, 0.1, ar1, 0.5)
sr.ub.a.1.5[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 5, 0, 0.1, ar1, 0.5)
sr.ub.a.1.5[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 5, 0.2, 0.1, ar1, 0.5)
sr.ub.a.1.5[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 5, 0.5, 0.1, ar1, 0.5)
sr.ub.a.1.5[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 5, 0, 0.1, ar1, 0.5)
sr.ub.a.1.5[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 5, 0.2, 0.1, ar1, 0.5)
sr.ub.a.1.5[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 5, 0.5, 0.1, ar1, 0.5)





sr.ub.a.5.5 <- matrix(0, 3, 6)
sr.ub.a.5.5[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 5, 0, 0.5, ar1, 0.5)
sr.ub.a.5.5[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 5, 0.2, 0.5, ar1, 0.5)
sr.ub.a.5.5[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 5, 0.5, 0.5, ar1, 0.5)
sr.ub.a.5.5[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 5, 0, 0.5, ar1, 0.5)
sr.ub.a.5.5[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 5, 0.2, 0.5, ar1, 0.5)
sr.ub.a.5.5[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 5, 0.5, 0.5, ar1, 0.5)
sr.ub.a.5.5[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 5, 0, 0.5, ar1, 0.5)
sr.ub.a.5.5[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 5, 0.2, 0.5, ar1, 0.5)
sr.ub.a.5.5[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 5, 0.5, 0.5, ar1, 0.5)


sr.ub.a.9.5 <- matrix(0, 3, 6)
sr.ub.a.9.5[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 5, 0, 0.9, ar1, 0.5)
sr.ub.a.9.5[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 5, 0.2, 0.9, ar1, 0.5)
sr.ub.a.9.5[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 5, 0.5, 0.9, ar1, 0.5)
sr.ub.a.9.5[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 5, 0, 0.9, ar1, 0.5)
sr.ub.a.9.5[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 5, 0.2, 0.9, ar1, 0.5)
sr.ub.a.9.5[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 5, 0.5, 0.9, ar1, 0.5)
sr.ub.a.9.5[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 5, 0, 0.9, ar1, 0.5)
sr.ub.a.9.5[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 5, 0.2, 0.9, ar1, 0.5)
sr.ub.a.9.5[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 5, 0.5, 0.9, ar1, 0.5)










sr.ub.a.1.10 <- matrix(0, 3, 6)
sr.ub.a.1.10[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 10, 0, 0.1, ar1, 0.5)
sr.ub.a.1.10[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 10, 0.2, 0.1, ar1, 0.5)
sr.ub.a.1.10[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 10, 0.5, 0.1, ar1, 0.5)
sr.ub.a.1.10[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 10, 0, 0.1, ar1, 0.5)
sr.ub.a.1.10[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 10, 0.2, 0.1, ar1, 0.5)
sr.ub.a.1.10[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 10, 0.5, 0.1, ar1, 0.5)
sr.ub.a.1.10[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 10, 0, 0.1, ar1, 0.5)
sr.ub.a.1.10[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 10, 0.2, 0.1, ar1, 0.5)
sr.ub.a.1.10[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 10, 0.5, 0.1, ar1, 0.5)





sr.ub.a.5.10 <- matrix(0, 3, 6)
sr.ub.a.5.10[1, 1:2] <- simpower(4000, 0.05, TRUE, 20, 10, 0, 0.5, ar1, 0.5)
sr.ub.a.5.10[1, 3:4] <- simpower(4000, 0.05, TRUE, 20, 10, 0.2, 0.5, ar1, 0.5)
sr.ub.a.5.10[1, 5:6] <- simpower(4000, 0.05, TRUE, 20, 10, 0.5, 0.5, ar1, 0.5)
sr.ub.a.5.10[2, 1:2] <- simpower(4000, 0.05, TRUE, 50, 10, 0, 0.5, ar1, 0.5)
sr.ub.a.5.10[2, 3:4] <- simpower(4000, 0.05, TRUE, 50, 10, 0.2, 0.5, ar1, 0.5)
sr.ub.a.5.10[2, 5:6] <- simpower(4000, 0.05, TRUE, 50, 10, 0.5, 0.5, ar1, 0.5)
sr.ub.a.5.10[3, 1:2] <- simpower(4000, 0.05, TRUE, 100, 10, 0, 0.5, ar1, 0.5)
sr.ub.a.5.10[3, 3:4] <- simpower(4000, 0.05, TRUE, 100, 10, 0.2, 0.5, ar1, 0.5)
sr.ub.a.5.10[3, 5:6] <- simpower(4000, 0.05, TRUE, 100, 10, 0.5, 0.5, ar1, 0.5)


sr.ub.a.9.10 <- matrix(0, 3, 6)
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












library(xtable)

sum.tab.rs.e2 <- as.data.frame(sum.tab.rs.e / 4000)
sum.tab.rs.a2 <- as.data.frame(sum.tab.rs.a / 4000)
sum.tab.sr.e2 <- as.data.frame(sum.tab.sr.e / 4000)
sum.tab.sr.a2 <- as.data.frame(sum.tab.sr.a / 4000)

colnames(sum.tab.rs.e2) <- paste0(rep(c("RGL", "DS"), 3), c("(d=0)", "", "(d=0.2)", "", "(d=0.5)", ""))
sum.tab.rs.e3 <- cbind(N = c("20", "50", "100"), sum.tab.rs.e2)
sum.tab.rs.e3 <- cbind(rho = c("0.1, 0.1", "", "" ,"0.5, 0.5", "", "","$-$0.1, 0.9", "", ""), sum.tab.rs.e3)
ni <- c(2, rep("", 8), 5, rep("", 8), 10, rep("", 8))
sum.tab.rs.e3 <- cbind(ni = ni, sum.tab.rs.e3)
missingrate <- c(0, rep("", 17), 0.5, rep("", 8))
sum.tab.rs.e3 <- cbind(missingrate = missingrate, sum.tab.rs.e3)
grplev <- c("cluster", rep("", 26), "subunit", rep("", 26))
sum.tab.rs.e3 <- cbind(grplev = grplev, sum.tab.rs.e3)
print(xtable(sum.tab.rs.e3, caption = "Empirical rejection percentage of the RGL and the DS methods for rank-sum tests at nominal significance level 0.05 when intracluster correlation is exchangable. The results are based on 4000 datasets.", digits = 1), include.rownames = FALSE, sanitize.text.function = function(x) {x})







colnames(sum.tab.rs.a2) <- paste0(rep(c("RGL", "DS"), 3), c("(d=0)", "", "(d=0.2)", "", "(d=0.5)", ""))
sum.tab.rs.a3 <- cbind(N = c("20", "50", "100"), sum.tab.rs.a2)
sum.tab.rs.a3 <- cbind(rho = c("0.1, 0.1", "", "" ,"0.5, 0.5", "", "","$-$0.1, 0.9", "", ""), sum.tab.rs.a3)
ni <- c(2, rep("", 8), 5, rep("", 8), 10, rep("", 8))
sum.tab.rs.a3 <- cbind(ni = ni, sum.tab.rs.a3)
missingrate <- c(0, rep("", 17), 0.5, rep("", 8))
sum.tab.rs.a3 <- cbind(missingrate = missingrate, sum.tab.rs.a3)
grplev <- c("cluster", rep("", 26), "subunit", rep("", 26))
sum.tab.rs.a3 <- cbind(grplev = grplev, sum.tab.rs.a3)
print(xtable(sum.tab.rs.a3, caption = "Empirical rejection percentage of the RGL and the DS methods for rank-sum tests at nominal significance level 0.05 when intracluster correlation is AR1. The results are based on 4000 datasets.", digits = 1), include.rownames = FALSE, sanitize.text.function = function(x) {x})


colnames(sum.tab.sr.e2) <- paste0(rep(c("RGL", "DS"), 3), c("(d=0)", "", "(d=0.2)", "", "(d=0.5)", ""))
sum.tab.sr.e3 <- cbind(N = c("20", "50", "100"), sum.tab.sr.e2)
sum.tab.sr.e3 <- cbind(rho = c("0.1", "", "" ,"0.5", "", "","0.9", "", ""), sum.tab.sr.e3)
ni <- c(2, rep("", 8), 10, rep("", 8), 5, rep("", 8), 10, rep("", 8))
sum.tab.sr.e3 <- cbind(ni = ni, sum.tab.sr.e3)
missingrate <- c(0, rep("", 17), 0.5, rep("", 17))
sum.tab.sr.e3 <- cbind(missingrate = missingrate, sum.tab.sr.e3)
print(xtable(sum.tab.sr.e3, caption = "Empirical rejection percentage of the RGL and DS methods for signed-rank tests at nominal significance level 0.05 with exchangable intracluster correlation. The results are based on 4000 datasets.", digits = 1), include.rownames = FALSE)






colnames(sum.tab.sr.a2) <- paste0(rep(c("RGL", "DS"), 3), c("(d=0)", "", "(d=0.2)", "", "(d=0.5)", ""))
sum.tab.sr.a3 <- cbind(N = c("20", "50", "100"), sum.tab.sr.a2)
sum.tab.sr.a3 <- cbind(rho = c("0.1", "", "" ,"0.5", "", "","0.9", "", ""), sum.tab.sr.a3)
ni <- c(2, rep("", 8), 10, rep("", 8), 5, rep("", 8), 10, rep("", 8))
sum.tab.sr.a3 <- cbind(ni = ni, sum.tab.sr.a3)
missingrate <- c(0, rep("", 17), 0.5, rep("", 17))
sum.tab.sr.a3 <- cbind(missingrate = missingrate, sum.tab.sr.a3)
print(xtable(sum.tab.sr.a3, caption = "Empirical rejection percentage of the RGL and DS methods for signed-rank tests at nominal significance level 0.05 with AR1 intracluster correlation. The results are based on 4000 datasets."), include.rownames = FALSE)


