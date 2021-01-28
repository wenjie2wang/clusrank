## R code for simulation tables in 'clusrank.Rnw'
## need to source('clusrank.R') first
source("clusrank.R")

nrep <- 4000

####################################################################
### Rank sum, exchangable, rho = 0.1, 0.1, clus size = 2, delta = 0,
### 0.2, 0.5, grpsize = 20, 50, balanced data
####################################################################

rs.b.e.11.2 <- matrix(0, 2, 6)
rs.b.e.11.2[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 2, 0, c(0.1, 0.1), ex, 0)
rs.b.e.11.2[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 2, 0.2, c(0.1, 0.1), ex, 0)
rs.b.e.11.2[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 2, 0.5, c(0.1, 0.1), ex, 0)
rs.b.e.11.2[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 2, 0, c(0.1, 0.1), ex, 0)
rs.b.e.11.2[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 2, 0.2, c(0.1, 0.1), ex, 0)
rs.b.e.11.2[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 2, 0.5, c(0.1, 0.1), ex, 0)

####################################################################
### Rank sum, exchangable, rho = 0.5, 0.5, clus size = 2, delta = 0,
### 0.2, 0.5, grpsize = 20, 50, balanced data
####################################################################
rs.b.e.55.2 <- matrix(0, 2, 6)
rs.b.e.55.2[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 2, 0, c(0.5, 0.5), ex, 0)
rs.b.e.55.2[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 2, 0.2, c(0.5, 0.5), ex, 0)
rs.b.e.55.2[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 2, 0.5, c(0.5, 0.5), ex, 0)
rs.b.e.55.2[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 2, 0, c(0.5, 0.5), ex, 0)
rs.b.e.55.2[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 2, 0.2, c(0.5, 0.5), ex, 0)
rs.b.e.55.2[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 2, 0.5, c(0.5, 0.5), ex, 0)

####################################################################
### Rank sum, exchangable, rho = -0.1, 0.9, clus size = 2, delta = 0,
### 0.2, 0.5, grpsize = 20, 50, balanced data
####################################################################
rs.b.e.19.2 <- matrix(0, 2, 6)
rs.b.e.19.2[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 2, 0, c(-0.1, 0.9), ex, 0)
rs.b.e.19.2[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 2, 0.2, c(-0.1, 0.9), ex, 0)
rs.b.e.19.2[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 2, 0.5, c(-0.1, 0.9), ex, 0)
rs.b.e.19.2[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 2, 0, c(-0.1, 0.9), ex, 0)
rs.b.e.19.2[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 2, 0.2, c(-0.1, 0.9), ex, 0)
rs.b.e.19.2[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 2, 0.5, c(-0.1, 0.9), ex, 0)

####################################################################
### Rank sum, exchangable, rho = 0.1, 0.1, clus size = 5, delta = 0,
### 0.2, 0.5, grpsize = 20, 50, balanced data
####################################################################
rs.b.e.11.5 <- matrix(0, 2, 6)
rs.b.e.11.5[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 5, 0, c(0.1, 0.1), ex, 0)
rs.b.e.11.5[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 5, 0.2, c(0.1, 0.1), ex, 0)
rs.b.e.11.5[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 5, 0.5, c(0.1, 0.1), ex, 0)
rs.b.e.11.5[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 5, 0, c(0.1, 0.1), ex, 0)
rs.b.e.11.5[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 5, 0.2, c(0.1, 0.1), ex, 0)
rs.b.e.11.5[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 5, 0.5, c(0.1, 0.1), ex, 0)

####################################################################
### Rank sum, exchangable, rho = 0.5, 0.5, clus size = 5, delta = 0,
### 0.2, 0.5, grpsize = 20, 50, balanced data
####################################################################
rs.b.e.55.5 <- matrix(0, 2, 6)
rs.b.e.55.5[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 5, 0, c(0.5, 0.5), ex, 0)
rs.b.e.55.5[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 5, 0.2, c(0.5, 0.5), ex, 0)
rs.b.e.55.5[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 5, 0.5, c(0.5, 0.5), ex, 0)
rs.b.e.55.5[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 5, 0, c(0.5, 0.5), ex, 0)
rs.b.e.55.5[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 5, 0.2, c(0.5, 0.5), ex, 0)
rs.b.e.55.5[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 5, 0.5, c(0.5, 0.5), ex, 0)

####################################################################
### Rank sum, exchangable, rho = -0.1, 0.9, clus size = 5, delta = 0,
### 0.2, 0.5, grpsize = 20, 50, balanced data
####################################################################
rs.b.e.19.5 <- matrix(0, 2, 6)
rs.b.e.19.5[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 5, 0, c(-0.1, 0.9), ex, 0)
rs.b.e.19.5[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 5, 0.2, c(-0.1, 0.9), ex, 0)
rs.b.e.19.5[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 5, 0.5, c(-0.1, 0.9), ex, 0)
rs.b.e.19.5[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 5, 0, c(-0.1, 0.9), ex, 0)
rs.b.e.19.5[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 5, 0.2, c(-0.1, 0.9), ex, 0)
rs.b.e.19.5[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 5, 0.5, c(-0.1, 0.9), ex, 0)

####################################################################
### Rank sum, exchangable, rho = 0.1, 0.1, clus size = 10, delta = 0,
### 0.2, 0.5, grpsize = 20, 50, unbalanced data, missingrate = 0.5
####################################################################
rs.ub.e.11.10 <- matrix(0, 2, 6)
rs.ub.e.11.10[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 10, 0, c(0.1, 0.1), ex, 0.5)
rs.ub.e.11.10[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 10, 0.2, c(0.1, 0.1), ex, 0.5)
rs.ub.e.11.10[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 10, 0.5, c(0.1, 0.1), ex, 0.5)
rs.ub.e.11.10[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 10, 0, c(0.1, 0.1), ex, 0.5)
rs.ub.e.11.10[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 10, 0.2, c(0.1, 0.1), ex, 0.5)
rs.ub.e.11.10[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 10, 0.5, c(0.1, 0.1), ex, 0.5)

####################################################################
### Rank sum, exchangable, rho = 0.5, 0.5, clus size = 10, delta = 0,
### 0.2, 0.5, grpsize = 20, 50, unbalanced data, missingrate = 0.5
####################################################################
rs.ub.e.55.10 <- matrix(0, 2, 6)
rs.ub.e.55.10[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 10, 0, c(0.5, 0.5), ex, 0.5)
rs.ub.e.55.10[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 10, 0.2, c(0.5, 0.5), ex, 0.5)
rs.ub.e.55.10[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 10, 0.5, c(0.5, 0.5), ex, 0.5)
rs.ub.e.55.10[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 10, 0, c(0.5, 0.5), ex, 0.5)
rs.ub.e.55.10[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 10, 0.2, c(0.5, 0.5), ex, 0.5)
rs.ub.e.55.10[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 10, 0.5, c(0.5, 0.5), ex, 0.5)

####################################################################
### Rank sum, exchangable, rho = -0.1, 0.9, clus size = 10, delta = 0,
### 0.2, 0.5, grpsize = 20, 50, unbalanced data, missingrate = 0.5
####################################################################
rs.ub.e.19.10 <- matrix(0, 2, 6)
rs.ub.e.19.10[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 10, 0, c(-0.1, 0.9), ex, 0.5)
rs.ub.e.19.10[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 10, 0.2, c(-0.1, 0.9), ex, 0.5)
rs.ub.e.19.10[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 10, 0.5, c(-0.1, 0.9), ex, 0.5)
rs.ub.e.19.10[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 10, 0, c(-0.1, 0.9), ex, 0.5)
rs.ub.e.19.10[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 10, 0.2, c(-0.1, 0.9), ex, 0.5)
rs.ub.e.19.10[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 10, 0.5, c(-0.1, 0.9), ex, 0.5)

####################################################################
### Rank sum, sub-unit grouping exchangable, rho = 0.1, 0.1, clus size
### = 2, delta = 0, 0.2, 0.5, grpsize = 20, 50
####################################################################
rs.b.e.11.2.sub <- matrix(0, 2, 6)
rs.b.e.11.2.sub[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 2, 0, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.2.sub[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 2, 0.2, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.2.sub[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 2, 0.5, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.2.sub[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 2, 0, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.2.sub[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 2, 0.2, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.2.sub[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 2, 0.5, c(0.1, 0.1), ex, 0, FALSE)

####################################################################
### Rank sum, sub-unit grouping, exchangable, rho = 0.5, 0.5, clus size
### = 2, delta = 0, 0.2, 0.5, grpsize = 20, 50
####################################################################
rs.b.e.55.2.sub <- matrix(0, 2, 6)
rs.b.e.55.2.sub[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 2, 0, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.2.sub[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 2, 0.2, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.2.sub[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 2, 0.5, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.2.sub[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 2, 0, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.2.sub[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 2, 0.2, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.2.sub[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 2, 0.5, c(0.5, 0.5), ex, 0, FALSE)

####################################################################
### Rank sum, sub-unit grouping exchangable, rho = -0.1, 0.9, clus size
### = 2, delta = 0, 0.2, 0.5, grpsize = 20, 50
####################################################################
rs.b.e.19.2.sub <- matrix(0, 2, 6)
rs.b.e.19.2.sub[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 2, 0, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.2.sub[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 2, 0.2, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.2.sub[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 2, 0.5, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.2.sub[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 2, 0, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.2.sub[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 2, 0.2, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.2.sub[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 2, 0.5, c(-0.1, 0.9), ex, 0, FALSE)

####################################################################
### Rank sum, sub-unit grouping exchangable, rho = 0.1, 0.1, clus size
### = 5, delta = 0, 0.2, 0.5, grpsize = 20, 50
####################################################################
rs.b.e.11.5.sub <- matrix(0, 2, 6)
rs.b.e.11.5.sub[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 5, 0, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.5.sub[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 5, 0.2, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.5.sub[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 5, 0.5, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.5.sub[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 5, 0, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.5.sub[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 5, 0.2, c(0.1, 0.1), ex, 0, FALSE)
rs.b.e.11.5.sub[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 5, 0.5, c(0.1, 0.1), ex, 0, FALSE)

####################################################################
### Rank sum, sub-unit grouping exchangable, rho = 0.5, 0.5, clus size
### = 5, delta = 0, 0.2, 0.5, grpsize = 20, 50
####################################################################
rs.b.e.55.5.sub <- matrix(0, 2, 6)
rs.b.e.55.5.sub[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 5, 0, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.5.sub[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 5, 0.2, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.5.sub[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 5, 0.5, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.5.sub[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 5, 0, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.5.sub[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 5, 0.2, c(0.5, 0.5), ex, 0, FALSE)
rs.b.e.55.5.sub[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 5, 0.5, c(0.5, 0.5), ex, 0, FALSE)

####################################################################
### Rank sum, sub-unit grouping exchangable, rho = -0.1, 0.9, clus size
### = 5, delta = 0, 0.2, 0.5, grpsize = 20, 50
####################################################################
rs.b.e.19.5.sub <- matrix(0, 2, 6)
rs.b.e.19.5.sub[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 5, 0, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.5.sub[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 5, 0.2, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.5.sub[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 5, 0.5, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.5.sub[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 5, 0, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.5.sub[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 5, 0.2, c(-0.1, 0.9), ex, 0, FALSE)
rs.b.e.19.5.sub[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 5, 0.5, c(-0.1, 0.9), ex, 0, FALSE)

####################################################################
### Rank sum, sub-unit grouping exchangable, rho = 0.1, 0.1, clus size
### = 10, delta = 0, 0.2, 0.5, grpsize = 20, 50, unbalanced data,
### missingrate = 0.5
####################################################################
rs.ub.e.11.10.sub <- matrix(0, 2, 6)
rs.ub.e.11.10.sub[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 10, 0, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.e.11.10.sub[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 10, 0.2, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.e.11.10.sub[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 10, 0.5, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.e.11.10.sub[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 10, 0, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.e.11.10.sub[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 10, 0.2, c(0.1, 0.1), ex, 0.5, FALSE)
rs.ub.e.11.10.sub[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 10, 0.5, c(0.1, 0.1), ex, 0.5, FALSE)
####################################################################
### Rank sum, sub-unit grouping exchangable, rho = 0.5, 0.5, clus size
### = 10, delta = 0, 0.2, 0.5, grpsize = 20, 50, unbalanced data,
### missingrate = 0.5
####################################################################

rs.ub.e.55.10.sub <- matrix(0, 2, 6)
rs.ub.e.55.10.sub[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 10, 0, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.e.55.10.sub[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 10, 0.2, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.e.55.10.sub[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 10, 0.5, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.e.55.10.sub[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 10, 0, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.e.55.10.sub[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 10, 0.2, c(0.5, 0.5), ex, 0.5, FALSE)
rs.ub.e.55.10.sub[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 10, 0.5, c(0.5, 0.5), ex, 0.5, FALSE)
####################################################################
### Rank sum, sub-unit grouping exchangable, rho = -0.1, 0.9, clus size
### = 10, delta = 0, 0.2, 0.5, grpsize = 20, 50, unbalanced data,
### missingrate = 0.5
####################################################################

rs.ub.e.19.10.sub <- matrix(0, 2, 6)
rs.ub.e.19.10.sub[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 10, 0, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.e.19.10.sub[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 10, 0.2, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.e.19.10.sub[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 10, 0.5, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.e.19.10.sub[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 10, 0, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.e.19.10.sub[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 10, 0.2, c(-0.1, 0.9), ex, 0.5, FALSE)
rs.ub.e.19.10.sub[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 10, 0.5, c(-0.1, 0.9), ex, 0.5, FALSE)





####################################################################
### Rank sum, ar1, rho = 0.1, 0.1, clus size = 2, delta = 0,
### 0.2, 0.5, grpsize = 20, 50, balanced data
####################################################################

rs.b.a.11.2 <- matrix(0, 2, 6)
rs.b.a.11.2[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 2, 0, c(0.1, 0.1), ar1, 0)
rs.b.a.11.2[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 2, 0.2, c(0.1, 0.1), ar1, 0)
rs.b.a.11.2[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 2, 0.5, c(0.1, 0.1), ar1, 0)
rs.b.a.11.2[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 2, 0, c(0.1, 0.1), ar1, 0)
rs.b.a.11.2[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 2, 0.2, c(0.1, 0.1), ar1, 0)
rs.b.a.11.2[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 2, 0.5, c(0.1, 0.1), ar1, 0)

####################################################################
### Rank sum, ar1, rho = 0.5, 0.5, clus size = 2, delta = 0,
### 0.2, 0.5, grpsize = 20, 50, balanced data
####################################################################
rs.b.a.55.2 <- matrix(0, 2, 6)
rs.b.a.55.2[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 2, 0, c(0.5, 0.5), ar1, 0)
rs.b.a.55.2[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 2, 0.2, c(0.5, 0.5), ar1, 0)
rs.b.a.55.2[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 2, 0.5, c(0.5, 0.5), ar1, 0)
rs.b.a.55.2[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 2, 0, c(0.5, 0.5), ar1, 0)
rs.b.a.55.2[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 2, 0.2, c(0.5, 0.5), ar1, 0)
rs.b.a.55.2[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 2, 0.5, c(0.5, 0.5), ar1, 0)

####################################################################
### Rank sum, ar1, rho = -0.1, 0.9, clus size = 2, delta = 0,
### 0.2, 0.5, grpsize = 20, 50, balanced data
####################################################################
rs.b.a.19.2 <- matrix(0, 2, 6)
rs.b.a.19.2[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 2, 0, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.2[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 2, 0.2, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.2[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 2, 0.5, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.2[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 2, 0, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.2[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 2, 0.2, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.2[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 2, 0.5, c(-0.1, 0.9), ar1, 0)

####################################################################
### Rank sum, ar1, rho = 0.1, 0.1, clus size = 5, delta = 0,
### 0.2, 0.5, grpsize = 20, 50, balanced data
####################################################################
rs.b.a.11.5 <- matrix(0, 2, 6)
rs.b.a.11.5[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 5, 0, c(0.1, 0.1), ar1, 0)
rs.b.a.11.5[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 5, 0.2, c(0.1, 0.1), ar1, 0)
rs.b.a.11.5[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 5, 0.5, c(0.1, 0.1), ar1, 0)
rs.b.a.11.5[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 5, 0, c(0.1, 0.1), ar1, 0)
rs.b.a.11.5[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 5, 0.2, c(0.1, 0.1), ar1, 0)
rs.b.a.11.5[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 5, 0.5, c(0.1, 0.1), ar1, 0)

####################################################################
### Rank sum, ar1, rho = 0.5, 0.5, clus size = 5, delta = 0,
### 0.2, 0.5, grpsize = 20, 50, balanced data
####################################################################
rs.b.a.55.5 <- matrix(0, 2, 6)
rs.b.a.55.5[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 5, 0, c(0.5, 0.5), ar1, 0)
rs.b.a.55.5[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 5, 0.2, c(0.5, 0.5), ar1, 0)
rs.b.a.55.5[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 5, 0.5, c(0.5, 0.5), ar1, 0)
rs.b.a.55.5[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 5, 0, c(0.5, 0.5), ar1, 0)
rs.b.a.55.5[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 5, 0.2, c(0.5, 0.5), ar1, 0)
rs.b.a.55.5[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 5, 0.5, c(0.5, 0.5), ar1, 0)

####################################################################
### Rank sum, ar1, rho = -0.1, 0.9, clus size = 5, delta = 0,
### 0.2, 0.5, grpsize = 20, 50, balanced data
####################################################################
rs.b.a.19.5 <- matrix(0, 2, 6)
rs.b.a.19.5[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 5, 0, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.5[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 5, 0.2, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.5[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 5, 0.5, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.5[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 5, 0, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.5[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 5, 0.2, c(-0.1, 0.9), ar1, 0)
rs.b.a.19.5[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 5, 0.5, c(-0.1, 0.9), ar1, 0)

####################################################################
### Rank sum, ar1, rho = 0.1, 0.1, clus size = 10, delta = 0,
### 0.2, 0.5, grpsize = 20, 50, unbalanced data, missingrate = 0.5
####################################################################
rs.ub.a.11.10 <- matrix(0, 2, 6)
rs.ub.a.11.10[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 10, 0, c(0.1, 0.1), ar1, 0.5)
rs.ub.a.11.10[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 10, 0.2, c(0.1, 0.1), ar1, 0.5)
rs.ub.a.11.10[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 10, 0.5, c(0.1, 0.1), ar1, 0.5)
rs.ub.a.11.10[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 10, 0, c(0.1, 0.1), ar1, 0.5)
rs.ub.a.11.10[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 10, 0.2, c(0.1, 0.1), ar1, 0.5)
rs.ub.a.11.10[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 10, 0.5, c(0.1, 0.1), ar1, 0.5)

####################################################################
### Rank sum, ar1, rho = 0.5, 0.5, clus size = 10, delta = 0,
### 0.2, 0.5, grpsize = 20, 50, unbalanced data, missingrate = 0.5
####################################################################
rs.ub.a.55.10 <- matrix(0, 2, 6)
rs.ub.a.55.10[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 10, 0, c(0.5, 0.5), ar1, 0.5)
rs.ub.a.55.10[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 10, 0.2, c(0.5, 0.5), ar1, 0.5)
rs.ub.a.55.10[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 10, 0.5, c(0.5, 0.5), ar1, 0.5)
rs.ub.a.55.10[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 10, 0, c(0.5, 0.5), ar1, 0.5)
rs.ub.a.55.10[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 10, 0.2, c(0.5, 0.5), ar1, 0.5)
rs.ub.a.55.10[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 10, 0.5, c(0.5, 0.5), ar1, 0.5)

####################################################################
### Rank sum, ar1, rho = -0.1, 0.9, clus size = 10, delta = 0,
### 0.2, 0.5, grpsize = 20, 50, unbalanced data, missingrate = 0.5
####################################################################
rs.ub.a.19.10 <- matrix(0, 2, 6)
rs.ub.a.19.10[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 10, 0, c(-0.1, 0.9), ar1, 0.5)
rs.ub.a.19.10[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 10, 0.2, c(-0.1, 0.9), ar1, 0.5)
rs.ub.a.19.10[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 10, 0.5, c(-0.1, 0.9), ar1, 0.5)
rs.ub.a.19.10[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 10, 0, c(-0.1, 0.9), ar1, 0.5)
rs.ub.a.19.10[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 10, 0.2, c(-0.1, 0.9), ar1, 0.5)
rs.ub.a.19.10[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 10, 0.5, c(-0.1, 0.9), ar1, 0.5)

####################################################################
### Rank sum, sub-unit grouping ar1, rho = 0.1, 0.1, clus size
### = 2, delta = 0, 0.2, 0.5, grpsize = 20, 50
####################################################################
rs.b.a.11.2.sub <- matrix(0, 2, 6)
rs.b.a.11.2.sub[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 2, 0, c(0.1, 0.1), ar1, 0, FALSE)
rs.b.a.11.2.sub[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 2, 0.2, c(0.1, 0.1), ar1, 0, FALSE)
rs.b.a.11.2.sub[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 2, 0.5, c(0.1, 0.1), ar1, 0, FALSE)
rs.b.a.11.2.sub[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 2, 0, c(0.1, 0.1), ar1, 0, FALSE)
rs.b.a.11.2.sub[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 2, 0.2, c(0.1, 0.1), ar1, 0, FALSE)
rs.b.a.11.2.sub[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 2, 0.5, c(0.1, 0.1), ar1, 0, FALSE)

####################################################################
### Rank sum, sub-unit grouping, ar1, rho = 0.5, 0.5, clus size
### = 2, delta = 0, 0.2, 0.5, grpsize = 20, 50
####################################################################
rs.b.a.55.2.sub <- matrix(0, 2, 6)
rs.b.a.55.2.sub[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 2, 0, c(0.5, 0.5), ar1, 0, FALSE)
rs.b.a.55.2.sub[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 2, 0.2, c(0.5, 0.5), ar1, 0, FALSE)
rs.b.a.55.2.sub[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 2, 0.5, c(0.5, 0.5), ar1, 0, FALSE)
rs.b.a.55.2.sub[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 2, 0, c(0.5, 0.5), ar1, 0, FALSE)
rs.b.a.55.2.sub[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 2, 0.2, c(0.5, 0.5), ar1, 0, FALSE)
rs.b.a.55.2.sub[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 2, 0.5, c(0.5, 0.5), ar1, 0, FALSE)

####################################################################
### Rank sum, sub-unit grouping ar1, rho = -0.1, 0.9, clus size
### = 2, delta = 0, 0.2, 0.5, grpsize = 20, 50
####################################################################
rs.b.a.19.2.sub <- matrix(0, 2, 6)
rs.b.a.19.2.sub[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 2, 0, c(-0.1, 0.9), ar1, 0, FALSE)
rs.b.a.19.2.sub[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 2, 0.2, c(-0.1, 0.9), ar1, 0, FALSE)
rs.b.a.19.2.sub[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 2, 0.5, c(-0.1, 0.9), ar1, 0, FALSE)
rs.b.a.19.2.sub[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 2, 0, c(-0.1, 0.9), ar1, 0, FALSE)
rs.b.a.19.2.sub[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 2, 0.2, c(-0.1, 0.9), ar1, 0, FALSE)
rs.b.a.19.2.sub[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 2, 0.5, c(-0.1, 0.9), ar1, 0, FALSE)

####################################################################
### Rank sum, sub-unit grouping ar1, rho = 0.1, 0.1, clus size
### = 5, delta = 0, 0.2, 0.5, grpsize = 20, 50
####################################################################
rs.b.a.11.5.sub <- matrix(0, 2, 6)
rs.b.a.11.5.sub[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 5, 0, c(0.1, 0.1), ar1, 0, FALSE)
rs.b.a.11.5.sub[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 5, 0.2, c(0.1, 0.1), ar1, 0, FALSE)
rs.b.a.11.5.sub[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 5, 0.5, c(0.1, 0.1), ar1, 0, FALSE)
rs.b.a.11.5.sub[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 5, 0, c(0.1, 0.1), ar1, 0, FALSE)
rs.b.a.11.5.sub[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 5, 0.2, c(0.1, 0.1), ar1, 0, FALSE)
rs.b.a.11.5.sub[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 5, 0.5, c(0.1, 0.1), ar1, 0, FALSE)

####################################################################
### Rank sum, sub-unit grouping ar1, rho = 0.5, 0.5, clus size
### = 5, delta = 0, 0.2, 0.5, grpsize = 20, 50
####################################################################
rs.b.a.55.5.sub <- matrix(0, 2, 6)
rs.b.a.55.5.sub[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 5, 0, c(0.5, 0.5), ar1, 0, FALSE)
rs.b.a.55.5.sub[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 5, 0.2, c(0.5, 0.5), ar1, 0, FALSE)
rs.b.a.55.5.sub[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 5, 0.5, c(0.5, 0.5), ar1, 0, FALSE)
rs.b.a.55.5.sub[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 5, 0, c(0.5, 0.5), ar1, 0, FALSE)
rs.b.a.55.5.sub[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 5, 0.2, c(0.5, 0.5), ar1, 0, FALSE)
rs.b.a.55.5.sub[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 5, 0.5, c(0.5, 0.5), ar1, 0, FALSE)

####################################################################
### Rank sum, sub-unit grouping ar1, rho = -0.1, 0.9, clus size
### = 5, delta = 0, 0.2, 0.5, grpsize = 20, 50
####################################################################
rs.b.a.19.5.sub <- matrix(0, 2, 6)
rs.b.a.19.5.sub[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 5, 0, c(-0.1, 0.9), ar1, 0, FALSE)
rs.b.a.19.5.sub[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 5, 0.2, c(-0.1, 0.9), ar1, 0, FALSE)
rs.b.a.19.5.sub[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 5, 0.5, c(-0.1, 0.9), ar1, 0, FALSE)
rs.b.a.19.5.sub[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 5, 0, c(-0.1, 0.9), ar1, 0, FALSE)
rs.b.a.19.5.sub[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 5, 0.2, c(-0.1, 0.9), ar1, 0, FALSE)
rs.b.a.19.5.sub[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 5, 0.5, c(-0.1, 0.9), ar1, 0, FALSE)

####################################################################
### Rank sum, sub-unit grouping ar1, rho = 0.1, 0.1, clus size
### = 10, delta = 0, 0.2, 0.5, grpsize = 20, 50, unbalanced data,
### missingrate = 0.5
####################################################################
rs.ub.a.11.10.sub <- matrix(0, 2, 6)
rs.ub.a.11.10.sub[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 10, 0, c(0.1, 0.1), ar1, 0.5, FALSE)
rs.ub.a.11.10.sub[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 10, 0.2, c(0.1, 0.1), ar1, 0.5, FALSE)
rs.ub.a.11.10.sub[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 10, 0.5, c(0.1, 0.1), ar1, 0.5, FALSE)
rs.ub.a.11.10.sub[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 10, 0, c(0.1, 0.1), ar1, 0.5, FALSE)
rs.ub.a.11.10.sub[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 10, 0.2, c(0.1, 0.1), ar1, 0.5, FALSE)
rs.ub.a.11.10.sub[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 10, 0.5, c(0.1, 0.1), ar1, 0.5, FALSE)
####################################################################
### Rank sum, sub-unit grouping ar1, rho = 0.5, 0.5, clus size
### = 10, delta = 0, 0.2, 0.5, grpsize = 20, 50, unbalanced data,
### missingrate = 0.5
####################################################################

rs.ub.a.55.10.sub <- matrix(0, 2, 6)
rs.ub.a.55.10.sub[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 10, 0, c(0.5, 0.5), ar1, 0.5, FALSE)
rs.ub.a.55.10.sub[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 10, 0.2, c(0.5, 0.5), ar1, 0.5, FALSE)
rs.ub.a.55.10.sub[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 10, 0.5, c(0.5, 0.5), ar1, 0.5, FALSE)
rs.ub.a.55.10.sub[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 10, 0, c(0.5, 0.5), ar1, 0.5, FALSE)
rs.ub.a.55.10.sub[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 10, 0.2, c(0.5, 0.5), ar1, 0.5, FALSE)
rs.ub.a.55.10.sub[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 10, 0.5, c(0.5, 0.5), ar1, 0.5, FALSE)
####################################################################
### Rank sum, sub-unit grouping ar1, rho = -0.1, 0.9, clus size
### = 10, delta = 0, 0.2, 0.5, grpsize = 20, 50, unbalanced data,
### missingrate = 0.5
####################################################################

rs.ub.a.19.10.sub <- matrix(0, 2, 6)
rs.ub.a.19.10.sub[1, 1:2] <- simpower(nrep, 0.05, FALSE, 20, 10, 0, c(-0.1, 0.9), ar1, 0.5, FALSE)
rs.ub.a.19.10.sub[1, 3:4] <- simpower(nrep, 0.05, FALSE, 20, 10, 0.2, c(-0.1, 0.9), ar1, 0.5, FALSE)
rs.ub.a.19.10.sub[1, 5:6] <- simpower(nrep, 0.05, FALSE, 20, 10, 0.5, c(-0.1, 0.9), ar1, 0.5, FALSE)
rs.ub.a.19.10.sub[2, 1:2] <- simpower(nrep, 0.05, FALSE, 50, 10, 0, c(-0.1, 0.9), ar1, 0.5, FALSE)
rs.ub.a.19.10.sub[2, 3:4] <- simpower(nrep, 0.05, FALSE, 50, 10, 0.2, c(-0.1, 0.9), ar1, 0.5, FALSE)
rs.ub.a.19.10.sub[2, 5:6] <- simpower(nrep, 0.05, FALSE, 50, 10, 0.5, c(-0.1, 0.9), ar1, 0.5, FALSE)







####################################################################
### Signed-rank, exchangable, rho = 0.1, clus size
### = 2, delta = 0, 0.2, 0.5, grpsize = 20, 50
####################################################################
sr.b.e.1.2 <- matrix(0, 2, 6)
sr.b.e.1.2[1, 1:2] <- simpower(nrep, 0.05, TRUE, 20, 2, 0, 0.1, ex, 0)
sr.b.e.1.2[1, 3:4] <- simpower(nrep, 0.05, TRUE, 20, 2, 0.2, 0.1, ex, 0)
sr.b.e.1.2[1, 5:6] <- simpower(nrep, 0.05, TRUE, 20, 2, 0.5, 0.1, ex, 0)
sr.b.e.1.2[2, 1:2] <- simpower(nrep, 0.05, TRUE, 50, 2, 0, 0.1, ex, 0)
sr.b.e.1.2[2, 3:4] <- simpower(nrep, 0.05, TRUE, 50, 2, 0.2, 0.1, ex, 0)
sr.b.e.1.2[2, 5:6] <- simpower(nrep, 0.05, TRUE, 50, 2, 0.5, 0.1, ex, 0)

####################################################################
### Signed-rank, exchangable, rho = 0.5, clus size
### = 2, delta = 0, 0.2, 0.5, grpsize = 20, 50
####################################################################
sr.b.e.5.2 <- matrix(0, 2, 6)
sr.b.e.5.2[1, 1:2] <- simpower(nrep, 0.05, TRUE, 20, 2, 0, 0.5, ex, 0)
sr.b.e.5.2[1, 3:4] <- simpower(nrep, 0.05, TRUE, 20, 2, 0.2, 0.5, ex, 0)
sr.b.e.5.2[1, 5:6] <- simpower(nrep, 0.05, TRUE, 20, 2, 0.5, 0.5, ex, 0)
sr.b.e.5.2[2, 1:2] <- simpower(nrep, 0.05, TRUE, 50, 2, 0, 0.5, ex, 0)
sr.b.e.5.2[2, 3:4] <- simpower(nrep, 0.05, TRUE, 50, 2, 0.2, 0.5, ex, 0)
sr.b.e.5.2[2, 5:6] <- simpower(nrep, 0.05, TRUE, 50, 2, 0.5, 0.5, ex, 0)

####################################################################
### Signed-rank, exchangable, rho = 0.9, clus size
### = 2, delta = 0, 0.2, 0.5, grpsize = 20, 50
####################################################################
sr.b.e.9.2 <- matrix(0, 2, 6)
sr.b.e.9.2[1, 1:2] <- simpower(nrep, 0.05, TRUE, 20, 2, 0, 0.9, ex, 0)
sr.b.e.9.2[1, 3:4] <- simpower(nrep, 0.05, TRUE, 20, 2, 0.2, 0.9, ex, 0)
sr.b.e.9.2[1, 5:6] <- simpower(nrep, 0.05, TRUE, 20, 2, 0.5, 0.9, ex, 0)
sr.b.e.9.2[2, 1:2] <- simpower(nrep, 0.05, TRUE, 50, 2, 0, 0.9, ex, 0)
sr.b.e.9.2[2, 3:4] <- simpower(nrep, 0.05, TRUE, 50, 2, 0.2, 0.9, ex, 0)
sr.b.e.9.2[2, 5:6] <- simpower(nrep, 0.05, TRUE, 50, 2, 0.5, 0.9, ex, 0)

####################################################################
### Signed-rank, exchangable, rho = 0.1, clus size
### = 10, delta = 0, 0.2, 0.5, grpsize = 20, 50
####################################################################
sr.b.e.1.10 <- matrix(0, 2, 6)
sr.b.e.1.10[1, 1:2] <- simpower(nrep, 0.05, TRUE, 20, 10, 0, 0.1, ex, 0)
sr.b.e.1.10[1, 3:4] <- simpower(nrep, 0.05, TRUE, 20, 10, 0.2, 0.1, ex, 0)
sr.b.e.1.10[1, 5:6] <- simpower(nrep, 0.05, TRUE, 20, 10, 0.5, 0.1, ex, 0)
sr.b.e.1.10[2, 1:2] <- simpower(nrep, 0.05, TRUE, 50, 10, 0, 0.1, ex, 0)
sr.b.e.1.10[2, 3:4] <- simpower(nrep, 0.05, TRUE, 50, 10, 0.2, 0.1, ex, 0)
sr.b.e.1.10[2, 5:6] <- simpower(nrep, 0.05, TRUE, 50, 10, 0.5, 0.1, ex, 0)

####################################################################
### Signed-rank, exchangable, rho = 0.5, clus size
### = 10, delta = 0, 0.2, 0.5, grpsize = 20, 50
####################################################################
sr.b.e.5.10 <- matrix(0, 2, 6)
sr.b.e.5.10[1, 1:2] <- simpower(nrep, 0.05, TRUE, 20, 10, 0, 0.5, ex, 0)
sr.b.e.5.10[1, 3:4] <- simpower(nrep, 0.05, TRUE, 20, 10, 0.2, 0.5, ex, 0)
sr.b.e.5.10[1, 5:6] <- simpower(nrep, 0.05, TRUE, 20, 10, 0.5, 0.5, ex, 0)
sr.b.e.5.10[2, 1:2] <- simpower(nrep, 0.05, TRUE, 50, 10, 0, 0.5, ex, 0)
sr.b.e.5.10[2, 3:4] <- simpower(nrep, 0.05, TRUE, 50, 10, 0.2, 0.5, ex, 0)
sr.b.e.5.10[2, 5:6] <- simpower(nrep, 0.05, TRUE, 50, 10, 0.5, 0.5, ex, 0)
####################################################################
### Signed-rank, exchangable, rho = 0.9, clus size
### = 10, delta = 0, 0.2, 0.5, grpsize = 20, 50
####################################################################
sr.b.e.9.10 <- matrix(0, 2, 6)
sr.b.e.9.10[1, 1:2] <- simpower(nrep, 0.05, TRUE, 20, 10, 0, 0.9, ex, 0)
sr.b.e.9.10[1, 3:4] <- simpower(nrep, 0.05, TRUE, 20, 10, 0.2, 0.9, ex, 0)
sr.b.e.9.10[1, 5:6] <- simpower(nrep, 0.05, TRUE, 20, 10, 0.5, 0.9, ex, 0)
sr.b.e.9.10[2, 1:2] <- simpower(nrep, 0.05, TRUE, 50, 10, 0, 0.9, ex, 0)
sr.b.e.9.10[2, 3:4] <- simpower(nrep, 0.05, TRUE, 50, 10, 0.2, 0.9, ex, 0)
sr.b.e.9.10[2, 5:6] <- simpower(nrep, 0.05, TRUE, 50, 10, 0.5, 0.9, ex, 0)

####################################################################
### Signed-rank, exchangable, rho = 0.1, clus size
### = 5, delta = 0, 0.2, 0.5, grpsize = 20, 50, missingrate = 0.5
####################################################################
sr.ub.e.1.5 <- matrix(0, 2, 6)
sr.ub.e.1.5[1, 1:2] <- simpower(nrep, 0.05, TRUE, 20, 5, 0, 0.1, ex, 0.5)
sr.ub.e.1.5[1, 3:4] <- simpower(nrep, 0.05, TRUE, 20, 5, 0.2, 0.1, ex, 0.5)
sr.ub.e.1.5[1, 5:6] <- simpower(nrep, 0.05, TRUE, 20, 5, 0.5, 0.1, ex, 0.5)
sr.ub.e.1.5[2, 1:2] <- simpower(nrep, 0.05, TRUE, 50, 5, 0, 0.1, ex, 0.5)
sr.ub.e.1.5[2, 3:4] <- simpower(nrep, 0.05, TRUE, 50, 5, 0.2, 0.1, ex, 0.5)
sr.ub.e.1.5[2, 5:6] <- simpower(nrep, 0.05, TRUE, 50, 5, 0.5, 0.1, ex, 0.5)
####################################################################
### Signed-rank, exchangable, rho = 0.5, clus size
### = 5, delta = 0, 0.2, 0.5, grpsize = 20, 50, missingrate = 0.5
####################################################################
sr.ub.e.5.5 <- matrix(0, 2, 6)
sr.ub.e.5.5[1, 1:2] <- simpower(nrep, 0.05, TRUE, 20, 5, 0, 0.5, ex, 0.5)
sr.ub.e.5.5[1, 3:4] <- simpower(nrep, 0.05, TRUE, 20, 5, 0.2, 0.5, ex, 0.5)
sr.ub.e.5.5[1, 5:6] <- simpower(nrep, 0.05, TRUE, 20, 5, 0.5, 0.5, ex, 0.5)
sr.ub.e.5.5[2, 1:2] <- simpower(nrep, 0.05, TRUE, 50, 5, 0, 0.5, ex, 0.5)
sr.ub.e.5.5[2, 3:4] <- simpower(nrep, 0.05, TRUE, 50, 5, 0.2, 0.5, ex, 0.5)
sr.ub.e.5.5[2, 5:6] <- simpower(nrep, 0.05, TRUE, 50, 5, 0.5, 0.5, ex, 0.5)
####################################################################
### Signed-rank, exchangable, rho = 0.9, clus size
### = 5, delta = 0, 0.2, 0.5, grpsize = 20, 50, missingrate = 0.5
####################################################################
sr.ub.e.9.5 <- matrix(0, 2, 6)
sr.ub.e.9.5[1, 1:2] <- simpower(nrep, 0.05, TRUE, 20, 5, 0, 0.9, ex, 0.5)
sr.ub.e.9.5[1, 3:4] <- simpower(nrep, 0.05, TRUE, 20, 5, 0.2, 0.9, ex, 0.5)
sr.ub.e.9.5[1, 5:6] <- simpower(nrep, 0.05, TRUE, 20, 5, 0.5, 0.9, ex, 0.5)
sr.ub.e.9.5[2, 1:2] <- simpower(nrep, 0.05, TRUE, 50, 5, 0, 0.9, ex, 0.5)
sr.ub.e.9.5[2, 3:4] <- simpower(nrep, 0.05, TRUE, 50, 5, 0.2, 0.9, ex, 0.5)
sr.ub.e.9.5[2, 5:6] <- simpower(nrep, 0.05, TRUE, 50, 5, 0.5, 0.9, ex, 0.5)
####################################################################
### Signed-rank, exchangable, rho = 0.1, clus size
### = 10, delta = 0, 0.2, 0.5, grpsize = 20, 50, missingrate = 0.5
####################################################################
sr.ub.e.1.10 <- matrix(0, 2, 6)
sr.ub.e.1.10[1, 1:2] <- simpower(nrep, 0.05, TRUE, 20, 10, 0, 0.1, ex, 0.5)
sr.ub.e.1.10[1, 3:4] <- simpower(nrep, 0.05, TRUE, 20, 10, 0.2, 0.1, ex, 0.5)
sr.ub.e.1.10[1, 5:6] <- simpower(nrep, 0.05, TRUE, 20, 10, 0.5, 0.1, ex, 0.5)
sr.ub.e.1.10[2, 1:2] <- simpower(nrep, 0.05, TRUE, 50, 10, 0, 0.1, ex, 0.5)
sr.ub.e.1.10[2, 3:4] <- simpower(nrep, 0.05, TRUE, 50, 10, 0.2, 0.1, ex, 0.5)
sr.ub.e.1.10[2, 5:6] <- simpower(nrep, 0.05, TRUE, 50, 10, 0.5, 0.1, ex, 0.5)
####################################################################
### Signed-rank, exchangable, rho = 0.5, clus size
### = 10, delta = 0, 0.2, 0.5, grpsize = 20, 50, missingrate = 0.5
####################################################################
sr.ub.e.5.10 <- matrix(0, 2, 6)
sr.ub.e.5.10[1, 1:2] <- simpower(nrep, 0.05, TRUE, 20, 10, 0, 0.5, ex, 0.5)
sr.ub.e.5.10[1, 3:4] <- simpower(nrep, 0.05, TRUE, 20, 10, 0.2, 0.5, ex, 0.5)
sr.ub.e.5.10[1, 5:6] <- simpower(nrep, 0.05, TRUE, 20, 10, 0.5, 0.5, ex, 0.5)
sr.ub.e.5.10[2, 1:2] <- simpower(nrep, 0.05, TRUE, 50, 10, 0, 0.5, ex, 0.5)
sr.ub.e.5.10[2, 3:4] <- simpower(nrep, 0.05, TRUE, 50, 10, 0.2, 0.5, ex, 0.5)
sr.ub.e.5.10[2, 5:6] <- simpower(nrep, 0.05, TRUE, 50, 10, 0.5, 0.5, ex, 0.5)
####################################################################
### Signed-rank, exchangable, rho = 0.9, clus size
### = 10, delta = 0, 0.2, 0.5, grpsize = 20, 50, missingrate = 0.5
####################################################################
sr.ub.e.9.10 <- matrix(0, 2, 6)
sr.ub.e.9.10[1, 1:2] <- simpower(nrep, 0.05, TRUE, 20, 10, 0, 0.9, ex, 0.5)
sr.ub.e.9.10[1, 3:4] <- simpower(nrep, 0.05, TRUE, 20, 10, 0.2, 0.9, ex, 0.5)
sr.ub.e.9.10[1, 5:6] <- simpower(nrep, 0.05, TRUE, 20, 10, 0.5, 0.9, ex, 0.5)
sr.ub.e.9.10[2, 1:2] <- simpower(nrep, 0.05, TRUE, 50, 10, 0, 0.9, ex, 0.5)
sr.ub.e.9.10[2, 3:4] <- simpower(nrep, 0.05, TRUE, 50, 10, 0.2, 0.9, ex, 0.5)
sr.ub.e.9.10[2, 5:6] <- simpower(nrep, 0.05, TRUE, 50, 10, 0.5, 0.9, ex, 0.5)











####################################################################
### Signed-rank, AR1, rho = 0.1, clus size
### = 2, delta = 0, 0.2, 0.5, grpsize = 20, 50
####################################################################
sr.b.a.1.2 <- matrix(0, 2, 6)
sr.b.a.1.2[1, 1:2] <- simpower(nrep, 0.05, TRUE, 20, 2, 0, 0.1, ar1, 0)
sr.b.a.1.2[1, 3:4] <- simpower(nrep, 0.05, TRUE, 20, 2, 0.2, 0.1, ar1, 0)
sr.b.a.1.2[1, 5:6] <- simpower(nrep, 0.05, TRUE, 20, 2, 0.5, 0.1, ar1, 0)
sr.b.a.1.2[2, 1:2] <- simpower(nrep, 0.05, TRUE, 50, 2, 0, 0.1, ar1, 0)
sr.b.a.1.2[2, 3:4] <- simpower(nrep, 0.05, TRUE, 50, 2, 0.2, 0.1, ar1, 0)
sr.b.a.1.2[2, 5:6] <- simpower(nrep, 0.05, TRUE, 50, 2, 0.5, 0.1, ar1, 0)

####################################################################
### Signed-rank, AR1, rho = 0.5, clus size
### = 2, delta = 0, 0.2, 0.5, grpsize = 20, 50
####################################################################
sr.b.a.5.2 <- matrix(0, 2, 6)
sr.b.a.5.2[1, 1:2] <- simpower(nrep, 0.05, TRUE, 20, 2, 0, 0.5, ar1, 0)
sr.b.a.5.2[1, 3:4] <- simpower(nrep, 0.05, TRUE, 20, 2, 0.2, 0.5, ar1, 0)
sr.b.a.5.2[1, 5:6] <- simpower(nrep, 0.05, TRUE, 20, 2, 0.5, 0.5, ar1, 0)
sr.b.a.5.2[2, 1:2] <- simpower(nrep, 0.05, TRUE, 50, 2, 0, 0.5, ar1, 0)
sr.b.a.5.2[2, 3:4] <- simpower(nrep, 0.05, TRUE, 50, 2, 0.2, 0.5, ar1, 0)
sr.b.a.5.2[2, 5:6] <- simpower(nrep, 0.05, TRUE, 50, 2, 0.5, 0.5, ar1, 0)

####################################################################
### Signed-rank, AR1, rho = 0.9, clus size
### = 2, delta = 0, 0.2, 0.5, grpsize = 20, 50
####################################################################
sr.b.a.9.2 <- matrix(0, 2, 6)
sr.b.a.9.2[1, 1:2] <- simpower(nrep, 0.05, TRUE, 20, 2, 0, 0.9, ar1, 0)
sr.b.a.9.2[1, 3:4] <- simpower(nrep, 0.05, TRUE, 20, 2, 0.2, 0.9, ar1, 0)
sr.b.a.9.2[1, 5:6] <- simpower(nrep, 0.05, TRUE, 20, 2, 0.5, 0.9, ar1, 0)
sr.b.a.9.2[2, 1:2] <- simpower(nrep, 0.05, TRUE, 50, 2, 0, 0.9, ar1, 0)
sr.b.a.9.2[2, 3:4] <- simpower(nrep, 0.05, TRUE, 50, 2, 0.2, 0.9, ar1, 0)
sr.b.a.9.2[2, 5:6] <- simpower(nrep, 0.05, TRUE, 50, 2, 0.5, 0.9, ar1, 0)

####################################################################
### Signed-rank, AR1, rho = 0.1, clus size
### = 10, delta = 0, 0.2, 0.5, grpsize = 20, 50
####################################################################
sr.b.a.1.10 <- matrix(0, 2, 6)
sr.b.a.1.10[1, 1:2] <- simpower(nrep, 0.05, TRUE, 20, 10, 0, 0.1, ar1, 0)
sr.b.a.1.10[1, 3:4] <- simpower(nrep, 0.05, TRUE, 20, 10, 0.2, 0.1, ar1, 0)
sr.b.a.1.10[1, 5:6] <- simpower(nrep, 0.05, TRUE, 20, 10, 0.5, 0.1, ar1, 0)
sr.b.a.1.10[2, 1:2] <- simpower(nrep, 0.05, TRUE, 50, 10, 0, 0.1, ar1, 0)
sr.b.a.1.10[2, 3:4] <- simpower(nrep, 0.05, TRUE, 50, 10, 0.2, 0.1, ar1, 0)
sr.b.a.1.10[2, 5:6] <- simpower(nrep, 0.05, TRUE, 50, 10, 0.5, 0.1, ar1, 0)

####################################################################
### Signed-rank, AR1, rho = 0.5, clus size
### = 10, delta = 0, 0.2, 0.5, grpsize = 20, 50
####################################################################
sr.b.a.5.10 <- matrix(0, 2, 6)
sr.b.a.5.10[1, 1:2] <- simpower(nrep, 0.05, TRUE, 20, 10, 0, 0.5, ar1, 0)
sr.b.a.5.10[1, 3:4] <- simpower(nrep, 0.05, TRUE, 20, 10, 0.2, 0.5, ar1, 0)
sr.b.a.5.10[1, 5:6] <- simpower(nrep, 0.05, TRUE, 20, 10, 0.5, 0.5, ar1, 0)
sr.b.a.5.10[2, 1:2] <- simpower(nrep, 0.05, TRUE, 50, 10, 0, 0.5, ar1, 0)
sr.b.a.5.10[2, 3:4] <- simpower(nrep, 0.05, TRUE, 50, 10, 0.2, 0.5, ar1, 0)
sr.b.a.5.10[2, 5:6] <- simpower(nrep, 0.05, TRUE, 50, 10, 0.5, 0.5, ar1, 0)
####################################################################
### Signed-rank, AR1, rho = 0.9, clus size
### = 10, delta = 0, 0.2, 0.5, grpsize = 20, 50
####################################################################
sr.b.a.9.10 <- matrix(0, 2, 6)
sr.b.a.9.10[1, 1:2] <- simpower(nrep, 0.05, TRUE, 20, 10, 0, 0.9, ar1, 0)
sr.b.a.9.10[1, 3:4] <- simpower(nrep, 0.05, TRUE, 20, 10, 0.2, 0.9, ar1, 0)
sr.b.a.9.10[1, 5:6] <- simpower(nrep, 0.05, TRUE, 20, 10, 0.5, 0.9, ar1, 0)
sr.b.a.9.10[2, 1:2] <- simpower(nrep, 0.05, TRUE, 50, 10, 0, 0.9, ar1, 0)
sr.b.a.9.10[2, 3:4] <- simpower(nrep, 0.05, TRUE, 50, 10, 0.2, 0.9, ar1, 0)
sr.b.a.9.10[2, 5:6] <- simpower(nrep, 0.05, TRUE, 50, 10, 0.5, 0.9, ar1, 0)

####################################################################
### Signed-rank, AR1, rho = 0.1, clus size
### = 5, delta = 0, 0.2, 0.5, grpsize = 20, 50, missingrate = 0.5
####################################################################
sr.ub.a.1.5 <- matrix(0, 2, 6)
sr.ub.a.1.5[1, 1:2] <- simpower(nrep, 0.05, TRUE, 20, 5, 0, 0.1, ar1, 0.5)
sr.ub.a.1.5[1, 3:4] <- simpower(nrep, 0.05, TRUE, 20, 5, 0.2, 0.1, ar1, 0.5)
sr.ub.a.1.5[1, 5:6] <- simpower(nrep, 0.05, TRUE, 20, 5, 0.5, 0.1, ar1, 0.5)
sr.ub.a.1.5[2, 1:2] <- simpower(nrep, 0.05, TRUE, 50, 5, 0, 0.1, ar1, 0.5)
sr.ub.a.1.5[2, 3:4] <- simpower(nrep, 0.05, TRUE, 50, 5, 0.2, 0.1, ar1, 0.5)
sr.ub.a.1.5[2, 5:6] <- simpower(nrep, 0.05, TRUE, 50, 5, 0.5, 0.1, ar1, 0.5)
####################################################################
### Signed-rank, AR1, rho = 0.5, clus size
### = 5, delta = 0, 0.2, 0.5, grpsize = 20, 50, missingrate = 0.5
####################################################################
sr.ub.a.5.5 <- matrix(0, 2, 6)
sr.ub.a.5.5[1, 1:2] <- simpower(nrep, 0.05, TRUE, 20, 5, 0, 0.5, ar1, 0.5)
sr.ub.a.5.5[1, 3:4] <- simpower(nrep, 0.05, TRUE, 20, 5, 0.2, 0.5, ar1, 0.5)
sr.ub.a.5.5[1, 5:6] <- simpower(nrep, 0.05, TRUE, 20, 5, 0.5, 0.5, ar1, 0.5)
sr.ub.a.5.5[2, 1:2] <- simpower(nrep, 0.05, TRUE, 50, 5, 0, 0.5, ar1, 0.5)
sr.ub.a.5.5[2, 3:4] <- simpower(nrep, 0.05, TRUE, 50, 5, 0.2, 0.5, ar1, 0.5)
sr.ub.a.5.5[2, 5:6] <- simpower(nrep, 0.05, TRUE, 50, 5, 0.5, 0.5, ar1, 0.5)
####################################################################
### Signed-rank, AR1, rho = 0.9, clus size
### = 5, delta = 0, 0.2, 0.5, grpsize = 20, 50, missingrate = 0.5
####################################################################
sr.ub.a.9.5 <- matrix(0, 2, 6)
sr.ub.a.9.5[1, 1:2] <- simpower(nrep, 0.05, TRUE, 20, 5, 0, 0.9, ar1, 0.5)
sr.ub.a.9.5[1, 3:4] <- simpower(nrep, 0.05, TRUE, 20, 5, 0.2, 0.9, ar1, 0.5)
sr.ub.a.9.5[1, 5:6] <- simpower(nrep, 0.05, TRUE, 20, 5, 0.5, 0.9, ar1, 0.5)
sr.ub.a.9.5[2, 1:2] <- simpower(nrep, 0.05, TRUE, 50, 5, 0, 0.9, ar1, 0.5)
sr.ub.a.9.5[2, 3:4] <- simpower(nrep, 0.05, TRUE, 50, 5, 0.2, 0.9, ar1, 0.5)
sr.ub.a.9.5[2, 5:6] <- simpower(nrep, 0.05, TRUE, 50, 5, 0.5, 0.9, ar1, 0.5)
####################################################################
### Signed-rank, AR1, rho = 0.1, clus size
### = 10, delta = 0, 0.2, 0.5, grpsize = 20, 50, missingrate = 0.5
####################################################################
sr.ub.a.1.10 <- matrix(0, 2, 6)
sr.ub.a.1.10[1, 1:2] <- simpower(nrep, 0.05, TRUE, 20, 10, 0, 0.1, ar1, 0.5)
sr.ub.a.1.10[1, 3:4] <- simpower(nrep, 0.05, TRUE, 20, 10, 0.2, 0.1, ar1, 0.5)
sr.ub.a.1.10[1, 5:6] <- simpower(nrep, 0.05, TRUE, 20, 10, 0.5, 0.1, ar1, 0.5)
sr.ub.a.1.10[2, 1:2] <- simpower(nrep, 0.05, TRUE, 50, 10, 0, 0.1, ar1, 0.5)
sr.ub.a.1.10[2, 3:4] <- simpower(nrep, 0.05, TRUE, 50, 10, 0.2, 0.1, ar1, 0.5)
sr.ub.a.1.10[2, 5:6] <- simpower(nrep, 0.05, TRUE, 50, 10, 0.5, 0.1, ar1, 0.5)
####################################################################
### Signed-rank, AR1, rho = 0.5, clus size
### = 10, delta = 0, 0.2, 0.5, grpsize = 20, 50, missingrate = 0.5
####################################################################
sr.ub.a.5.10 <- matrix(0, 2, 6)
sr.ub.a.5.10[1, 1:2] <- simpower(nrep, 0.05, TRUE, 20, 10, 0, 0.5, ar1, 0.5)
sr.ub.a.5.10[1, 3:4] <- simpower(nrep, 0.05, TRUE, 20, 10, 0.2, 0.5, ar1, 0.5)
sr.ub.a.5.10[1, 5:6] <- simpower(nrep, 0.05, TRUE, 20, 10, 0.5, 0.5, ar1, 0.5)
sr.ub.a.5.10[2, 1:2] <- simpower(nrep, 0.05, TRUE, 50, 10, 0, 0.5, ar1, 0.5)
sr.ub.a.5.10[2, 3:4] <- simpower(nrep, 0.05, TRUE, 50, 10, 0.2, 0.5, ar1, 0.5)
sr.ub.a.5.10[2, 5:6] <- simpower(nrep, 0.05, TRUE, 50, 10, 0.5, 0.5, ar1, 0.5)
####################################################################
### Signed-rank, AR1, rho = 0.9, clus size
### = 10, delta = 0, 0.2, 0.5, grpsize = 20, 50, missingrate = 0.5
####################################################################
sr.ub.a.9.10 <- matrix(0, 2, 6)
sr.ub.a.9.10[1, 1:2] <- simpower(nrep, 0.05, TRUE, 20, 10, 0, 0.9, ar1, 0.5)
sr.ub.a.9.10[1, 3:4] <- simpower(nrep, 0.05, TRUE, 20, 10, 0.2, 0.9, ar1, 0.5)
sr.ub.a.9.10[1, 5:6] <- simpower(nrep, 0.05, TRUE, 20, 10, 0.5, 0.9, ar1, 0.5)
sr.ub.a.9.10[2, 1:2] <- simpower(nrep, 0.05, TRUE, 50, 10, 0, 0.9, ar1, 0.5)
sr.ub.a.9.10[2, 3:4] <- simpower(nrep, 0.05, TRUE, 50, 10, 0.2, 0.9, ar1, 0.5)
sr.ub.a.9.10[2, 5:6] <- simpower(nrep, 0.05, TRUE, 50, 10, 0.5, 0.9, ar1, 0.5)









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


save.image("sim.RData")

