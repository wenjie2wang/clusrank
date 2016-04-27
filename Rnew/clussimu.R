library(MASS)
clus.simu.2 <- function(M, rho, n, pn, delta, pc) {
    ## Treatments are assigned at cluster level.
    ##
    ## M: A numerical vector with two entries, indicating the number of clusters under each
    ## treatment.
    ## rho: A numerical vector, the intracluster correlation
    ## for cluster under each treatment, the larger the strong the correlation is.
    ## n: A numerical vector, the cluster sizes.
    ## pn: A numerical vector, the distribution of cluster size.
    ## delta: The location shift between two groups.
    ## pc: A numerical vector, the probability distribution of a cluster of a specific size
    ##    to be assigned treatment 0.

    if(length(n) == (length(pn) + 1)) {
        pn <- c(pn, 1 - sum(pn))
    }
    if( length(n) != length(pn)) {
        stop("n and pn are not of equal length")
    }
    if(length(M) == (length(pc) + 1)) {
        pc <- c(pc, 1 - sum(pc))
    }
    if( length(M) != length(pc)) {
        stop(" M and pc are not of equal length")
    }

    n.c.uniq <- length(n)
    cov.list <- vector("list", n.c.uniq)
    for( i in 1 : n.c.uniq) {
        cov.list[[i]] <- (1 - rho[i]) * diag(n[i]) + rho[i] * matrix(1, n[i], n[i])
    }
    
    n.dist <- rmultinom(1, sum(M), pn)
    c.id <- c(1 : sum(M)) ## cluster id
    c.size <- rep.int(n, times = n.dist) ## size of each cluster
    c.id <- rep.int(c.id, times = c.size)

    n.obs <- sum(c.size)

    group <- rep(0, sum(M))
    count <- 0
    count0 <- 1
    for( i in 1 : n.c.uniq) {
        count <- count +  n.dist[i]
        group[count0 : count] <- rbinom(n.dist[i], 1, c(pc, 1 - pc))
        count0 <- count + 1
    }

    group <- rep.int(group, times = c.size) ## group id for each observation

    Y.list <- vector("list", n.c.uniq) ## Generate multinormal data
    X.list <- vector("list", n.c.uniq) ## STore generate data
    for( i in 1 : n.c.uniq) {
        Y.list[[i]] <- mvrnorm(n.dist[i], numeric(c.size[i]), cov.list[[i]])
        X.list[[i]] <- exp(c(t(Y.list[[i]])))
    }
    X <- unlist(X.list)
    X <- X + delta * group
    X
}
