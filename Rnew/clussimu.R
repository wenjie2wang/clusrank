library(MASS)
clus.simu.2 <- function(M, rho, n, delta, pc) {
    ## Treatments are assigned at cluster level.
    ##
    ## M: A numerical vector with two entries, indicating the number of clusters under each
    ## treatment.
    ## rho: A numerical vector, the intracluster correlation, specific for each cluster size
    ##      If only one value is supplied, then it is used for all cluster size
    ## for cluster under each treatment, the larger the strong the correlation is.
    ## n: A numerical vector, the cluster sizes.
    ## delta: The location shift between two groups.
    ## pc: A numerical vector or a matrix, the probability distribution
    ## treatment. If it is a matrix, each row represent a treatment group
    ##     
    ## 0 to be assigned to a cluster with a specific size
    ## Example:  clus.simu.2(c(2, 3), c(0.2, 0.5), c(2, 3), 0.2, 1, 0.2)

    n.grp <- length(M)
    grp.uniq <- c(0 : (n.grp - 1))
    if(!is.vector(rho)) {
        rho <- rep(rho, length(n))
    }
    if(length(rho) != length(n)) {
        stop("rho is not of the same length as n")
    }
    if(is.matrix(pc)) {
        n.in.pc <- nrow(pc)
    } else if(is.vector(pc)) {
        n.in.pc <- length(pc)
        pc <- replicate(length(M), pc)
        pc <- as.matrix(pc)
    }
    if(length(M) == (n.in.pc + 1)) {
        pc <- cbind(pc, 1 - rowSums(pc))
        n.in.pc <- n.in.pc + 1
    }
    
    if( length(M) != n.in.pc) {
        stop(" M and pc are not of equal length")
    }

    n.c.uniq <- length(n)
    cov.list <- vector("list", n.c.uniq)
    for( i in 1 : n.c.uniq) {
        cov.list[[i]] <- (1 - rho[i]) * diag(n[i]) + rho[i] * matrix(1, n[i], n[i])
    }

    
    c.id <- c(1 : sum(M)) ## cluster id

    n.dist.l <- c.size.l <-  vector("list", length(M))
    for( i in 1 : length(M)) {
        n.dist.l[[i]] <- rmultinom(1, M[i], pc[i, ])
        c.size.l[[i]] <- rep.int(n, times = n.dist.l[[i]])
        
    }
    c.size <- unlist(c.size.l)
    c.id <- rep.int(c.id, times = c.size)
    group <- rep.int(rep.int(grp.uniq, times = M), times = c.size)
    

    Y.list <- vector("list", n.grp) ## Generate multinormal data
    X.list <- vector("list", n.grp) ## STore generate data
    x.list <- vector("list", n.grp)
    for( i in 1 : n.grp) {
        counter <- 1
        for( j in 1 : length(n)) {
            if(n.dist.l[[i]][[j]] == 0) {
               next
            } else {      
                Y.list[[i]][[counter]] <- mvrnorm(n.dist.l[[i]][j],
                                            numeric(n[j]),
                                            cov.list[[j]])
                X.list[[i]][[counter]] <- exp(c(t(Y.list[[i]][[counter]]))) +
                    grp.uniq[i] * delta
                counter <- counter + 1
            }
        }
        x.list[[i]] <- unlist(X.list[[i]])
    }
            
    X <- unlist(x.list)
    data.frame(X, c.id, group)
}



clus.simu.m <- function(M, rho, n, delta, pc) {
    ## Treatments are assigned at cluster level.
    ##
    ## M: A numerical vector with two entries, indicating the number of clusters under each
    ## treatment.
    ## rho: A numerical vector, the intracluster correlation
    ## for cluster under each treatment, the larger the strong the correlation is.
    ## n: A numerical vector, the cluster sizes.
    ## delta: A vector, the location shift from 0 of the treatment groups, if length is n.grp - 1,
    ##        the first element is filled with 0.
    ## pc: A numerical vector or a matrix, the probability distribution
    ## treatment. If it is a matrix, each row represent a treatment group
    ##     
    ## 0 to be assigned to a cluster with a specific size
    ## Example:  clus.simu.2(c(2, 3), c(0.2, 0.5), c(2, 3), 0.2, 1, 0.2)

    n.grp <- length(M)
    grp.uniq <- c(1 : (n.grp))
    if(length(delta) == (n.grp- 1)) {
        delta <- c(0, delta)
    }
    if(length(delta) != n.grp) {
        stop("delta should have the same length as n")
    }
    
      if(!is.vector(rho)) {
        rho <- rep(rho, length(n))
    }
    if(length(rho) != length(n)) {
        stop("rho is not of the same length as n")
    }
    if(is.matrix(pc)) {
        n.in.pc <- nrow(pc)
    } else if(is.vector(pc)) {
        n.in.pc <- length(M)
        pc <- replicate(length(M), pc)
        pc <- t(as.matrix(pc))
    }
    if(length(n) == (ncol(pc) - 1)) {
        pc <- cbind(pc, 1 - colSums(pc))
    } 
    
    if( length(M) != n.in.pc) {
        stop(" M and pc are not of equal length")
    }

    n.c.uniq <- length(n)
    cov.list <- vector("list", n.c.uniq)
    for( i in 1 : n.c.uniq) {
        cov.list[[i]] <- (1 - rho[i]) * diag(n[i]) + rho[i] * matrix(1, n[i], n[i])
    }

    
    c.id <- c(1 : sum(M)) ## cluster id

    n.dist.l <- c.size.l <-  vector("list", length(M))
    for( i in 1 : length(M)) {
        n.dist.l[[i]] <- rmultinom(1, M[i], pc[i, ])
        c.size.l[[i]] <- rep.int(n, times = n.dist.l[[i]])
    }
    c.size <- unlist(c.size.l)
    c.id <- rep.int(c.id, times = c.size)
    group <- rep.int(rep.int(grp.uniq, times = M), times = c.size)    

    Y.list <- vector("list", n.grp) ## Generate multinormal data
    X.list <- vector("list", n.grp) ## STore generate data
    x.list <- vector("list", n.grp)
    for( i in 1 : n.grp) {
        counter <- 1
        for( j in 1 : length(n)) {
            if(n.dist.l[[i]][[j]] == 0) {
               next
            } else {      
                Y.list[[i]][[counter]] <- mvrnorm(n.dist.l[[i]][j],
                                            numeric(n[j]),
                                            cov.list[[j]])
                X.list[[i]][[counter]] <- exp(c(t(Y.list[[i]][[counter]]))) +
                    delta[i]
                counter <- counter + 1
            }
        }
        x.list[[i]] <- unlist(X.list[[i]])
    }
            
    X <- unlist(x.list)
    data.frame(X, c.id, group)
}


clus.simu.sr1 <- function(M, n, delta, rho) {
    ## Simplest data generating model in DS paper signed rank test
    ## Xij = sign(Hij)exp(|Hij|)
    ## Hij = Wi + eij, Wi ~N(delta, rho_i), eij ~ N(0, 1 - rho_i)
    ## Arguments:
    ##  M: a vector or a numeric, contains the number of clusters for each size
    ##  n: a vector or a numeric, possible cluster size
    ## rho: rho is suppposed to have the length of hte number of clusters
    if(length(rho) == 1) {
        rho <- rep(rho, sum(M))
    }
    if(length(rho) == length(M)) {
        rho <- rho[rep.int(1 : length(M), times = M)]
    }
    clus.size <- rep.int(n, times = M)
    cluster <- rep.int(1 : sum(M), times = clus.size)

    W <- rnorm(sum(M))
    W <- W[rep.int(1 : length(W), times = clus.size)]

    W <- W * sqrt(rho)[rep.int(1 : sum(M), times = clus.size)]
    W <- delta + W
    n.obs <- sum(M * n)
    e <- rnorm(n.obs)
    e <- e * sqrt(1 - rho)[rep.int( 1 : sum(M), times = clus.size)]
    H <- W + e
    X <- sign(H) * exp(abs(H))
    data.frame(X, cluster)
}


clus.simu.sr2 <- function(M, n, delta, rho, phi, gamma) {
    ## Simplest data generating model in DS paper signed rank test
    ## intracluster ecorrelation is not exchangable
    ## Xij = sign(Hij)exp(|Hij|)
    ## Hij = Wi + B_phij + eij, Wi ~N(delta, rho_i), eij ~ N(0, 1 - rho_i, gamma j)
    ## Bj ~ N(0, gamma j)
    ## A single cluster size.
    ## Arguments:
    ##  M: a vector or a numeric, contains the number of clusters for each size
    ##  n: a vector or a numeric, possible cluster size
    ## rho: rho is suppposed to have the length of hte number of clusters
    ## phi: code for membrs in ea h cluster, to create non0interchangeable correlation structure
    ## gamma: variance of Bj coressponding to each code

    if(length(rho) == 1) {
        rho <- rep(rho, sum(M))
    }
    if(length(rho) == length(M)) {
        rho <- rho[rep.int(1 : length(M), times = M)]
    }
    clus.size <- rep.int(n, times = M)
    cluster <- rep.int(1 : sum(M), times = clus.size)
    
    if(length(phi) != n) {
        stop("phi should have the length as the cluster size")
    }
    if(length(gamma) != length(unique(gamma))) {
        stop("gamma should have the same length as the unique code in phi")
    }
    W <- rnorm(sum(M))
    W <- W[rep.int(1 : length(W), times = clus.size)]

    W <- W * sqrt(rho)[rep.int(1 : sum(M), times = clus.size)]
    W <- delta + W

    phi.tab <- table(phi)
    phi.temp <- rep.int(1 : length(unique(phi)), times = phi.tab)
    
    B<- rnorm( M * length(gamma)) * sqrt(gamma)
    B.l <- length(B)
    phi.temp <- (rep( (1 : (M)), each = n) - 1) * length(unique(phi)) + phi.temp
    B <- B[phi.temp]


    gamma <- gamma[rep.int(1 : length(unique(phi)), times = table(phi))]


    n.obs <- sum(M * n)
    e <- rnorm(n.obs)
    e <- e * sqrt(1 - rho - gamma)[rep.int( 1 : sum(M), times = clus.size)]
    H <- W + e + B
    X <- sign(H) * exp(abs(H))
    data.frame(X, cluster)
}



clus.simu.sr.info <- function(M, n, beta, delta, rho) {
    ## Simplest data generating model in DS paper signed rank test
    ## informative cluster size (for the sign of the cluster)
    ## Xij = sign(Hij)exp(|Hij|)
    ## Hij = Wi + eij, Wi ~N(delta, rho_i), eij ~ N(0, 1 - rho_i)
    ## Arguments:
    ##  M:  a numeric,  total number of clusters 
    ##  n: a vector or a numeric, possible cluster size
    ##  beta: a numeric, the parameter of the beta distribution for
    ##      determine the sample size and sign of the observation
    ##  rho: rho is suppposed to have the length of hte number of clusters
    ran <- rbeta(M, 1, beta)
    cluster.size <- ifelse(ran < 0.5, n[1], n[2])
    sn <- ifelse(ran < 0.5, 1, -1)

    M <- table(sort(cluster.size))
    n <- sort(n)

    if(length(rho) == 1) {
        rho <- rep(rho, sum(M))
    }
    if(length(rho) == length(M)) {
        rho <- rho[rep.int(1 : length(M), times = M)]
    }
    
    clus.size <- rep.int(n, times = M)
    cluster <- rep.int(1 : sum(M), times = clus.size)

    W <- rnorm(sum(M))
    W <- W[rep.int(1 : length(W), times = clus.size)]

    W <- W * sqrt(rho)[rep.int(1 : sum(M), times = clus.size)]
    W <- delta + W
    n.obs <- sum(M * n)
    e <- rnorm(n.obs)
    e <- e * sqrt(1 - rho)[rep.int( 1 : sum(M), times = clus.size)]
    H <- W + e
    H <- H + e
    X <- sn * exp(abs(H))
    data.frame(X, cluster)
}


do1.sr1 <- function(M, n, delta, rho, method = c("rgl", "ds"), paired = FALSE) {
    dat <- clus.simu.sr1(M, n, delta, rho)
    method <- match.arg(method)
    if(method == "rgl") {
        pval <- cluswilcox.test(X, cluster = cluster, data = dat, paired) $ p.value
    }
    if(method == "ds") {
        pval <- cluswilcox.test(X, cluster = cluster, data = dat, paired, method = "ds") $ p.value
    }
    return(pval)
}

do.sr1 <- function(simn, M, n, delta, rho, method = c("rgl", "ds"), paired = FALSE) {
    res <- replicate(simn, do1.sr1(M, n, delta, rho, method, paired))
    return(res)
}

power.sr1 <- function(th, simn, M, n, delta, rho, method = c("rgl", "ds"), paired = FALSE) {
    res <- do.sr1(simn, M, n, delta, rho, method, paired)
    power <- mean(res < th)
    return(power)
}
    
