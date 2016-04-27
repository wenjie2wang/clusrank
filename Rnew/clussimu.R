library(MASS)
clus.simu.2 <- function(M, rho, n, delta, pc) {
    ## Treatments are assigned at cluster level.
    ##
    ## M: A numerical vector with two entries, indicating the number of clusters under each
    ## treatment.
    ## rho: A numerical vector, the intracluster correlation
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
    X
}
