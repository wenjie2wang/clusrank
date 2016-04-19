cluswilcox.test.ranksum.rgl <- function(x, cluster, group, stratum,
                                         alternative, exact, mu,
                                         DNAME = NULL, METHOD = NULL) {
    l.clus.grp <- length(unlist(lapply(split(cluster, group), unique)))
    l.clus <- length(unique(cluster))
     arglist <- setNames(list(x, cluster, group, stratum, alternative,
                                 mu, DNAME, METHOD, exact),
                            c("x", "cluster", "group", "stratum",
                              "alternative", "mu", "DNAME", "METHOD",
                              "exact"))
    if(l.clus != l.clus.grp) {
            result <- do.call("cluswilcox.test.ranksum.rgl.sub", c(arglist))
         return(result)
    } else {
        result <- do.call("cluswilcox.test.ranksum.rgl.clus", c(arglist))
        return(result)
    }
}

cluswilcox.test.ranksum.rgl.clus <- function(x, cluster, group,
                                             stratum, alternative, exact,
                                             mu, DNAME = NULL, METHOD = NULL) {
## The input data should be already in the increasing order with cluster.
    n.obs <- length(x)
    one <- rep(1, n.obs)
    x[which(group == 1)] <- x[which(group == 1)] - mu
    xrank <- rank(x)
    rksum <-stats::aggregate(xrank ~ cluster, FUN = sum)[, 2]
    csize <-stats::aggregate(one ~ cluster, FUN = sum)[, 2]
    grp <- stats::aggregate(group ~ cluster, FUN = mean)[, 2]
    strt <- stats::aggregate(stratum ~ cluster, FUN = mean)[, 2]
    clus <- unique(cluster)
    dat <- data.frame(clus, strt, grp, csize, rksum)
    bal <- (!(length(table(cluster)) != 1L))

    csize.uniq <- unique(csize)
    strt.uniq <- unique(strt)
    l.csu <- length(csize.uniq)
    l.stu <- length(strt.uniq)

    dat.l <- split(dat, strt) ## Split the rank data by stratumum
    csize.split <- function(dat) {
        split(dat, dat$"csize")
    }
    dat.l <- lapply(dat.l, csize.split)
    if(is.null(exact)) {
        exact <- FALSE
    }
    if(exact == TRUE) {
        ## Small sample permutation test
        perm.l <- vector("list", l.csu * l.stu)
        counter <- 1
        perm.len <- rep(0, l.csu * l.stu)
        for( i in csize.uniq) {
            for( j in strt.uniq) {
                j.ch <- as.character(j)
                i.ch <- as.character(i)
                temp <- dat.l[[j.ch]][[i.ch]] ## data with csize i and stratum j
                temp.rksum <- temp[, "rksum"]
                temp.grp <- temp[, "grp"]
                mgv <- length(temp.grp[temp.grp == 1])
                perm.l[[counter]] <- utils::combn(temp.rksum, mgv, FUN = sum)
                perm.len[counter] <- length(perm.l[[counter]])
                counter <- counter + 1
             }
         }
        all.perm <- base::expand.grid(perm.l)
        all.Wc <- rowSums(all.perm)
        ##all.Wc is all possible values of Wc
        Wc <- sum(dat[dat$grp == 1, "rksum"])
        ecdf.Wc <- ecdf(all.Wc)
        pval<- switch(alternative,
                      less = ecdf.Wc(abs(Wc)),
                      greater = 1 - ecdf.wc(abs(Wc)),
                      two.sided = 2 * min(ecdf.wc(abs(Wc)),
                                          1 - ecdf.wc(abs(Wc)),
                                          0.5))

        names(mu) <- "location"

        names(Wc) <- "Rank sum statistic"

        result <- list(rstatistic = Wc, p.value = pval,
                       null.value = mu, alternative = alternative,
                       data.name = DNAME, method = METHOD,
                       balance = balance)
        class(result) <- "ctest"
        return(result)


        
    } else {
        mgv <- ngv <- Ngv <- Rsumgv <- VRgv <- numeric(l.csu * l.stu)
    ## mgv: number of clusters under trt X with csize g in stratumum v
    ## ngv: number of clusters under trt Y with csize g in stratumum v
    ## Ngv: number of clusters with csize g in stratumum v
    ## Rsumgv: rank sum of clusters with size g in stratumum v
    ## VRgv: variance of rank of clusters with size g in stratumum v
    Csize <- rep(csize.uniq, l.stu)
    Str <- rep(strt.uniq, each = l.csu)
    summary.mat <- data.frame(mgv, ngv, Ngv, Rsumgv, VRgv, Csize, Str)
    for( i in csize.uniq) {
        for( j in strt.uniq) {
            j.ch <- as.character(j)
            i.ch <- as.character(i)
            temp <- dat.l[[j.ch]][[i.ch]] ## data with csize i and stratumum j
            temp.clus <- temp[, "clus"]
            temp.rksum <- temp[, "rksum"]
            temp.grp <- temp[, "grp"]
            summary.mat[which(Csize == i & Str == j), "mgv"] <- length(temp.grp[temp.grp == 1])
            summary.mat[which(Csize == i & Str == j), "Ngv"] <- length(temp.clus)
            summary.mat[which(Csize == i & Str == j), "Rsumgv"] <- sum(temp.rksum)
            summary.mat[which(Csize == i & Str == j), "VRgv"] <- var(temp.rksum) * (length(temp.clus) - 1)
        }
    }
        summary.mat[, "ngv"] <- summary.mat[, "Ngv"] - summary.mat[, "mgv"]
    ## Wc is the ranksum stat
        Wc <- sum(dat[dat$grp == 1, "rksum"])
        Ngv <- summary.mat[, "Ngv"]
        mgv <- summary.mat[, "mgv"]
        ngv <- summary.mat[, "ngv"]
        VRgv <- summary.mat[, "VRgv"]
        Rsumgv <- summary.mat[, "Rsumgv"]
        EWc <- sum(mgv * Rsumgv / Ngv)
        VarWc <- sum(mgv * ngv / (Ngv * (Ngv - 1)) * VRgv)
        Zc <- (Wc - EWc) / sqrt(VarWc)

        pval <- switch(alternative, less = pnorm(abs(Zc)),
                       greater = pnorm(abs(Zc), lower.tail = FALSE),
                       two.sided = 2 * min(pnorm(abs(Zc)),
                                           pnorm(abs(Zc), lower.tail = FALSE)))
  
        names(Wc) <- "Rank sum statistic"
        names(EWc) <- "Expected value of rank sum statistic"
        names(VarWc) <- "Variance of rank sum statistic"
        names(Zc) <- "Test statistic"
        names(mu) <- "difference in locations"
        result <- list(rstatistic = Wc, erstatistic = EWc,
                 vrstatistic = VarWc,
                 statistic = Zc, p.value = pval,
                 alternative = alternative,
                 null.value = mu,
                 data.name = DNAME,
                 method = METHOD,
                 balance = bal)
        class(result) <- "ctest"
        return(result)  
    }

}



cluswilcox.test.ranksum.rgl.sub <- function(x, cluster, group, alternative,
                                       exact, mu, DNAME = NULL, METHOD = NULL) {
    ## The input data should be already arranged
    ## check balance of data.
    
    bal <- (length(table(table(cluster))) == 1)
    clus <- unique(cluster)
    n.clus <- length(clus)
    n.obs <- length(x)
    x[which(group == 1)] <- x[which(group == 1)] - mu
    if(bal == TRUE) {
        xrank <- rank(x)
        rksum <- stats::aggregate(xrank ~ cluster, FUN = sum)[, 2]
        csize <- table(table(cluster))
        count.grp <- function(x) length(x[ x==1])
        ## q: obs under trt x in each cluster
        q <-stats::aggregate(group ~ cluster, FUN = count.grp)[, 2]
        ones <- rep(1, n.clus)
        nq.mat <- stats::aggregate(ones ~ q, FUN = sum)
        vrk <- stats::aggregate(xrank ~ cluster, FUN = var)[, 2]
        q.uniq <- nq.mat[, 1]
        nq <- nq.mat[, 2]
        varQ <- (sum(q.uniq^2 * nq) - (sum(q.uniq * nq)) ^ 2 / n.clus) / n.clus
        sb.2 <- sum((rksum - csize * (csize * n.clus + 1) / 2) ^ 2) / n.clus
        sw.2 <- sum(vrk) / n.clus
        Wc <- sum(as.numeric(group == 1) * xrank)
        EWc <- (n.clus * csize + 1) / 2 * sum(q.uniq * nq)
        VarWc <- n.clus *
            ((n.clus / (n.clus - 1)) * varQ * sb.2 / csize ^ 2) +
            sum(q.uniq * (csize - q.uniq) * nq) * sw.2 / csize
        Zc <- (Wc - EWc) / sqrt(VarWc)
         pval <- switch(alternative,
                 less = pnorm(abs(Zc)),
                 greater = pnorm(abs(Zc), lower.tail = FALSE),
                 two.sided = 2 * min(pnorm(abs(Zc)),
                                     pnorm(abs(Zc), lower.tail = FALSE)))
        names(Wc) <- "Rank sum statistic"
        names(EWc) <- "Expected value of rank sum statistic"
        names(VarWc) <- "Variance of rank sum statistic"
        names(Zc) <- "Test statistic"
        names(mu) <- "difference in locations"
        result <- list(rstatistic = Wc, erstatistic = EWc,
                 vrstatistic = VarWc,
                 statistic = Zc, p.value = pval,
                 alternative = alternative,
                 null.value = mu,
                 data.name = DNAME,
                 method = METHOD,
                 balance = bal)
        class(result) <- "ctest"
        
    } else {
        dat <- data.frame(x, cluster, group)
        ones <- rep(1, n.obs)
        csize <- stats::aggregate(ones ~ cluster, FUN = sum)
        n.csize <- length(table(csize[, "ones"]))
        colnames(csize)[2] <- "csize"
        dat <- merge(dat, csize, by = "cluster")
        dat.l <- split(dat, dat$csize)
        getWc <- function(dat) {
            xrank <- rank(dat[, "x"])
            sum(xrank * as.numeric(dat[, "group"] == 1))
        }
        Wc.l <- lapply(dat.l, getWc)
        Wc <- unlist(Wc.l)
        csize.grp <- as.numeric(names(dat.l))
        CountQplus <- function(dat) {
            sum(length(which(dat[, "group"] == 1 )))
        }
        Qplus <- lapply(dat.l, CountQplus)
        Qplus <- unlist(Qplus)

        CountN <- function(dat) {
            length(table(dat[, "cluster"]))
        }

        N.clus <- lapply(dat.l, CountN)
        N.clus <- unlist(N.clus)

        getRank <- function(dat) {
            rank(dat[, "x"])
        }
        rank.l <- lapply(dat.l, getRank)
        
        varR <- lapply(rank.l, var)
        varR <- unlist(varR)
        getGrp <- function(dat) {
            dat[, "group"]
        }
        grp.l <- lapply(dat.l, getGrp)
   

        ## theta: estimated prob of an obs under trt1 has a larger rank than an obs from trt 2
        ## deno: the denominator used in the estimation of intercluster correlation
        theta <- deno <- numeric(n.csize)
        trans.rk <- vector("list", n.csize)
        for( i in 1 : n.csize) {
            theta[i] <-( Wc[i] - Qplus[i] * (Qplus[i] + 1) / 2 ) / (Qplus[i] * (csize.grp[i] * N.clus[i] - Qplus[i]))
            deno[i] <- csize.grp[i] * N.clus[i] + 1
            trans.rk[[i]] <- rank.l[[i]] / deno[i]
        }
        pscore <- unlist(lapply(trans.rk, qnorm))
        aovout <- stats::aov(pscore ~ as.factor(unlist(grp.l)))
        aovout <- summary(aovout)

         errorms <- aovout[[1]][["Mean Sq"]][2]
        modelms <- aovout[[1]][["Mean Sq"]][1]
        sumn <- n.obs
        sumnsq <- sum(unlist(N.clus) * csize.grp ^ 2)

        g0 <- (sumn - sumnsq / sumn) / (n.obs - 1)
        sigsqa <- (modelms - errorms) / g0
        rhoh <- sigsqa / (sigsqa + errorms)
        if(rhoh < 0) {
      rhoh <- 0
    }
    rhoh <- rhoh * ( 1 + ( 1 - rhoh ^ 2) / (2 * (n.obs - 4)))
        rhor <- (6 / pi) * asin(rhoh / 2)
        varR <- unlist(varR)
        ssb <- csize.grp * varR * ( 1 + (csize.grp - 1) * rhoh)
        ssw <- varR * ( 1 - rhor)

        getEQgQ <- function(dat) {
            tab.clus <- table(dat[, "cluster"])
            csize <- tab.clus[1]
            n.clus <- length(tab.clus)
            summary.q <- table(dat[, "cluster"] * as.numeric(dat[, "group"] == 1))
            sum.n <- names(summary.q)
            summary.q <- summary.q[sum.n[-which(sum.n == '0')]]
            sum.q <- table(summary.q)
            q.in <- as.numeric(names(sum.q))
            E <- sum(q.in * sum.q * (csize - q.in)) / n.clus
            return(E)
        }

        EQgQ <- lapply(dat.l, getEQgQ)
        EQgQ <- unlist(EQgQ)

        getVarQ <- function(dat) {
            tab.clus <- table(dat[, "cluster"])
            n.clus <- length(tab.clus)
            summary.q <- table(dat[, "cluster"] * as.numeric(dat[, "group"] == 1))
            sum.n <- names(summary.q)
            summary.q <- summary.q[sum.n[-which(sum.n == '0')]]
            sum.q <- table(summary.q)
            q.in <- as.numeric(names(sum.q))
            V <- (sum(q.in ^ 2 * sum.q) - sum(q.in * sum.q) ^ 2 / n.clus) / n.clus
            return(V)
        }

        VarQ <- lapply(dat.l, getVarQ)
        VarQ <- unlist(VarQ)

        VarTheta <- N.clus * (N.clus / (N.clus - 1) * as.numeric(N.clus > 1) *
                              VarQ * ssb / csize.grp ^ 2 + EQgQ * ssw / csize.grp) /
            (Qplus * (csize.grp * N.clus - Qplus)) ^ 2
        VarTheta[which(N.clus < 2)] <- 0
        wt <- rep(0, length(csize.grp))
        wt[which(VarTheta > 0)] <- VarTheta[which(VarTheta > 0)]

        sumtheta <- sum(theta * wt)
        sumwt <- sum(wt)

        theta <- sumtheta / sumwt
        sdtheta <- sqrt(1 / sumwt)
        Zc <- (theta - 1/2) / sdtheta
        pval <- switch(alternative, less = pnorm(abs(Zc)),
                       greater = pnorm(abs(zc), lower.tail = FALSE),
                       two.sided = 2 * min(pnorm(abs(Zc)),
                                           pnorm(abs(Zc), lower.tail = FALSE)))
 
        names(Zc) <- "Test statistic"
        
        names(mu) <- "difference in locations"
        result <- list( statistic = Zc, p.value = pval,
                       alternative = alternative, null.value = mu,
                       data.name = DNAME, method = METHOD,
                       balance = bal)
        class(result) <- "ctest"
        return(result)
    }
}


cluswilcox.test.ranksum.ds <- function(x, cluster, group, 
                                       alternative,
                                       mu,
                                       DNAME, METHOD) {
    group.uniq <- length(unique(group))
    if(group.uniq == 1) {
        stop("invalid group variable, should contain at least 2 groups")
    }
 
  if(group.uniq == 2) {
    #####calculate quantity 2 (using the pooled estimate of F)
      n<-length(x)
      x[which(group == 1)] <- x[which(group == 1)] - mu
      F.hat<-numeric(n)
      for (i in 1:n){
      F.hat[i] <- (sum(x <= x[i]) + sum( x < x[i])) / (2 * n)
    }
    #####calculate quantity 1 (using ECD-F for each cluster)
    #### M is No. of clusters, n is No. of observations
    M<-length(unique(cluster)) 
    n.i <- table(cluster)
    F.prop <-  numeric(n)
    for(ii in 1:n){
      F.j<-numeric(M)
      for (i in 1:M){
        F.j[i]<-(sum(x[cluster==i]<x[ii])+0.5*sum(x[cluster==i]==x[ii]))/(n.i[i])
      } 
      F.prop[ii]<-sum(F.j[-cluster[ii]])
    }   
    
    ###########calculate S=E(W*|x,g)
    a<-numeric(M)
    b<-1+F.prop
    for (i in 1:M){
      a[i]<-sum(((group[cluster==i] - 1) *b[cluster==i])/(n.i[i]))
    }    
    c<-1/(M+1)
    S<-c*sum(a)
    ########note: for m groups maybe can use group[cluster==i&group=m]
    
#########Calculate E(S)=E(W*)
     group.temp <- as.numeric(group == 2)
      n.i1 <- aggregate(group.temp ~ cluster, FUN = sum)
      n.i1 <- n.i1[, 2]
      d<-n.i1/n.i
    E.S<-(1/2)*sum(d)
    
    #######Calculate estimate of variance of S
    W.hat<-numeric(M)        #####first calculate W.hat for each cluster
    a<-n.i1/n.i
      for (i in 1:M){
          b<-1/(n.i[i]*(M+1))
          c<-(group[cluster==i] - 1) * (M-1)
          d<-sum(a[-i])
      W.hat[i]<-b*sum((c-d)*F.hat[cluster==i])
    }
      a<-n.i1/n.i
      E.W<-(M/(2*(M+1)))*(a-sum(a)/M)    ##second, calculate E(W)
    
      var.s<-sum((W.hat-E.W)^2) #calculate var(s)
      Zc<-(S-E.S)/sqrt(var.s)   #calculate the test statistic
      pval <- switch(alternative,
                 less = pnorm(abs(Zc)),
                 greater = pnorm(abs(Zc), lower.tail = FALSE),
                 two.sided = 2 * min(pnorm(abs(Zc)),
                                     pnorm(abs(Zc), lower.tail = FALSE)))
      names(Zc) <- "test statistic"
      names(mu) <- "mu"
      result <- list(statistic = Zc, p.value = pval,
                     alternative = alternative, null.value = mu,
                     data.name = DNAME, method = METHOD)
                
      class(result) <- "ctest"
      result

       
  } else {
#####calculate quantity 2 (using the pooled estimate of F)
      if(!is.null(mu) & mu != 0) {
          warning("comparison between m (m > 2) groups cannot set location shift parameter")
      }
      mu <- 0
      n<-length(x)
      F.hat<-numeric(n)
      for (i in 1:n){
          F.hat[i]<-(sum(x<=x[i])+sum(x<x[i]))/(2*n)
    }
    #####calculate quantity 1 (using ECD-F for each cluster)
    #### M is No. of clusters, n is No. of observations
      M<-length(unique(cluster)) 
      n.i<-table(cluster)
      F.prop<-numeric(n)
      for(ii in 1:n){
          F.j<-numeric(M)
      for (i in 1:M){
        F.j[i]<-(sum(x[cluster==i]<x[ii])+0.5*sum(x[cluster==i]==x[ii]))/(n.i[i])
      } 
      F.prop[ii]<-sum(F.j[-cluster[ii]])
    }   
    
    ###########calculate S(j)=E(W*|x,g=j), where m is the number of groups
    m<-length(unique(group))
    a<-matrix(0,m,M)
    b<-1+F.prop
    for(j in 1:m){
      for (i in 1:M){
        gik.j<-ifelse(group==j,1,0)
        a[j,i]<-sum((gik.j[cluster==i]*b[cluster==i])/(n.i[i]))
      } 
    }   
    c<-1/(M+1)
    S.j<-c*(apply(a,1,sum))
    
    #########Calculate E(S)=E(W*)
    n.ij<-matrix(0,m,M)
    for (i in 1:m){
      n.ij[i,]<-table(cluster[group==i])           
    }
    d<-apply(n.ij,1,FUN=function(x){x/n.i})
    E.S.j<-(1/2)*(apply(d,2,sum))
    
    #######Calculate estimate of variance of S
    W.hat<-matrix(0,m,M) 
    a<-t(d)       #####first calculate W.hat for each cluster
    for (i in 1:M){
      for (j in 1:m){
        gik.j<-ifelse(group[cluster==i]==j,1,0)
        b<-1/(n.i[i]*(M+1))
        c<-(gik.j)*(M-1)
        d<-sum(a[j,-i])
        W.hat[j,i]<-b*sum((c-d)*F.hat[cluster==i])
      }
    }
    E.W<-matrix(0,m,M)
    for (j in 1:m){
      E.W[j,]<-(M/(2*(M+1)))*(a[j,]-sum(a[j,])/M)
    }                            ##second, calculate E(W)
    
    ##########calculate sample variance
    dev.W<-W.hat-E.W
    term.old<-matrix(0,m,m)
    for (i in 1:M){
      term<-dev.W[,i]%*%t(dev.W[,i])
      term.old<-term+term.old
    }
    V.hat<-(1/M)*term.old
    
    ######calculate the test statistic
    T<-(t(S.j-E.S.j) %*%MASS::ginv(V.hat) %*% (S.j-E.S.j))*(1/M)
    pval <- pchisq(T, df=(m-1),lower.tail=F)
      names(T) <- "test statistic"
      ngrp <- group.uniq
      df <- ngrp - 1
      names(ngrp) <- "Number of groups: "
      names(df) <- "Degree of freedom: "
      METHOD <- paste(METHOD, "using Chisq test")
                                        #calculate the test statistic
      result <- list(statistic = T, p.value = pval, n.group = ngrp,
                     df = df,
                     data.name = DNAME, method = METHOD)
      class(result) <- "ctest"
      result
  }    
}
