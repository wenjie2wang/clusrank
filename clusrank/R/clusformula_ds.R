cluswilcox.test.ranksum.ds <- function(X, cluster, grp, 
                                       alternative,
                                       mu,
                                       DNAME, METHOD) {
  
  group.n <- length(unique(grp))
  if(group.n == 1) {
    stop("invalid group variable, should contain at least 2 groups")
  }
  if(group.n == 2) {
    #####calculate quantity 2 (using the pooled estimate of F)
    n<-length(X)
    F.hat<-numeric(n)
    for (i in 1:n){
      F.hat[i]<-(sum(X<=X[i])+sum(X<X[i]))/(2*n)
    }
    #####calculate quantity 1 (using ECD-F for each cluster)
    #### M is No. of clusters, n is No. of observations
    M<-length(unique(cluster)) 
    n.i<-table(cluster)
    F.prop<-numeric(n)
    for(ii in 1:n){
      F.j<-numeric(M)
      for (i in 1:M){
        F.j[i]<-(sum(X[cluster==i]<X[ii])+0.5*sum(X[cluster==i]==X[ii]))/(n.i[i])
      } 
      F.prop[ii]<-sum(F.j[-cluster[ii]])
    }   
    
    ###########calculate S=E(W*|X,g)
    a<-numeric(M)
    b<-1+F.prop
    for (i in 1:M){
      a[i]<-sum((grp[cluster==i]*b[cluster==i])/(n.i[i]))
    }    
    c<-1/(M+1)
    S<-c*sum(a)
    ########note: for m groups maybe can use grp[cluster==i&grp=m]
    
    #########Calculate E(S)=E(W*)
    n.i1<-table(cluster[grp==1])
    d<-n.i1/n.i
    E.S<-(1/2)*sum(d)
    
    #######Calculate estimate of variance of S
    W.hat<-numeric(M)        #####first calculate W.hat for each cluster
    a<-n.i1/n.i
    for (i in 1:M){
      b<-1/(n.i[i]*(M+1))
      c<-(grp[cluster==i])*(M-1)
      d<-sum(a[-i])
      W.hat[i]<-b*sum((c-d)*F.hat[cluster==i])
    }
    a<-n.i1/n.i
    E.W<-(M/(2*(M+1)))*(a-sum(a)/M)    ##second, calculate E(W)
    
    var.s<-sum((W.hat-E.W)^2) #calculate var(s)
    stat<-(S-E.S)/sqrt(var.s)   #calculate the test statistic
    p.value<-2*pnorm(stat,lower.tail=F)
    return(list(S=S,E.S=E.S,Var.S=var.s,z.stat=stat,p.value=p.value))
  } else {
    #####calculate quantity 2 (using the pooled estimate of F)
    n<-length(X)
    F.hat<-numeric(n)
    for (i in 1:n){
      F.hat[i]<-(sum(X<=X[i])+sum(X<X[i]))/(2*n)
    }
    #####calculate quantity 1 (using ECD-F for each cluster)
    #### M is No. of clusters, n is No. of observations
    M<-length(unique(cluster)) 
    n.i<-table(cluster)
    F.prop<-numeric(n)
    for(ii in 1:n){
      F.j<-numeric(M)
      for (i in 1:M){
        F.j[i]<-(sum(X[cluster==i]<X[ii])+0.5*sum(X[cluster==i]==X[ii]))/(n.i[i])
      } 
      F.prop[ii]<-sum(F.j[-cluster[ii]])
    }   
    
    ###########calculate S(j)=E(W*|X,g=j), where m is the number of groups
    m<-length(unique(grp))
    a<-matrix(0,m,M)
    b<-1+F.prop
    for(j in 1:m){
      for (i in 1:M){
        gik.j<-ifelse(grp==j,1,0)
        a[j,i]<-sum((gik.j[cluster==i]*b[cluster==i])/(n.i[i]))
      } 
    }   
    c<-1/(M+1)
    S.j<-c*(apply(a,1,sum))
    
    #########Calculate E(S)=E(W*)
    n.ij<-matrix(0,m,M)
    for (i in 1:m){
      n.ij[i,]<-table(cluster[grp==i])           
    }
    d<-apply(n.ij,1,FUN=function(x){x/n.i})
    E.S.j<-(1/2)*(apply(d,2,sum))
    
    #######Calculate estimate of variance of S
    W.hat<-matrix(0,m,M) 
    a<-t(d)       #####first calculate W.hat for each cluster
    for (i in 1:M){
      for (j in 1:m){
        gik.j<-ifelse(grp[cluster==i]==j,1,0)
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
    library(MASS)
    T<-(t(S.j-E.S.j)%*%ginv(V.hat)%*%(S.j-E.S.j))*(1/M)
    p.value<-p.value<-pchisq(T, df=(m-1),lower.tail=F)
    
    #calculate the test statistic
    return(list(S=S.j,E.S=E.S.j,Chisq.stat=T,df=m-1,p.value=p.value))
  }
  
  
}
