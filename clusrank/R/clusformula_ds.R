################################################################################
##
##   R package clusrank by Mei-Ling Ting Lee, Jun Yan, and Yujing Jiang
##   Copyright (C) 2016
##
##   This file is part of the R package clusrank.
##
##   The R package clusrank is free software: you can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package clusrank is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package clusrank. If not, see <http://www.gnu.org/licenses/>.
##
################################################################################
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
      a[i]<-sum((group[cluster==i]*b[cluster==i])/(n.i[i]))
    }    
    c<-1/(M+1)
    S<-c*sum(a)
    ########note: for m groups maybe can use group[cluster==i&group=m]
    
    #########Calculate E(S)=E(W*)
    n.i1<-table(cluster[group==1])
    d<-n.i1/n.i
    E.S<-(1/2)*sum(d)
    
    #######Calculate estimate of variance of S
    W.hat<-numeric(M)        #####first calculate W.hat for each cluster
    a<-n.i1/n.i
    for (i in 1:M){
      b<-1/(n.i[i]*(M+1))
      c<-(group[cluster==i])*(M-1)
      d<-sum(a[-i])
      W.hat[i]<-b*sum((c-d)*F.hat[cluster==i])
    }
    a<-n.i1/n.i
    E.W<-(M/(2*(M+1)))*(a-sum(a)/M)    ##second, calculate E(W)
    
    var.s<-sum((W.hat-E.W)^2) #calculate var(s)
    stat<-(S-E.S)/sqrt(var.s)   #calculate the test statistic
    pvalue<-2*pnorm(stat,lower.tail=F)
    return(list(S=S,E.S=E.S,Var.S=var.s,z.stat=stat,pvalue=pvalue))
  } else {
    #####calculate quantity 2 (using the pooled estimate of F)
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
    library(MASS)
    T<-(t(S.j-E.S.j)%*%ginv(V.hat)%*%(S.j-E.S.j))*(1/M)
    pvalue<-pvalue<-pchisq(T, df=(m-1),lower.tail=F)
    
    #calculate the test statistic
    return(list(S=S.j,E.S=E.S.j,Chisq.stat=T,df=m-1,pvalue=pvalue))
  }
  
  
}
