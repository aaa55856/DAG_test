library(dHSIC)
library(mgcv)
library(Matrix)
library(MASS)
library(parallel)
library(graph)
library(RBGL)
library(DAAG)
library(splitstackshape)
library(gRbase)
library(Rgraphviz)
library(energy)
library(mvtnorm)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("graph")
# BiocManager::install("RBGL")
# BiocManager::install("Rgraphviz")


## For d univariate random variables 
## Input x is a n*d data matrix, n being the sample size


v.center <- function(x){
  if (is.matrix(x)) {n=dim(x)[1]
  if (isSymmetric(x)) A=x else A=as.matrix(dist(x))} else {n=length(x); A=as.matrix(dist(x))}
  R=rowSums(A)
  C=colSums(A)
  T=sum(A)                             
  r=matrix(rep(R,n),n,n)/n
  c=t(matrix(rep(C,n),n,n))/n
  t=matrix(T/n^2,n,n)
  UA=-(A-r-c+t)
  return(UA)
}

Jdcov.sq.V <- function(x,c){    ## Computes the V-Statistic type estimator of JdCov
  n=dim(x)[1]
  d=dim(x)[2]
  A <- v.center(x[,1]) + c 
  for(i in 2:d) A <- A*( v.center(x[,i]) + c )
  return( sum(A)/n^2 - c^d )
}

u.center <- function(x){
  if (is.matrix(x)) {n=dim(x)[1]
  if (isSymmetric(x)) A=x else A=as.matrix(dist(x))} else {n=length(x); A=as.matrix(dist(x))}
  R=rowSums(A)
  C=colSums(A)
  T=sum(A)
  r=matrix(rep(R,n),n,n)/(n-2)
  c=t(matrix(rep(C,n),n,n))/(n-2)
  t=matrix(T/(n-1)/(n-2),n,n)
  UA=-(A-r-c+t)
  diag(UA)=0
  return(UA)
}

Jdcov.sq.U <- function(x,c){    ## Computes the U-Statistic type estimator of JdCov
  n=dim(x)[1]
  d=dim(x)[2]
  A <- u.center(x[,1]) + c 
  for(i in 2:d) A <- A*( u.center(x[,i]) + c )
  return( sum(A)/n/(n-3)-c^d*n/(n-3) )
}

dCov.sq.U <- function(x,n) {
  A <- u.center(x)
  return( sum(A^2)/n/(n-3) )
}

Jdcov.sq.US <- function(x,c){   ## Computes the bias-corrected estimator of JdCov_S
  n=dim(x)[1]
  d=dim(x)[2]
  A <- u.center(x[,1])/sqrt( dCov.sq.U(x[,1],n) ) + c         
  for(i in 2:d) A <- A*( u.center(x[,i])/sqrt( dCov.sq.U(x[,i],n) ) + c )
  return( sum(A)/n/(n-3)-c^d*n/(n-3) )
}


## For d random vectors of arbitrary dimensions p_1, .. , p_d


v.center.list <- function(x){      ## Input x is a n*p matrix, p=p_i for some i in 1:d
  n=nrow(x)
  A=as.matrix(dist(x))
  R=rowSums(A)
  C=colSums(A)
  T=sum(A)                             
  r=matrix(rep(R,n),n,n)/n
  c=t(matrix(rep(C,n),n,n))/n
  t=matrix(T/n^2,n,n)
  UA=-(A-r-c+t)
  return(UA)
}

## The following computes the V-Statistic type estimator of JdCov

Jdcov.sq.V.list <- function(x,c){         ## input x is a list of d elements, each element being a matrix
  ## with n rows & p_i columns, i=1,..,d
  n.sample=as.matrix(sapply(x,dim))[1,]
  if(length(unique(n.sample))!=1) stop("Unequal Sample Sizes")
  n=n.sample[1]
  d=length(x)
  p=as.matrix(sapply(x,dim))[2,]
  
  A <- v.center.list(x[[1]]) + c 
  for(i in 2:d) A <- A*( v.center.list(x[[i]]) + c )
  return( sum(A)/n^2 - c^d )
}



u.center.list <- function(x){            ## Input x is a n*p matrix, p=p_i for some i in 1:d
  n=nrow(x)
  A=as.matrix(dist(x))
  R=rowSums(A)
  C=colSums(A)
  T=sum(A)
  r=matrix(rep(R,n),n,n)/(n-2)
  c=t(matrix(rep(C,n),n,n))/(n-2)
  t=matrix(T/(n-1)/(n-2),n,n)
  UA=-(A-r-c+t)
  diag(UA)=0
  return(UA)
}


## The following computes the U-Statistic type estimator of JdCov

Jdcov.sq.U.list <- function(x,c){             ## input x is a list of d elements, each element being a matrix
  ## with n rows & p_i columns, i=1,..,d
  n.sample=as.matrix(sapply(x,dim))[1,]
  if(length(unique(n.sample))!=1) stop("Unequal Sample Sizes")
  n=n.sample[1]
  d=length(x)
  p=as.matrix(sapply(x,dim))[2,]
  A <- u.center.list(x[[1]]) + c 
  for(i in 2:d) A <- A*( u.center.list(x[[i]]) + c )
  return( sum(A)/n/(n-3)-c^d*n/(n-3) )
}



dCov.sq.U.list <- function(x,n) {          ## Input x is a n*p matrix, p=p_i for some i in 1:d
  A <- u.center.list(x)
  return( sum(A^2)/n/(n-3) )
}


## The following computes the bias-corrected estimator of JdCov

Jdcov.sq.US.list <- function(x,c){         ## input x is a list of d elements, each element being a matrix
  ## with n rows & p_i columns, i=1,..,d
  n.sample=as.matrix(sapply(x,dim))[1,]
  if(length(unique(n.sample))!=1) stop("Unequal Sample Sizes")
  n=n.sample[1]
  d=length(x)
  p=as.matrix(sapply(x,dim))[2,]
  A <- u.center.list(x[[1]])/sqrt( dCov.sq.U.list(x[[1]],n) ) + c         
  for(i in 2:d) A <- A*( u.center.list(x[[i]])/sqrt( dCov.sq.U.list(x[[i]],n) ) + c )
  return( sum(A)/n/(n-3)-c^d*n/(n-3) )
}



Sim.size.list <- function(n, d=3, p=5, cc=1, index=1, type, Mp=500, reptim=1000, Ustat=FALSE)
{
  set.seed(100)
  rej  <- rej.r <- rej.s <- rej3 <- rej4 <- matrix(NA, reptim, 2)
  if(Ustat) F <- Jdcov.sq.U.list else F <- Jdcov.sq.V.list
  for(jj in 1:reptim)
  {
    if(type==3){
      x=list() 
      x[[1]]=x[[2]]=x[[3]]=matrix(0,n,p)
      set.seed(jj*12)
      x[[1]]=rmvnorm(n,rep(0,p),diag(p))
      set.seed(jj*20)
      x[[2]]=rmvnorm(n,rep(0,p),diag(p))
      set.seed(jj*30)
      W=rexp(n,1/sqrt(2))
      x[[3]][,1]=(sign(x[[1]][,1] * x[[2]][,1])) * W
      set.seed(jj*40)                                         
      x[[3]][,2:p]=rmvnorm(n,rep(0,(p-1)),diag(p-1))
    }
    if(type==4){
      x=list()
      x[[1]]=x[[2]]=x[[3]]=matrix(0,n,p)
      set.seed(jj*100)
      x[[1]]=rmvnorm(n,rep(0,p),diag(p))
      set.seed(jj*110)
      x[[2]]=rmvnorm(n,rep(0,p),diag(p))
      set.seed(jj*120)
      e=runif(1,-1,1)
      set=1:3
      set.seed(jj*130)
      set.pick=sample(set,1)
      if(set.pick==1) x[[3]][,1] = (x[[1]][,1])^2 + e
      if(set.pick==2) x[[3]][,1] = (x[[2]][,1])^2 + e
      if(set.pick==3) x[[3]][,1] = (x[[1]][,1])*(x[[2]][,1]) + e
      set.seed(jj*140)                                        
      x[[3]][,2:p]=rmvnorm(n,rep(0,(p-1)),diag(p-1))
    }
    # X <- split(x, rep(1:ncol(x), each = nrow(x)))
    x.r <- x
    for(j in 1:d) {
      for(l in 1:ncol(x[[j]])){
        f.cdf <- ecdf(x[[j]][,l]); x.r[[j]][,l] <- f.cdf(x[[j]][,l])
      }
    }
    
    stat <- F(x,cc)       
    stat.r <- F(x.r,cc)
    stat.s <- Jdcov.sq.US.list(x,cc)   
    stat4 <- 0
    for(j in 1:(d-1)) stat4 <- stat4 + dcovU(x[[j]],do.call(cbind,x[(j+1):d]))^2  
    
    stat.p <- stat.pr <- stat.ps  <- stat.p4 <- rep(0, Mp)
    
    for(i in 1:Mp)
    {
      x.pr <- x.p <- x
      for(j in 1:d) x.p[[j]] <- x[[j]][sample(1:n,n,replace=TRUE),]
      # X.p <- split(x.p, rep(1:ncol(x.p), each = nrow(x.p)))
      for(j in 1:d) {
        for(l in 1:ncol(x[[j]])){
          f.cdf <- ecdf(x.p[[j]][,l]); x.pr[[j]][,l] <- f.cdf(x.p[[j]][,l])
        }
      }
      stat.p[i] <- F(x.p,cc)    ## okay
      stat.pr[i] <- F(x.pr,cc)
      stat.ps[i] <- Jdcov.sq.US.list(x.p,cc)  ## okay
      
      for(j in 1:(d-1)) stat.p4[i] <- stat.p4[i] + dcovU(x.p[[j]],do.call(cbind,x.p[(j+1):d]))^2
    }
    
    crit <- quantile(stat.p, c(0.90,0.95))
    crit.r <- quantile(stat.pr, c(0.90,0.95))
    crit.s <- quantile(stat.ps, c(0.90,0.95))
    crit4 <- quantile(stat.p4, c(0.90,0.95))
    
    rej[jj,] <- stat > crit
    rej.r[jj,] <- stat.r > crit.r
    rej.s[jj,] <- stat.s > crit.s
    rej4[jj,] <- stat4 > crit4
    
    out <- dhsic.test(x, alpha = 0.05, method = "bootstrap", kernel = "gaussian", B = Mp, pairwise = FALSE)
    rej3[jj,] <- c(out$p.value < 0.10, out$p.value < 0.05)  
    
    #if(jj%%50==0) print(jj)       
  }
  
  result <- rbind( c( apply(rej,2,mean), apply(rej.s,2,mean), apply(rej.r,2,mean), apply(rej3,2,mean), apply(rej4,2,mean) ) )
  return(result)
}