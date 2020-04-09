time.start <- Sys.time()

##### install.packages("dHSIC")
##### install.packages("energy")

library(dHSIC)
library(Matrix)
library(energy)
library(mvtnorm)



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


## ------------------------------------------------------------------------------------------------------ ##



##### Ex 5.1

Sim.size <- function(n, d, cc=1, index=1, type, Mp=500, reptim=1000, Ustat=FALSE)
{
  set.seed(100)
  rej  <- rej.r <- rej.s <- rej3 <- rej4 <- matrix(NA, reptim, 2)
  mu=rep(0,d) ; sigma=diag(d)
  if(Ustat) F <- Jdcov.sq.U else F <- Jdcov.sq.V
  for(jj in 1:reptim)
  {
    set.seed(jj)
    z=rmvnorm(n,mu,sigma)
    if(type==0) x <- z
    if(type==1) x <- z^3
    if(type==2) x <- sign(z)*abs(z)^(1/3)
    
    X <- split(x, rep(1:ncol(x), each = nrow(x)))
    x.r <- x
    for(j in 1:d) { f.cdf <- ecdf(x[,j]); x.r[,j] <- f.cdf(x[,j]) }
    
    stat <- F(x,cc)
    stat.r <- F(x.r,cc)
    stat.s <- Jdcov.sq.US(x,cc)
    stat4 <- 0
    for(j in 1:(d-1)) stat4 <- stat4 + dcovU(x[,j],x[,(j+1):d])^2
    
    stat.p <- stat.pr <- stat.ps  <- stat.p4 <- rep(0, Mp)
    
    for(i in 1:Mp)
    {
      x.pr <- x.p <- x
      for(j in 1:d) x.p[,j] <- x[sample(1:n,n,replace=TRUE),j]
      X.p <- split(x.p, rep(1:ncol(x.p), each = nrow(x.p)))
      for(j in 1:d) { f.cdf <- ecdf(x.p[,j]); x.pr[,j] <- f.cdf(x.p[,j]) }
      stat.p[i] <- F(x.p,cc)
      stat.pr[i] <- F(x.pr,cc)
      stat.ps[i] <- Jdcov.sq.US(x.p,cc)
      
      for(j in 1:(d-1)) stat.p4[i] <- stat.p4[i] + dcovU(x.p[,j],x.p[,(j+1):d])^2
    }
    
    crit <- quantile(stat.p, c(0.90,0.95))
    crit.r <- quantile(stat.pr, c(0.90,0.95))
    crit.s <- quantile(stat.ps, c(0.90,0.95))
    crit4 <- quantile(stat.p4, c(0.90,0.95))
    
    rej[jj,] <- stat > crit
    rej.r[jj,] <- stat.r > crit.r
    rej.s[jj,] <- stat.s > crit.s
    rej4[jj,] <- stat4 > crit4
    
    out <- dhsic.test(X, alpha = 0.05, method = "bootstrap", kernel = "gaussian", B = Mp, pairwise = FALSE)
    rej3[jj,] <- c(out$p.value < 0.10, out$p.value < 0.05)  
    
  }
  
  result <- rbind( c( apply(rej,2,mean), apply(rej.s,2,mean), apply(rej.r,2,mean), apply(rej3,2,mean), apply(rej4,2,mean) ) )
  return(result)
}



## n=100, d=5, Mp=500, reptim=1000

set.seed(801)
Eg1.1a = Sim.size(n=100, d=5, cc=1, index=1, type=0, Mp=500, reptim=1000, Ustat=TRUE)
set.seed(802)
Eg1.2a = Sim.size(n=100, d=5, cc=1, index=1, type=2, Mp=500, reptim=1000, Ustat=TRUE)
set.seed(803)
Eg1.3a = Sim.size(n=100, d=5, cc=1, index=1, type=1, Mp=500, reptim=1000, Ustat=TRUE)


## n=100, d=10, Mp=500, reptim=1000

set.seed(810)
Eg1.1b = Sim.size(n=100, d=10, cc=1, index=1, type=0, Mp=500, reptim=1000, Ustat=TRUE)
set.seed(811)
Eg1.2b = Sim.size(n=100, d=10, cc=1, index=1, type=2, Mp=500, reptim=1000, Ustat=TRUE)
set.seed(812)
Eg1.3b = Sim.size(n=100, d=10, cc=1, index=1, type=1, Mp=500, reptim=1000, Ustat=TRUE)



## ------------------------------------------------------------------------------------------------------ ##



## Example 5.2


Sim.size8.2 <- function(n, d, rho, cc=1, index=1, type, Mp=500, reptim=1000, Ustat=FALSE)
{
  set.seed(100)
  rej  <- rej.r <- rej.s <- rej3 <- rej4 <- matrix(NA, reptim, 2)
  mu=rep(0,d) ; sigma=matrix(0,d,d)
  if(Ustat) F <- Jdcov.sq.U else F <- Jdcov.sq.V
  
  if(type==1) sigma=toeplitz(rho^(0:(d-1)))
  
  if(type==2) {
    for(i in 1:d){
      for(j in 1:d){
        if(i==j) sigma[i,j]=1
        if(abs(i-j)>=1 & abs(i-j)<=2) sigma[i,j]=rho
      }
    }
  }
  
  if(type==3) {
    sig.block=toeplitz(rho^(c(0,rep(1,5-1))))
    sigma=kronecker(diag(d/5),sig.block)
  }
  
  for(jj in 1:reptim)
  {
    set.seed(jj)
    x=rmvnorm(n,mu,sigma)
    
    X <- split(x, rep(1:ncol(x), each = nrow(x)))
    x.r <- x
    for(j in 1:d) { f.cdf <- ecdf(x[,j]); x.r[,j] <- f.cdf(x[,j]) }
    
    stat <- F(x,cc)
    stat.r <- F(x.r,cc)
    stat.s <- Jdcov.sq.US(x,cc)
    stat4 <- 0
    for(j in 1:(d-1)) stat4 <- stat4 + dcovU(x[,j],x[,(j+1):d])^2
    
    stat.p <- stat.pr <- stat.ps  <- stat.p4 <- rep(0, Mp)
    
    for(i in 1:Mp)
    {
      x.pr <- x.p <- x
      for(j in 1:d) x.p[,j] <- x[sample(1:n,n,replace=TRUE),j]
      X.p <- split(x.p, rep(1:ncol(x.p), each = nrow(x.p)))
      for(j in 1:d) { f.cdf <- ecdf(x.p[,j]); x.pr[,j] <- f.cdf(x.p[,j]) }
      stat.p[i] <- F(x.p,cc)
      stat.pr[i] <- F(x.pr,cc)
      stat.ps[i] <- Jdcov.sq.US(x.p,cc)
      
      for(j in 1:(d-1)) stat.p4[i] <- stat.p4[i] + dcovU(x.p[,j],x.p[,(j+1):d])^2
    }
    
    crit <- quantile(stat.p, c(0.90,0.95))
    crit.r <- quantile(stat.pr, c(0.90,0.95))
    crit.s <- quantile(stat.ps, c(0.90,0.95))
    crit4 <- quantile(stat.p4, c(0.90,0.95))
    
    rej[jj,] <- stat > crit
    rej.r[jj,] <- stat.r > crit.r
    rej.s[jj,] <- stat.s > crit.s
    rej4[jj,] <- stat4 > crit4
    
    out <- dhsic.test(X, alpha = 0.05, method = "bootstrap", kernel = "gaussian", B = Mp, pairwise = FALSE)
    rej3[jj,] <- c(out$p.value < 0.10, out$p.value < 0.05)  
    
  }
  
  result <- rbind( c( apply(rej,2,mean), apply(rej.s,2,mean), apply(rej.r,2,mean), apply(rej3,2,mean), apply(rej4,2,mean) ) )
  return(result)
}



## n=100, d=5, Mp=500, reptim=1000, rho=0.25

set.seed(500)
Eg2.1a = Sim.size8.2(n=100, d=5, rho=0.25, cc=1, index=1, type=1, Mp=500, reptim=1000, Ustat=TRUE)
set.seed(501)
Eg2.2a = Sim.size8.2(n=100, d=5, rho=0.25, cc=1, index=1, type=2, Mp=500, reptim=1000, Ustat=TRUE)
set.seed(502)
Eg2.3a = Sim.size8.2(n=100, d=5, rho=0.25, cc=1, index=1, type=3, Mp=500, reptim=1000, Ustat=TRUE)


## n=100, d=10, Mp=500, reptim=1000, rho=0.25

set.seed(600)
Eg2.1b = Sim.size8.2(n=100, d=10, rho=0.25, cc=1, index=1, type=1, Mp=500, reptim=1000, Ustat=TRUE)
set.seed(601)
Eg2.2b = Sim.size8.2(n=100, d=10, rho=0.25, cc=1, index=1, type=2, Mp=500, reptim=1000, Ustat=TRUE)
set.seed(602)
Eg2.3b = Sim.size8.2(n=100, d=10, rho=0.25, cc=1, index=1, type=3, Mp=500, reptim=1000, Ustat=TRUE)




## ------------------------------------------------------------------------------------------------------ ##



#### Example 5.3

Sim.size8.4 <- function(n, d=3, cc=1, index=1, type=1, Mp=500, reptim=1000, Ustat=FALSE)
{
  set.seed(100)
  rej  <- rej.r <- rej.s <- rej3 <- rej4 <- matrix(NA, reptim, 2)
  if(Ustat) F <- Jdcov.sq.U else F <- Jdcov.sq.V
  for(jj in 1:reptim)
  {
    if(type==1){
      set.seed(jj)
      z=rmvnorm(n,rep(0,2),diag(2))
      set.seed(jj*101)
      x=cbind(z,( sign(z[,1]*z[,2]) * rexp(n,1/sqrt(2)) ) )
    }
    
    if(type==2){
      set.seed(jj)
      z1 <- rbinom(n, 1, 0.5)
      set.seed(jj*100)
      z2 <- rbinom(n, 1, 0.5)
      z3 <- 1-abs(z1-z2)
      x <- matrix(cbind(z1,z2,z3), n, d)
    }
    
    X <- split(x, rep(1:ncol(x), each = nrow(x)))
    x.r <- x
    for(j in 1:d) { f.cdf <- ecdf(x[,j]); x.r[,j] <- f.cdf(x[,j]) }
    
    stat <- F(x,cc)
    stat.r <- F(x.r,cc)
    stat.s <- Jdcov.sq.US(x,cc)
    stat4 <- 0
    for(j in 1:(d-1)) stat4 <- stat4 + dcovU(x[,j],x[,(j+1):d])^2
    
    stat.p <- stat.pr <- stat.ps  <- stat.p4 <- rep(0, Mp)
    
    for(i in 1:Mp)
    {
      x.pr <- x.p <- x
      for(j in 1:d) x.p[,j] <- x[sample(1:n,n,replace=TRUE),j]
      X.p <- split(x.p, rep(1:ncol(x.p), each = nrow(x.p)))
      for(j in 1:d) { f.cdf <- ecdf(x.p[,j]); x.pr[,j] <- f.cdf(x.p[,j]) }
      stat.p[i] <- F(x.p,cc)
      stat.pr[i] <- F(x.pr,cc)
      stat.ps[i] <- Jdcov.sq.US(x.p,cc)
      
      for(j in 1:(d-1)) stat.p4[i] <- stat.p4[i] + dcovU(x.p[,j],x.p[,(j+1):d])^2
    }
    
    crit <- quantile(stat.p, c(0.90,0.95))
    crit.r <- quantile(stat.pr, c(0.90,0.95))
    crit.s <- quantile(stat.ps, c(0.90,0.95))
    crit4 <- quantile(stat.p4, c(0.90,0.95))
    
    rej[jj,] <- stat > crit
    rej.r[jj,] <- stat.r > crit.r
    rej.s[jj,] <- stat.s > crit.s
    rej4[jj,] <- stat4 > crit4
    
    out <- dhsic.test(X, alpha = 0.05, method = "bootstrap", kernel = "gaussian", B = Mp, pairwise = FALSE)
    rej3[jj,] <- c(out$p.value < 0.10, out$p.value < 0.05)  
    
  }
  
  result <- rbind( c( apply(rej,2,mean), apply(rej.s,2,mean), apply(rej.r,2,mean), apply(rej3,2,mean), apply(rej4,2,mean) ) )
  return(result)
}

## n=100, d=3, Mp=500, reptim=1000

set.seed(5000)
Ex3i = Sim.size8.4(n=100, d=3, cc=1, index=1, type=1, Mp=500, reptim=1000, Ustat=TRUE)

set.seed(5001)
Exbinom = Sim.size8.4(n=100, d=3, cc=1, index=1, type=2, Mp=500, reptim=1000, Ustat=TRUE)


Sys.time() - time.start