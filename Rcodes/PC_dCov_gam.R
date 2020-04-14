## correlation to txt for TETRAD

rm(list = ls())
graphics.off()
setwd('D:/Code/Econometrics/DAG/working/Rcodes_JdCov/Rcodes')
source('./JdCov.R')

# df <- read.csv('./KAM2016-2018.csv')
# df <- df[df$year == 2018, 6:36]
# # 6:36
# write.csv(df, './KAM2018.csv', row.names = FALSE)

df <- read.csv('./KAM2018.csv')

# df <- read.csv(file='./pima.csv')

rn <- names(df)
n <- nrow(df)
d <- ncol(df)

df.cor <- cor(df)

# df.cor[upper.tri(df.cor)] <- NA
# write.table(nrow(df), './KAM2018_cor.txt', append=FALSE, 
#             row.names=FALSE, col.names=FALSE)
# cat(names(df), file='./KAM2018_cor.txt', sep=',', append=TRUE)
# cat('\n', file='./KAM2018_cor.txt', sep=',', append=TRUE)
# apply(df.cor, 1, function(x){
#   x <- x[!is.na(x)]
#   cat(x, file='./KAM2018_cor.txt', sep=',', append=TRUE)
#   cat('\n', file='./KAM2018_cor.txt', sep=',', append=TRUE)
# })


## PC
df.fit <- pc(suffStat = list(C = df.cor, n = n),
            indepTest = gaussCItest, 
            alpha=0.1, labels = rn, verbose = FALSE)
df.mat <- t(matrix(as.numeric(as(df.fit, "amat")), nrow = d))

## undirected
mat.bd <- df.mat==1 & t(df.mat)==1  ## bidirected edges matrix
df.mat[mat.bd] <- 0
mat.bd[upper.tri(mat.bd)] <- FALSE  ## only leave lower part
ind.bd <- which(mat.bd == TRUE)  ## bidirected edges index (lower matrix index)

allDags <- list()
for(i in 1:length(ind.bd)) {
  tmp <- combn(1:length(ind.bd), i)
  ind <- matrix(0, nrow=nrow(mat.bd), ncol=ncol(mat.bd))
  x <- ind.bd[tmp]
  ind[x] <- 1

  allDags <- cbind(allDags, list(as.numeric(df.mat + ind)), list(as.numeric(df.mat + t(ind))))
}

## PC
df.fit <- pc(suffStat = list(C = df.cor, n = n),
             indepTest = gaussCItest, skel.method = 'stable',
             alpha=0.1, labels = rn, verbose = FALSE)
df.mat <- t(matrix(as.numeric(as(df.fit, "amat")), nrow = d))

## undirected
mat.bd <- df.mat==1 & t(df.mat)==1  ## bidirected edges matrix
df.mat[mat.bd] <- 0
mat.bd[upper.tri(mat.bd)] <- FALSE  ## only leave lower part
ind.bd <- which(mat.bd == TRUE)  ## bidirected edges index (lower matrix index)

for(i in 1:length(ind.bd)) {
  tmp <- combn(1:length(ind.bd), i)
  ind <- matrix(0, nrow=nrow(mat.bd), ncol=ncol(mat.bd))
  x <- ind.bd[tmp]
  ind[x] <- 1
  
  allDags <- cbind(allDags, list(as.numeric(df.mat + ind)), list(as.numeric(df.mat + t(ind))))
}

allDags <- lapply(allDags, 
                       function(x) matrix(x, d, d, dimnames = list(rn, rn)))
numDAGs <- length(allDags)

pdf('./KAM2018.pdf')
for(i in 1:2) {
  plot(graphAM(adjMat = allDags[[i]], edgemode="directed"))
}
dev.off()

# ## pretreatment
# df.means = colMeans(df)
# pima_final_c=sweep(df,2,df.means) 
# weights=apply(pima_final_c,2,function(x) (sqrt(sum(x^2)))/sqrt(nrow(pima_final_c)))
# df=sweep(pima_final_c,2,weights,FUN = "/") 

X <- lapply(df, as.matrix)
names(X) <- paste0('X', 1:d) 

## initialize parameter
B=100 ; cc=1

# initalize pvalue vector
pval_dhsic_p <- vector("numeric",numDAGs)
pval_Jdcov_p <- vector("numeric",numDAGs)
pval_JdcovUS_p <- vector("numeric",numDAGs)
pval_JdcovR_p <- vector("numeric",numDAGs)
pval_MS16_p <- vector("numeric",numDAGs)

# compute the p-value for each DAG model

for(i in 1:numDAGs){
  print(i)
  
  A <- allDags[[i]] 
  # initialize residual vector
  # res <- list()
  # res[[1]]=res[[2]]=res[[3]]=res[[4]]=res[[5]]=matrix(0,nrow=n,ncol=1)
  # fitgam <- list()
  # fitgam[[1]]=fitgam[[2]]=fitgam[[3]]=fitgam[[4]]=fitgam[[5]]=matrix(0,nrow=n,ncol=1)
  res <- rep(list(matrix(0, nrow=n, ncol=1)), d)
  fitgam <- rep(list(matrix(0, nrow=n, ncol=1)), d)
  
  # regress each node on its variables and collect the residuals
  for(j in 1:d){
    if(sum(A[,j])==0){
      res[[j]][,1] <- X[[j]]
    }
    else{
      formula <- paste("X",toString(j),"~",sep="")
      for(k in 1:d){
        if(A[k,j]==1){
          if((sum(A[,j])-sum(A[1:k,j]))==0){
            formula <- paste(formula,"s(X",toString(k),")",sep="")
          }
          else{
            formula <- paste(formula,"s(X",toString(k),")+",sep="")
          }
        }
      }
      #print(formula)
      # print(formula)
      g <- gam(as.formula(formula), data=X)
      res[[j]][,1] <- g$residuals
      fitgam[[j]][,1] <- g$fitted.values
    }
  }
  
  x=res
  ## computing values of the original statistics
  stat.dhsic <- dhsic(X=x,kernel="gaussian")$dHSIC
  stat <- Jdcov.sq.U.list(x,cc) 
  stat.s <- Jdcov.sq.US.list(x,cc)
  
  x.r <- x
  for(j in 1:d) {
    for(l in 1:ncol(x[[j]])){
      f.cdf <- ecdf(x[[j]][,l]); x.r[[j]][,l] <- f.cdf(x[[j]][,l])
    }
  }
  
  stat.r <- Jdcov.sq.U.list(x.r,cc)
  stat4 <- 0
  for(b in 1:(d-1)) stat4 <- stat4 + dcovU(x[[b]],do.call(cbind,x[(b+1):d]))^2
  
  ## creating storage spaces for bootstrap resamples
  stat.pdhsic <- stat.p <- stat.pr <- stat.ps  <- stat.p4 <- rep(0, B)
  
  
  for(ii in 1:B) {
    
    # res.B <- list() ; res.new <- list()
    # res.B[[1]]=res.B[[2]]=res.B[[3]]=res.B[[4]]=res.B[[5]]=matrix(0,nrow=n,ncol=1)
    # res.new[[1]]=res.new[[2]]=res.new[[3]]=res.new[[4]]=res.new[[5]]=matrix(0,nrow=n,ncol=1)
    res.B <- rep(list(matrix(0,nrow=n,ncol=1)), d)
    res.new <- rep(list(matrix(0,nrow=n,ncol=1)), d)
    Y.B = list()
    
    for(j in 1:d){
      res.B[[j]]=as.matrix(res[[j]][sample(1:n,n,replace=TRUE),])
      Y.B[[j]]=fitgam[[j]]+res.B[[j]]
    }
    # XB=list(XB1=Y.B[[1]],XB2=Y.B[[2]],XB3=Y.B[[3]],XB4=Y.B[[4]],XB5=Y.B[[5]])
    XB <- Y.B
    names(XB) <- paste0('XB', 1:d)
    for(j in 1:d){
      if(sum(A[,j])==0){
        res.new[[j]][,1] <- XB[[j]]
      }
      else{
        formula <- paste("XB",toString(j),"~",sep="")
        for(k in 1:d){
          if(A[k,j]==TRUE){
            if((sum(A[,j])-sum(A[1:k,j]))==0){
              formula <- paste(formula,"s(XB",toString(k),")",sep="")
            }
            else{
              formula <- paste(formula,"s(XB",toString(k),")+",sep="")
            }
          }
        }
        # print(formula)
        gB=gam(as.formula(formula),data=XB)
        res.new[[j]][,1] <- gB$residuals
      }
    }
    
    x=res.new
    stat.pdhsic[ii] <- dhsic(X=x,kernel="gaussian")$dHSIC
    stat.p[ii] <- Jdcov.sq.U.list(x,cc) 
    stat.ps[ii] <- Jdcov.sq.US.list(x,cc)
    
    x.r <- x
    for(j in 1:d) {
      for(l in 1:ncol(x[[j]])){
        f.cdf <- ecdf(x[[j]][,l]); x.r[[j]][,l] <- f.cdf(x[[j]][,l])
      }
    }
    
    stat.pr[ii] <- Jdcov.sq.U.list(x.r,cc)
    stat.p4[ii] <- 0
    for(b in 1:(d-1)) stat.p4[ii] <- stat.p4[ii] + dcovU(x[[b]],do.call(cbind,x[(b+1):d]))^2
    
  }    
  
  pval_dhsic_p[i] = length(stat.pdhsic[stat.pdhsic > stat.dhsic]) /B 
  pval_Jdcov_p[i] = length(stat.p[stat.p > stat]) /B
  pval_JdcovUS_p[i] = length(stat.ps[stat.ps > stat.s]) /B
  pval_JdcovR_p[i] = length(stat.pr[stat.pr > stat.r]) /B
  pval_MS16_p[i] = length(stat.p4[stat.p4 > stat4]) /B
  
}

pval_dhsic_p
pval_Jdcov_p
pval_JdcovUS_p
pval_JdcovR_p
pval_MS16_p

max(pval_dhsic_p) ; which.max(pval_dhsic_p)      ## returning the maximum p-values and the corresponding models
max(pval_Jdcov_p) ; which.max(pval_Jdcov_p)
max(pval_JdcovUS_p) ; which.max(pval_JdcovUS_p)
max(pval_JdcovR_p) ; which.max(pval_JdcovR_p)
max(pval_MS16_p) ; which.max(pval_MS16_p)

max.idx <- c(which.max(which.max(pval_dhsic_p)),
             which.max(pval_Jdcov_p),
             which.max(pval_JdcovUS_p),
             which.max(pval_JdcovR_p),
             which.max(pval_MS16_p))
cri.idx <- c('dHIC', 'JdCov', 'JdCovUS', 'JdCovR', 'MT')
par(mfrow = c(3, 2))
for(i in 1:length(max.idx)) {
  plot(graphAM(adjMat = allDags[[max.idx[i]]], edgemode="directed"),
       main = cri.idx[i])
}
par(mfrow = c(1, 1))
