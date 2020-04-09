# library(dHSIC)
# library(mgcv)
# library(Matrix)
# library(MASS)
# library(parallel)
# library(graph)
# library(RBGL)
# library(DAAG)
# library(splitstackshape)
# library(gRbase)
# library(Rgraphviz)


## ------------------------------ Pima Indians Diabetes -------------------------------- ##


## Reading the data

rm(list = ls())
setwd('D:/Code/Econometrics/DAG/working/Rcodes_JdCov/R codes')
# pima=read.table("D:/Distance Correlation  Research/pima-indians-diabetes.data.txt",sep="," ,header = FALSE)
source('./JdCov.R')

pima = read.table("./pima-indians-diabetes.data.txt",
                  sep=",", header = FALSE)
pima=as.data.frame(pima)
rn=c("npreg","glu","DBP","skin","SI","BMI","ped","age","Class")
colnames(pima)=rn
head(pima); dim(pima)
pima_new=pima[,c("age","BMI","SI","glu","DBP")]
head(pima_new)

## We do a centering a scaling

row_sub = apply(pima_new, 1, function(row) all(row !=0 ))
pima_final=pima_new[row_sub,]
head(pima_final) ; dim(pima_final)
pima_means=colMeans(pima_final)
pima_final_c=sweep(pima_final,2,pima_means) 
weights=apply(pima_final_c,2,function(x) (sqrt(sum(x^2)))/sqrt(nrow(pima_final_c)))
pima_final_cs=sweep(pima_final_c,2,weights,FUN = "/") 
head(pima_final_cs)


## Generating 48 candidate DAG models for the five variables

allDags_pima=rep(list(matrix(0, 5, 5, dimnames = list(c("age","BMI","SI","glu","DBP"),c("age","BMI","SI","glu","DBP")))), 48)
allDags_pima[[2]][1,2]=1
allDags_pima[[3]][1,3]=1
allDags_pima[[4]][1,4]=1
allDags_pima[[5]][1,5]=1
allDags_pima[[6]][1,c(2,3)]=1
allDags_pima[[7]][1,c(2,4)]=1
allDags_pima[[8]][1,c(2,5)]=1
allDags_pima[[9]][1,c(3,4)]=1
allDags_pima[[10]][1,c(3,5)]=1
allDags_pima[[11]][1,c(4,5)]=1
allDags_pima[[12]][1,c(2,3,4)]=1
allDags_pima[[13]][1,c(2,3,5)]=1
allDags_pima[[14]][1,c(3,4,5)]=1
allDags_pima[[15]][1,c(2,3,4,5)]=1
allDags_pima[[16]]=t(allDags_pima[[2]])
allDags_pima[[17]]=t(allDags_pima[[3]])
allDags_pima[[18]]=t(allDags_pima[[4]])
allDags_pima[[19]]=t(allDags_pima[[5]])
allDags_pima[[20]]=t(allDags_pima[[6]])
allDags_pima[[21]]=t(allDags_pima[[7]])
allDags_pima[[22]]=t(allDags_pima[[8]])
allDags_pima[[23]]=t(allDags_pima[[9]])
allDags_pima[[24]]=t(allDags_pima[[10]])
allDags_pima[[25]]=t(allDags_pima[[11]])
allDags_pima[[26]]=t(allDags_pima[[12]])
allDags_pima[[27]]=t(allDags_pima[[13]])
allDags_pima[[28]]=t(allDags_pima[[14]])
allDags_pima[[29]]=t(allDags_pima[[15]])
allDags_pima[[30]][1,c(2,3)]=1 ; allDags_pima[[30]][3,2]=1
allDags_pima[[31]][1,c(2,3)]=1 ; allDags_pima[[31]][2,3]=1
allDags_pima[[32]][1,c(2,3)]=1 ; allDags_pima[[32]][3,2]=1 ; allDags_pima[[32]][2,3]=1
allDags_pima[[33]][1,c(4,5)]=1 ; allDags_pima[[33]][4,5]=1
allDags_pima[[34]][1,c(4,5)]=1 ; allDags_pima[[34]][5,4]=1
allDags_pima[[35]][1,c(3,4)]=1 ; allDags_pima[[35]][3,4]=1
allDags_pima[[36]][1,c(3,4)]=1 ; allDags_pima[[36]][4,3]=1
allDags_pima[[37]][1,c(3,4)]=1 ; allDags_pima[[37]][4,3]=1 ; allDags_pima[[37]][3,4]=1
allDags_pima[[38]][1,c(3,4,5)]=1 ; allDags_pima[[38]][4,3]=1 ; allDags_pima[[38]][4,5]=1
allDags_pima[[39]][1,c(3,4,5)]=1 ; allDags_pima[[39]][3,4]=1 ; allDags_pima[[39]][4,5]=1
allDags_pima[[40]][1,c(3,4,5)]=1 ; allDags_pima[[40]][3,4]=1 ; allDags_pima[[40]][4,5]=1 ; allDags_pima[[40]][4,3]=1
allDags_pima[[41]][1,c(2,3,4,5)]=1 ; allDags_pima[[41]][3,2]=1 ; allDags_pima[[41]][3,4]=1 ; allDags_pima[[41]][4,5]=1
allDags_pima[[42]][1,c(2,3,4,5)]=1 ; allDags_pima[[42]][3,2]=1 ; allDags_pima[[42]][4,3]=1 ; allDags_pima[[42]][4,5]=1
allDags_pima[[43]][1,c(2,3,4,5)]=1 ; allDags_pima[[43]][3,2]=1 ; allDags_pima[[43]][3,4]=1 ; allDags_pima[[43]][4,3]=1 ; allDags_pima[[43]][4,5]=1
allDags_pima[[44]][1,c(2,3,5)]=1 ; allDags_pima[[44]][3,2]=1 ; allDags_pima[[44]][3,4]=1 ; allDags_pima[[44]][4,5]=1
allDags_pima[[45]][1,c(2,3,5)]=1 ; allDags_pima[[45]][3,2]=1 ; allDags_pima[[45]][4,3]=1 ; allDags_pima[[45]][4,5]=1
allDags_pima[[46]][1,c(2,3,5)]=1 ; allDags_pima[[46]][3,2]=1 ; allDags_pima[[46]][3,4]=1 ; allDags_pima[[46]][4,3]=1 ; allDags_pima[[46]][4,5]=1
allDags_pima[[47]]=matrix(c(0,1,1,1,1,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0),5,5,byrow=T,dimnames=list(c("age","BMI","SI","glu","DBP"),c("age","BMI","SI","glu","DBP")))
allDags_pima[[48]]=matrix(c(0,1,1,1,1,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0),5,5,byrow=T,dimnames=list(c("age","BMI","SI","glu","DBP"),c("age","BMI","SI","glu","DBP")))


X <- list(X1=as.matrix(pima_final_cs[,1]),X2=as.matrix(pima_final_cs[,2]),X3=as.matrix(pima_final_cs[,3]),X4=as.matrix(pima_final_cs[,4]),X5=as.matrix(pima_final_cs[,5]))

d <- 5
numDAGs_pima <- length(allDags_pima)

# initalize pvalue vector
pval_dhsic_p <- vector("numeric",numDAGs_pima)
pval_Jdcov_p <- vector("numeric",numDAGs_pima)
pval_JdcovUS_p <- vector("numeric",numDAGs_pima)
pval_JdcovR_p <- vector("numeric",numDAGs_pima)
pval_MS16_p <- vector("numeric",numDAGs_pima)


B=100 ; cc=1 ; n=392


# compute the p-value for each DAG model

for(i in 1:48){
  t.start <- Sys.time()
  
  A <- allDags_pima[[i]] 
  # initialize residual vector
  res <- list()
  res[[1]]=res[[2]]=res[[3]]=res[[4]]=res[[5]]=matrix(0,nrow=n,ncol=1)
  fitgam <- list()
  fitgam[[1]]=fitgam[[2]]=fitgam[[3]]=fitgam[[4]]=fitgam[[5]]=matrix(0,nrow=n,ncol=1)
  
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
      g=gam(as.formula(formula),data=X)
      res[[j]][,1] = g$residuals
      fitgam[[j]][,1] = g$fitted.values
    }
  }
  
  x=res
  ## computing values of the original statistics
  F <- Jdcov.sq.U.list
  stat.dhsic <- dhsic(X=x,kernel="gaussian")$dHSIC
  stat <- F(x,cc) 
  stat.s <- Jdcov.sq.US.list(x,cc)
  
  x.r <- x
  for(j in 1:d) {
    for(l in 1:ncol(x[[j]])){
      f.cdf <- ecdf(x[[j]][,l]); x.r[[j]][,l] <- f.cdf(x[[j]][,l])
    }
  }
  
  stat.r <- F(x.r,cc)
  stat4 <- 0
  for(b in 1:(d-1)) stat4 <- stat4 + dcovU(x[[b]],do.call(cbind,x[(b+1):d]))^2
  
  ## creating storage spaces for bootstrap resamples
  stat.pdhsic <- stat.p <- stat.pr <- stat.ps  <- stat.p4 <- rep(0, B)
  
  
  for(ii in 1:B)
  {
    res.B <- list() ; res.new <- list()
    res.B[[1]]=res.B[[2]]=res.B[[3]]=res.B[[4]]=res.B[[5]]=matrix(0,nrow=n,ncol=1)
    res.new[[1]]=res.new[[2]]=res.new[[3]]=res.new[[4]]=res.new[[5]]=matrix(0,nrow=n,ncol=1)
    Y.B = list()
    
    for(j in 1:d){
      res.B[[j]]=as.matrix(res[[j]][sample(1:n,n,replace=TRUE),])
      Y.B[[j]]=fitgam[[j]]+res.B[[j]]
    }
    XB=list(XB1=Y.B[[1]],XB2=Y.B[[2]],XB3=Y.B[[3]],XB4=Y.B[[4]],XB5=Y.B[[5]])
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
        #print(formula)
        gB=gam(as.formula(formula),data=XB)
        res.new[[j]][,1] <- gB$residuals
      }
    }
    
    x=res.new
    stat.pdhsic[ii] <- dhsic(X=x,kernel="gaussian")$dHSIC
    stat.p[ii] <- F(x,cc) 
    stat.ps[ii] <- Jdcov.sq.US.list(x,cc)
    
    x.r <- x
    for(j in 1:d) {
      for(l in 1:ncol(x[[j]])){
        f.cdf <- ecdf(x[[j]][,l]); x.r[[j]][,l] <- f.cdf(x[[j]][,l])
      }
    }
    
    stat.pr[ii] <- F(x.r,cc)
    stat.p4[ii] <- 0
    for(b in 1:(d-1)) stat.p4[ii] <- stat.p4[ii] + dcovU(x[[b]],do.call(cbind,x[(b+1):d]))^2
    
  }    
  
  pval_dhsic_p[i] = length(stat.pdhsic[stat.pdhsic > stat.dhsic]) /B 
  pval_Jdcov_p[i] = length(stat.p[stat.p > stat]) /B
  pval_JdcovUS_p[i] = length(stat.ps[stat.ps > stat.s]) /B
  pval_JdcovR_p[i] = length(stat.pr[stat.pr > stat.r]) /B
  pval_MS16_p[i] = length(stat.p4[stat.p4 > stat4]) /B
  
  print(Sys.time() - t.start)
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


M41p=as.matrix(allDags_pima[[41]])                               ## DAG model selected by dHSIC                     
DAG41p<-new("graphAM", adjMat=M41p, edgemode="directed")
plot(DAG41p)

M47p=as.matrix(allDags_pima[[47]])                               ## DAG model selected by Jdcov, JdcovUS,                       
DAG47p<-new("graphAM", adjMat=M47p, edgemode="directed")         ## JdcovR and T_MT
plot(DAG47p)

