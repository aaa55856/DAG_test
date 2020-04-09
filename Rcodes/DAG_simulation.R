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
# 
# # if (!requireNamespace("BiocManager", quietly = TRUE))
# #   install.packages("BiocManager")
# # BiocManager::install("graph")
# # BiocManager::install("RBGL")
# # BiocManager::install("Rgraphviz")


## Reading the data

rm(list = ls())
setwd('D:/Code/Econometrics/DAG/working/Rcodes_JdCov/R codes')
# pima=read.table("D:/Distance Correlation  Research/pima-indians-diabetes.data.txt",sep="," ,header = FALSE)
source('./JdCov.R')

pima = read.table("./pima-indians-diabetes.data.txt",
                  sep=",", header = FALSE)

pima = as.data.frame(pima)
rn = c("npreg","glu","DBP","skin","SI","BMI","ped","age","Class")
colnames(pima) = rn
head(pima); dim(pima)
pima_new=pima[,c("age","BMI","SI","glu","DBP")]
head(pima_new)

## we do a centering and scaling of the original data

row_sub = apply(pima_new, 1, function(row) all(row !=0 ))
pima_final=pima_new[row_sub,]
head(pima_final) ; dim(pima_final)
pima_means=colMeans(pima_final)
pima_final_c=sweep(pima_final,2,pima_means) 
weights=apply(pima_final_c,2,function(x) (sqrt(sum(x^2)))/sqrt(nrow(pima_final_c)))
pima_final_cs=sweep(pima_final_c,2,weights,FUN = "/") 
head(pima_final_cs)
X <- list(X1=as.matrix(pima_final_cs[,1]),
          X2=as.matrix(pima_final_cs[,4]),
          X3=as.matrix(pima_final_cs[,5]))

# load("D:/Distance Correlation  Research/supplementary_code/causality_example/allDagsWith3Nodes.RData")
load("./allDagsWith3Nodes.RData")

d <- 3

## Generating 27 candidate DAG models for the three variables

allDags_num=1*allDags
r26=c(0,0,0,1,0,1,0,1,0)
r27=c(0,0,0,1,0,1,1,1,0)
allDags_num=rbind(allDags_num,r26,r27)
numDAGs <- nrow(allDags_num)

# initalize pvalue vector
pval_dhsic.sim <- vector("numeric",numDAGs)
pval_Jdcov.sim <- vector("numeric",numDAGs)
pval_JdcovUS.sim <- vector("numeric",numDAGs)
pval_JdcovR.sim <- vector("numeric",numDAGs)
pval_MS16.sim <- vector("numeric",numDAGs)

B=100 ; cc=1 ; n=392 ; m=100


Sim_val=list()
Sim_val[[1]]=Sim_val[[2]]=Sim_val[[3]]=matrix(0,nrow=n,ncol=1)
Actual_adj = matrix(allDags_num[15,], d, d,
                    dimnames=list(c("age","glu","DBP"),
                                  c("age","glu","DBP")))
maximalp_dhsic=numeric(m)
maximalp_Jdcov=numeric(m)
maximalp_JdcovUS=numeric(m)
maximalp_JdcovR=numeric(m)
maximalp_MS16=numeric(m)


## Repeating the model-diagnostic checking 100 times

for(jj in 1:m){
  
  
  t.start <- Sys.time()
  
  ## generate data for parent, i.e. Age
  j=1 
  set.seed(jj*7)
  Sim_val[[j]][,1] <- rnorm(n,0,1)        ## var_parent=1
  ## sd(X1) ==== 1.001278  (This is why we take var_parent=1) 
  
  
  ## generate data for the child DBP
  j=3
  X1=as.matrix(pima_final_cs[,1])
  X2=as.matrix(pima_final_cs[,4])
  X3=as.matrix(pima_final_cs[,5])

  X <- list(X1,X2,X3)
  formula <- paste("X",toString(j),"~",sep="")
  for(k in 1:d){
    if(Actual_adj[k,j]==1){
      if((sum(Actual_adj[,j])-sum(Actual_adj[1:k,j]))==0){
        formula <- paste(formula,"s(X",toString(k),")",sep="")
      }
      else{
        formula <- paste(formula,"s(X",toString(k),")+",sep="")
      }
    }
  }
  g=gam(as.formula(formula),data=X)
  
  ## sd(g$residuals) ==== 0.95 (This is why we take var_child=0.95 for DBP)
  
  X1=Sim_val[[1]][,1]
  set.seed(jj*8)
  Sim_val[[j]][,1]=predict.gam(g, data.frame(X1,X2=rep(0,n),X3=rep(0,n))) + rnorm(n,0,0.95)  ## var_child=0.95
                                                                                            ## for DBP
  
  ## generate data for the child glu
  
  j=2
  X1=as.matrix(pima_final_cs[,1])
  X2=as.matrix(pima_final_cs[,4])
  X3=as.matrix(pima_final_cs[,5])
  
  X <- list(X1,X2,X3)
  formula <- paste("X",toString(j),"~",sep="")
  for(k in 1:d){
    if(Actual_adj[k,j]==1){
      if((sum(Actual_adj[,j])-sum(Actual_adj[1:k,j]))==0){
        formula <- paste(formula,"s(X",toString(k),")",sep="")
      }
      else{
        formula <- paste(formula,"s(X",toString(k),")+",sep="")
      }
    }
  }
  g=gam(as.formula(formula),data=X)
  
  ## sd(g$residuals) ==== 0.918 (This is why we choose var_child=0.918 for glu)

  X1=Sim_val[[1]][,1]
  #X3=Sim_val[[3]][,1]
  set.seed(jj*9)
  Sim_val[[j]][,1]=predict.gam(g,data.frame(X1,X2=rep(0,n),X3=rep(0,n))) + rnorm(n,0,0.918)  ## var_child=0.918
                                                                                             ## for glu
  X=list(X1=Sim_val[[1]],X2=Sim_val[[2]],X3=Sim_val[[3]])
  
  ## we now have simulated data for a particular jj in 1:m
  ## now do model diagnostic checking for 27 candidate DAG models
  
  for(i in 1:numDAGs){
    print(i)
    
    A <- matrix(allDags_num[i,],d,d,
                dimnames=list(c("age","glu","DBP"),c("age","glu","DBP")))
    # initialize residual vector
    res <- list()
    res[[1]]=res[[2]]=res[[3]]=matrix(0,nrow=n,ncol=1)
    fitgam <- list()
    fitgam[[1]]=fitgam[[2]]=fitgam[[3]]=matrix(0,nrow=n,ncol=1)
    
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
      res.B[[1]]=res.B[[2]]=res.B[[3]]=matrix(0,nrow=n,ncol=1)
      res.new[[1]]=res.new[[2]]=res.new[[3]]=matrix(0,nrow=n,ncol=1)
      Y.B = list()
      
      for(j in 1:d){
        res.B[[j]]=as.matrix(res[[j]][sample(1:n,n,replace=TRUE),])
        Y.B[[j]]=fitgam[[j]]+res.B[[j]]
      }
      XB=list(XB1=Y.B[[1]],XB2=Y.B[[2]],XB3=Y.B[[3]])
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
    
    pval_dhsic.sim[i] = length(stat.pdhsic[stat.pdhsic > stat.dhsic]) /B 
    pval_Jdcov.sim[i] = length(stat.p[stat.p > stat]) /B
    pval_JdcovUS.sim[i] = length(stat.ps[stat.ps > stat.s]) /B
    pval_JdcovR.sim[i] = length(stat.pr[stat.pr > stat.r]) /B
    pval_MS16.sim[i] = length(stat.p4[stat.p4 > stat4]) /B
    
    
  }    
  
  maximalp_dhsic[jj]=which.max(pval_dhsic.sim)
  maximalp_Jdcov[jj]=which.max(pval_Jdcov.sim)
  maximalp_JdcovUS[jj]=which.max(pval_JdcovUS.sim)
  maximalp_JdcovR[jj]=which.max(pval_JdcovR.sim)
  maximalp_MS16[jj]=which.max(pval_MS16.sim)      
  
  print(Sys.time() - t.start)
}


# result_final = read.table("D:/Distance Correlation  Research/DAG_sim_final.txt",sep="",header=TRUE)
result_final = read.table("./DAG_sim_final.txt", 
                          sep="", header=TRUE)
result_final
dim(result_final)
names(which.max(table(result_final[,1])))           
names(which.max(table(result_final[,2])))
names(which.max(table(result_final[,3])))
names(which.max(table(result_final[,4])))
names(which.max(table(result_final[,5])))

sum(result_final[,1]==15)
sum(result_final[,2]==15)
sum(result_final[,3]==15)
sum(result_final[,4]==15)
sum(result_final[,5]==15) 

