distmat_f <- function(x,y){
  if(!is.matrix(x)){
    x <- t(as.matrix(x))
  }
  if(!is.matrix(y)){
    y <- t(as.matrix(y))
  }
  n <- nrow(x)
  m <- nrow(y)
  xmat <- apply(x, 1, crossprod)
  ymat <- apply(y,1,crossprod)
  mat1 <- matrix(xmat, nrow=n, ncol=m)
  mat2 <- matrix(ymat,nrow =n,ncol = m,byrow = T )
  mat3 <- tcrossprod(x,y)
  mat4 <- mat1 + mat2 - 2*mat3
  mat5 <- sqrt(mat4)
  return(mat5)
}

k_means <- function(x,k,c){
  n=nrow(x)
  m=ncol(x)
  ce=as.matrix(aggregate(x,list(c),mean)[,-1])
  ic=matrix(nrow=k,ncol=ncol(x))
  while(identical(ic,ce)==F){
    ic= ce
    dm=distmat_f(x,ce)
    c=apply(dm,1,which.min)
    ce=as.matrix(aggregate(x,list(c),mean)[,-1])
  }
  return(c)
}

kmeans_star=function(x,k,steps=20){
  n=nrow(x)
  m=ncol(x)
  x_new=matrix(nrow=k,ncol=m)
  if(m>2){
    x_new[,-m]=colMeans(x[,-m])
  }else if(m!=1){
    x_new[,-m]=mean(x[,-m])
  }
  trim=quantile(x[,m],c(.05,.95))
  x_new[,m]=seq(trim[1],trim[2],length=k)
  c=sample(1:k,n,replace = T)
  x_transform=matrix(nrow=n,ncol=m)
  for(i in 1:steps){
    for(j in 1:n){
      x_transform[j,]=x_new[c[j],]+i*(x[j,]-x_new[c[j],])/steps
    }
    c=k_means(x_transform,k,c)
  }
  c
}

library(tidyverse)
library(MASS)

##1
gau1=mvrnorm(100,Sigma = matrix(c(1,0,0,1),2,2),mu = c(0,0))
gau2=mvrnorm(100,Sigma = matrix(c(1,1,2,2),2,2),mu = c(5,10))
gau3=mvrnorm(200,Sigma = matrix(c(3,0,0,3),2,2),mu = c(-3,9))
gau=rbind(gau1,gau2,gau3)
colnames(gau)= c('x','y')
gau_tb=as.tibble(gau)
c=kmeans_star(gau,3,20)
ggplot(gau_tb,aes(x,y))+geom_point(aes(color=as.character(c)))
