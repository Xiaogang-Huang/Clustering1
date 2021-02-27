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

findcluster=function(d){
  cl=NULL
  k=0
  n=nrow(d)
  check=rep(0,n)
  for(i in 1:n){
    s=NULL
    m=1
    if(check[i]==0){
      check[i]=1
      for(j in 1:n){
        if(all(d[j,]==d[i,])){
          s[m]=j
          check[j]=1
          m=m+1
        }
      }
      if(m>1){
        k=k+1
        cl[c(i,s)]=k
      }
    }
  }
  return(cl)
}

dynamicclustering=function(x,eps){
  n=nrow(x)
  t=0
  r_c=NULL
  r_c_temp=NULL
  converge=0
  x_new=x
  while(converge==0){
    dm=distmat_f(x,x)
    for(i in 1:n){
      neighbors=which(dm[i,]<eps)
      m=length(neighbors)
      if(m>1){
        x_new[i,]=x[i,]+colMeans(sin(x[neighbors,]-x[i,]))
        r_c_temp[i]=mean(exp(-distmat_f(x[neighbors,],x[i,])))
      }else{
        x_new[i,]=x[i,]+sin(x[neighbors,]-x[i,])
        r_c_temp[i]=mean(exp(-distmat_f(x[neighbors,],x[i,])))
      }
    }
    t=t+1
    r_c[t]=sum(r_c_temp)/n
    if(t>1){
      if(r_c[t]-r_c[t-1]<10^-5){
        converge=1
      }
    }
    x=x_new
  }
  findcluster(x)
}


library(tidyverse)
library(MASS)

##1
gau1=mvrnorm(100,Sigma = matrix(c(.1,0,0,.1),2,2),mu = c(0,0))
gau2=mvrnorm(100,Sigma = matrix(c(.1,0,0,.1),2,2),mu = c(5,10))
gau3=mvrnorm(100,Sigma = matrix(c(.1,0,0,.1),2,2),mu = c(-3,9))
gau=rbind(gau1,gau2,gau3)
colnames(gau)= c('x','y')
gau_tb=as.tibble(gau)
ggplot(gau_tb,aes(x,y))+geom_point()
cl=dynamicclustering(gau,1)
