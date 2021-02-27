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
reach_dist=function(p,o,k_dist,dm){
  return(max(k_dist[o],dm[p,o]))
}
LDBSCAN=function(data,LOFUB,pct,minpts){
  n=nrow(data)
  dm=distmat_f(data,data)
  diag(dm)=1/0
  k_dist=NULL
  k_dist_neig=list()
  for(i in 1:n){
    a=dm[i,]
    k_dist[i]=a[rank(a)==minpts]
    k_dist_neig[[i]]=which(a<=k_dist[i])
  }
  lrd=NULL
  lof=NULL
  for(i in 1:n){
    n_i=k_dist_neig[[i]]
    s=0
    for(j in n_i){
      s=s+reach_dist(i,j,k_dist,dm)
    }
    lrd[i]=length(n_i)/s
  }
  for(i in 1:n){
    n_i=k_dist_neig[[i]]
    lof[i]=sum(lrd[n_i])/(length(n_i)*lrd[i])
  }
  k=0
  c=rep(0,n)
  visited <- rep(0,n)
  noise <- rep(0,n)
  for(i in 1:n){
    if(visited[i]==0){
      n_i=k_dist_neig[[i]]
      m=length(n_i)
      visited[i]=1
      if(lof[i]<=LOFUB){
        k=k+1
        c[i]=k
        t=1
        ###############################
        ##        ÀàµÄÀ©ÕÅ           ##
        ###############################
        while(!all(visited[n_i]==1)){
          for(j in t:m){
            p <- n_i[j]
            if(visited[p]==0){
              visited[p] <- 1
              n_i_new <- k_dist_neig[[p]]
              m_new <- length(n_i_new)
              if(lof[p]<=LOFUB){
                n_i <- union(n_i,n_i_new)
                if(c[p]==0){
                  c[p] <- k
                }
              }
            }
          }
          t <- m+1
          m <- length(n_i)
        }
      }
    }
  }
  c[which(c==0)]=k+1
  return(c)
}
library(factoextra)
data("multishapes")
df<-multishapes[,1:2]
df1 <- as.matrix(df)
c1=LDBSCAN(df1,1.5,0.3,6)
plot(df1,col=c1)
plot(df1,col=multishapes[,3])
