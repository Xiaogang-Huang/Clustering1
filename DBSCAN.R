#¼ÆËã¾àÀë¾ØÕó
distmat_fast <- function(mat){
  smat <- apply(mat, 1, crossprod)
  mat1 <- matrix(smat, nrow=nrow(mat), ncol=nrow(mat))
  mat3 <- tcrossprod(mat)
  mat4 <- mat1 + t(mat1) - 2*mat3
  diag(mat4) <- 0
  mat5 <- sqrt(mat4)
  return(mat5)
}
#dbscanËã·¨
DBSCAN <- function(x,eps,minpts){
  n <- nrow(x)
  c <- rep(0,n)
  k <- 0
  visited <- rep(0,n)
  noise <- rep(0,n)
  dm <- distmat_fast(x)
  for(i in 1:n){
    if(visited[i]==0){
      visited[i] <- 1
      neighborpts <- which(dm[i,]<eps)
      m <- length(neighborpts)
      if(m<minpts){
        noise[i] <- 1
      }else{
        k <- k+1
        c[i] <- k
        t <- 1
        while(!all(visited[neighborpts]==1)){
          for(j in t:m){
            p <- neighborpts[j]
            if(visited[p]==0){
              visited[p] <- 1
              neighborpts_new <- which(dm[p,]<eps)
              m_new <- length(neighborpts_new)
              if(m_new>=minpts){
                neighborpts <- union(neighborpts,neighborpts_new)
                if(c[p]==0){
                  c[p] <- k
                }
              }else{
                noise[p] <- 1
                c[p] <- k
              }
            }
          }
          t <- m+1
          m <- length(neighborpts)
        }
      }
    }
  }
  c[which(noise==1&c==0)] <- k+1
  return(c)
}

library(factoextra)
data("multishapes")
df<-multishapes[,1:2]
df <- as.matrix(df)
cl <- DBSCAN(df,0.15,5)
plot(df,col=cl)
