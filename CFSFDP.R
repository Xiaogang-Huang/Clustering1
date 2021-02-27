distvec <- function(x,y){
  a <- as.vector(crossprod(x))
  b <- apply(y,1,crossprod)
  c <- as.vector(tcrossprod(x,y))
  return(sqrt(a+b-2*c))
}
distmat_fast <- function(mat){
  smat <- apply(mat, 1, crossprod)
  mat1 <- matrix(smat, nrow=nrow(mat), ncol=nrow(mat))
  mat3 <- tcrossprod(mat)
  mat4 <- mat1 + t(mat1) - 2*mat3
  diag(mat4) <- 0
  mat5 <- sqrt(mat4)
  return(mat5)
}
delta <- function(data_set,r){
  dm <- distmat_fast(data_set)
  midu=colSums(dm<r)
  delta=NULL
  for (i in 1:nrow(data_set)) {
    d <- dm[i,]
    if(midu[i]==max(midu)){
      delta[i]=max(d)
    }
    else{
      delta[i]=min(d[which(midu>midu[i])]) 
    }
  }
  return(list(delta,midu,dm))
}
find_center <- function(data_set,midu,delta,k,theta){
  a=delta
  a[midu<=2]=0
  centers_index=NULL
  centers_index <- which(a==sort(a,decreasing = T)[1])
  i <- 2
  while(length(centers_index)<k){
    centers_candidate <- which(delta==sort(delta,decreasing = T)[i])
    if(which.max(delta) %in% centers_candidate){centers_candidate <- centers_candidate[-1]}
    for(j in 1:length(centers_candidate)){
      if(length(centers_index)==1){t=t(as.matrix(data_set[centers_index,]))}else{
        t <- data_set[centers_index,]
      }
      cen_dist <- distvec(data_set[centers_candidate[j],],t)
      if(any(cen_dist<r)) {}else{
        centers_index <- c(centers_index ,centers_candidate[j])
      }
    }
    i <- i+length(centers_candidate)
  }
  if(length(centers_index)>k){centers_index <- centers_index[1:k]}
  return(centers_index)
}
assigned <- function(data_set,midu,centers_index,dm){
  assignment <- rep(0,nrow(data_set))
  assignment[centers_index] <- 1:length(centers_index)
  m <- midu
  outliers_index <- which(m==1)
  m[centers_index] <- 0
  i <- 1
  t <- 1
  while(i<nrow(data_set)){
    tth_midu <- sort(m,decreasing = T)[t]
    num <- which(m==tth_midu)
    for(j in 1:length(num)){
      if(num[j] %in% centers_index || num[j] %in% outliers_index){}else{
        dist <- dm[num[j],]
        dist[num[j]] <- 1/0
        dist[midu<tth_midu] <- 1/0
        x <- cbind(dist,assignment)
        v <- x[order(x[,1]),][,2]
        assignment[num[j]] <- v[which(v!=0)[1]]
      }
    }
    i <- i+length(num)
    t <- t+length(num)
  }
  assignment[outliers_index] <- length(centers_index)+1
  return(assignment)
}

#例子
library(rgl)
library(factoextra)
data("multishapes")
df<-multishapes[,1:2]
df=as.matrix(df)
mdelta <- delta(df,r=0.15)
mdelta <- delta(df,r=0.2)
md <- mdelta[[2]]
del <- mdelta[[1]]
plot(md,del)
abline(h=0.2)
dm <- mdelta[[3]]
centers_index <- find_center(df,md,del,5,theta=1)
fclust <- assigned(df,midu=md,centers_index = centers_index,dm)
plot(df,col=fclust)
data_set_midu <- as.data.frame(cbind(df,mmidu)) 


#spirals数据集
library(kernlab)
data(spirals)
mdelta <- delta(spirals,r=0.2)
md <- mdelta[[2]]
del <- mdelta[[1]]
plot(md,del)
plot3d(spirals[,1],spirals[,2],md)
dm <- mdelta[[3]]
centers_index <- find_center(spirals,md,del,2,theta=0.2)
fclust <- assigned(spirals,midu=md,centers_index = centers_index,dm)
plot(spirals,col=fclust)
