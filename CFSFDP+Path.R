distmat_fast <- function(mat){
  smat <- apply(mat, 1, crossprod)
  mat1 <- matrix(smat, nrow=nrow(mat), ncol=nrow(mat))
  mat3 <- tcrossprod(mat)
  mat4 <- mat1 + t(mat1) - 2*mat3
  diag(mat4) <- 0
  mat5 <- sqrt(mat4)
  return(mat5)
}
euclid <- function(x,y){
  a <- sqrt(sum((x-y)^2))
  return(a)
}
distvec <- function(x,t){
  a <- apply(t,1,euclid,x=x)
  return(a)
}
mst_path_dist <- function(dm){
  n <- nrow(dm)
  visited <- rep(0,n)
  visited[1] <- 1 #记录数据点是否已进入最小生成树
  d <- dm[1,] #记录点到已进入生成树点的距离
  p <- 1  #记录点进入顺序
  pd <- diag(0,nrow = n)  #路径距离矩阵
  pe <- diag(0,nrow = n)  #决定路径距离的边
  while(sum(visited)<n){
    min=Inf
    for(i in 1:n){
      if(visited[i]==0&d[i]<min){
        min <- d[i]
        t <- i
      }
    } #找到离已经找过的点距离最小值
    p_visited <- which(visited==1)
    t2 <- which(dm[,t]==min & visited==1)[1]
    for(k in 1:sum(visited)){
      t1 <- p_visited[k]
      if(pd[t1,t2]<dm[t2,t]){
        pd[t1,t] <- dm[t2,t]
        pd[t,t1] <- pd[t1,t]
        pe[t1,t] <- paste(t2,t,sep = ',')
        pe[t,t1] <- pe[t1,t]
      }else{
        pd[t1,t] <- pd[t1,t2]
        pd[t,t1] <- pd[t1,t]
        pe[t1,t] <- pe[t1,t2]
        pe[t,t1] <- pe[t1,t]
      }
    }
    visited[t] <- 1
    p <- c(p,t)
    for(j in 1:n){
      if(visited[j]==0&d[j]>dm[t,j]){
        d[j] <- dm[t,j]
      }
    } #更新d向量
  }
  return(list(pd,pe))
}
cfsfdp_path <- function(data_set,k,r){
  n=nrow(data_set)
  dm <- distmat_fast(data_set)
  midu=colSums(dm<r)
  d <- mst_path_dist(dm)
  dm <- d[[1]]
  edge <- d[[2]]
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
  delta[midu<=floor(n/100)]=0
  centers_index=NULL
  centers_index <- which(delta==sort(delta,decreasing = T)[1])[1]
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
  a <- apply(dm[,centers_index],1,min)
  c <- apply(dm[,centers_index],1,which.min)
  for(i in 1:nrow(dm)){
    if(sum(dm[i,centers_index]==a[i])>1&a[i]!=0){
      e <- as.numeric(unlist(strsplit(edge[i,centers_index[c[i]]],','))) 
      if(min(dm[e[1],centers_index])<a[i]){
        e1 <- as.numeric(unlist(strsplit(edge[e[1],centers_index[c[e[1]]]],',')))
        if(min(dm[e1[1],centers_index])<min(dm[e[1],centers_index])){
          c[i] <- c[e1[1]]
        }else{
          c[i] <- c[e1[2]]
        }
      }else{
        e2 <- as.numeric(unlist(strsplit(edge[e[2],centers_index[c[e[2]]]],',')))
        if(min(dm[e2[1],centers_index])<min(dm[e[1],centers_index])){
          c[i] <- c[e2[1]]
        }else{
          c[i] <- c[e2[2]]
        }
      }
    }
  }
  return(list(centers_index,c,midu,delta))
}


library(factoextra)
data("multishapes")
df<-multishapes[,1:2]
df <- as.matrix(df)
c1 <- cfsfdp_path(df,5,0.19)
library(rgl)
plot3d(df[,1],df[,2],c1[[3]])
plot(df,col=c1[[2]])


library(kernlab)
data(spirals)
c2 <- cfsfdp_path(spirals,2,0.2)
plot(spirals,col=c2[[2]])


data(iris)
iris1 <- as.matrix(iris[,1:4])
c3 <- cfsfdp_path(iris1,3,0.55)
c3


library(MASS)
means <- c(0,0)
vars <- matrix(c(1,0,0,1),nrow = 2)
a_normal <- mvrnorm(10000,means,vars)
up_data <- a_normal[0.6<rowSums(a_normal^2)&rowSums(a_normal^2)<1&a_normal[,2]>0,]
a_normal_1 <- a_normal[,1]-1
a_normal_1 <- cbind(a_normal_1,a_normal[,2]-0.5)
low_data <- a_normal[0.8<rowSums(a_normal_1^2)&rowSums(a_normal_1^2)<1.2&a_normal[,2]<0.5,]
up_data_1 <- up_data[sample(1:nrow(up_data),80),]
low_data_1 <- low_data[sample(1:nrow(low_data),80),]
moon_data <- rbind(up_data_1,low_data_1)
plot(moon_data)
c4 <- cfsfdp_path(moon_data,2,0.3)
plot(moon_data,col=c4[[2]])
