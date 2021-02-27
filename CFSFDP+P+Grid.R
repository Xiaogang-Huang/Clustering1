#计算距离矩阵
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
#计算距离向量
distvec <- function(x,y){
  a <- as.vector(crossprod(x))
  b <- apply(y,1,crossprod)
  c <- as.vector(tcrossprod(x,y))
  return(sqrt(a+b-2*c))
}
#x是二维数据
grid_divide = function(x,k){
  n <- nrow(x)
  x_axis_low <-  min(x[,1])
  x_axis_up <-  max(x[,1])
  y_axis_low <- min(x[,2])
  y_axis_up <- max(x[,2])
  d_x <- (x_axis_up-x_axis_low)/k
  d_y <- (y_axis_up-y_axis_low)/k
  grid_index <- NULL
  for(i in 1:n){
    dot <- x[i,]
    grid_index[i] <- floor((dot[1]-x_axis_low)/d_x)+floor((dot[2]-y_axis_low)/d_y)*k+1
  }
  return(grid_index)
}
mst_path_dist <- function(dm){
  n <- nrow(dm)
  visited <- rep(0,n)
  visited[1] <- 1 #记录数据点是否已进入最小生成树
  d <- dm[1,] #记录点到已进入生成树点的距离
  p <- 1  #记录点进入顺序
  e <- matrix(nrow=n-1,ncol = 2) #按进入顺序记录生成树边
  pd <- diag(0,nrow = n)  #路径距离矩阵
  pe <- diag(0,nrow = n)  #决定路径距离的边
  di <- NULL #记录路径长度
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
    e[sum(visited),2] <- t
    e[sum(visited),1] <- sort(p)[which(dm[visited==1,t]==min)[1]]
    di <- c(di,min)
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
  return(list(pd,pe,e,di,p))
}
#网格减少时间复杂度的密度峰算法
dp_grid <- function(x,r,k,kn){
  n <- nrow(x)
  c <- grid_divide(x,kn)
  ce <- as.matrix(aggregate(x,list(c),mean)[,-1])
  dm <- distmat_f(x,ce)
  c1 <- apply(dm,1,which.min)
  h <- nrow(ce)
  extracts <- NULL
  for(t in 1:h){
    if(sum(c1==t)==0){
      extracts <- c(extracts,t)
    }
  }
  if(length(extracts)>0){
    dm <- dm[,-extracts]
    ce <- ce[-extracts,]
  }
  c <- apply(dm,1,which.min)
  k1 <- nrow(ce)
  md <- colSums(dm<r)
  dm_cen <- diag(0,k1)
  for(i in 1:(k1-1)){
    cq_to_zq <- dm[c==i,i]
    for(j in (i+1):k1){
      cq_to_zl <- dm[c==i,j]
      cl_to_zl <- dm[c==j,j]
      cl_to_zq <- dm[c==j,i]
      cq_ind <- which(rank(abs(cq_to_zq-cq_to_zl))<3)
      cl_ind <- which(rank(abs(cl_to_zq-cl_to_zl))<3)
      cq_point <- x[which(c==i)[cq_ind],]
      cl_point <- x[which(c==j)[cl_ind],]
      d1 <- distmat_f(cq_point,cl_point)
      dm_cen[i,j] <- min(d1)+max(d1)
    }
  }
  dm_cen <- dm_cen+t(dm_cen)
  dm_path <- mst_path_dist(dm_cen)
  dm_ce <- dm_path[[1]]
  edge <- dm_path[[2]]
  tr <- dm_path[[3]]
  di <- dm_path[[4]]
  p <- dm_path[[5]]
  delta <- NULL
  for(i in 1:k1){
    d <- dm_ce[i,]
    if(md[i]==max(md)){
      delta[i]=max(d)
    }else{
      delta[i]=min(d[which(md>md[i])])
    }
  }
  outliers <- which(table(c)<2)
  if(length(outliers)>0){
    delta[outliers]=0
  }
  delta[md<=12]=0
  centers_index=NULL
  centers_index <- which(delta==sort(delta,decreasing = T)[1])[1]
  i <- 2
  while(length(centers_index)<k){
    centers_candidate <- which(delta==sort(delta,decreasing = T)[i])
    if(which.max(delta) %in% centers_candidate){centers_candidate <- centers_candidate[-1]}
    for(j in 1:length(centers_candidate)){
      if(length(centers_index)==1){
        t=t(as.matrix(ce[centers_index,]))
        cen_dist <- sqrt(sum((ce[centers_candidate[j],]-t)^2))
      }else{
        t <- ce[centers_index,]
        cen_dist <- distvec(ce[centers_candidate[j],],t)
      }
      if(any(cen_dist<r)) {}else{
        centers_index <- c(centers_index ,centers_candidate[j])
      }
    }
    i <- i+length(centers_candidate)
  }
  if(length(centers_index)>k){centers_index <- centers_index[1:k]}
  a <- apply(dm_ce[,centers_index],1,min)
  c_prot <- apply(dm_ce[,centers_index],1,which.min)
  for(i in 1:k1){
    if(sum(dm_ce[i,centers_index]==a[i])>1&a[i]!=0&!(i %in% outliers)){
      e <- as.numeric(unlist(strsplit(edge[i,centers_index[c_prot[i]]],','))) 
      if(min(dm_ce[e[1],centers_index])<a[i]){
        e1 <- as.numeric(unlist(strsplit(edge[e[1],centers_index[c_prot[e[1]]]],',')))
        if(min(dm_ce[e1[1],centers_index])<min(dm_ce[e[1],centers_index])){
          c_prot[i] <- c_prot[e1[1]]
        }else{
          c_prot[i] <- c_prot[e1[2]]
        }
      }else{
        e2 <- as.numeric(unlist(strsplit(edge[e[2],centers_index[c_prot[e[2]]]],',')))
        if(min(dm_ce[e2[1],centers_index])<min(dm_ce[e[1],centers_index])){
          c_prot[i] <- c_prot[e2[1]]
        }else{
          c_prot[i] <- c_prot[e2[2]]
        }
      }
    }
  }
  c_final <- NULL
  if(length(outliers)>0){
    c_prot[outliers] <- k+1
    c_final[c %in% which(c_prot==k+1)] <- k+1
  }
  for(i in 1:(k)){
    c_final[c %in% which(c_prot==i)] <- i
  }
  return(list(centers_index,c_final,tr,di,p,md))
}

library(ggplot2)
library(kernlab)
data("spirals")
x <- as.data.frame(spirals)
cl <- dp_grid(spirals,r=0.2,k=2,kn=11)
ggplot(x,aes(V1,V2,color=cl[[2]]))+geom_point()+scale_color_drsimonj(discrete = FALSE, palette = 'mixed')

library(factoextra)
data("multishapes")
df<-multishapes[,1:2]
df1 <- as.matrix(df)
cl1 <- dp_grid(df1,r=0.2,k=5,kn=16)
ggplot(df,aes(x,y,color=cl1[[2]]))+geom_point()+scale_color_drsimonj(discrete = FALSE, palette = 'mixed')

y <- read_csv("F:/聚类/聚类展示/DPC+K+P/Gaussiandata.csv")
y <- as.matrix(y[,2:3])
cl2 <- dp_grid(y,r=0.1,k=13,kn=22)
ggplot(as.data.frame(y),aes(V1,V2,color=cl2[[2]]))+geom_point()+scale_color_drsimonj(discrete = FALSE, palette = 'mixed')
