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
distmat_fast <- function(mat){
  smat <- apply(mat, 1, crossprod)
  mat1 <- matrix(smat, nrow=nrow(mat), ncol=nrow(mat))
  mat3 <- tcrossprod(mat)
  mat4 <- mat1 + t(mat1) - 2*mat3
  diag(mat4) <- 0
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
#近似计算类之间的最短距离
dist_prot <- function(cq,cl){
  zq <- colMeans(cq)
  zl <- colMeans(cl)
  cq_to_zq <- distvec(zq,cq)
  cl_to_zq <- distvec(zq,cl)
  cl_to_zl <- distvec(zl,cl)
  cq_to_zl <- distvec(zl,cq)
  cq_ind <- which(rank(abs(cq_to_zq-cq_to_zl))<4)
  cl_ind <- which(rank(abs(cl_to_zq-cl_to_zl))<4)
  dm <- distmat_f(cq[cq_ind,],cl[cl_ind,])
  cq_cl_dist <- (mean(apply(dm,1,min))+mean(apply(dm,2,min)))/2
  return(cq_cl_dist)
}
#k个类两两间距离
dist_f <- function(x,dm,c,ce){
  k <- nrow(ce)
  d <- diag(0,k)
  for(i in 1:(k-1)){
    for(j in 2:k){
      cq_to_zq <- dm[c==i,i]
      cq_to_zl <- dm[c==i,j]
      cl_to_zl <- dm[c==j,j]
      cl_to_zq <- dm[c==j,i]
      cq_ind <- which(rank(abs(cq_to_zq-cq_to_zl))<4)
      cl_ind <- which(rank(abs(cl_to_zq-cl_to_zl))<4)
      d1 <- distmat_f(x[which(c==i)[cq_ind],],x[which(c==j)[cl_ind],])
      d[i,j] <- (mean(apply(d1,1,min))+mean(apply(d1,2,min)))/2
      d2 <- distmat_fast(x[which(c==i)[cq_ind],])
      d3 <- distmat_fast(x[which(c==j)[cl_ind],])
      for(t in 1:length(cq_ind)){
        if(max(d2[t,])>max(d1[t,])){
          c[cq_ind[t]] <- j
        }
      }
      for(l in 1:length(cl_ind)){
        if(max(d3[l,])>max(d1[,l])){
          c[cl_ind[l]] <- i
        }
      }
    }
  }
  return(d,c)
}
#k_means算法
k_means <- function(x,k){
  n <- nrow(x)
  ic <- x[sample(n,k),]
  dm <- distmat_f(x,ic)
  c <- apply(dm,1,which.min)
  ce <- as.matrix(aggregate(x,list(c),mean)[,-1])
  while(identical(ic,ce)==F){
    ic <- ce
    dm <- distmat_f(x,ce)
    c <- apply(dm,1,which.min)
    ce <- as.matrix(aggregate(x,list(c),mean)[,-1])
  }
  return(list(c,dm,ce))
}#c是数据点分类情况，dm是点到类中心的距离矩阵，ce是中心
#将不纯的类分为两个子类
split_subcluster <- function(x){
  n <- nrow(x)
  center <- colMeans(x)
  dv <- distmat_f(x,center)
  dmin <- min(dv[dv!=min(dv)])
  dmin_ind <- which(dv==dmin)
  num <- sum(dv<2*dmin)
  if(num>=floor(n/3)){
    split_cluster <- k_means(x,floor(log(n,2)))
    clu <-  split_cluster[[1]]
    d_to_centers <-  split_cluster[[2]]
    centers_new <-  split_cluster[[3]]
  }else{
    clu <- NULL
    d_to_centers <- dv
    centers_new <- center
  }
  return(list(clu,d_to_centers,centers_new))
}
#将数据集分为K个小类，且这些类中的点的真实类别尽量相同
k_means_split <- function(x,k){
  first_k_means <- k_means(x,k)
  c <- first_k_means[[1]]
  dm <- first_k_means[[2]]
  ce <- first_k_means[[3]]
  t <- 0
  c_new <- NULL
  dm_new <- NULL
  centers_new <- NULL
  for(i in 1:k){
    m <- sum(c==i)
    center <- ce[i,]
    dv <- dm[c==i,i]
    dmin <- min(dv[dv!=min(dv)])
    dmin_ind <- which(dv==dmin)
    num <- sum(dv<2*dmin)
    if(num>=floor(m/4)){
      split_cluster <- k_means(x[c==i,],ceiling(log(m,2)))
      c_new[c==i] <- split_cluster[[1]]+t
      centers_new <-  rbind(centers_new,split_cluster[[3]])
      dm_new <- cbind(dm_new,distmat_f(x,split_cluster[[3]]))
      t <- t+ceiling(log(m,2))
    }else{
      c_new[c==i] <- t+1
      t <- t+1
      centers_new <- rbind(centers_new,center)
      dm_new <- cbind(dm_new,dm[,i])
    }
  }
  c_final <- apply(dm_new,1,which.min)
  while(identical(c_final,c_new)==F){
    c_final <- c_new
    dm <- distmat_f(x,centers_new)
    c_new <- apply(dm,1,which.min)
    centers_new <- as.matrix(aggregate(x,list(c_new),mean)[,-1])
  }
  return(list(c_final,dm_new,centers_new))
}
#prim算法求最小瓶颈路
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
#k_means减少时间复杂度的密度峰算法
dp_k <- function(x,r,k){
  n <- nrow(x)
  kc <- k_means_split(x,15)
  c <- kc[[1]]
  dm <- kc[[2]]
  ce <- kc[[3]]
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
      dm_cen[i,j] <- mean(apply(d1,1,min))+mean(apply(d1,2,min))
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
  delta[md<=10]=0
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
    if(sum(dm_ce[i,centers_index]==a[i])>1&a[i]!=0){
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
  for(i in 1:k){
    c_final[c%in%which(c_prot==i)] <- i
  }
  return(list(centers_index,c_final,tr,di,p,md,ce))
}






library(kernlab)
library(ggplot2)
data(spirals)
clr <- dp_k(spirals,r=0.3,k=2)
plot(spirals,col=clr[[2]])
tr <- clr[[3]]
di <- clr[[4]]
p <- clr[[5]]
ce <- clr[[7]]
pat <- matrix(nrow =2*nrow(tr),ncol=3)
for(i in 1:nrow(tr)){
  pat[2*i-1,] <- c(ce[tr[i,1],],z=i)
  pat[2*i,] <- c(ce[tr[i,2],],z=i)
}
pat <- as.data.frame(pat)
names(pat) <- c('x','y','z')
ggplot(pat,aes(x,y,group=z))+geom_path()+geom_point()
di <- round(c(0,di),2)
ggplot(data=as.data.frame(ce[p,]),aes(V1,V2))+geom_text(aes(label=di))
spirals1 <- as.data.frame(spirals)
ggplot(spirals1,aes(V1,V2))+geom_point(color=c)+geom_path(aes(x,y,group=z),data=pat)
ggplot(spirals1,aes(V1,V2))+geom_point(color=clr[[2]])


library(factoextra)
data("multishapes")
df<-multishapes[,1:2]
df <- as.matrix(df)
c1 <- dp_k(df,r=0.9,k=5)
plot(df,col=c1[[2]])
tr <- c1[[3]]
di <- round(c1[[4]],3)
p <- c1[[5]]
ce <- c1[[7]]
pat <- matrix(nrow =2*nrow(tr),ncol=3)
for(i in 1:nrow(tr)){
  pat[2*i-1,] <- c(ce[tr[i,1],],z=i)
  pat[2*i,] <- c(ce[tr[i,2],],z=i)
}
pat <- as.data.frame(pat)
names(pat) <- c('x','y','z')
di <- c(0,di)
de <- as.data.frame(ce[p,])
ggplot(pat,aes(x,y,group=z))+geom_path()+geom_point()
ggplot(data=as.data.frame(ce[p,]),aes(x,y))+geom_text(aes(label=di))
df1 <- as.data.frame(df)
ggplot(df1,aes(x,y))+geom_point(color=c1[[2]])  
ggplot(df1,aes(x,y))+geom_point(color=c1[[2]])+geom_path(aes(x,y,group=z),data=pat)

