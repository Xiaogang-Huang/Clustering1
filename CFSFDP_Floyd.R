#计算任意两点路径距离floyd算法
floyd_path_distance <- function(dm){
  n <- nrow(dm)
  for(k in 1:n){
    for(i in 1:n){
      for(j in 1:n){
        if(max(dm[i,k],dm[k,j])<dm[i,j]){
          dm[i,j] <- max(dm[i,k],dm[k,j])
        }
      }
    }
  }
  return(dm)
}
#计算密度、距离和匹配向量
delta <- function(data_set,r){
  dm <- as.matrix(dist(data_set,diag = TRUE))
  midu=colSums(dm<r)
  dm <- floyd_path_distance(dm)
  delta=NULL
  match_vector <- rep(0,nrow(data_set))
  for (i in 1:nrow(data_set)) {
    d <- dm[i,]
    if(midu[i]==max(midu)){
      delta[i]=max(d)
      match_vector[i] <- which.max(d)
    }
    else{
      delta[i]=min(d[which(midu>midu[i])]) 
      d[which(midu<=midu[i])]=1/0
      match_vector[i] <- which.min(d)
    }
  }
  return(list(delta,midu,match_vector,floyd_matrix=dm))
}
#距离计算函数
euclid <- function(x,y){
  a <- sqrt(sum((x-y)^2))
  return(a)
}
distvec <- function(x,t){
  a <- apply(t,1,euclid,x=x)
  return(a)
}
distmat <- function(s,t){
  a <- apply(s,1,distvec,t=t)
  return(a)
}
#寻找类中心
find_center <- function(data_set,midu,delta,k,theta){
  a=delta
  a[midu<=8]=0
  centers_index=NULL
  centers_index <- which(a==sort(a,decreasing = T)[1])[1]
  i <- 2
  while(length(centers_index)<k){
    centers_candidate <- which(a==sort(a,decreasing = T)[i])
    if(which.max(a) %in% centers_candidate){centers_candidate <- centers_candidate[-1]}
    for(j in 1:length(centers_candidate)){
      if(length(centers_index)==1){t=t(as.matrix(data_set[centers_index,]))}else{
        t <- data_set[centers_index,]
      }
      cen_dist <- distvec(data_set[centers_candidate[j],],t)
      if(any(cen_dist<theta)) {}else{
        centers_index <- c(centers_index ,centers_candidate[j])
      }
    }
    i <- i+length(centers_candidate)
  }
  if(length(centers_index)>k){centers_index <- centers_index[1:k]}
  return(centers_index)
}
#数据点分类
assigned <- function(data_set,midu,centers_index,dm,outliers=T){
  assignment <- rep(0,nrow(data_set))
  assignment[centers_index] <- 1:length(centers_index)
  m <- midu
  if(outliers==T){
    outliers_index <- which(m<=3)
  }else{outliers_index <- NULL}
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



#多图形数据集聚类效果
library(factoextra)
data("multishapes")
df<-multishapes[,1:2]
set.seed(1)
df <- as.matrix(df[sample(1100,500),])
mdelta <- delta(df,r=0.2)
md <- mdelta[[2]]
del <- mdelta[[1]]
dm <- mdelta$floyd_matrix
plot(md,del)
library(rgl)
plot3d(df[,1],df[,2],md)
centers_index <- find_center(df,md,del,5,theta=0.2)
col=rep(1,500)
col[centers_index] <- 2
plot(df,col=col)
fclust <- assigned(df,midu=md,centers_index = centers_index,dm,outliers = F)
plot(df,col=fclust)

fclust1=apply(dm[,centers_index],1,which.min)
plot(df,col=fclust1)

heatmap(dm)

#spirals数据集
library(kernlab)
data(spirals)
mdelta <- delta(spirals,r=0.2)
md <- mdelta[[2]]
del <- mdelta[[1]]
dm <- mdelta$floyd_matrix
plot(md,del)
centers_index <- find_center(spirals,md,del,2,theta=0.2)
col=rep(1,300)
col[centers_index] <- 2
plot(spirals,col=col)
fclust <- assigned(spirals,midu=md,centers_index = centers_index,dm,outliers = F)
plot(spirals,col=fclust)




#月牙型数据聚类效果
##生成数据集
library(MASS)
means <- c(0,0)
vars <- matrix(c(1,0,0,1),nrow = 2)
a_normal <- mvrnorm(10000,means,vars)
up_data <- a_normal[0.6<rowSums(a_normal^2)&rowSums(a_normal^2)<1&a_normal[,2]>0,]
a_normal_1 <- a_normal[,1]-1
a_normal_1 <- cbind(a_normal_1,a_normal[,2]-0.5)
low_data <- a_normal[0.8<rowSums(a_normal_1^2)&rowSums(a_normal_1^2)<1.2&a_normal[,2]<0.5,]
up_data_1 <- up_data[sample(1:nrow(up_data),50),]
low_data_1 <- low_data[sample(1:nrow(low_data),50),]
moon_data <- rbind(up_data_1,low_data_1)
plot(moon_data)
##聚类
mdelta1 <- delta(moon_data,r=0.2)
md1 <- mdelta1[[2]]
del1 <- mdelta1[[1]]
dm1 <- mdelta1$floyd_matrix
plot(md1,del1)
library(rgl)
plot3d(moon_data[,1],moon_data[,2],md1)
centers_index <- find_center(moon_data,md1,del1,3,theta=.1)
col=rep(1,100)
col[centers_index] <- 2
plot(moon_data,col=col)
fclust <- assigned(moon_data,midu=md1,centers_index = centers_index,dm1,outliers = F)
plot(moon_data,col=fclust)
which(fclust==2&moon_data[,1]>1.5)
dm1[which(fclust==2&moon_data[,1]>1.5),]
fclust1=apply(dm1[,centers_index],1,which.min)
plot(moon_data,col=fclust1)

data(iris)
d <- iris[,3:4]
mdelta2 <- delta(d,r=0.3)
md2 <- mdelta2[[2]]
del2 <- mdelta2[[1]]
dm2 <- mdelta2$floyd_matrix
library(rgl)
plot3d(d[,1],d[,2],md2)
centers_index <- find_center(d,md2,del2,3,theta=.3)
fclust21 <- assigned(d,midu=md2,centers_index = centers_index,dm2,outliers = F)
fclust21
fclust22 <- apply(dm2[,centers_index],1,which.min)
fclust22
table(fclust21)
