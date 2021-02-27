library(kernlab)
library(ggplot2)
#引入数据
data(spirals)
#将数据设置为数据框格式
df<-as.data.frame(spirals)
#重新命名
names(df)<-c("x1","x2")
#查看原始数据
ggplot(df,aes(x=x1,y=x2))+geom_point()

#构造相似矩阵
rbf <- function(x,y,r){
  a=exp(-r*sum((x-y)^2))
  return(a)
}
rbfvec <- function(x,t,r){
  a <- apply(t,1,rbf,x=x,r=r)
  return(a)
}
rbfmat <- function(s,t,r){
  a <- apply(s,1,rbfvec,t=t,r=r)
  return(a)
}
epsilon_matrix <- function(data_set,epsilon){
  dm <- as.matrix(dist(data_set))
  dm[dm>epsilon] <- 0
  dm[dm<=epsilon] <- epsilon
  return(dm)
}
kn_matrix <- function(data_set,k,r){
  n <- nrow(data_set)
  dm <- as.matrix(dist(data_set))
  rbf_m <- rbfmat(data_set,data_set,r=r)
  d <- apply(dm,1,rank,ties.method='min')
  d1 <- d
  d1[d<=t(d)] <- d[d<=t(d)]
  d1[d>t(d)] <- t(d)[d>t(d)]
  d <- d1
  rbf_m[d>k+1] <- 0
  return(rbf_m)
}
mutual_kn_matrix <- function(data_set,k){
  dm <- as.matrix(dist(data_set))
  d <- apply(dm,1,rank,ties.method='min')
  rbf_m <- rbfmat(data_set,data_set,r=r)
  d1 <- d
  d1[d>=t(d)] <- t(d)[d>=t(d)]
  d1[d<t(d)] <- d[d<t(d)]
  d <- d1
  rbf_m[d>k+1] <- 0
  return(rbf_m)
}
spec_clu <- function(data_set,k,lei,r=1){
  w <- kn_matrix(data_set,k,r)
  d <- diag(rowSums(w))
  l <- d-w
  z <- eigen(l)
  eva <- z$values
  eve <- z$vectors
  v <- eve[,rank(eva)<=lei]
  return(kmeans(v,centers=lei))
}
t <- spec_clu(spirals,k=3,lei=2)
plot(spirals,col=t$cluster)
lei <- 2
r <- 2
for(i in 1:4){
  k <- 2+i
  w <- kn_matrix(spirals,k,r)
  d <- diag(rowSums(w))
  l <- d-w
  z <- eigen(l)
  eva <- z$values
  eve <- z$vectors
  v <- eve[,rank(eva)<=lei]
  plot(v,col=t$cluster)
}
library(MASS)
means <- c(0,0)
vars <- matrix(c(1,0,0,1),nrow = 2)
a_normal <- mvrnorm(10000,means,vars)
up_data <- a_normal[0.6<rowSums(a_normal^2)&rowSums(a_normal^2)<1&a_normal[,2]>0,]
a_normal_1 <- a_normal[,1]-1
a_normal_1 <- cbind(a_normal_1,a_normal[,2]-0.5)
low_data <- a_normal[0.8<rowSums(a_normal_1^2)&rowSums(a_normal_1^2)<1.2&a_normal[,2]<0.5,]
up_data_1 <- up_data[sample(1:nrow(up_data),200),]
low_data_1 <- low_data[sample(1:nrow(low_data),200),]
moon_data <- rbind(up_data_1,low_data_1)
plot(moon_data)
t <- spec_clu(moon_data,k=5,lei=2)
plot(moon_data,col=t$cluster)
for(i in 1:4){
  k <- 50+i
  w <- kn_matrix(moon_data,k,r=0.2)
  d <- diag(rowSums(w))
  l <- d-w
  z <- eigen(l)
  eva <- z$values
  eve <- z$vectors
  v <- eve[,rank(eva)<=lei]
  plot(v,col=t$cluster)
}
