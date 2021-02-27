library(MASS)
sigma1 <- matrix(c(1,0,0,1),2)
x <- mvrnorm(n=100,rep(0,2),sigma1)
sigma2 <- matrix(c(1,0,0,1),2)
y <- mvrnorm(100,c(5,4),sigma2)
data_set <- rbind(x,y)
plot(data_set)


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
k <- 2
r <- 1
kkm <- function(dataset,k,r){
  n <- nrow(dataset)
  initcenter <- dataset[sample(n,k),]
  initcluster <- NULL
  for (i in 1:n) {
    initcluster[i] <- which.max(rbfvec(dataset[i,],initcenter,r))
  }
  fcluster <- initcluster
  xcluster <- rep(0,n)
  mat <- rbfmat(dataset,dataset,r)
  while(!all(fcluster==xcluster)){
    fcluster <- xcluster
    for (i in 1:n) {
      a <- NULL
      x <- dataset[i,]
      for (j in 1:k) {
        c <- dataset[initcluster==j,]
        m <- nrow(c)
        f <- -2*sum(mat[i,initcluster==j])/m
        g <- sum(mat[initcluster==j,initcluster==j])/m^2
        a[j] <- 1+g+f
      }  
      xcluster[i] <- which.min(a)
    }
  }
  return(xcluster)
}
clu <- kkm(data_set,2,r=1)
plot(data_set,col=clu)  #仍然是局部最优
library(kernlab)
a <- kkmeans(data_set,kernel="rbfdot",kpar=list(sigma=0.01),centers=2)
plot(data_set,col=a)


#月牙型数据
#kernlab的kkmeans函数聚类效果
a <- kkmeans(moon_data,kernel="rbfdot",kpar=list(sigma=1),centers=2)
plot(moon_data,col=a)
#kmeans
b <- kmeans(moon_data,centers = moon_data[c(200,201),])
plot(moon_data,col=b$cluster)
#kkm
a <- kkm(moon_data,k=2,r = 0.1)
plot(moon_data,col=a)

#错误分类的分析
sum(a[1:200]==1)
sum(a[201:400]==1)
lei2_1 <- moon_data[a[1:200]==1,]#类2的点被分到类1
lei1_2 <- moon_data[which(a[201:400]==2)+200,]#类1的点被分到类2
lei2 <- moon_data[a==2,]##图中红色的点
lei1 <- moon_data[a==1,]##图中黑色的点
#计算点到类距离的函数
pnt_clu <- function(x,c,r){
  mat <-rbfmat(c,c,r)
  dv <- rbfvec(x,t=c,r)
  m <- nrow(c)
  f <- -2*sum(dv)/m
  g <- sum(mat)/m^2
  return(1+g+f)
}


m <- which.max(lei2_1[,2])
points(lei2_1[m,][1],lei2_1[m,][2],col=3)
rv=seq(0.1,20,length.out = 100)
d <- matrix(nrow = 100,ncol = 2)
for (i in 1:100) {
  d[i,2] <- pnt_clu(lei2_1[19,],lei2,r=rv[i])
  d[i,1] <- pnt_clu(lei2_1[19,],lei1,r=rv[i])
}
d
apply(d,1,diff)


n <- which.min(lei1_2[,2])
points(lei1_2[n,][1],lei1_2[n,][2],col=4)
d1 <- matrix(nrow = 20,ncol = 2)
for (i in 1:20) {
  d1[i,2] <- pnt_clu(lei1_2[n,],lei2,r=rv[i+80])
  d1[i,1] <- pnt_clu(lei1_2[n,],lei1,r=rv[i+80])
}
d1
apply(d1,1,diff)


library(factoextra)
data("multishapes")
df<-multishapes[,1:2]
plot(df)
df <- as.matrix(df)
a <- kkmeans(df,kernel="rbfdot",kpar=list(sigma=1),centers=5)
plot(df,col=a)


clu <- kkm(df,5,r=2)
plot(df,col=clu)



#用核函数计算密度的函数
midu <- function(data_set,r){
  dm <- rbfmat(data_set,data_set,r)
  diag(dm) <- 0
  midu=colSums(dm)
  return(midu)
}
kkm_1 <- function(dataset,k,mmidu){
  n <- nrow(dataset)
  initcenter <- dataset[sample(n,k),]
  initcluster <- NULL
  for (i in 1:n) {
    initcluster[i] <- which.max(rbfvec(dataset[i,],initcenter,r))
  }
  fcluster <- initcluster
  xcluster <- rep(0,n)
  mat <- diag(0,n)
  r <- 10/mmidu
  for (i in 1:(n-1)) {
    for(j in (i+1):n){
      mat[i,j] <- rbf(dataset[i,],dataset[j,],r=(r[i]+r[j])/2)
    }
  }
  mat <- mat+t(mat)+diag(1,n)
  while(!all(fcluster==xcluster)){
    fcluster <- xcluster
    for (i in 1:n) {
      a <- NULL
      x <- dataset[i,]
      for (j in 1:k) {
        c <- dataset[initcluster==j,]
        m <- nrow(c)
        f <- -2*sum(mat[i,initcluster==j])/m
        g <- sum(mat[initcluster==j,initcluster==j])/m^2
        a[j] <- 1+g+f
      }  
      xcluster[i] <- which.min(a)
    }
  }
  return(xcluster)
}
r_1 <- seq(0.01,10,length.out = 20)
for(i in 1:20){
  mmidu=midu(moon_data,r=r_1[i])
  a <- kkm_1(moon_data,k=2,mmidu = mmidu)
  plot(moon_data,col=a)
}



a <- kkm(moon_data,k=2,r = 0.1)
plot(moon_data,col=a)



mmidu<- midu(moon_data,10)
cl <- kkm_2(moon_data,2,1,mmidu)
plot(moon_data,col=cl)


