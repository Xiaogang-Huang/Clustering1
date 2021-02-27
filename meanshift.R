#T和S的距离矩阵
euclid <- function(x,y){
  a <- sqrt(sum((x-y)^2))
  return(a)
}
distvec <- function(x,t){
  a <- apply(t,1,euclid,x=x)
  return(a)
}
#flat kernel meanshif
mx <- function(v,s,lamda){
  a <- distvec(x=v,t=s)
  b <- sum(a<lamda)
  if(b >1){
     mx <- apply(s[a < lamda,],2,mean)
  } else if(b==0) {
    mx <- v
    } else {mx <- s[a < lamda,]}
  return(mx)
}
mo <- function(x){return(sqrt(sum(x^2)))}

meanshift <- function(s,t,lamda,tao){
  b <- s
  subject <- 1:dim(s)[1]
  t1 <- apply(t,1,mx,s=s,lamda=lamda)
  m_t <- t(t1)
  b <- rbind(s,m_t)
  i <- 2
  while(sum(apply(m_t-t,1,mo)>tao)>0){
    t <- m_t
    t1 <- apply(t,1,mx,s=s,lamda=lamda)
    m_t <- t(t1)
    b <- rbind(b,m_t)
    i=i+1
  }
  subject <- rep(subject,i)
  b <- cbind(subject,b)
  l <- list(b=b,m_t=m_t,iteration=i-1)
  return(l)
}

blurmeanshift <- function(s,t,lamda,tao){
  b <- s
  subject <- 1:dim(s)[1]
  t1 <- apply(t,1,mx,s=s,lamda=lamda)
  m_t <- t(t1)
  s <- m_t
  b <- rbind(s,m_t)
  i <- 2
  while(sum(apply(m_t-t,1,mo)>tao)>0){
    t <- m_t
    t1 <- apply(t,1,mx,s=s,lamda=lamda)
    m_t <- t(t1)
    s <- m_t
    b <- rbind(b,m_t)
    i=i+1
  }
  subject <- rep(subject,i)
  b <- cbind(subject,b)
  l <- list(b=b,m_t=m_t,iteration=i-1)
  return(l)
}

#案例
library(MASS) 
library(ggplot2)
Sigma <- matrix(c(1,0,0,1),2,2) 
x=mvrnorm(n=50, rep(0, 2), Sigma)
y=mvrnorm(n=50, c(4,4), Sigma)
a <- rbind(x,y)
{ms <- meanshift(s=a,t=a,lamda =1.3,tao = 10^-4)
ms$m_t
ms$iteration
b <- as.data.frame(ms$b)
head(b)
colnames(b) <- c("subject","x","y")
ggplot(b,aes(x,y,group=subject))+geom_path()+geom_point()}
x <- runif(100,0,1)
y <- runif(100,0,1)
a <- cbind(x,y)
un <- as.data.frame(a)
head(un)
ggplot(data = un,aes(x,y))+geom_point()
ms <- meanshift(s=a,t=a,lamda = 0.5,tao = 10^-4)
b <- as.data.frame(ms$b)
head(b)
colnames(b) <- c("subject","x","y")
ggplot(b,aes(x,y,group=subject))+geom_line()+geom_point()



#月牙型数据
means <- c(0,0)
vars <- matrix(c(1,0,0,1),nrow = 2)
a_normal <- mvrnorm(10000,means,vars)
up_data <- a_normal[0.6<rowSums(a_normal^2)&rowSums(a_normal^2)<1&a_normal[,2]>0,]
a_normal_1 <- a_normal[,1]-1
a_normal_1 <- cbind(a_normal_1,a_normal[,2]-0.5)
low_data <- a_normal[0.8<rowSums(a_normal_1^2)&rowSums(a_normal_1^2)<1.2&a_normal[,2]<0.5,]
up_data_1 <- up_data[sample(1:679,150),]
low_data_1 <- low_data[sample(1:519,150),]
moon_data <- rbind(up_data_1,low_data_1)
y <- meanshift(s=moon_data,t=moon_data,lamda =0.4 ,tao = 10^-4)
b <- as.data.frame(y$b)
colnames(b) <- c("subject","x","y")
ggplot(b,aes(x,y,group=subject))+geom_path()+geom_point()
m_t <- y$m_t
unique(m_t)


library(factoextra)
data("multishapes")
df<-multishapes[,1:2]
plot(df)
df <- as.matrix(df)
df1 <- meanshift(df,df,0.3,10^-2)
b <- as.data.frame(df1$b)
colnames(b) <- c("subject","x","y")
ggplot(b,aes(x,y,group=subject))+geom_path()+geom_point()
m_t <- df1$m_t
unique(m_t)
