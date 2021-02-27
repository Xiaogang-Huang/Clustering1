library(readxl)
library(data.table)
library(ggplot2)
library(reshape2)
x <- read_xls("F:/FUHE/2008.xls",sheet = 1,na = "NULL")
dt <- as.data.table(x)
setnames(dt,c("日期"),c("riqi"))
dt <- dt[riqi=="2008-06-01 00:00:00"]
df <- data.frame(dt)
df[,8:103] <- apply(df[,8:103],2,as.numeric)
apply(is.na(df[,8:103]),1,sum)
df <- df[complete.cases(df[,8:103]),]
shuju <- df[,8:103]
y <- 1:96
yi <- as.data.frame(cbind(y=y,i=t(shuju[1,])))
msp <- lm(yi[,2] ~ splines::bs(y,df=19,degree = 3))
yi$fitted <- fitted(msp)
names(yi) <- c("x","1","xx")
fit <- melt(yi,id.vars="x")
ggplot(fit,aes(x,value))+geom_line(aes(color=variable))
knot_fit <- attr(splines::bs(y,df=19,degree = 3),"knots")
names(knot_fit) <- NULL
knot_fit <- c(1,knot_fit,96)
coe <- msp$coefficients
names(coe) <- NULL
x <- seq(1,96,length.out=96)
fit1 <- NULL          ###自编b样条拟合曲线
for (t in 1:length(x)) {
  bsp <- NULL
  for (j in 1:16) {
    bsp[j] <- basic_spline(knot_fit,x[t],i = j,p=3)
  }
  fit1 <- c(fit1,coe[1]+coe[2]*x[i]+coe[3]*x[i]^2+coe[4]*x[2]^3+sum(coe[5:20]*bsp))
}
length(fit1)
fit1 <- cbind(x,fit1)
yi <- cbind(yi,fit1)
names(yi)
fit <- melt(yi,id.vars="y")
names(fit)
ggplot(fit,aes(y,value))+geom_line(aes(color=variable))

#meanshift
ms <- meanshift(s=shuju,t=shuju,lamda = 500,tao = 10^-5)
m_t <- ms$m_t
y <- 1:dim(m_t)[1]
julei <- cbind(m_t,y)
ronghe <- melt(julei,id.vars="y")
ggplot(ronghe,aes(Var2,value,group=Var1))+geom_line()
zhongxin <- unique(as.data.frame(m_t))
lei <- 1:dim(zhongxin)[1]
zxyl <- cbind(zhongxin,lei)
final1 <- merge(julei,zxyl,by=colnames(m_t))
fenlei <- cbind(shuju[final1$y,],y=final1$y,lei=final1$lei)
d <- melt(as.data.table(fenlei),id.vars=c("lei","y"))
ggplot(d,aes(variable,value,group=y))+geom_line(aes(color=lei))+scale_colour_gradientn(colours=rainbow(length(lei)))
{ j <- as.data.table(julei)[KW1==zhongxin[265,1]]
  if(dim(j)[1]>1){
    j <- cbind(shuju[j$y,],y=j$y)
    j<- melt(j,id.vars="y")
    ggplot(j,aes(variable,value,group=y))+geom_line()} else {
      j <- cbind(t(shuju[j$y,]),y=j$y)
      colnames(j) <- c("value","y")
      j <- cbind(j,x=1:96)
      j <- as.data.frame(j)
      ggplot(j,aes(x,value))+geom_line()
    }
}

gulidian <- function(x,maxi){
  if(sum(x>maxi)>0) return(1)
  else return(0)
}
a <- apply(shuju,1,gulidian,maxi=50000)
shuju <- shuju[!a,]
dim(shuju)
junzhi <- apply(shuju,1,mean)
ggplot(data = NULL, mapping = aes(x = junzhi)) + geom_histogram(binwidth =200)
a <- biaozhuncha/junzhi #变异系数
complete.cases(a)
a <- a[complete.cases(a)]
ggplot(data = NULL, mapping = aes(x = a)) + geom_histogram(binwidth =1)

###特征重要性排序
sim <- function(a,x,y){
  s <- exp(-a*sqrt(sum((x-y)^2)))
  return(s)
}
exy <- function(a,x,y){
  s <- sim(a,x,y)
  if(s==0 |s==1){
    return(0)
  }else return(-s*log(s)-(1-s)*log(1-s))
}
evec <- function(a,x,m){
  evec <- apply(m,1,exy,a=a,x=x)
  return(evec)
}
Entr <- function(a,m){
  if(dim(m)[2]==1){
    emat <- exy(a,m,m)
  }else{
    emat <- apply(m,1,evec,m=m,a=1)
  }
  e <- sum(emat)
  return(e)
}
fr <- function(a,m){
  e <- Entr(a,m)
  p <- c()
  for(i in 1:dim(m)[2]){
    n <- as.matrix(m[,-i])
    p[i] <- Entr(a,m=n) 
  }
  names(p) <- colnames(m)
  return(rank(p))
} 
m <- iris[,-5]
plot(iris[,2:4])
head(m)
m <- m/(max(m)-min(m))
fr(a,m)

###数据极值标准化
bianhuan <- function(x){
  d <- min(x)
  u <- max(x)
  if(u==d){
    return(x-d)
  }else{
    return((x-d)/(u-d))
  }
}
bzh <- t(apply(shuju,1,bianhuan))
dim(bzh)
y <- 1:96
yi1 <- as.data.frame(cbind(y=y,i=bzh[1,]))
msp <- lm(yi1[,2] ~ splines::bs(y,df=19,degree = 3))
yi1$fitted <- fitted(msp)
names(yi1) <- c("x","1","xx")
fit <- melt(yi1,id.vars="x")
ggplot(fit,aes(x,value))+geom_line(aes(color=variable))

###原始数据聚类并画图
shuju <- bzh
ms <- meanshift(s=shuju,t=shuju,lamda =3,tao = 10^-5)
m_t <- ms$m_t
y <- 1:dim(m_t)[1]
julei <- cbind(m_t,y)
ronghe <- melt(julei,id.vars="y")
zhongxin <- unique(as.data.frame(m_t))
lei <- 1:dim(zhongxin)[1]
zxyl <- cbind(zhongxin,lei)
final1 <- merge(julei,zxyl,by=colnames(m_t))
fenlei <- cbind(shuju[final1$y,],y=final1$y,lei=final1$lei)
d <- melt(as.data.table(fenlei),id.vars=c("lei","y"))
ggplot(d,aes(variable,value,group=y))+geom_line(aes(color=lei))+scale_colour_gradientn(colours=rainbow(length(lei)))
j <- as.data.table(julei)[KW1==zhongxin[1,1]]
if(dim(j)[1]>1){
  j <- as.data.frame(cbind(shuju[j$y,],y=j$y))
  j <- melt(j,id.vars="y")
  ggplot(j,aes(variable,value,group=y))+geom_line()} else {
    j <- cbind(t(shuju[j$y,]),y=j$y)
    colnames(j) <- c("value","y")
    j <- cbind(j,x=1:96)
    j <- as.data.frame(j)
    ggplot(j,aes(x,value))+geom_line()
  }

###高斯核版meanshift
gaosi_ker <- function(beta,x,y){
  return(exp(-beta*sum((x-y)^2)))
}
gaosi_dv <- function(x,t,beta){
  a <- apply(t,1,gaosi_ker,x=x,beta=beta)
  return(a)
}
gaosi_mx <- function(v,s,lamda,beta){
  a <- gaosi_dv(x=v,t=s,beta)
  b <- sum(a<lamda)
  a[a>lamda] <- 0
  if(b >1){
    weight <- a/sum(a)
    mx <- apply(s*a,2,mean)
  } else if(b==0) {
    mx <- v
  } else {mx <- s[a < lamda,]}
  return(mx)
}
gaosi_meanshift <- function(s,t,lamda,tao,beta){
  b <- s
  subject <- 1:dim(s)[1]
  t1 <- apply(t,1,gaosi_mx,s=s,lamda=lamda,beta=beta)
  m_t <- t(t1)
  b <- rbind(s,m_t)
  i <- 2
  while(sum(apply(m_t-t,1,mo)>tao)>0){
    t <- m_t
    t1 <- apply(t,1,gaosi_mx,s=s,lamda=lamda,beta=beta)
    m_t <- t(t1)
    b <- rbind(b,m_t)
    i=i+1
  }
  subject <- rep(subject,i)
  b <- cbind(subject,b)
  l <- list(b=b,m_t=m_t,iteration=i-1)
  return(l)
}


###数据降维
y <- 1:96
bspline <- function(x,y,weidu,ci){
  bsp  <- lm(x ~ splines::bs(y,df=weidu,degree = ci))
  return(bsp$coefficients)
}
jiangwei <- function(sj,y,weidu,ci){
  return(t(apply(sj,1,bspline,y=y,weidu=weidu,ci=ci)))
}
jw_shuju <- jiangwei(sj=shuju,y=y,weidu=19,ci=3)
jw_shuju <- t(apply(jw_shuju,1,bianhuan))

###高斯meanshift聚类（数据降维但未标准化）
ms <- gaosi_meanshift(s=jw_shuju,t=jw_shuju,lamda =2,tao = 10^-3,beta =0.05)
m_t <- ms$m_t
m_t <- round(m_t,5)
y <- 1:dim(m_t)[1]
julei <- cbind(m_t,y)
ronghe <- melt(julei,id.vars="y")
zhongxin <- unique(as.data.frame(m_t))
lei <- 1:dim(zhongxin)[1]
zxyl <- cbind(zhongxin,lei)
final1 <- merge(julei,zxyl,by=colnames(m_t))
fenlei <- cbind(jw_shuju[final1$y,],y=final1$y,lei=final1$lei)
d <- melt(as.data.table(fenlei),id.vars=c("lei","y"))
ggplot(d,aes(variable,value,group=y))+geom_line(aes(color=lei))+scale_colour_gradientn(colours=rainbow(length(lei)))
j <- as.data.table(julei)[x1==zhongxin[1,1]]
if(dim(j)[1]>1){
  j <- cbind(shuju[j$x21,],y=j$x21)
  i <- j*as.numeric(1-apply(j,1,gulidian,maxi=1000))
  i <- melt(j,id.vars="y")
  ggplot(i,aes(variable,value,group=y))+geom_line()} else {
    j <- cbind(t(shuju[j$y,]),y=j$y)
    colnames(j) <- c("value","y")
    j <- cbind(j,x=1:96)
    j <- as.data.frame(j)
    ggplot(j,aes(x,value))+geom_line()
  }






       