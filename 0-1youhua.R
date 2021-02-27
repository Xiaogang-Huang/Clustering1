library(lpSolve)
library(Rdonlp2)
miu1 <- function(Q,r,a,b){
  lamda <-max(eigen(Q)$values) 
  Q1 <- Q-diag(lamda,dim(Q)[1])
  const.rhs <- c(b,c(rep(1,dim(Q)[1]),rep(0,dim(Q)[1])))
  f.l <- lp(objective.in=r+lamda+apply(Q1*(Q1<0),1,sum),
          const.mat=rbind(a,diag(1,dim(Q)[1]),diag(-1,dim(Q)[1])),
          const.rhs=const.rhs,
          const.dir=rep("<=",length(const.rhs)),
          direction="min")$objval
  if(is.null(a)) miu1 <- lamda else{
    miu1 <- lamda+2*(sum(Q[Q>0])+sum(r>0) -f.l)*max(abs(a))
  }
  return(miu1)
}
jie <- function(){
  fn <- function(x){
    t(x)%*%Q%*%x+sum(r*x)+miu*sum(x*(1-x))
  } 
  fn1 <- function(x){
    x%*%Q%*%x+sum(r*x)
  } 
  weidu <- length(r)
  p <- matrix(runif(weidu*10,-1,1),nrow = 10,ncol = weidu)
  candidate <- matrix(nrow = 10,ncol=weidu+1)
  colnames(candidate) <- c(paste("x",1:weidu,sep=""),"fmin")
  m <- c()
  for(i in 1:10){
  fin <- donlp2(p[i,],fn, par.u=u, par.l=l,a,lin.l=ll,lin.u=lu)
  candidate[i,] <- append(fin$par,fn1(fin$par))
  m[i] <- fin$fx
  }
  f_u <- sum(Q[Q>0])+sum(r>0)
  ifinteger <- function(x){
    all(x %% 1==0|x==0)
  }
  if(min(m)>f_u){
    lt <- list(fn_min=min(m),f_u=f_u,"原问题无解")
    print(lt)
  } else{
    candidate <- round(candidate,2)
    zsd <- apply(candidate[,1:weidu],1,ifinteger)
    candidate <- candidate[zsd,]
    candidate <- candidate[candidate[,weidu+1]==min(candidate[,weidu+1]),]
    candidate <- unique(candidate)
    print(candidate)
  }
}
#例子实现
Q <- matrix(c(0,-1,-1,0),2,byrow = T)
r <- c(0,-1)
l <- c(0,0)
u <- c(1,1)
b <- c(1,0)
ll <- c(-Inf,-Inf)
lu <- c(1,0)
a <- matrix(c(1,1,-1,1),2,byrow = T)
miu <- miu1(Q,r,a,b)
jieguo <- jie()
#无解的例子
a <- matrix(c(4,3,-4,3),2,byrow = T)
b <- c(3,-1)
lu <- b
miu <- miu1(Q,r,a,b)
jie()
#四维例子
Q <- matrix(c(0,26,44,-73,26,0,-45,11,44,-45,0,84,-73,11,84,0),nrow = 4,byrow = T)
r <- c(-119,27,-187,-2)
b <- c()
a <- c()
miu <- miu1(Q,r,a,b)
miu <- miu+1
ll <- rep(-Inf,4)
lu <- rep(Inf,4)
u <- c(1,1,1,1)
l <- c(0,0,0,0)
jieguo <- jie()
