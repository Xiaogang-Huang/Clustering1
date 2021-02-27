library(rgl)
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
library(factoextra)
data("multishapes")
df<-multishapes[,1:2]
dm <- rbfmat(df,df,r=0.001)
midu=colSums(dm)
plot3d(df[,1],df[,2],midu)
library(kernlab)
data(spirals)
dm <- rbfmat(spirals,spirals,r=1)
midu=colSums(dm)
plot3d(spirals[,1],spirals[,2],midu)
