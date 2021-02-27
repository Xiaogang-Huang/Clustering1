###basic spline function
basic <- function(knots,x,i){
  if(knots[i]<=x&x<knots[i+1]){
    return(1)
  }else {
    return(0)
  }
}
basic_spline <- function(knots,x,i,p){
  if(p==0){
    return(basic(knots,x,i))
  }else{
    if(knots[i+p]-knots[i]==0){a <- 0}else{
      a <- (x-knots[i])*basic_spline(knots,x,i,p-1)/(knots[i+p]-knots[i])
    }
    if(knots[i+p+1]-knots[i+1]==0){b <- 0}else{
      b <- (knots[i+p+1]-x)*basic_spline(knots,x,i+1,p-1)/(knots[i+p+1]-knots[i+1])
    }
    return(a+b)
  }
}
t <- c(0, 0, 0, 0.3, 0.5, 0.5, 0.6, 1, 1, 1 )
basic_spline(t,x=0.8,i=7,p=2)
x <- seq(0,1,length.out=100)
y <- matrix(nrow = 100,ncol = 7)
for (j in 1:7) {
  for (i in 1:100) {
    y[i,j] <- basic_spline(t,x[i],j,2)
  }
}
sx <- as.data.frame(cbind(t(y),i=1:7))
library(reshape2)
library(ggplot2)
colnames(sx) <- c(1:100,"i")
sx_melt <- melt(sx,id.vars = "i")
ggplot(sx_melt,aes(variable,value,group=i))+geom_line(aes(color=i))

##一次样条的具体函数
one_degree_bspline <- function(knots,x,i){
  bi_zero <- basic(knots,x,i)
  biplus_zero <- basic(knots,x,i+1)
  a <- bi_zero/(knots[i+1]-knots[i])-biplus_zero/(knots[i+2]-knots[i+1])
  b <- -bi_zero*knots[i]/(knots[i+1]-knots[i])+biplus_zero*knots[i+2]/(knots[i+2]-knots[i+1])
  xishu <- c(a,b)
  return(xishu)
}
t <- c( 0, 0.25, 0.5, 0.75, 1 )
one_degree_bspline(t,x=0.25,i=2)

##三次样条的二阶导的函数表达式
sec_order_derivative <- function(knots,x,i){
  a1 <- knots[i+3]-knots[i]
  a2 <- knots[i+2]-knots[i]
  a3 <- knots[i+3]-knots[i+1]
  a4 <- knots[i+4]-knots[i+1]
  a5 <- knots[i+4]-knots[i+2]     
  b1 <- one_degree_bspline(knots,x,i)
  b2 <- one_degree_bspline(knots,x,i+1)
  b3 <- one_degree_bspline(knots,x,i+2)
  return(6*(b1/(a1*a2)-b2/(a3*a1)-b2/(a4*a3)+b3/(a4*a5)))
}
sec_order_derivative(t,x=0.875,i=1)



