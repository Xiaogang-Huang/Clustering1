#生成月牙型数据
library(MASS)
#定义点到类的距离
CBF_kernel <- function(x,y,r){
  return(exp(-r*sum((x-y)^2)))
}
kernel_distance <- function(x,c,r){
  n <- nrow(c)
  g <- diag(1,nrow = n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      g[i,j] <- CBF_kernel(c[i,],c[j,],r)
    }
  }
  g <- g+t(g)-diag(1,nrow = n)
  dist <- 1+sum(g)/n^2-2*sum(apply(c,MARGIN = 1,FUN = CBF_kernel,x=x,r=r))/n
  return(dist)
}
means <- c(0,0)
vars <- matrix(c(1,0,0,1),nrow = 2)
a_normal <- mvrnorm(10000,means,vars)
plot(a_normal)
up_data <- a_normal[0.6<rowSums(a_normal^2)&rowSums(a_normal^2)<1&a_normal[,2]>0,]
plot(up_data)
a_normal_1 <- a_normal[,1]-1
a_normal_1 <- cbind(a_normal_1,a_normal[,2]-0.5)
low_data <- a_normal[0.8<rowSums(a_normal_1^2)&rowSums(a_normal_1^2)<1.2&a_normal[,2]<0.5,]
plot(low_data)
dim(up_data)
dim(low_data)
up_data_1 <- up_data[sample(1:679,200),]
low_data_1 <- low_data[sample(1:519,200),]
point_analyze <- low_data_1[which.max(low_data_1[,2]),]
moon_data <- rbind(up_data_1,low_data_1)
co <- rep(1,nrow(moon_data))
co[which.max(low_data_1[,2])+dim(up_data_1)[1]] <- 2
plot(moon_data,col=co)
#r=1时
point_to_up <- kernel_distance(point_analyze,up_data_1,r=1)
point_to_low <- kernel_distance(point_analyze,low_data_1[-which.max(low_data_1[,2]),],r=1)
r_distance <- matrix(nrow = 10,ncol = 2)
for(i in 1:10){
  r_distance[i,] <- c(kernel_distance(point_analyze,up_data_1,r=i),
                      kernel_distance(point_analyze,low_data_1[-which.max(low_data_1[,2]),],r=i))
}
colnames(r_distance) <- c('到上月牙的距离','到下月牙的距离')
rownames(r_distance) <- paste('r=',1:10,sep = '')
r_distance








