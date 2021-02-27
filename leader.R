#leader
distvec <- function(x,t){
  m <- nrow(t)
  a <- NULL
  for(i in 1:m){
    a[i] <- sqrt(sum((x-t[i,])^2))
  }
  return(a)
}
leader <- function(data_set,tao){
  n <- nrow(data_set)
  count <- NULL
  l <- t(as.matrix(data_set[1,]))
  count[1] <- 1 
  f <- list(1)
  j <- 1
  for(i in 2:n){
    d <- distvec(data_set[i,],l)
    candidates <- which(d<tao)
    if(length(candidates)<1){
      l <- rbind(l,data_set[i,])
      count <- c(count,1)
      j <- j+1
      f[[j]] <- candidates
    }else{
      candidates <- candidates[[1]]
      f[[candidates]] <- c(f[[candidates]],i)
      count[candidates] <- count[candidates]+1
    }
  }
  return(list(l,f,count))
}


library(factoextra)
data("multishapes")
df<-multishapes[,1:2]
df <- as.matrix(df)
plot(df)
df_l <- leader(df,0.08)
plot(df_l[[1]])
dim(df_l[[1]])
c1 <- cfsfdp_path(df_l[[1]],5,0.3)
plot(df_l[[1]],col=c1[[2]])

library(kernlab)
data(spirals)
spirals_l <- leader(spirals,0.08)
dim(spirals_l[[1]])
plot(spirals_l[[1]])
c2 <- cfsfdp_path(spirals_l[[1]],2,0.3)
plot(spirals_l[[1]],col=c2[[2]])


iris_l <- leader(iris1,0.2)
dim(iris_l[[1]])
c3 <- cfsfdp_path(iris_l[[1]],3,0.55)
c3

