ap <- function(x,lamda){
  #similarity:s;responsibility:r;avalability:a
  n <- nrow(x)
  s <- -dist(x)
  med <- median(s)
  s <- as.matrix(s)
  diag(s) <- med
  a <- diag(0,n)
  r <- diag(0,n)
  r_new <- a
  a_new <- r
  ir <- 0
  while(ir<=100){
    for(i in 1:n){
      for(k in 1:n){
        r_new[i,k] <- s[i,k]-max(a[i,-k]+s[i,-k])
        z <- r[,k]
        if(i==k){
          a_new[i,k] <- sum(z[-k][z[-k]>0])
        }else{
          a_new[i,k] <- min(0,r[k,k]+sum(z[-c(i,k)][z[-c(i,k)]>0]))
        }
      }
    }
    r <- lamda*r+(1-lamda)*r_new
    a <- lamda*a+(1-lamda)*a_new
    ir <- ir+1
  }
  c <- apply(a+r,1,which.max)
  return(c)
}

x <- cbind(c(0,1)+rnorm(50,sd=0.1),c(0,2)+rnorm(50,sd=0.1))
plot(x)
c <- ap(x,0.9)
plot(x,col=c)
