#¼ÆËã¾àÀë¾ØÕó
distmat_f <- function(x,y){
  if(!is.matrix(x)){
    x <- t(as.matrix(x))
  }
  if(!is.matrix(y)){
    y <- t(as.matrix(y))
  }
  n <- nrow(x)
  m <- nrow(y)
  xmat <- apply(x, 1, crossprod)
  ymat <- apply(y,1,crossprod)
  mat1 <- matrix(xmat, nrow=n, ncol=m)
  mat2 <- matrix(ymat,nrow =n,ncol = m,byrow = T )
  mat3 <- tcrossprod(x,y)
  mat4 <- mat1 + mat2 - 2*mat3
  mat5 <- sqrt(mat4)
  return(mat5)
}
#k_meansËã·¨
k_means <- function(x,k){
  n <- nrow(x)
  ic <- x[sample(n,k),]
  dm <- distmat_f(x,ic)
  c <- apply(dm,1,which.min)
  ce <- as.matrix(aggregate(x,list(c),mean)[,-1])
  while(identical(ic,ce)==F){
    ic <- ce
    dm <- distmat_f(x,ce)
    c <- apply(dm,1,which.min)
    ce <- as.matrix(aggregate(x,list(c),mean)[,-1])
  }
  return(list(c,dm,ce))
}
