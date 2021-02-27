#julijuzhen
euclid <- function(x,y){
  a <- sqrt(sum((x-y)^2))
  return(a)
}
distvec <- function(x,t){
  a <- apply(t,1,euclid,x=x)
  return(a)
}
distmat <- function(s,t){
  a <- apply(s,1,distvec,t=t)
  return(a)
}
distmat_fast <- function(mat){
  smat <- apply(mat, 1, crossprod)
  mat1 <- matrix(smat, nrow=2000, ncol=2000)
  mat3 <- tcrossprod(mat)
  mat4 <- mat1 + t(mat1) - 2*mat3
  diag(mat4) <- 0
  mat5 <- sqrt(mat4)
}
