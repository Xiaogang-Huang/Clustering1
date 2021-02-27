#判断一个向量是否与矩阵的某行（列）相等
vec_in_matrix <- function(x,m,marg){
  a <- apply(m, marg, identical,y=x)
  if(sum(a)==0){
    return(FALSE)
  }else{
    return(TRUE)
  }
}