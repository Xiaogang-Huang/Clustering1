grid_size=function(x){
  n=nrow(x)
  s1=(apply(x,2,function(y) diff(range(y)))/n)^(1/n)
  grid_length=1
  for(i in seq_along(s1)){
    grid_length=grid_length*s1[i] 
  }
  return(s2)
}