my_dtw=function(x,y){
  v=vector('list',2)
  n=length(x)
  m=length(y)
  cost=matrix(nrow=n,ncol=m)
  path=array(dim=c(n,m,2))
  cost[1,1]=abs(x[1]-y[1])
  path[1,1,]=c(0,0)
  for(i in 2:n){
    cost[i,1]=cost[i-1,1]+abs(x[i]-y[1])
    path[i,1,]=c(i-1,1)
  }
  for(j in 2:m){
    cost[1,j]=cost[1,j-1]+abs(x[1]-y[j])
    path[1,j,]=c(1,j-1)
  }
  for(i in 2:n){
    for(j in 2:m){
      minmum=cost[i-1,j]
      p_temp=c(i-1,j)
      if(minmum>cost[i-1,j-1]){
        if(cost[i-1,j-1]>cost[i,j-1]){
          minmum=cost[i,j-1]
          p_temp=c(i,j-1)
        }else{
          minmum=cost[i-1,j-1]
          p_temp=c(i-1,j-1)
        }
      }else if(minmum>cost[i,j-1]){
        minmum=cost[i,j-1]
        p_temp=c(i,j-1)
      }
      cost[i,j]=minmum+abs(x[i]-y[j])
      path[i,j,]=p_temp
    }
  }
  v[[1]]=cost
  v[[2]]=path
  return(v)
}

dba=function(c,s){
  lc=length(c)
  ls=ncol(s)
  ns=nrow(s)
  assoctab=vector('list',lc)
  for(i in 1:ns){
    v=my_dtw(c,s[i,])
    k=lc
    h=ls
    while(k>=1&&h>=1){
      assoctab[[k]]=c(assoctab[[k]],s[i,h])
      t=v[[2]][k,h,]
      k=t[1]
      h=t[2]
    }
  }
  for(i in 1:lc){
    c[i]=mean(assoctab[[i]])
  }
  return(c)
}
library(tidyverse)
x=rnorm(10)
y=rnorm(10)
df=tibble(x,y,t=1:10)
gather(df,x,y,key=line,value=value)%>%ggplot(aes(t,value,color=line))+geom_line()
c=rep(0,10)
c=dba(c,rbind(x,y))
df=tibble(x,y,c,t=1:10)
gather(df,x,y,c,key=line,value=value)%>%ggplot(aes(t,value,color=line))+geom_line()


  