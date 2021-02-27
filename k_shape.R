shape_based_distance=function(x,y){
  len=length(x)
  len_y=length(y)
  t=log2(2*len-1)
  if(ceiling(t)!=t){
    len_2=2^ceiling(t)
    x[(len+1):len_2]=0
    y[(len+1):len_2]=0
  }
  ncc=Re(fft(fft(x)*Conj(fft(y)),inverse = T)/len_2)
  den=drop(sqrt(x%*%x)*sqrt(y%*%y))
  if(den==0){
    den=Inf
  }
  ncc=ncc/den
  ncc=ncc[c((len_2+2-len):len_2,1:len)]
  index=which.max(ncc)
  distance=1-ncc[index]
  shift=index-len
  if(shift>=0){
    y=c(rep(0,shift),y[1:(len-shift)])
  }else{
    y=c(y[(1-shift):len],rep(0,-shift))
  }
  return(list(dist=distance,y=y))
}  

shapeextraction=function(x,c){
  n=nrow(x)
  m=ncol(x)
  x_prime=matrix(nrow=n,ncol=m)
  for(i in 1:n){
    if(all(c==0)){
      c=x[i,]
      x_prime[i,]=c
    }else{
      sbd=shape_based_distance(c,x[i,])
      x_prime[i,]=sbd$y
    }
  }
  x_prime=t(apply(x_prime,1,scale))
  x_prime=x_prime[complete.cases(x_prime),]
  if(!is.matrix(x_prime)){
    x_prime=t(as.matrix(x_prime))
  }
  s=t(x_prime)%*%x_prime
  q=diag(1,m)-matrix(1,m,m)/m
  mt=t(q)%*%s%*%q
  eig=eigen(mt)
  centroid=scale(eig$vectors[,which.max(eig$value)])
  finddistance1 = sqrt(sum(x_prime[1,] - centroid)^2)
  finddistance2 = sqrt(sum(x_prime[1,] + centroid)^2)
  if(finddistance1 >= finddistance2){
    centroid =-centroid
  }
  return(centroid)
}

kshape=function(x,k){
  n=nrow(x)
  m=ncol(x)
  iter=0
  c=x[sample(n,k),]
  initcenters=c
  idx=NULL
  for(i in 1:n){
    minidist=1/0
    for(j in 1:k){
      sbd=shape_based_distance(c[j,],x[i,])
      if(sbd$dist<minidist){
        minidist=sbd$dist
        idx[i]=j
      }
    }
  }
  idx_prime=rep(0,n)
  while(!all(idx_prime==idx)&iter<100){
    idx_prime=idx
    for(j in 1:k){
      x_prime=NULL
      for(i in 1:n){
        if(idx[i]==j){
          x_prime=rbind(x_prime,x[i,])
        }
      }
      if(!is.null(x_prime)){
        c[j,]=shapeextraction(x_prime,c[j,])
      }
    }
    for(i in 1:n){
      minidist=1/0
      for(j in 1:k){
        sbd=shape_based_distance(c[j,],x[i,])
        if(sbd$dist<minidist){
          minidist=sbd$dist
          idx[i]=j
        }
      }
    }
    iter=iter+1
  }
  return(list(centers=c,idx=idx,initcenters))
}

#中心用相应维度坐标均值计算
sbd_kmeans=function(x,k){
  n=nrow(x)
  m=ncol(x)
  iter=0
  c=x[sample(n,k),]
  idx=NULL
  for(i in 1:n){
    minidist=1/0
    for(j in 1:k){
      sbd=shape_based_distance(c[j,],x[i,])
      if(sbd$dist<minidist){
        minidist=sbd$dist
        idx[i]=j
      }
    }
  }
  idx_prime=rep(0,n)
  while(!all(idx_prime==idx)&iter<100){
    idx_prime=idx
    c <- as.matrix(aggregate(x,list(idx),mean)[,-1])
    for(i in 1:n){
      minidist=1/0
      for(j in 1:k){
        sbd=shape_based_distance(c[j,],x[i,])
        if(sbd$dist<minidist){
          minidist=sbd$dist
          idx[i]=j
        }
      }
    }
    iter=iter+1
  }
  return(list(centers=c,idx=idx))
}

library(tidyverse)
x=matrix(c(1:4,0:3,-1,1,-1,1,1,2,2,3),nrow=4,byrow = T)
x_scale=t(apply(x,1,scale))
ks=kshape(x_scale,2)
ks$idx
ks=ks$centers
ks_tb=as.tibble(ks)
ks_tb$t=1:2
ks_tb=gather(ks_tb,key=n,value,-t)
x_tb=as.tibble(x)
x_tb$t=1:4
x_tb%>%gather(V1:V4,key=n,value=value)%>%ggplot(aes(n,value,group=t))+
  geom_line()+geom_line(data=ks_tb,color='red')


