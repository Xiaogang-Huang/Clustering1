kd_tree=function(data,kd,index){
  n_data=nrow(data)
  n=length(index)
  k=ncol(data)
  if(n>0){
    split=0
    s=-1/0
    for(i in 1:k){
      a=my_var(data[index,i])
      if(a>s){
        s=a
        split=i
      }
    }
    node=index[rank(data[index,split],ties.method = 'first')==floor(n/2+1)]
    a=data[node,split]
    left=NULL
    right=NULL
    for(i in index[!index%in%node]){
      if(data[i,split]>a){
        right=c(right,i)
      }else{
        left=c(left,i)
      }
    }
    kd[1,1]=node+1         ##根节点位置
    kd[node+1,c(1,2,5)]=c(split,a,1)
    if(length(left)>0){
      kd_left=kd_tree(data,kd,left)
      kd_left_root=kd_left[1,1]
      kd[left+1,1:5]=kd_left[left+1,1:5]
      kd[kd_left_root,5]=node+1
      kd[node+1,3]=kd_left_root
    }
    if(length(right)>0){
      kd_right=kd_tree(data,kd,right)
      kd_right_root=kd_right[1,1]
      kd[right+1,1:5]=kd_right[right+1,1:5] 
      kd[kd_right_root,5]=node+1
      kd[node+1,4]=kd_right_root
    }
  }
  return(kd)
}


fix_radius_neighbor=function(data,t,root,x,r){
  if(is.na(root)){
    return(0)
  }
  p=root
  split=t[p,1]
  a=t[p,2]
  m=0
  while(!all(is.na(t[p,3:4]))){
    split=t[p,1]
    if(is.na(t[p,3])){
      p=t[p,4]
    }else if(x[split]<=data[p-1,split]|is.na(t[p,4])){
      p=t[p,3]
    }else{
      p=t[p,4]
    }
  }
  d=sqrt(sum((x-data[p-1,])^2))
  if(d<r){
    m=m+1
  }
  q=p
  tmp=q
  while(q!=root){
    q=t[q,5]
    d=sqrt(sum((x-data[q-1,])^2))
    if(d<r){
      m=m+1
    }
    split=t[q,1]
    x_to_split=abs(x[split]-t[q,2])
    if(x_to_split<r){
      if(tmp==t[q,3]){
        node_tmp=t[q,4]
        tmpresult=fix_radius_neighbor(data,t,node_tmp,x,r)
      }else if(tmp==t[q,4]){
        node_tmp=t[q,3]
        tmpresult=fix_radius_neighbor(data,t,node_tmp,x,r)
      }
      m=m+tmpresult
    }
    tmp=q
  }
  return(m)
}
library(kernlab)
data(spirals)

Rprof()
kd=matrix(nrow=301,ncol=5)
kd=kd_tree(spirals,kd,1:300)
midu=list()
for(i in 1:300){
  midu[i]=nrow(fix_radius_neighbor(spirals,kd,kd[1,1],spirals[i,],r))
}
Rprof(NULL)
summaryRprof()

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
Rprof()
midu1=colSums(distmat_f(spirals,spirals)<r)
Rprof(NULL)
summaryRprof()


library(factoextra)
data("multishapes")
df=as.matrix(multishapes[,1:2])
Rprof()
kd=matrix(nrow=1101,ncol=5)
kd=kd_tree(df,kd,1:1100)
midu=list()
for(i in 1:1100){
  midu[[i]]=nrow(fix_radius_neighbor(df,kd,kd[1,1],df[i,],r))
}
Rprof(NULL)
summaryRprof()

Rprof()
kd=matrix(nrow=1101,ncol=5)
kd=kd_tree(df,kd,1:1100)
midu=NULL
for(i in 1:1100){
  midu[i]=fix_radius_neighbor(df,kd,kd[1,1],df[i,],r)
}
Rprof(NULL)
summaryRprof()


Rprof()
midu1=colSums(distmat_f(df,df)<r)
Rprof(NULL)
summaryRprof()

r=0.2
x=read_csv("F:/聚类/聚类展示/DPC+K+P/Gaussiandata.csv")
x=as.matrix(x[sample(15000,6000,replace = F),])
Rprof()
kd=matrix(nrow=6001,ncol=5)
kd=kd_tree(x,kd,1:6000)
midu=NULL
for(i in 1:6000){
  midu[i]=fix_radius_neighbor(x,kd,kd[1,1],x[i,],r)
}
Rprof(NULL)
summaryRprof()


Rprof()
midu1=rep(0,6000)
for(i in 1:6000){
  for(j in 1:6000){
    midu1[i]=midu1[i]+(sqrt(sum((x[i,]-x[j,])^2))<r)
  }
}
Rprof(NULL)
summaryRprof()


