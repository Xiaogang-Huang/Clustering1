r=0.2
library(rgl)
library(kernlab)
data(spirals)
kd=kd_tree(spirals,1:300)
r_neighbors=list()
for(i in 1:300){
  r_neighbors[[i]]=fix_radius_neighbor(kd,kd[1,1],spirals[i,],r)
}
midu=NULL
for(i in 1:300){
  midu[i]=sum(exp(-r_neighbors[[i]][,2]^2))-1
}
c=NULL
for(i in 1:nrow(spirals)){
  nei=r_neighbors[[i]][,1]
  nei1=nei[which(r_neighbors[[i]][,2]!=0)]
  if(all(midu[nei1]<midu[i])&midu[i]>0){
    c[i]=1
  }else{
    c[i]=0
  }
}
cl=rep(1,nrow(spirals))
cl[which(c==1)]=2
plot(spirals,col=cl)
plot3d(spirals[,1],spirals[,2],midu,col = cl)

library(factoextra)
data("multishapes")
df<-multishapes[,1:2]
df <- as.matrix(df)
kd=kd_tree(df,1:1100)
r_neighbors=list()
r=0.2
for(i in 1:1100){
  r_neighbors[[i]]=fix_radius_neighbor(kd,kd[1,1],df[i,],r)
}
midu=NULL
for(i in 1:1100){
  midu[i]=sum(exp(-r_neighbors[[i]][,2]^2))-1
}
c=NULL
for(i in 1:1100){
  nei=r_neighbors[[i]][,1]
  nei1=nei[which(r_neighbors[[i]][,2]!=0)]
  if(all(midu[nei1]<midu[i])&midu[i]>0){
    c[i]=1
  }else{
    c[i]=0
  }
}
cl=rep(1,1100)
cl[which(c==1)]=2
plot(df,col=cl)
plot3d(df[,1],df[,2],midu,col = cl)

p=which(c==1)
delta=NULL
for(i in 1:1100){
  if(!(i %in% p)){
    nei=r_neighbors[[i]]
    nei=nei[nei[,2]!=0,]
    if(is.matrix(nei)){
      nei=nei[midu[nei[,1]]>midu[i],]
    }
    if(is.matrix(nei)){
      if(nrow(nei)>0){
        delta[i]=nei[which.min(nei[,2]),1]
      }
    }else{
      delta[i]=nei[1]
    }
  }
}
c=rep(0,1100)
k=length(p)
for(i in 1:k){
  t=p[i]
  c[t]=i
  q=which(delta==t)
  repeat{
    c[q]=i
    q=which(delta%in%q)
    if(length(q)==0){
      break
    }
  }
}


plot(df,col=c)
c1=c
for(i in 1:1100){
  nei=r_neighbors[[i]]
  if(length(unique(c1[nei]))>1){
    if(diff(table(c1[nei]))/min(table(c1[nei]))<2){
      c1[which(c1%in%unique(c1[nei]))]=min(unique(c1[nei]))
    }
  }
}
u=unique(c1)
for(i in 1:length(u)){
  c1[which(c1==u[i])]=i
}
table(c1)
plot(df,col=c1)



r=0.6
ywh=iris[,1:4]
kd=kd_tree(ywh,1:150)
r_neighbors=list()
for(i in 1:150){
  r_neighbors[[i]]=fix_radius_neighbor(kd,kd[1,1],ywh[i,],r)
}
midu=NULL
for(i in 1:150){
  midu[i]=sum(exp(-r_neighbors[[i]][,2]^2))-1
}
c=NULL
for(i in 1:150){
  nei=r_neighbors[[i]][,1]
  nei1=nei[which(r_neighbors[[i]][,2]!=0)]
  if(all(midu[nei1]<midu[i])&midu[i]>0){
    c[i]=1
  }else{
    c[i]=0
  }
}
which(c==1)
cl=rep(1,150)
cl[which(c==1)]=2
plot(ywh[,1:4],col=cl)
midu[which(c==1)]




