#距离计算函数
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
#dijkstra算法,输入数据集
dijkstra <- function(data_set,i){
  dist_vec <- function(data_set,v){
    vec <- NULL
    for(i in 1:(length(v)-1)){
      vec[i] <- euclid(data_set[v[i],],data_set[v[i+1],])
    }
    return(vec)
  }
  n <- nrow(data_set)
  dm <- as.matrix(dist(data_set,diag = TRUE))
  juli <- NULL
  path <- as.list(rep(i,n))
  visit <- rep(0,n)
  visit[i] <- 1
  juli = dm[i,]
  juli[i] <- 0
  max_path <- matrix(nrow = n,ncol = 2)
  for(j in 1:n){
    min_cost <- 1/0
    for(t in 1:n){
      if(visit[t]==0&juli[t]<min_cost){
        min_cost <-  juli[t]
        min_cost_index  <-  t
      }
    }
    visit[min_cost_index] <- 1
    for(k in 1:n){
      if (dm[min_cost_index,k]!=1/0&max(dm[min_cost_index,k],min_cost) < juli[k]){
        juli[k] = max(dm[min_cost_index,k],min_cost)
        path[[k]] <- c(path[[min_cost_index]],min_cost_index)
      }
    }
  }
  for(l in 1:n){
    path[[l]] <- c(path[[l]],l)
  }
  #数据点向量对应的距离向量
  for(j in 1:n){
    ve <- path[[j]]
    d <- which.max(dist_vec(data_set,ve))
    max_path[j,] <- ve[d:(d+1)]
  }
  return(list(juli=juli,path=path,max_path=max_path))
}
apsp_dijkstra <- function(data_set){
  n <- nrow(data_set)
  ma <- matrix(nrow = n,ncol = n)
  for(i in 1:n){
    ma[i,] <- dijkstra(data_set,i)[[1]]
  }
  return(ma)
}
#最小生成树
prim <- function(data_set){
  n <- nrow(data_set)
  dm <- as.matrix(dist(data_set,diag = T))
  visited <- rep(0,n)
  visited[1] <- 1
  d <- dm[1,]
  s <- 0
  p <- 1
  q <- 1
  e <- matrix(nrow=n-1,ncol = 2)
  c <- 1
  while(c<n){
    min=Inf
    for(i in 1:n){
      if(visited[i]==0&d[i]<min){
        min <- d[i]
        t <- i
      }
    }#找到离已经找过的点距离最小值
    e[c,2] <- t
    e[c,1] <- q[which(dm[visited==1,t]==min)[1]]
    visited[t] <- 1
    s <- s+min
    p <- c(p,t)
    q <- sort(p)
    for(j in 1:n){
      if(visited[j]==0&d[j]>dm[t,j]){
        d[j] <- dm[t,j]
      }
    }#更新d向量
    c <- c+1
  }
  return(list(s,p,e))
}
#计算密度、距离和匹配向量
delta <- function(data_set,r){
  dm <- as.matrix(dist(data_set,diag = TRUE))
  midu=colSums(dm<r)
  dm <- floyd_path_distance(dm)
  delta=NULL
  match_vector <- rep(0,nrow(data_set))
  for (i in 1:nrow(data_set)) {
    d <- dm[i,]
    if(midu[i]==max(midu)){
      delta[i]=max(d)
      match_vector[i] <- which.max(d)
    }
    else{
      delta[i]=min(d[which(midu>midu[i])]) 
      d[which(midu<=midu[i])]=1/0
      match_vector[i] <- which.min(d)
    }
  }
  return(list(delta,midu,match_vector,floyd_matrix=dm))
}
#寻找类中心
find_center <- function(data_set,midu,delta,k,theta){
  a=delta
  a[midu<=8]=0
  centers_index=NULL
  centers_index <- which(a==sort(a,decreasing = T)[1])[1]
  i <- 2
  while(length(centers_index)<k){
    centers_candidate <- which(a==sort(a,decreasing = T)[i])
    if(which.max(a) %in% centers_candidate){centers_candidate <- centers_candidate[-1]}
    for(j in 1:length(centers_candidate)){
      if(length(centers_index)==1){t=t(as.matrix(data_set[centers_index,]))}else{
        t <- data_set[centers_index,]
      }
      cen_dist <- distvec(data_set[centers_candidate[j],],t)
      if(any(cen_dist<theta)) {}else{
        centers_index <- c(centers_index ,centers_candidate[j])
      }
    }
    i <- i+length(centers_candidate)
  }
  if(length(centers_index)>k){centers_index <- centers_index[1:k]}
  return(centers_index)
}
#数据点分类
assigned <- function(data_set,midu,centers_index,dm,outliers=T){
  assignment <- rep(0,nrow(data_set))
  assignment[centers_index] <- 1:length(centers_index)
  m <- midu
  if(outliers==T){
    outliers_index <- which(m<=3)
  }else{outliers_index <- NULL}
  m[centers_index] <- 0
  i <- 1
  t <- 1
  while(i<nrow(data_set)){
    tth_midu <- sort(m,decreasing = T)[t]
    num <- which(m==tth_midu)
    for(j in 1:length(num)){
      if(num[j] %in% centers_index || num[j] %in% outliers_index){}else{
        dist <- dm[num[j],]
        dist[num[j]] <- 1/0
        dist[midu<tth_midu] <- 1/0
        x <- cbind(dist,assignment)
        v <- x[order(x[,1]),][,2]
        assignment[num[j]] <- v[which(v!=0)[1]]
      }
    }
    i <- i+length(num)
    t <- t+length(num)
  }
  assignment[outliers_index] <- length(centers_index)+1
  return(assignment)
}
library(factoextra)
library(ggplot2)
data("multishapes")
dt<-multishapes[,1:2]
set.seed(1)
df <- as.data.frame(dt[sort(sample(1100,500)),])
d <- NULL
max_p <- NULL
for(j in 1:nrow(df)){
  a <- dijkstra(df,j)
  max_path <- a$max_path
  lujing <- a$path
  for(i in 1:nrow(df)){
    #d <- rbind(d,cbind(df[lujing[[i]],],i))
    max_p <- rbind(max_p,cbind(df[max_path[i,],],(j-1)*nrow(df)+i))
  }
}
colnames(max_p)[3] <- 'i'
p1 <- ggplot(df,aes(x,y))+geom_point()
p1+geom_path(data=d,aes(x,y,group=i))+
  geom_line(data=max_p,aes(x,y,group=i),colour='blue')
p1+geom_line(data=max_p,aes(x,y,group=i),colour='blue')
es <- prim(df)[[3]]
eds <- NULL
for(i in 1:nrow(es)){
  eds <- rbind(eds,cbind(df[es[i,],],i))
}
p1+geom_line(data=eds,aes(x,y,group=i),colour='blue')


t <- 150
a <- dijkstra(df,t)
lujing <- a$path
d <- NULL
for(i in 1:nrow(df)){
  d <- rbind(d,cbind(df[lujing[[i]],],i))
  }
p1+geom_path(data=d,aes(x,y,group=i))+
  geom_line(data=max_p[(t-1)*nrow(df)<max_p$i&max_p$i<=t*nrow(df),],aes(x,y,group=i),colour='red')+
  geom_text(data=df[t,],aes(x,y-0.1,label=t))











