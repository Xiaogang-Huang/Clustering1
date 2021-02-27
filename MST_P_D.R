mst_path_dist <- function(dm){
  n <- nrow(dm)
  visited <- rep(0,n)
  visited[1] <- 1 #记录数据点是否已进入最小生成树
  d <- dm[1,] #记录点到已进入生成树点的距离
  s <- 0 #记录最短路径
  p <- 1 #记录点进入顺序
  e <- matrix(nrow=n-1,ncol = 2) #按进入顺序记录生成树边
  pd <- diag(0,nrow = n) #路径距离矩阵
  while(sum(visited)<n){
    min=Inf
    for(i in 1:n){
      if(visited[i]==0&d[i]<min){
        min <- d[i]
        t <- i
      }
    }#找到离已经找过的点距离最小值
    e[sum(visited),2] <- t
    e[sum(visited),1] <- sort(p)[which(dm[visited==1,t]==min)[1]]
    for(k in 1:sum(visited)){
      pd[which(visited==1)[k],t] <- pd[which(visited==1)[k],e[sum(visited),1]]+
        dm[e[sum(visited),1],t]
      pd[t,which(visited==1)[k]] <- pd[which(visited==1)[k],e[sum(visited),1]]+
        dm[e[sum(visited),1],t]
    }
    visited[t] <- 1
    s <- s+min
    p <- c(p,t)
    for(j in 1:n){
      if(visited[j]==0&d[j]>dm[t,j]){
        d[j] <- dm[t,j]
      }
    }#更新d向量
  }
  return(list(s,p,e,pd))
}
