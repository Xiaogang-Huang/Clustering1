mst_path_dist <- function(dm){
  n <- nrow(dm)
  visited <- rep(0,n)
  visited[1] <- 1 #��¼���ݵ��Ƿ��ѽ�����С������
  d <- dm[1,] #��¼�㵽�ѽ�����������ľ���
  s <- 0 #��¼���·��
  p <- 1 #��¼�����˳��
  e <- matrix(nrow=n-1,ncol = 2) #������˳���¼��������
  pd <- diag(0,nrow = n) #·���������
  while(sum(visited)<n){
    min=Inf
    for(i in 1:n){
      if(visited[i]==0&d[i]<min){
        min <- d[i]
        t <- i
      }
    }#�ҵ����Ѿ��ҹ��ĵ������Сֵ
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
    }#����d����
  }
  return(list(s,p,e,pd))
}