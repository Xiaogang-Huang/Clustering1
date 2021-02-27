library(readxl)
library(data.table)
library(ggplot2)
library(reshape2)
library(dtw)
library(grid)
shapedtw_d <- function(x,y,l,g){
  d <- shapedtw(x,y,l,g)$shapedtw_dist
  return(d)
}
distvec <- function(x,t,l,g){
  a <- apply(t,1,shapedtw_d,x=x,l=l,g=g)
  return(a)
}
distmat <- function(s,t,l,g){
  a <- apply(s,1,distvec,t=t,l=l,g=g)
  return(a)
}
delta <- function(data_set,r,l,g){
  dm <- distmat(data_set,data_set,l,g)
  midu=colSums(dm<r)
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
  return(list(midu,delta,match_vector))
}
find_center <- function(data_set,midu,delta,k,l,g,theta){
  a=delta
  a[midu<=2]=0
  centers_index=NULL
  centers_index <- which(a==sort(a,decreasing = T)[1])
  i <- 2
  while(length(centers_index)<k){
    centers_candidate <- which(a==sort(a,decreasing = T)[i])
    cen_dist <- distvec(data_set[centers_candidate,],data_set[centers_index,],l,g)
    if(any(cen_dist<theta)) {}else{
      centers_index <- c(centers_index ,centers_candidate)
    }
    i <- i+1
  }
  return(centers_index)
}
assigned <- function(data_set,midu,centers_index,l,g){
  assignment <- rep(0,nrow(data_set))
  assignment[centers_index] <- 1:length(centers_index)
  m <- midu
  outliers_index <- which(m==1)
  m[centers_index] <- 0
  i <- 1
  dm <- distmat(data_set,data_set,l,g)
  t <- 1
  first_non_zero <- function(v){
    return(v[which(v!=0)[1] ])
  }
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
        assignment[num[j]] <- first_non_zero(v)
      }
    }
    i <- i+length(num)
    t <- t+length(num)
  }
  return(rowSums(assignment))
}
#数据
fuhe <- read_xls("F:/FUHE/2008.xls",sheet = 1,na = "NULL")
dt <- as.data.table(fuhe)
setnames(dt,c("日期"),c("riqi"))
dt <- dt[riqi=="2008-06-01 00:00:00"]
df <- data.frame(dt)
df[,8:103] <- apply(df[,8:103],2,as.numeric)
df <- df[complete.cases(df[,8:103]),]
shuju <- as.matrix(df[1:20,8:103])
mmidu <- midu(shuju,r=1000,l=48,g=24)
mdelta <- delta(shuju,r=1000,48,24)
mdelta1 <- mdelta[[1]]
match_vector<- mdelta[[2]]
plot(mmidu,mdelta1)
outliers <- shuju[which(midu==0),]
centers_index <- find_center(shuju,mmidu,mdelta1,3)
centers <- shuju[centers_index,]
fclust <- assigned(shuju,midu=mmidu,centers_index = centers_index,48,24)
fclust[which(mdelta[[2]]==1)] <- 0
data_set_clus <- cbind(shuju,fclust)


dm <- distmat(shuju,shuju,48,24)

fclust <- assigned1(shuju,midu=mmidu,centers_index = centers_index,48,24,dm)

assigned1 <- function(data_set,midu,centers_index,l,g,dm){
  assignment <- matrix(0,nrow=nrow(data_set),ncol = length(centers_index))
  assignment[centers_index,] <- diag(1:length(centers_index),nrow = length(centers_index))
  m <- midu
  outliers_index <- which(m==0)
  m[centers_index] <- 0
  i <- 1
  t <- 1
  while(i<nrow(data_set)){
    tth_midu <- sort(m,decreasing = T)[t]
    num <- which(m==tth_midu)
    for(j in 1:length(num)){
      if(num[j] %in% centers_index || num[j] %in% outliers_index){ i <- i+1}else{
        dist <- dm[j,]*(midu>tth_midu)
        dist[midu<=tth_midu] <- 1/0
        assignment[num[j],] <- assignment[which.min(dist),]
      }
    }
    i <- i+length(num)
    t <- t+length(num)
  }
  return(rowSums(assignment))
}
fclust <- assigned1(shuju,midu=mmidu,centers_index = centers_index,48,24,dm)
fc <- fclust[!fclust==0]
shuju_clus <- cbind(shuju[-as.vector(outliers_index),],y=1:nrow(shuju[-as.vector(outliers_index),]))
shuju_melt <- melt(shuju_clus,id.vars='y')
ggplot(shuju_melt,aes(Var2,value,group=Var1))+geom_line(aes(color=rep(fc,97)))+scale_colour_gradientn(colours=rainbow(3)) 








