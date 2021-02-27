grid_divide = function(x,gamma){
  n <- nrow(x)
  x_axis_low <-  min(x[,1])
  x_axis_up <-  max(x[,1])
  y_axis_low <- min(x[,2])
  y_axis_up <- max(x[,2])
  k_x=ceiling((x_axis_up-x_axis_low)/gamma)
  k_y=ceiling((y_axis_up-y_axis_low)/gamma)
  grid_index <- NULL
  for(i in 1:n){
    dot <- x[i,]
    grid_index[i] <- as.numeric(floor((dot[1]-x_axis_low)/gamma)+floor((dot[2]-y_axis_low)/gamma)*k_x+1)
  }
  return(list(grid_index,k_x,k_y))
}
findneighbors=function(i,k_x,k_y,grid_indexs){
  neis=NULL
  if(i%%k_x==0){
    l=c(i-k_x,i-k_x-1,i-1,i+k_x-1,i+k_x)
    neis=l[l%in%grid_indexs]
  }else if(i%%k_x==1){
    l=c(i-k_x,i-k_x+1,i+1,i+k_x,i+k_x+1)
    neis=l[l%in%grid_indexs]
  }else{
    l=c(i-k_x-1,i-k_x,i-k_x+1,i-1,i+1,i+k_x-1,i+k_x,i+k_x+1)
    neis=l[l%in%grid_indexs]
  }
  return(neis)
}
Lin_DBSCAN <- function(x,gamma,minpts){
  grids=grid_divide(x,gamma)
  k_x=grids[[2]]
  k_y=grids[[3]]
  grids=grids[[1]]
  grid_indexs=sort(unique(grids))
  grid_density=table(grids)
  gn=length(grid_indexs)
  n <- nrow(x)
  c <- rep(0,gn)
  noise <- rep(0,gn)
  k <- 0
  visited <- rep(0,gn)
  for(i in 1:gn){
    if(visited[i]==0){
      visited[i] <- 1
      cell=grid_indexs[i]
      neighborpts=findneighbors(cell,k_x,k_y,grid_indexs)
      neighborpts_index=which(grid_indexs%in%neighborpts)
      m=length(neighborpts)
      grid_pts=grid_density[i]
      if(grid_pts<minpts){
        noise[i] <- 1
      }else{
        k <- k+1
        c[i] <- k
        t <- 1
        while(!all(visited[neighborpts_index]==1)){
          for(j in t:m){
            p <- neighborpts_index[j]
            cell_p=neighborpts[j]
            if(visited[p]==0){
              visited[p] <- 1
              neighborpts_new <- findneighbors(cell_p,k_x,k_y,grid_indexs)
              neighborpts_index_new=which(grid_indexs%in%neighborpts_new)
              m_new <- length(neighborpts_new)
              grid_pts_new=grid_density[p]
              if(grid_pts_new>=minpts){
                neighborpts <- union(neighborpts,neighborpts_new)
                neighborpts_index <- union(neighborpts_index,neighborpts_index_new)
                if(c[p]==0){
                  c[p] <- k
                }
              }else{
                noise[p] <- 1
                c[p] <- k
              }
            }
          }
          t <- m+1
          m <- length(neighborpts)
        }
      }
    }
  }
  c[which(noise==1&c==0)] <- k+1
  clust=NULL
  for(i in 1:(k+1)){
    grid_i=grid_indexs[which(c==i)]
    points_i=which(grids%in%grid_i)
    clust[points_i]=i
  }
  return(clust)
}


library(factoextra)
data("multishapes")
df<-multishapes[,1:2]
eps=0.15
minpts=1
co=Lin_DBSCAN(df,0.6,5)
plot(df,col=co)
