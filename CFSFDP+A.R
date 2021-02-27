#求点的密度和密度高于该点的点与该点的最近距离
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
fast_m <- function(x,r,k){
  n <- nrow(x)
  a <- kmeans(x,centers = k)
  cen <- a$centers
  clu <- a$cluster
  dm <- t(distmat(x,cen))
  md <- rep(0,n)
  delta <- rep(0,n)
  for(i in 1:n){
    for(l in 1:k){
      h <- clu[i]
      b_upper <- which(abs(dm[,l]-dm[i,l])<r & abs(dm[i,h]-dm[,h])<r & clu==l)
      b_lower <- which((dm[,l]+dm[i,l]<r | dm[i,h]+dm[,h]<r) & clu==l)
      b_n <- b_upper[!b_upper %in% b_lower]
      md[i] <- md[i]+length(b_lower)
      if(length(b_n)>=1){
        for(j in 1:length(b_n)){
          if(euclid(x[b_n[j],],x[i,])<r){
            md[i] <- md[i]+1
          }
        }
      }
    }
  }
  for(i in 1:n){
    if(md[i]==max(md)){
      delta[i]=max(distvec(x[i,],x))
    }else{
      t <- sample(which(md>md[i]),1)
      r <- clu[t]
      h <- clu[i]
      d <- min(dm[i,r]+dm[t,r],dm[i,h]+dm[t,h])
      for(l in 1:k){
        g_upper <- which(abs(dm[,l]-dm[i,l])<d & abs(dm[i,h]-dm[,h])<d & md[i]<md & clu==l)
        g_lower <- which((dm[,l]+dm[i,l]<d | dm[i,h]+dm[,h]<d) & md[i]<md & clu==l)
        if(length(g_lower)>0){
          mi <- apply(cbind(dm[i,l]+dm[,l],dm[i,h]+dm[,h]),1,min)
          ma <- apply(cbind(abs(dm[,l]-dm[i,l]),abs(dm[i,h]-dm[,h])),1,max)
          d_upper <- g_upper[which(mi[g_upper]<=min(ma[g_upper]))]
          d_lower <- g_lower[which(mi[g_lower]<=min(ma[g_lower]))]
          if(length(d_lower)>0){
            d <- min(dm[i,l]+dm[d_lower,l],dm[i,h]+dm[d_lower,h])
          }
          if(length(d_upper)>0){
            for(j in 1:length(d_upper)){
              if(min(dm[i,l]+dm[d_upper[j],l],dm[i,h]+dm[d_upper[j],h])){
                d_z <- euclid(x[i,],x[d_upper[j],])
                if(d_z<d){
                  d <- d_z
                  delta[i] <- d
                  }
                }
              }
            }
          }
        }
      }
    }
  return(list(md,delta))
}
md <- fast_m(spirals,0.2,k=4)

