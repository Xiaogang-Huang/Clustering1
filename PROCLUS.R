manhattan_segmental_distance <- function(x,y,d){
  d <- sum(abs(x[d]-y[d]))/length(d)
  return(d)
}
greedy <- function(s,k){
  n <- nrow(s)
  dist_matrix <- diag(1,nrow = n)
  for (i in 1:(n-1)) {
    for(j in (i+1):n){
      dist_matrix[i,j] <- sqrt(sum((s[i,]-s[j,])^2))
    }
  }
  dist_matrix <- dist_matrix+t(dist_matrix)-diag(1,nrow = n)
  m <- sample(1:n,1)
  dist_x <- dist_matrix[m,]
  for (i in 2:k) {
    m <- c(m,which.max(dist_x))
    dist_x <- apply(dist_matrix[,m], 1, min)
  }
  return(s[m,])
}
FindDimensions <- function(s,medoid,l,big_l,d){
  k <- nrow(medoid)
  y <- NULL
  x <- matrix(nrow = k,ncol = d)
  z <- matrix(nrow = k,ncol = d)
  D <- matrix(0,nrow = k,ncol = d)
  delta <- NULL
  for(i in 1:k){
    x[i,] <- colMeans(abs(s[big_l==i,]-matrix(medoid[i,], nrow =sum(big_l==i) ,ncol =d, byrow = T))) 
    y[i] <- sum(x[i,])/d
    delta[i] <- sqrt((sum(x[i,]-y[i])^2)/(d-1))
    for(j in 1:d){
      z[i,j] <- (x[i,j]-y[i])/delta[i]
    }
    t <- sort(unique(z[i,]))
    u <- 1
    while(u<=k*l||sum(D[i,])<2){
      D[i,which(z[i,]==t[1])] <- 1
      u <- u+1 
    }
  }
  return(D)
}
assignpoints <- function(D,medoid,s){
  n <- nrow(s)
  k <- nrow(medoid)
  assignments <- NULL
  p <- matrix(nrow = n,ncol = k)
  for(i in 1:n){
    for(j in 1:k){
      p[i,j] <- manhattan_segmental_distance(s[i,],medoid[j,],d=D[j,])
    }
  }
  assignments <- apply(p,1,which.min)
  return(assignments)
}
evaluatecluster <- function(assignments,medoid,D,s){
  N <- nrow(s)
  k <- nrow(D)
  d <- ncol(D)
  y <- matrix(0,nrow = k,ncol = d)
  w <- NULL
  for(i in 1:k){
    c <- s[assignments==i,]
    n <- nrow(c)
    for(j in which(D[i,]==1)){
      y[i,j] <- colMeans(abs(c-matrix(medoid[i,], nrow = n ,ncol =d, byrow = T)))[j]
    }
    w[i] <- n*sum(y[i,])/sum(D[i,]==1)
  }
  return(sum(w)/N)
}
dist_matrix <- function(s,t){
  n <- nrow(s)
  m <- nrow(t)
  dist_matrix <- matrix(nrow = n,ncol = m)
  for (i in 1:n) {
    for (j in 1:m) {
      dist_matrix[i,j] <- sqrt(sum((s[i,]-t[j,])^2))
    }
  }
  return(dist_matrix)
}
proclus <- function(k,l,s,A,B,mindeviation){
  n <- nrow(s)
  d <- ncol(s)
  candidate_medoids <- s[sample(1:n,A*k),]
  M <- greedy(candidate_medoids,B*k)
  bag <- 1:nrow(M)
  bestobjective <- 1/0
  sam <- sample(bag,k)
  M_current <- M[sam,]
  while(length(th)==0){#µü´úÍ£Ö¹Ìõ¼ş
    big_l <- NULL
    distmat <- dist_matrix(M_current,M_current)
    for (i in 1:k) {
      delta_m <- apply(distmat, 1, min)
      distmat_current <- dist_matrix(s,M_current)
      big_l[distmat_current[,i]<delta_m[i]] <- i
    }
    D_current <- FindDimensions(s,M_current,l,big_l,d)
    assignment <- assignpoints(D_current,M_current,s)
    objectivefunction <- evaluatecluster(assignment,M_current,D,s)
    if(objectivefunction<bestobjective){
      bestobjective <- objectivefunction
      M_best <- M_current
      a <- as.data.frame(table(assignment))
      bad_medoids <- as.numeric(a$assignment[a$Freq<n*mindeviation/k])
    }
    bag <- bag[-sam[bad_medoids]]
    th <- sample(bag,length(bad_medoids))
    sam[bad_medoids] <- th
    M_current <- rbind(M_best[-bad_medoids,],M[th,])
  }
  big_l <- assignment
  D <- FindDimensions(s,M_best,l,big_l,d)
  assignment <- assignpoints(D,M_best,s)
  return(list(medoids=M_best,D=D,assignment=assignment))
}