c_step_x <- function(x,y,i,j,c){
  if((x[i-1]<=x[i]&&x[i]<=y[j]) || (x[i-1]>=x[i]&&x[i]>=y[j])) {
    c_cost_X <- c
    }else{
      c_cost_X <- c+min(abs(x[i]-x[i-1]),abs(x[i]-y[j]))
    }
  return(c_cost_X)
}
c_step_y <- function(x,y,i,j,c){
  if((y[j-1]<=y[j]&&y[j]<=x[i]) || (y[j-1]>=y[j]&&y[j]>=x[i])){c_cost_y <- c}else{
    c_cost_y <- c+min(abs(y[j]-y[j-1]),abs(y[j]-x[i]))
  }
  return(c_cost_y)
}
msm_dist <- function(x,y,c){
  n <- length(x)
  m <- length(y)
  cost <- matrix(nrow =n,ncol = m )
  cost[1,1] <- abs(x[1]-y[1])
  for(i in 2:n){
    cost[i,1] <- cost[i-1,1]+c_step_x(x,y,i,j=1,c)
  }
  for(j in 2:m){
    cost[1,j] <- cost[1,j-1]+c_step_y(x,y,1,j,c)
  }
  for(i in 2:m){
    for(j in 2:n){
      cost[i,j] <- min(cost[i-1,j-1]+abs(x[i]-y[j]),cost[i,j-1]+c_step_y(x,y,i,j,c),cost[i-1,j]+c_step_x(x,y,i,j,c))
    }
  }
  return(cost)
}
q <- as.vector(as.matrix(b[1,]))
p <- as.vector(as.matrix(b[9,]))
c1 <- msm_dist(p,q,c=0.1)
dist_19 <- c1[96,96]
q <- as.vector(as.matrix(b[1,]))
p <- as.vector(as.matrix(b[3,]))
c2 <- msm_dist(p,q,c=0.1)
dist_13 <- c2[96,96]
msmdist <- c(dist_13,dist_19)

#×÷Í¼
a_melt1 <- melt(a[c(1,3),],id.vars = 'y')
p1 <- ggplot(a_melt1,aes(variable,value,group=y))+geom_line()
a_melt2 <- melt(a[c(1,9),],id.vars = 'y')
p2 <- ggplot(a_melt2,aes(variable,value,group=y))+geom_line()
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
vplayout = function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(p1, vp = vplayout(1, 1))
print(p2, vp = vplayout(1, 2))