library(readxl)
library(data.table)
library(ggplot2)
library(reshape2)
library(dtw)
library(grid)
fuhe <- read_xls("F:/FUHE/2008.xls",sheet = 1,na = "NULL")
dt <- as.data.table(fuhe)
setnames(dt,c("ÈÕÆÚ"),c("riqi"))
dt <- dt[riqi=="2008-06-01 00:00:00"]
df <- data.frame(dt)
df[,8:103] <- apply(df[,8:103],2,as.numeric)
df <- df[complete.cases(df[,8:103]),]
shuju <- df[,8:103]
b <- shuju[1:10,]
a <- as.data.frame(cbind(b,y=1:10))
a_melt <- melt(a[c(1,3),],id.vars = 'y')
p1 <- ggplot(a_melt,aes(variable,value,group=y))+geom_line()
a_melt <- melt(a[c(1,9),],id.vars = 'y')
p2 <- ggplot(a_melt,aes(variable,value,group=y))+geom_line()
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
vplayout = function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(p1, vp = vplayout(1, 1))
print(p2, vp = vplayout(1, 2))

align <- dtw(as.vector(as.matrix(b[1,])), as.vector(as.matrix(b[3,])), step=asymmetricP1, keep=T)
dtwPlotTwoWay(align)

#·Ö¶Î»Ø¹é,hÄÜ½«ÐòÁÐaµÈ·Ö
slope_seg <- function(a,b,h){
  slope_vec <- NULL
  for (i in 1:(length(a)/h)) {
    z <- h*(i-1)+1
    v <- h*i
    reg <- lm(b[z:v]~a[z:v])
    slope_vec[i] <- reg$coefficients[2]
  }
  return(slope_vec)
}
#shapedtwÖ÷º¯Êý,ÆäÖÐlÄÜ±»gÕû³ý
shapedtw <- function(x,y,l,g){
  n <- length(x)
  m <- length(y)
  s <- matrix(nrow = n,ncol = l)
  t <- matrix(nrow = m,ncol = l)
  k <- floor(l/2)
  x_c <- c(x[(n-k+1):n],x,x[1:k])
  for (i in 1:n) {
    s[i,] <- x_c[i:(i+l-1)]
  }
  dx <- t(apply(s,1,slope_seg,a=1:l,h=g))
  y_c <- c(y[(m-k+1):m],y,y[1:k])
  for (i in 1:m) {
    t[i,] <- y_c[i:(i+l-1)]
  }
  dy <- t(apply(t,1,slope_seg,a=1:l,h=g))
  #Ê±¼äÐòÁÐ¾àÀë¾ØÕó
  julizhen <- matrix(nrow = n,ncol = m)
  for (i in 1:n) {
    for(j in 1:m){
      julizhen[i,j] <- sum(abs(dx[i,]-dy[j,]))
    }
  }
  #ÀÛ»ý¾àÀë¾ØÕó
  align <- array(dim = c(n,m,2))
  align[1,1,] <- 1
  cummat <- matrix(1/0,nrow = n, ncol = m)
  cummat[1,1] <- julizhen[1,1]
  for (i in 2:n){
    cummat[i,1] <- cummat[i-1,1]+julizhen[i,1]
    align[i,1,] <- c(i-1,1)
  }
  for(j in 2:m){
    cummat[1,j] <- cummat[1,j-1]+julizhen[1,j]
    align[1,j,] <- c(1,j-1)
  }
  for (i in 2:n) {
    for(j in 2:m){
      cummat[i,j] <- min(cummat[i-1,j],cummat[i-1,j-1],cummat[i,j-1])+julizhen[i,j]
      before_x <- c(i-1,i-1,i)
      before_y <- c(j,j-1,j-1)
      align[i,j,1] <- before_x[which.min(c(cummat[i-1,j],cummat[i-1,j-1],cummat[i,j-1]))]
      align[i,j,2] <- before_y[which.min(c(cummat[i-1,j],cummat[i-1,j-1],cummat[i,j-1]))]
    }
  }
  d <- align[n,m,]
  lujing_x <- NULL
  lujing_y <- NULL
  while(d[1]!=1&&d[2]!=1) {
    lujing_x <- c(d[1],lujing_x)
    lujing_y <- c(d[2],lujing_y)
    d <- align[d[1],d[2],]
  }
  lujing <- cbind(x=c(1,lujing_x),y=c(1,lujing_y))
  return(list(shapedtw_dist=cummat[n,m],lujing=lujing))
}
a1 <- shapedtw(as.vector(as.matrix(b[1,])),as.vector(as.matrix(b[2,])),48,24)$shapedtw_dist
a2 <- shapedtw(as.vector(as.matrix(b[1,])),as.vector(as.matrix(b[9,])),48,24)$shapedtw_dist
a11 <- dist(b[c(1,3),],method = 'dtw')
a22 <- dist(b[c(1,9),],method = 'dtw')
#Æ¥Åä×÷Í¼
alignroad <- function(i,j,l,g){
  lujing <- shapedtw(as.vector(as.matrix(b[i,])),as.vector(as.matrix(b[j,])),l,g)$lujing
  x_y <- a[c(i,j),]
  colnames(x_y) <- c(1:96,'y')
  x_y$y <- c('a','b')
  x_y_melt <- melt(x_y,id.vars = 'y')
  alignments <- matrix(nrow=2*nrow(lujing),ncol = 3)
  for (i in 1:nrow(lujing)) {
    alignments[2*i-1,] <- c(i,lujing[i,][1],x_y[1,lujing[i,][1]])
    alignments[2*i,] <- c(i,lujing[i,][2],x_y[2,lujing[i,][2]])
  }
  alignments <- as.data.frame(alignments)
  colnames(alignments ) <- c('y','variable','value')
  x_y_melt <- rbind(x_y_melt,alignments)
  y <- x_y_melt$y
  y[1:(length(y)-192)] <-5
  y[(length(y)-192+1):length(y)] <- 1
  y1 <- as.numeric(y)
  return(ggplot(x_y_melt,aes(variable,value,group=y))+geom_line(linetype=y1))
}
p1 <- alignroad(1,2,48,24)
p2 <- alignroad(1,3,48,24)
p3 <- alignroad(1,4,48,24)
p4 <- alignroad(2,3,48,24)
p5 <- alignroad(2,4,48,24)
p6 <- alignroad(3,4,48,24)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 3)))
vplayout = function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
d12 <- shapedtw(as.vector(as.matrix(b[1,])),as.vector(as.matrix(b[2,])),48,24)$shapedtw_dist
d13 <- shapedtw(as.vector(as.matrix(b[1,])),as.vector(as.matrix(b[3,])),48,24)$shapedtw_dist
d14 <- shapedtw(as.vector(as.matrix(b[1,])),as.vector(as.matrix(b[4,])),48,24)$shapedtw_dist
d23 <- shapedtw(as.vector(as.matrix(b[2,])),as.vector(as.matrix(b[3,])),48,24)$shapedtw_dist
d24 <- shapedtw(as.vector(as.matrix(b[4,])),as.vector(as.matrix(b[2,])),48,24)$shapedtw_dist
d34 <- shapedtw(as.vector(as.matrix(b[4,])),as.vector(as.matrix(b[3,])),48,24)$shapedtw_dist
print(p1+labs(title =paste("ÇúÏßÆ¥ÅäÂ·¾¶,ÇÒ¾àÀë:",d12,sep="") ), vp = vplayout(1, 1))
print(p2+labs(title =paste("ÇúÏßÆ¥ÅäÂ·¾¶,ÇÒ¾àÀë:",d13,sep="") ), vp = vplayout(1, 2))
print(p3+labs(title =paste("ÇúÏßÆ¥ÅäÂ·¾¶,ÇÒ¾àÀë:",d14,sep="") ), vp = vplayout(1, 3))
print(p4+labs(title =paste("ÇúÏßÆ¥ÅäÂ·¾¶,ÇÒ¾àÀë:",d23,sep="") ), vp = vplayout(2,1))
print(p5+labs(title =paste("ÇúÏßÆ¥ÅäÂ·¾¶,ÇÒ¾àÀë:",d24,sep="") ), vp = vplayout(2, 2))
print(p6+labs(title =paste("ÇúÏßÆ¥ÅäÂ·¾¶,ÇÒ¾àÀë:",d34,sep="") ), vp = vplayout(2,3))

align <- dtw(as.vector(as.matrix(b[1,])), as.vector(as.matrix(b[2,])), step=asymmetricP1, keep=T)
dtwPlotTwoWay(align)


#dtw¾àÀë
mydtw <- function(x,y){
  n <- length(x)
  m <- length(y)
  julizhen <- matrix(nrow = n,ncol = m)
  for (i in 1:n) {
    for(j in 1:m){
      julizhen[i,j] <- abs(x[i]-y[j])
    }
  }
  #ÀÛ»ý¾àÀë¾ØÕó
  align <- array(dim = c(n,m,2))
  align[1,1,] <- 1
  cummat <- matrix(1/0,nrow = n, ncol = m)
  cummat[1,1] <- julizhen[1,1]
  for (i in 2:n){
    cummat[i,1] <- cummat[i-1,1]+julizhen[i,1]
    align[i,1,] <- c(i-1,1)
  }
  for(j in 2:m){
    cummat[1,j] <- cummat[1,j-1]+julizhen[1,j]
    align[1,j,] <- c(1,j-1)
  }
  for (i in 2:n) {
    for(j in 2:m){
      cummat[i,j] <- min(cummat[i-1,j],cummat[i-1,j-1],cummat[i,j-1])+julizhen[i,j]
      before_x <- c(i-1,i-1,i)
      before_y <- c(j,j-1,j-1)
      align[i,j,1] <- before_x[which.min(c(cummat[i-1,j],cummat[i-1,j-1],cummat[i,j-1]))]
      align[i,j,2] <- before_y[which.min(c(cummat[i-1,j],cummat[i-1,j-1],cummat[i,j-1]))]
    }
  }
  d <- align[n,m,]
  lujing_x <- NULL
  lujing_y <- NULL
  while(d[1]!=1&&d[2]!=1) {
    lujing_x <- c(d[1],lujing_x)
    lujing_y <- c(d[2],lujing_y)
    d <- align[d[1],d[2],]
  }
  lujing <- cbind(x=c(1,lujing_x),y=c(1,lujing_y))
  return(list(dtwdist=cummat[n,m],lujing=lujing))
}
a1_mdtw <- mydtw(as.vector(as.matrix(b[1,])),as.vector(as.matrix(b[3,])))$dtwdist
a2_mdtw <- mydtw(as.vector(as.matrix(b[1,])),as.vector(as.matrix(b[9,])))$dtwdist
dtwdist <- data.frame(shapedtw=c(a1,a2),dtw=c(a11,a22),mydtw=c(a1_mdtw,a2_mdtw),row.names = c("1 3","1 9"))
dtwdist
#Æ¥Åä×÷Í¼
lujing <- mydtw(as.vector(as.matrix(b[1,])),as.vector(as.matrix(b[3,])))$lujing
x_y <- a[c(1,3),]
colnames(x_y) <- c(1:96,'y')
x_y$y <- c('a','b')
x_y_melt <- melt(x_y,id.vars = 'y')
alignments <- matrix(nrow=2*nrow(lujing),ncol = 3)
for (i in 1:nrow(lujing)) {
  alignments[2*i-1,] <- c(i,lujing[i,][1],x_y[1,lujing[i,][1]])
  alignments[2*i,] <- c(i,lujing[i,][2],x_y[2,lujing[i,][2]])
}
alignments <- as.data.frame(alignments)
colnames(alignments ) <- c('y','variable','value')
x_y_melt <- rbind(x_y_melt,alignments)
ggplot(x_y_melt,aes(variable,value,group=y))+geom_line(color=)
