library(igraph)

#Function for a sngle  Spatial network
#Spatial graph of x,y coordinate join in space
spatial <- function(n,r) {
  v <- NULL
  for (i in 1:n) {
    x <- runif(1)
    y <- runif(1)
    v <- rbind(v,c(i,x,y))
    coords<<-rbind(coords, c(x,y))
  }
  
  E <- NULL
  R2 <- r*r;
  for (i in seq(1,n-1)) {
    for (j in seq(i+1,n)) {
      if (((v[i,2]-v[j,2])^2 + (v[i,3]-v[j,3])^2) < R2) {
        E <- as.data.frame(rbind(E,c(i,j)))
      }
    }
  }
  return(E)
}

coords=NULL

E = spatial(10,0.3)

G=graph_from_data_frame(E, directed=FALSE)

V(G)$color<-"grey"
G$layout=coords

##--Test case---
plot(G)

## (2) Replicating the Spatial network N times
set.seed(1)
spatial_net=list()
SPATIAL<-function(n,r,N){
  for (i in 1:N){
    spatial_net[[i]]=spatial(n,r)
  }
  return(spatial_net)
}
# h1=SPATIAL(20,.4,5)
 #h2=graph_from_data_frame(h1[[5]], directed=FALSE)
#plot(h2)
