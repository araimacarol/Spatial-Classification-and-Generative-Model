library(igraph)

# (1) Function to generate a single small world graph
sm_world<-function(d,s,n,p){
  sm_graph=sample_smallworld(d, s, n, p, loops = F, multiple = F)
}


##---Test case----
#f2=sm_world(1, 100, 5, 0.05)
#plot(f2)

# (2) Function for Replicating the small world network N times

set.seed(1)
sm_graph=list()
SM_World<-function(d,s,n,p,N){
  for (i in 1:N){
    sm_graph[[i]]=sample_smallworld(d, s, n, p, loops = F, multiple = F)
  }
  return(sm_graph)
}

##----Test case-----
# f2=sm_world(1, 10, 5, 0.05)
# plot(f2[[1]])