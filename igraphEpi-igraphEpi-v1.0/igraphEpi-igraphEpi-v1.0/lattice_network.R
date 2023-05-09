library(igraph)

# (1) Function to generate a single lattice graph
lattice_net<-function(d,l){
  lattice_graph=make_lattice(c(dim=d,length=l), directed= F, mutual= T, circular = F)
}



# (2) Function for Replicating the lattice network N times

lattice_graph=list()
Lattice_net<-function(d,l,N){
  for (i in 1:N){
    lattice_graph[[i]]=make_lattice(c(dim=d,length=l), directed= F, mutual= T, circular = F)
  }
  return(lattice_graph)
}

##---Test case----
#f3=Lattice(2,5)
#plot(f3)