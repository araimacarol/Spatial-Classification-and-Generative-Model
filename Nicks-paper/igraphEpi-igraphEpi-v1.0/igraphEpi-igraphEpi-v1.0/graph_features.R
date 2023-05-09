library(igraph)

GraphFeatures <- function(Graphs=NULL,nrep=NULL){

  
  Global_summary <- NULL
  
  network_size=vector(mode = "list", length(Graphs))
  FiedlerValue=vector(mode = "list", length(Graphs))
  Adj_val=vector(mode = "list", length(Graphs))
  Adj=vector(mode = "list", length(Graphs))
  network_size=vector(mode = "list", length(Graphs))
  most_infected_node=vector(mode = "list", length(Graphs))
  mod=vector(mode = "list", length(Graphs))
  deg=vector(mode = "list", length(Graphs))
  cent=vector(mode = "list", length(Graphs))
  trans=vector(mode = "list", length(Graphs))
  di=vector(mode = "list", length(Graphs))
  M=vector(mode = "list", length(Graphs));spec=vector(mode = "list", length(Graphs)); 
  Eigen_central=vector(mode = "list", length(Graphs));
  most_infected_nod=vector(mode = "list", length(Graphs))
  Rnot=vector(mode = "list", length(Graphs));wc=vector(mode = "list", length(Graphs));
  lc=vector(mode = "list", length(Graphs));Adjacency_spectrum=vector(mode = "list", length(Graphs))
  
for (g in 1:length(Graphs)){
  for(rep in 1:nrep){   
    network_size[[g]][[rep]]<-vcount(Graphs[[g]])
    #Adjacency matrix
    Adj[[g]][[rep]]<-as.matrix(Graphs[[g]][])
    #Adjacency spectrum
    Adjacency_spectrum[[g]][[rep]]<-eigen(Adj[[g]][[rep]])
      
    
    M[[g]][[rep]]<- laplacian_matrix(Graphs[[g]],sparse = FALSE)
    
    spec[[g]][[rep]]<- eigen(M[[g]][[rep]])

    FiedlerValue[[g]][[rep]]<- spec[[g]][[rep]]$values[length(spec[[g]][[rep]]$value)-1]
      
    
    Adj_val[[g]][[rep]]<-Adjacency_spectrum[[g]][[rep]]$values[1]
    
    Eigen_central[[g]][[rep]]<- Adjacency_spectrum[[g]][[rep]]$vectors[,ncol(Adjacency_spectrum[[g]][[rep]]$vectors)[1]]
  
    most_infected_node[[g]][[rep]]<- match(max(Eigen_central[[g]][[rep]]),Eigen_central[[g]][[rep]])# highest node degree

    Rnot[[g]][[rep]]<-(1/Adj_val[[g]][[rep]]) # spectral of adjacency for undirected network controls the epidemic threshold

    wc[[g]][[rep]]<- cluster_walktrap(Graphs[[g]])  
    lc[[g]][[rep]]<- cluster_louvain(Graphs[[g]]) 
    
    mod[[g]][[rep]]<- modularity(Graphs[[g]], membership(lc[[g]][[rep]]))
     
    deg[[g]][[rep]]<- mean(degree(Graphs[[g]], mode="in"))
      
   # cent[[g]][[rep]]<-  eigen_centrality(G_list[[g]], directed=T, weights=NA)
   # cent[[g]][[rep]] <- cent[[g]]$value
    
  #  trans[[g]]<- as.data.frame(transitivity(G_list[[g]], type="global"))
     
    di[[g]][[rep]]<- diameter(Graphs[[g]], directed=F, weights=NA)
    
    
    Global_summary[[g]]<- as.data.frame(cbind(FiedlerValue[[g]], Adj_val[[g]], network_size[[g]], most_infected_node[[g]], deg[[g]], di[[g]])) 
  }
} 
    Global_summary <- do.call(rbind,  Global_summary)
#    Global_summary  <- as.data.frame(Global_summary)
  
    colnames(Global_summary) <- c('Fiedler','Adj_val','network_size','most_infected_node', 'Mean_degree','diameter')

   return(Global_summary)
}


