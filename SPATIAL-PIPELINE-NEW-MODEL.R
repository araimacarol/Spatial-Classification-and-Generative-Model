setwd("C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi-New")
library(Matrix)
library(stats)
#library(spatstat)
library(igraph)
library(tidymodels)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Data generation for the machine learning model
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###---Simulate--spatial--networks---##
# simulate.spatial<-function(N=50,radius=0.4,nsim=100){
#   spatial.graph=NULL
#   for (i in 1:nsim){
#     spatial.graph[[i]]=makeSpatialGraphs(node.size=N,Radius=radius)
#   }
#   return(spatial.graph)
# }

##########################################################################################################
# Input an igraph object from a file, treating it as a static graph
##########################################################################################################
getGraphFromFile <- function(file, simplify=TRUE, useBiggestComponent=TRUE, asUndirected=TRUE) {
  
  dat <- read.table(file) # just read static graph: ignoring third column
  G <- graph_from_data_frame(dat)
  if (asUndirected==TRUE) {
    G <- as.undirected(G, "collapse")
  }
  
  #  g_names <- gsub(".edges","",networks[i]) # one edge for each pair of connect vertices (not sure what this is for)
  
  if (useBiggestComponent==TRUE) {
    netclust <- clusters(G) #look for subgraphs
    gcc <- V(G)[netclust$membership == which.max(netclust$csize)]#select vertices from the largest sub-graph
    G <- induced.subgraph(G, gcc) #make it a igraph object.
  }
  if (simplify==TRUE) {
    G <- igraph::simplify(G, remove.multiple = TRUE, remove.loops = TRUE)
  }
  return(G)
}

###--Normalized Laplacina function--##
normalized_laplacian=function(Graphs){
  laplacian_matrix(Graphs,normalized = T)
}



calcGraphFeatures <- function(Graphs=NULL) {
  
  features <- c(
    "order",                    # number of vertices
    "edges",                     # number of edges
    "connected",                # True / False
    "max_component",            # maximum component size (=order iff the graph is connected)
    "minDegree",                # minimum degree of any vertex
    "maxDegree",                # maximum degree of any vertex
    "mean_degree",                # average degree of any vertex
    "minCut",                   # minimum cut weight of the graph (might take a while to compute)
    "FiedlerValue",             # second-highest eigenvalue of the Laplacian matrix
    "Normalized_FiedlerValue",   # second-highest eigenvalue of the Normaized Laplacian matrix
    "closeness_centr",                # average inverse of distance between any pair of vertices
    "modularity",               # DEFINITION REQUIRED
    "diameter",                 # maximum distance between any two vertices (NAN if not connected)
    "betw_centr",              # max_{v} proportion of shortest paths going through vertex v
    "transitivity",             # aka Clustering Coefficient, is proportion of connected triples that form triangles: e.g., (a--b--c--a) when (a--b--c) is present.
    "threshold",                 # 1/max(eigen value of A)
    "spectral_radius"         # max (eigen value of A)
    
  )
  
  df <- as.data.frame(matrix(ncol=length(features),nrow=length(Graphs)))
  colnames(df)=features
  
  # Stuff that is simple to apply and needs no interim components:
  
  df$order = base::as.numeric(lapply(Graphs, gorder))
  df$edges = base::as.numeric(lapply(Graphs, gsize))
  df$connected = base::as.numeric(lapply(Graphs, is_connected))
  df$minCut = base::as.numeric(lapply(Graphs, min_cut))
  df$diameter = base::as.numeric(lapply(Graphs, diameter))
  df$transitivity = base::as.numeric(lapply(Graphs, transitivity))
  
  # stuff that needs interim things:
  degrees = lapply(Graphs,igraph::degree )
  df$minDegree = base::as.numeric(lapply(degrees, min))
  df$maxDegree = base::as.numeric(lapply(degrees, max))
  df$mean_degree = base::as.numeric(lapply(degrees, mean))
  
  # stuff that lapply doesn't like so has to be done in a loop:
  communities <-lapply(Graphs, cluster_walktrap) #lapply(Graphs, cluster_leading_eigen)
  Adj <- lapply(Graphs, as_adjacency_matrix)
  L <- lapply(Graphs, laplacian_matrix)
  Norm_Lap<-lapply(Graphs, normalized_laplacian)
  Fiedler.value=NULL
  norm.fiedler.value=NULL
  
  for (i in 1:length(Graphs)) {
    if (is.null(Graphs[[i]]$type)) { Graphs[[i]]$type = "untyped" }
    df$modularity[i] <- modularity(communities[[i]])
    df$spectral_radius[i] <- eigen(Adj[[i]], symmetric=TRUE, only.values=TRUE)$values[1]
    
    Fiedler.value[[i]]=eigen(L[[i]], symmetric=TRUE, only.values=TRUE)$values
    
    df$FiedlerValue[i] <- Fiedler.value[[i]][length(Fiedler.value[[i]])-1]
    
    norm.fiedler.value[[i]]=eigen(Norm_Lap[[i]], symmetric=TRUE, only.values=TRUE)$values
    
    df$Normalized_FiedlerValue[i] <- norm.fiedler.value[[i]][length(norm.fiedler.value[[i]])-1]
    
    df$eigen_centr[i] <- centr_eigen(Graphs[[i]])$centralization
    df$deg_centr[i] <- centr_degree(Graphs[[i]])$centralization
    
    df$betw_centr[i] <- centr_betw(Graphs[[i]])$centralization
    
    df$max_component[i] <- max(components(Graphs[[i]])$csize)
    df$mean_eccentr[i]<-mean(eccentricity(Graphs[[i]]))
    df$radius[i]<-radius(Graphs[[i]])
    df$mean_path_length[i]<-average.path.length(Graphs[[i]])
    #df$trace[i]<-sum(diag(Adj[[i]]))
    df$graph_energy[i]<-sum(abs(eigen(Adj[[i]], symmetric=TRUE, only.values=TRUE)$values))
    df$min_triangle[i]= min(count_triangles(Graphs[[i]]))
    df$mean_triangle[i]= mean(count_triangles(Graphs[[i]]))
    df$sd_triangle[i]= sd(count_triangles(Graphs[[i]]))
    df$max_triangle[i]= max(count_triangles(Graphs[[i]]))
    df$num_triangle[i]= sum(count_triangles(Graphs[[i]]))
    df$deg_assort_coef[i]=assortativity_degree(Graphs[[i]])
    
    df$threshold[i] <- 1/(df$spectral_radius[i])
    
    if (df$connected[i]==TRUE) {
      df$closeness_centr[i] = mean(closeness(Graphs[[i]]))
    } else { # handle the case where G isn't connected
      df$closeness_centr[i] = -1
    }
  }
  return (df)
}

#calcGraphFeatures(x)
##----------- Graph Features----------#####
## Used to perform runs on multiple simulated graphs on any given network
RunSimOnGraphFeatures<-function(Graphs, nreps=nreps,output_file=NULL, seed=-1) {
  set.seed(1)
  # ### Definition and initialization of parameters for graphfeatures
  graphProperties=list()
  # ### Definition and initialization of Graph Prefix
  graphid=list();graphreplicate=list(); graphname=list(); GraphPrefix=list(); analysis=list()
  for (g in 1:length(Graphs)){
    for (reps in 1:nreps) {
      ### Calculate the graph features for each simulated graph of all the synthetic networks
      print(paste("Calculating graph features on", Graphs[[g]]$name))
      #graphProperties[[reps]] <- calcGraphFeatures(Graphs[g])
      graphProperties[[reps]] <- calcGraphFeatures(Graphs[g])
    }
    graphname[[g]]=Graphs[[g]]$type
    graphid[[g]]=Graphs[[g]]$id
    graphreplicate[[g]]=c(1:nreps)
    GraphPrefix=cbind(graphname[[g]],graphid[[g]],graphreplicate[[g]])
    colnames(GraphPrefix)=c("GraphName","GraphID","GraphReplicate")
    analysis[[g]]=as.data.frame(cbind(GraphPrefix,graphProperties[[reps]]))
    row.names(analysis[[g]])=1:nreps
  }
  All_results=do.call(rbind,analysis)
 # write.csv(All_results, file=output_file)
  return( All_results)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [1] Spatial Expander Propagation Graph with weighted edges, 
# node attributes, scale free degree distribution, small world effect and community structures
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#library(quadtree)
# Function to generate a spatial scale-free expander graph
# with the given parameters
#
# Arguments:
#   N: number of nodes
#   r: cutoff distance for adjacency matrix
#   beta: power law parameter for distance-dependent probability
#   m: number of edges added for each new node
#   alpha: parameter controlling strength of degree effect
#   mu: parameter controlling preference when connecting nodes. Thus, preference to care about 
#       spatial distance between nodes when connecting them, vs caring only about vertex degree. 0:care aboutdegree
#1: care about distance
#   prewire: rewiring probability for small world effect
#   node_attrs: list of node attributes to add to graph
#   edge_weights: logical indicating whether to add edge weights
#
# Returns:
#   A graph object of class igraph
spatial_scale_free_expander_noncomm <- function(N=50, beta=2, m=2, prewire = 0.1,
                                                mu=0.2,add_edge_weight = FALSE,
                                                add_node_attr = FALSE,r=0.1,alpha=0.2) {
  
  
  ##--Generate spatial points on a unit square
  points <- matrix(runif(N*2), ncol=2)
  # calculate distance matrix between all pairs of points
  dist_mat <- as.matrix(dist(points))
  # Set up graph ogject
  graph <- list()
  graph$nodes <- 1:m
  # initialize the edge_list and adjacency matrix
  graph$edge_list <- data.frame(from=NULL, to=NULL)
  graph$adj_mat <- matrix(0, nrow = N, ncol = N)
  # initialize graph with an initial node and no edges
  graph$adj_mat[1]<-1
  # initialize the degree vector
  #graph$degrees <- numeric(N)
  graph$degrees=rep(m,N)
  
  
  # Grow the network with new nodes added one at a time
  for (i in 3:N) {
    # preferential attachment effect with scale free technique
    deg_probs<- graph$degrees[1:(i-1)]^alpha / sum(graph$degrees[1:(i-1)]^alpha)
    #--spatial distance effect incorporating short spatial distances (cut-off distance)
    spat_dist_effect <- ifelse(dist_mat[1:(i-1), i] <= r, 1/(1+(dist_mat[1:(i-1), i])^beta), 1/N)
    spatial_probs <- spat_dist_effect /sum(spat_dist_effect )#normalise
    # total probability of attachment which is normalised
    TotalprobsOfAttach <-(mu*spatial_probs)+((1-mu)*deg_probs)
    #Add m edges to existing nodes with probability proportional to their degree and distance
    targetnodes<- sample(1:(i-1), size = m, replace = TRUE, prob =  TotalprobsOfAttach)
    
    # Attach nodes
    graph$nodes <- c(graph$nodes, i)
    graph$adj_mat[i,  targetnodes] <- 1
    graph$adj_mat[targetnodes, i] <- 1
    
    # Update degree vectors
    graph$degrees[targetnodes] <- graph$degrees[targetnodes]+1
    graph$degrees[i] <- m 
    
    #update edge_list
    new_edges <- data.frame(from=i, to=targetnodes)
    graph$edge_list<- rbind(graph$edge_list, new_edges)
    
    # Small-world rewiring with probability p
    if (runif(1) < prewire) {
      # Select random node to rewire within cutoff distance
      neighbors <- which(graph$adj_mat[i,] == 1)
      if (length(neighbors) > 0) {
        new_neighbor <- sample(neighbors, 1)
        graph$adj_mat[i, new_neighbor] <- 0
        graph$adj_mat[new_neighbor, i] <- 0
        available_nodes <- setdiff(1:N, c(i, neighbors))
        new_neighbor <- sample(available_nodes, 1)
        graph$adj_mat[i, new_neighbor] <- 1
        graph$adj_mat[new_neighbor, i] <- 1
      }
    }
  }
  
  ### Graph object
  GraphModel <- graph.adjacency(as.matrix(graph$adj_mat), mode="undirected")
  GraphModel=igraph::simplify(GraphModel,remove.loops = TRUE,remove.multiple = TRUE)
  graph$GraphObject <- GraphModel
  
  # Compute graph Laplacian for the expansion property
  D <- diag(rowSums(graph$adj_mat))
  L <- D - graph$adj_mat
  # Compute eigenvalues and eigenvectors of Laplacian matrix
  eig <- eigen(Matrix(L))
  eig_values <- Re(eig$values)
  eig_vectors <- eig$vectors
  
  # Compute the expansion ratio
  lambda_2 <- eig_values[2]
  lambda_n <- eig_values[nrow(L)]
  expansion_ratio <- lambda_2 / lambda_n
  # Check for strong expansion properties
  #spectral_gap=eigen(Matrix::t(L) %*% L, symmetric=TRUE, only.values=TRUE)$values[n-1]
  if(!is.connected(GraphModel)||expansion_ratio <= 0){
    warning("Graph may not exhibit strong expansion properties or Graph is disconnected")
  }
  
  graph$ExpansionRatio=expansion_ratio
  
  if (add_edge_weight) {
    graph$adj_mat[graph$adj_mat == 1] <- runif(sum(graph$adj_mat == 1))
    graph$Weights=graph$adj_mat
  }
  
  # Create node attributes if requested
  if (add_node_attr) {
    node_attrs <- data.frame(x = points[,1], y = points[,2])
    graph$NodeAttributes=node_attrs
  }
  return(graph)
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2] Spatial Expander Propagation Graph with weighted edges, 
# node attributes, scale free degree distribution, small world effect and community structures
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++ Code in R to grow a spatial scale-free expander graph with the following properties. 
#+This code should be built from scratch and not contain any igraph packages etc.

# 1) Nodes should be generated on a unit square where the initial set of nodes are distributed spatially with an underlying spatial structure, e.g. by not permitting nodes to land too near each other, or not too far from each other using an in built 'quadtree' algorithm built from scratch
# 2) Use the distance matrix to create an adjacency matrix for the graph.o create an adjacency matrix based on the distance matrix, you can define a cutoff distance 'r' and set the adjacency matrix element A_{ij} to 1 if dist_{ij} <= 'r' and 0 otherwise . This parameter 'r' controls the strength of the spatial distance effect
# 3) This network should be grown and nodes should be organized into 
#non-overlapping communities ( or blocks) with nodes connecting with each 
#other within the same community or between communities with an attachment probability 
#favouring short spatial distances and higher degree nodes that follows scale-free sdegree distribution. 
#Include 'alpha': the parameter that controls the strength of the degree effect
# 4)Add a rewiring probability 'prewire' for small world effect to the resulting graph to enhance the connectivity to long range nodes. Thus you can randomly rewire each edge with probability 'prewire' to a random node within a certain distance
# 5) Add arguments to create edge weights, node attributes etc
# 6)Check the expansion properties of the graph by computing the Laplacian matrix. The Laplacian matrix represents the graph's connectivity and its ability to spread information. A good spatial graph should have a small eigenvalue gap, indicating strong expansion properties.
# 7)This graph should be a single function where we can tun', N, 'alpha' , 'prewire', 'm', to generate different graph models

# count nodes for each clusters
#NumOfNodesPerCluster <- table(graph$clusters)
# num_nodes_per_comm <- rep(1, numofComm)


spatial_scale_free_expander_comm <- function(N=50, beta=2, m=1, prewire = 0.1, numofComm=4, 
                                             mu=0.2,add_edge_weight = FALSE,
                                             Prob_Betw_comm=0.1,add_node_attr = FALSE,r=0.1,alpha=0.2) {
  
  
  
  ##--Generate spatial points on a unit square
  points <- matrix(runif(N*2), ncol=2)
  # calculate distance matrix between all pairs of points
  dist_mat <- as.matrix(dist(points))
  # Set up graph ogject
  graph <- list()
  graph$nodes <- 1:m
  # initialize the edge_list and adjacency matrix
  graph$edge_list <- data.frame(from=NULL, to=NULL)
  graph$adj_mat <- matrix(0, nrow = N, ncol = N)
  # initialize graph with an initial node and no edges
  graph$adj_mat[1]<-1
  # initialize the degree vector
  #graph$degrees <- numeric(N)
  graph$degrees=rep(m,N)
  
  # initialize the edge_list and community membership
  graph$edge_list <- as.data.frame(matrix(ncol=2, nrow=0))
  
  # create initial community for first node. Assin first node to the first community
  graph$clusters <- rep(1, N)
  
  # Step 2: Create initial community/block probability matrix
  P_ij <- matrix(0, nrow = numofComm, ncol = numofComm)
  
  Prob_Within_comm=1-Prob_Betw_comm
  
  if(is.numeric(Prob_Betw_comm)){
    diag(P_ij) <- Prob_Within_comm
    P_ij[lower.tri(P_ij)] <- Prob_Betw_comm
    P_ij[upper.tri(P_ij)] <- t(P_ij)[upper.tri(P_ij)]
    
  }else {
    P_ij <- matrix(runif(numofComm^2), numofComm, numofComm)
  }
  
  # Step 3: Add nodes to the network one by one
  # Grow the network by adding new nodes and edges
  for (i in 3:N) {
    # Step 3i: Create arbitrary number of communities/blocks
    comm_probs<-rep(1/numofComm, numofComm)#sample.int(numofComm, (i-1), replace = TRUE)
    graph$clusters[i] <- sample(numofComm, 1, prob = comm_probs)
    
    # Step 3ii: Define community probability matrix
    P_comm <- P_ij[1:numofComm, 1:numofComm]
    
    ##--Step 3iii: Define probability of attachment function for within and between-community connections
    # P_unique_within <- rep(0, (i-1))
    # P_unique_between <- rep(0, (i-1))
    
    spatial_probs <- ifelse(dist_mat[1:(i-1), i] <= r, 1/(1+(dist_mat[1:(i-1), i])^beta), 1/N)
    # graph$degrees<- rowSums(graph$adj_mat[1:(i-1), ])
    deg_probs <- (graph$degrees[1:(i-1)]^alpha) / sum( graph$degrees[1:(i-1)]^alpha)
    
    # Step 3iv: Add edges within and between communities based on community and attachment probabilty matrix 
    P_noncomm <- (mu*spatial_probs)+((1-mu)*deg_probs) 
    P_noncomm<-P_noncomm/sum(P_noncomm)
    if (numofComm==1){
      P_join <-P_comm[graph$clusters[(i-1)]] * P_noncomm
    }else{
      P_join <-P_comm[graph$clusters[1:(i-1)], graph$clusters[i]] * P_noncomm
    }
    
    #Connect nodes based on attachment probability
    targetnodes <- sample(1:(i-1), m, replace = TRUE, prob = P_join)
    graph$adj_mat[i,targetnodes] <- 1
    graph$adj_mat[targetnodes, i] <- 1
    graph$degrees[targetnodes] <- graph$degrees[targetnodes] + 1
    graph$degrees[i] <- m
    
    #update edge_list
    new_edges <- data.frame(from=i, to=targetnodes)
    graph$edge_list<- rbind(graph$edge_list, new_edges)
    
    # Step 4: Small-world rewiring with probability p
    if (runif(1) < prewire) {
      # Select random node to rewire within cutoff distance
      neighbors <- which(graph$adj_mat[i,] == 1)
      if (length(neighbors) > 0) {
        new_neighbor <- sample(neighbors, 1)
        graph$adj_mat[i, new_neighbor] <- 0
        graph$adj_mat[new_neighbor, i] <- 0
        available_nodes <- setdiff(1:N, c(i, neighbors))
        new_neighbor <- sample(available_nodes, 1)
        graph$adj_mat[i, new_neighbor] <- 1
        graph$adj_mat[new_neighbor, i] <- 1
      }
    }
  }
  
  ### Graph object
  GraphModel <- graph.adjacency(as.matrix(graph$adj_mat), mode="undirected")
  GraphModel<-igraph::simplify(GraphModel,remove.loops = TRUE,remove.multiple = TRUE)
  graph$GraphObject <- GraphModel
  
  # Compute graph Laplacian for the expansion property
  D <- diag(rowSums(graph$adj_mat))
  L <- D - graph$adj_mat
  # Compute eigenvalues and eigenvectors of Laplacian matrix
  eig <- eigen(Matrix(L))
  eig_values <- Re(eig$values)
  eig_vectors <- eig$vectors
  
  # Compute the expansion ratio
  lambda_2 <- eig_values[2]
  lambda_n <- eig_values[nrow(L)]
  expansion_ratio <- lambda_2 / lambda_n
  # Check for strong expansion properties
  #spectral_gap=eigen(Matrix::t(L) %*% L, symmetric=TRUE, only.values=TRUE)$values[n-1]
  if(!is.connected(GraphModel)||expansion_ratio <= 0){  
    warning("Graph may not exhibit strong expansion properties or Graph is disconnected")
  }
  
  graph$ExpansionRatio=expansion_ratio
  
  if (add_edge_weight) {
    graph$adj_mat[graph$adj_mat == 1] <- runif(sum(graph$adj_mat == 1))
    graph$Weights=graph$adj_mat
  }
  
  # Create node attributes if requested
  if (add_node_attr) {
    node_attrs <- data.frame(x =   points[,1], y =   points[,2])
    graph$NodeAttributes=node_attrs
  } 
  return(graph)
}      


##########################################################################################################
# Create a list of graphs for simulation experiment
##########################################################################################################
GenSpatialGraphNonComm<- function(N=200, beta=2, m=10, prewire = 0.01, 
                                  mu=0.2,add_edge_weight = FALSE,
                                  add_node_attr = T,r=0.1,alpha=0.2,nSamples=2) {
  Graphs = list()
  print("Creating Spatial scalefree networks without community")
  for (i in 1:nSamples) {
    Graphs[[i]] = spatial_scale_free_expander_noncomm (N=N, beta=beta, m=m, prewire = prewire, 
                                          mu=mu,add_edge_weight = add_edge_weight, 
                                          add_node_attr = add_node_attr,r=r,alpha=alpha)$GraphObject
    Graphs[[i]]$type = "Spatial_noncomm"
    Graphs[[i]]$name="Spatial_noncomm"
    Graphs[[i]]$id <- i
  }
  return(Graphs)
}

G=GenSpatialGraphNonComm(N=200,beta=2, m=10, prewire = 0.01, 
                         mu=0.2,add_edge_weight = FALSE,
                         add_node_attr = T,r=0.1,alpha=0.2,nSamples=1) 

#G

Sim_SpatialGraphNonComm=function(mu=c(0.3,0.4,0.6,0.05),N=100){
  net=NULL;data=NULL
  k=1
  for (i in 1:length(mu)){
    net[i]= GenSpatialGraphNonComm(N=N, beta=2, m=10, prewire = 0.1, 
                  mu=mu[i],add_edge_weight = FALSE,add_node_attr = T,r=Inf,alpha=0.2,nSamples=1)
    data=RunSimOnGraphFeatures(net,nreps = 1)
    k=k+1
  }
  df=cbind(mu,data)
  return(df)
}

### 50 nodes
mu0.50=Sim_SpatialGraphNonComm(mu=rep(0,1000),N=50)
write.csv(mu0.50,"mu0.50.csv")
 mu1.50=Sim_SpatialGraphNonComm(mu=rep(0.15,1000),N=50)
 write.csv(mu1.50,"mu1.50.csv")
 mu2.50=Sim_SpatialGraphNonComm(mu=rep(0.2,1000),N=50)
 write.csv(mu2.50,"mu2.50.csv")
 mu3.50=Sim_SpatialGraphNonComm(mu=rep(0.3,1000),N=50)
 write.csv(mu3.50,"mu3.50.csv")
 mu4.50=Sim_SpatialGraphNonComm(mu=rep(0.4,1000),N=50)
 write.csv(mu4.50,"mu4.50.csv")
 mu5.50=Sim_SpatialGraphNonComm(mu=rep(0.5,1000),N=50)
 write.csv(mu5.50,"mu5.50.csv")
 mu6.50=Sim_SpatialGraphNonComm(mu=rep(0.6,1000),N=50)
 write.csv(mu6.50,"mu6.50.csv")
 mu7.50=Sim_SpatialGraphNonComm(mu=rep(0.7,1000),N=50)
 write.csv(mu7.50,"mu7.50.csv")
 mu8.50=Sim_SpatialGraphNonComm(mu=rep(0.8,1000),N=50)
 write.csv(mu8.50,"mu8.50.csv")
 mu9.50=Sim_SpatialGraphNonComm(mu=rep(0.9,1000),N=50)
 write.csv(mu9.50,"mu9.50.csv")
 mu10.50=Sim_SpatialGraphNonComm(mu=rep(1,1000),N=50)
 write.csv(mu10.50,"mu10.50.csv")
 
 mu0_50=read.csv("mu0.50.csv");mu1_50=read.csv("mu1.50.csv");mu2_50=read.csv("mu2.50.csv");
 mu3_50=read.csv("mu3.50.csv");mu4_50=read.csv("mu4.50.csv");mu5_50=read.csv("mu5.50.csv");
 mu6_50=read.csv("mu6.50.csv");mu7_50=read.csv("mu7.50.csv");mu8_50=read.csv("mu8.50.csv");
 mu9_50=read.csv("mu9.50.csv");mu10_50=read.csv("mu10.50.csv")

 df_mu_50=rbind(mu0_50,mu1_50,mu2_50,mu3_50,mu4_50,mu5_50,mu6_50,
                mu7_50,mu8_50,mu9_50,mu10_50)
 write.csv( df_mu_50," df_mu_50.csv")
 dim(df_mu_50)
 node.num=50
 df_mu_50<-df_mu_50 %>%dplyr::filter(order==50)
 dim(df_mu_50)
write.csv(df_mu_50,"df_mu_50.csv")

### 100 nodes 

### 100 nodes
mu0.100=Sim_SpatialGraphNonComm(mu=rep(0,1000),N=100)
write.csv(mu0.100,"mu0.100.csv")
mu1.100=Sim_SpatialGraphNonComm(mu=rep(0.15,1000),N=100)
write.csv(mu1.100,"mu1.100.csv")
mu2.100=Sim_SpatialGraphNonComm(mu=rep(0.2,1000),N=100)
write.csv(mu2.100,"mu2.100.csv")
mu3.100=Sim_SpatialGraphNonComm(mu=rep(0.3,1000),N=100)
write.csv(mu3.100,"mu3.100.csv")
mu4.100=Sim_SpatialGraphNonComm(mu=rep(0.4,1000),N=100)
write.csv(mu4.100,"mu4.100.csv")
mu5.100=Sim_SpatialGraphNonComm(mu=rep(0.5,1000),N=100)
write.csv(mu5.100,"mu5.100.csv")
mu6.100=Sim_SpatialGraphNonComm(mu=rep(0.6,1000),N=100)
write.csv(mu6.100,"mu6.100.csv")
mu7.100=Sim_SpatialGraphNonComm(mu=rep(0.7,1000),N=100)
write.csv(mu7.100,"mu7.100.csv")
mu8.100=Sim_SpatialGraphNonComm(mu=rep(0.8,1000),N=100)
write.csv(mu8.100,"mu8.100.csv")
mu9.100=Sim_SpatialGraphNonComm(mu=rep(0.9,1000),N=100)
write.csv(mu9.100,"mu9.100.csv")
mu10.100=Sim_SpatialGraphNonComm(mu=rep(1,1000),N=100)
write.csv(mu10.100,"mu10.100.csv")

mu0_100=read.csv("mu0.100.csv");mu1_100=read.csv("mu1.100.csv");mu2_100=read.csv("mu2.100.csv");
mu3_100=read.csv("mu3.100.csv");mu4_100=read.csv("mu4.100.csv");mu5_100=read.csv("mu5.100.csv");
mu6_100=read.csv("mu6.100.csv");mu7_100=read.csv("mu7.100.csv");mu8_100=read.csv("mu8.100.csv");
mu9_100=read.csv("mu9.100.csv");mu10_100=read.csv("mu10.100.csv")

df_mu_100=rbind(mu0_100,mu1_100,mu2_100,mu3_100,mu4_100,mu5_100,mu6_100,
                mu7_100,mu8_100,mu9_100,mu10_100)
dim(df_mu_100)
node.num=100
df_mu_100<-df_mu_100 %>%dplyr::filter(order==100)
dim(df_mu_100)
write.csv(df_mu_100,"df_mu_100.csv")

### 150 nodes
mu0.150=Sim_SpatialGraphNonComm(mu=rep(0,1000),N=150)
write.csv(mu0.150,"mu0.150.csv")
mu1.150=Sim_SpatialGraphNonComm(mu=rep(0.15,1000),N=150)
write.csv(mu1.150,"mu1.150.csv")
mu2.150=Sim_SpatialGraphNonComm(mu=rep(0.2,1000),N=150)
write.csv(mu2.150,"mu2.150.csv")
mu3.150=Sim_SpatialGraphNonComm(mu=rep(0.3,1000),N=150)
write.csv(mu3.150,"mu3.150.csv")
mu4.150=Sim_SpatialGraphNonComm(mu=rep(0.4,1000),N=150)
write.csv(mu4.150,"mu4.150.csv")
mu5.150=Sim_SpatialGraphNonComm(mu=rep(0.5,1000),N=150)
write.csv(mu5.150,"mu5.150.csv")
mu6.150=Sim_SpatialGraphNonComm(mu=rep(0.6,1000),N=150)
write.csv(mu6.150,"mu6.150.csv")
mu7.150=Sim_SpatialGraphNonComm(mu=rep(0.7,1000),N=150)
write.csv(mu7.150,"mu7.150.csv")
mu8.150=Sim_SpatialGraphNonComm(mu=rep(0.8,1000),N=150)
write.csv(mu8.150,"mu8.150.csv")
mu9.150=Sim_SpatialGraphNonComm(mu=rep(0.9,1000),N=150)
write.csv(mu9.150,"mu9.150.csv")
mu10.150=Sim_SpatialGraphNonComm(mu=rep(1,1000),N=150)
write.csv(mu10.150,"mu10.150.csv")

mu0_150=read.csv("mu0.150.csv");mu1_150=read.csv("mu1.150.csv");mu2_150=read.csv("mu2.150.csv");
mu3_150=read.csv("mu3.150.csv");mu4_150=read.csv("mu4.150.csv");mu5_150=read.csv("mu5.150.csv");
mu6_150=read.csv("mu6.150.csv");mu7_150=read.csv("mu7.150.csv");mu8_150=read.csv("mu8.150.csv");
mu9_150=read.csv("mu9.150.csv");mu10_150=read.csv("mu10.150.csv")

df_mu_150=rbind(mu0_150,mu1_150,mu2_150,mu3_150,mu4_150,mu5_150,mu6_150,
                mu7_150,mu8_150,mu9_150,mu10_150)
dim(df_mu_150)
node.num=150
df_mu_150<-df_mu_150 %>%dplyr::filter(order==150)
dim(df_mu_150)
write.csv(df_mu_150,"df_mu_150.csv")


###-----200 nodes---#####
mu0.200=Sim_SpatialGraphNonComm(mu=rep(0,1000),N=200)
write.csv(mu0.200,"mu0.200.csv")
mu1.200=Sim_SpatialGraphNonComm(mu=rep(0.15,1000),N=200)
write.csv(mu1.200,"mu1.200.csv")
mu2.200=Sim_SpatialGraphNonComm(mu=rep(0.2,1000),N=200)
write.csv(mu2.200,"mu2.200.csv")
mu3.200=Sim_SpatialGraphNonComm(mu=rep(0.3,1000),N=200)
write.csv(mu3.200,"mu3.200.csv")
mu4.200=Sim_SpatialGraphNonComm(mu=rep(0.4,1000),N=200)
write.csv(mu4.200,"mu4.200.csv")
mu5.200=Sim_SpatialGraphNonComm(mu=rep(0.5,1000),N=200)
write.csv(mu5.200,"mu5.200.csv")
mu6.200=Sim_SpatialGraphNonComm(mu=rep(0.6,1000),N=200)
write.csv(mu6.200,"mu6.200.csv")
mu7.200=Sim_SpatialGraphNonComm(mu=rep(0.7,1000),N=200)
write.csv(mu7.200,"mu7.200.csv")
mu8.200=Sim_SpatialGraphNonComm(mu=rep(0.8,1000),N=200)
write.csv(mu8.200,"mu8.200.csv")
mu9.200=Sim_SpatialGraphNonComm(mu=rep(0.9,1000),N=200)
write.csv(mu9.200,"mu9.200.csv")
mu10.200=Sim_SpatialGraphNonComm(mu=rep(1,1000),N=200)
write.csv(mu10.200,"mu10.200.csv")

mu0_200=read.csv("mu0.200.csv");mu1_200=read.csv("mu1.200.csv");mu2_200=read.csv("mu2.200.csv");
mu3_200=read.csv("mu3.200.csv");mu4_200=read.csv("mu4.200.csv");mu5_200=read.csv("mu5.200.csv");
mu6_200=read.csv("mu6.200.csv");mu7_200=read.csv("mu7.200.csv");mu8_200=read.csv("mu8.200.csv");
mu9_200=read.csv("mu9.200.csv");mu10_200=read.csv("mu10.200.csv")

df_mu_200=rbind(mu0_200,mu1_200,mu2_200,mu3_200,mu4_200,mu5_200,mu6_200,
                mu7_200,mu8_200,mu9_200,mu10_200)
dim(df_mu_200)
node.num=200
df_mu_200<-df_mu_200 %>%dplyr::filter(order==200)
dim(df_mu_200)
write.csv(df_mu_200,"df_mu_200.csv")


###---250 nodes---#### 

### 250 nodes
mu0.250=Sim_SpatialGraphNonComm(mu=rep(0,1000),N=250)
write.csv(mu0.250,"mu0.250.csv")
mu1.250=Sim_SpatialGraphNonComm(mu=rep(0.15,1000),N=250)
write.csv(mu1.250,"mu1.250.csv")
mu2.250=Sim_SpatialGraphNonComm(mu=rep(0.2,1000),N=250)
write.csv(mu2.250,"mu2.250.csv")
mu3.250=Sim_SpatialGraphNonComm(mu=rep(0.3,1000),N=250)
write.csv(mu3.250,"mu3.250.csv")
mu4.250=Sim_SpatialGraphNonComm(mu=rep(0.4,1000),N=250)
write.csv(mu4.250,"mu4.250.csv")
mu5.250=Sim_SpatialGraphNonComm(mu=rep(0.5,1000),N=250)
write.csv(mu5.250,"mu5.250.csv")
mu6.250=Sim_SpatialGraphNonComm(mu=rep(0.6,1000),N=250)
write.csv(mu6.250,"mu6.250.csv")
mu7.250=Sim_SpatialGraphNonComm(mu=rep(0.7,1000),N=250)
write.csv(mu7.250,"mu7.250.csv")
mu8.250=Sim_SpatialGraphNonComm(mu=rep(0.8,1000),N=250)
write.csv(mu8.250,"mu8.250.csv")
mu9.250=Sim_SpatialGraphNonComm(mu=rep(0.9,1000),N=250)
write.csv(mu9.250,"mu9.250.csv")
mu10.250=Sim_SpatialGraphNonComm(mu=rep(1,1000),N=250)
write.csv(mu10.250,"mu10.250.csv")

mu0_250=read.csv("mu0.250.csv");mu1_250=read.csv("mu1.250.csv");mu2_250=read.csv("mu2.250.csv");
mu3_250=read.csv("mu3.250.csv");mu4_250=read.csv("mu4.250.csv");mu5_250=read.csv("mu5.250.csv");
mu6_250=read.csv("mu6.250.csv");mu7_250=read.csv("mu7.250.csv");mu8_250=read.csv("mu8.250.csv");
mu9_250=read.csv("mu9.250.csv");mu10_250=read.csv("mu10.250.csv")

df_mu_250=rbind(mu0_250,mu1_250,mu2_250,mu3_250,mu4_250,mu5_250,mu6_250,
                mu7_250,mu8_250,mu9_250,mu10_250)
dim(df_mu_250)
node.num=250
df_mu_250<-df_mu_250 %>%dplyr::filter(order==250)
dim(df_mu_250)
write.csv(df_mu_250,"df_mu_250.csv")


###-----300 nodes----#### 
mu1.300=Sim_SpatialGraphNonComm(mu=rep(0.15,200),N=300)
write.csv(mu1.300,"mu1.300.csv")

mu2.300=Sim_SpatialGraphNonComm(mu=rep(0.2,200),N=300)
write.csv(mu2.300,"mu2.300.csv")

mu3.300=Sim_SpatialGraphNonComm(mu=rep(0.3,200),N=300)
write.csv(mu3.300,"mu3.300.csv")

mu4.300=Sim_SpatialGraphNonComm(mu=rep(0.4,200),N=300)
write.csv(mu4.300,"mu4.300.csv")

mu5.300=Sim_SpatialGraphNonComm(mu=rep(0.5,200),N=300)
write.csv(mu5.300,"mu5.300.csv")

mu6.300=Sim_SpatialGraphNonComm(mu=rep(0.6,200),N=300)
write.csv(mu6.300,"mu6.300.csv")

mu7.300=Sim_SpatialGraphNonComm(mu=rep(0.7,200),N=300)
write.csv(mu7.300,"mu7.300.csv")

mu8.300=Sim_SpatialGraphNonComm(mu=rep(0.8,200),N=300)
write.csv(mu8.300,"mu8.300.csv")

mu9.300=Sim_SpatialGraphNonComm(mu=rep(0.9,200),N=300)
write.csv(mu9.300,"mu9.300.csv")

mu10.300=Sim_SpatialGraphNonComm(mu=rep(1,200),N=300)
write.csv(mu10.300,"mu10.300.csv")

######## Saving all 300 nodes data to csv
# mu1=read.csv("mu1.300.csv");mu2=read.csv("mu2.300.csv");mu3=read.csv("mu3.300.csv")
# mu4=read.csv("mu4.300.csv");mu5=read.csv("mu5.300.csv");mu6=read.csv("mu6.300.csv")
# mu7=read.csv("mu7.300.csv");mu8=read.csv("mu8.300.csv");mu9=read.csv("mu9.300.csv")
# mu10=read.csv("mu10.300.csv")

df_mu_300=rbind(mu1.300,mu2.300,mu3.300,mu4.300,mu5.300,mu6.300,
                mu7.300,mu8.300,mu9.300,mu10.300)
dim(df_mu_300)
node.num=300
df_mu_300<-df_mu_300 %>%dplyr::filter(order==300)
dim(df_mu_300)
write.csv(df_mu_300,"df_mu_300.csv")

### 350 nodes 
mu1.350=Sim_SpatialGraphNonComm(mu=rep(0.15,200),N=350)
write.csv(mu1.350,"mu1.350.csv")

mu2.350=Sim_SpatialGraphNonComm(mu=rep(0.2,200),N=350)
write.csv(mu2.350,"mu2.350.csv")

mu3.350=Sim_SpatialGraphNonComm(mu=rep(0.3,200),N=350)
write.csv(mu3.350,"mu3.350.csv")

mu4.350=Sim_SpatialGraphNonComm(mu=rep(0.4,200),N=350)
write.csv(mu4.350,"mu4.350.csv")

mu5.350=Sim_SpatialGraphNonComm(mu=rep(0.5,200),N=350)
write.csv(mu5.350,"mu5.350.csv")

mu6.350=Sim_SpatialGraphNonComm(mu=rep(0.6,200),N=350)
write.csv(mu6.350,"mu6.350.csv")

mu7.350=Sim_SpatialGraphNonComm(mu=rep(0.7,200),N=350)
write.csv(mu7.350,"mu7.350.csv")

mu8.350=Sim_SpatialGraphNonComm(mu=rep(0.8,200),N=350)
write.csv(mu8.350,"mu8.350.csv")

mu9.350=Sim_SpatialGraphNonComm(mu=rep(0.9,200),N=350)
write.csv(mu9.350,"mu9.350.csv")

mu10.350=Sim_SpatialGraphNonComm(mu=rep(1,200),N=350)
write.csv(mu10.350,"mu10.350.csv")

######## Saving all 350 nodes data to csv
# mu1=read.csv("mu1.350.csv");mu2=read.csv("mu2.350.csv");mu3=read.csv("mu3.350.csv")
# mu4=read.csv("mu4.350.csv");mu5=read.csv("mu5.350.csv");mu6=read.csv("mu6.350.csv")
# mu7=read.csv("mu7.350.csv");mu8=read.csv("mu8.350.csv");mu9=read.csv("mu9.350.csv")
# mu10=read.csv("mu10.350.csv")

df_mu_350=rbind(mu1.350,mu2.350,mu3.350,mu4.350,mu5.350,mu6.350,
                mu7.350,mu8.350,mu9.350,mu10.350)
dim(df_mu_350)
node.num=350
df_mu_350<-df_mu_350 %>%dplyr::filter(order==350)
dim(df_mu_350)
write.csv(df_mu_350,"df_mu_350.csv")

### 400 nodes 
mu1.400=Sim_SpatialGraphNonComm(mu=rep(0.15,200),N=400)
write.csv(mu1.400,"mu1.400.csv")

mu2.400=Sim_SpatialGraphNonComm(mu=rep(0.2,200),N=400)
write.csv(mu2.400,"mu2.400.csv")

mu3.400=Sim_SpatialGraphNonComm(mu=rep(0.3,200),N=400)
write.csv(mu3.400,"mu3.400.csv")

mu4.400=Sim_SpatialGraphNonComm(mu=rep(0.4,200),N=400)
write.csv(mu4.400,"mu4.400.csv")

mu5.400=Sim_SpatialGraphNonComm(mu=rep(0.5,200),N=400)
write.csv(mu5.400,"mu5.400.csv")

mu6.400=Sim_SpatialGraphNonComm(mu=rep(0.6,200),N=400)
write.csv(mu6.400,"mu6.400.csv")

mu7.400=Sim_SpatialGraphNonComm(mu=rep(0.7,200),N=400)
write.csv(mu7.400,"mu7.400.csv")

mu8.400=Sim_SpatialGraphNonComm(mu=rep(0.8,200),N=400)
write.csv(mu8.400,"mu8.400.csv")

mu9.400=Sim_SpatialGraphNonComm(mu=rep(0.9,200),N=400)
write.csv(mu9.400,"mu9.400.csv")

mu10.400=Sim_SpatialGraphNonComm(mu=rep(1,200),N=400)
write.csv(mu10.400,"mu10.400.csv")

######## Saving all 400 nodes data to csv
# mu1=read.csv("mu1.400.csv");mu2=read.csv("mu2.400.csv");mu3=read.csv("mu3.400.csv")
# mu4=read.csv("mu4.400.csv");mu5=read.csv("mu5.400.csv");mu6=read.csv("mu6.400.csv")
# mu7=read.csv("mu7.400.csv");mu8=read.csv("mu8.400.csv");mu9=read.csv("mu9.400.csv")
# mu10=read.csv("mu10.400.csv")

df_mu_400=rbind(mu1.400,mu2.400,mu3.400,mu4.400,mu5.400,mu6.400,
                mu7.400,mu8.400,mu9.400,mu10.400)
dim(df_mu_400)
node.num=400
df_mu_400<-df_mu_400 %>%dplyr::filter(order==400)
dim(df_mu_400)
write.csv(df_mu_400,"df_mu_400.csv")

### 450 nodes 
mu1.450=Sim_SpatialGraphNonComm(mu=rep(0.15,200),N=450)
write.csv(mu1.450,"mu1.450.csv")

mu2.450=Sim_SpatialGraphNonComm(mu=rep(0.2,200),N=450)
write.csv(mu2.450,"mu2.450.csv")

mu3.450=Sim_SpatialGraphNonComm(mu=rep(0.3,200),N=450)
write.csv(mu3.450,"mu3.450.csv")

mu4.450=Sim_SpatialGraphNonComm(mu=rep(0.4,200),N=450)
write.csv(mu4.450,"mu4.450.csv")

mu5.450=Sim_SpatialGraphNonComm(mu=rep(0.5,200),N=450)
write.csv(mu5.450,"mu5.450.csv")

mu6.450=Sim_SpatialGraphNonComm(mu=rep(0.6,200),N=450)
write.csv(mu6.450,"mu6.450.csv")

mu7.450=Sim_SpatialGraphNonComm(mu=rep(0.7,200),N=450)
write.csv(mu7.450,"mu7.450.csv")

mu8.450=Sim_SpatialGraphNonComm(mu=rep(0.8,200),N=450)
write.csv(mu8.450,"mu8.450.csv")

mu9.450=Sim_SpatialGraphNonComm(mu=rep(0.9,200),N=450)
write.csv(mu9.450,"mu9.450.csv")

mu10.450=Sim_SpatialGraphNonComm(mu=rep(1,200),N=450)
write.csv(mu10.450,"mu10.450.csv")

######## Saving all 450 nodes data to csv
# mu1=read.csv("mu1.450.csv");mu2=read.csv("mu2.450.csv");mu3=read.csv("mu3.450.csv")
# mu4=read.csv("mu4.450.csv");mu5=read.csv("mu5.450.csv");mu6=read.csv("mu6.450.csv")
# mu7=read.csv("mu7.450.csv");mu8=read.csv("mu8.450.csv");mu9=read.csv("mu9.450.csv")
# mu10=read.csv("mu10.450.csv")

df_mu_450=rbind(mu1.450,mu2.450,mu3.450,mu4.450,mu5.450,mu6.450,
                mu7.450,mu8.450,mu9.450,mu10.450)
dim(df_mu_450)
node.num=450
df_mu_450<-df_mu_450 %>%dplyr::filter(order==450)
dim(df_mu_450)
write.csv(df_mu_450,"df_mu_450.csv")


### 500 nodes 
mu1.500=Sim_SpatialGraphNonComm(mu=rep(0.15,200),N=500)
write.csv(mu1.500,"mu1.500.csv")

mu2.500=Sim_SpatialGraphNonComm(mu=rep(0.2,200),N=500)
write.csv(mu2.500,"mu2.500.csv")

mu3.500=Sim_SpatialGraphNonComm(mu=rep(0.3,200),N=500)
write.csv(mu3.500,"mu3.500.csv")

mu4.500=Sim_SpatialGraphNonComm(mu=rep(0.4,200),N=500)
write.csv(mu4.500,"mu4.500.csv")

mu5.500=Sim_SpatialGraphNonComm(mu=rep(0.5,200),N=500)
write.csv(mu5.500,"mu5.500.csv")

mu6.500=Sim_SpatialGraphNonComm(mu=rep(0.6,200),N=500)
write.csv(mu6.500,"mu6.500.csv")

mu7.500=Sim_SpatialGraphNonComm(mu=rep(0.7,200),N=500)
write.csv(mu7.500,"mu7.500.csv")

mu8.500=Sim_SpatialGraphNonComm(mu=rep(0.8,200),N=500)
write.csv(mu8.500,"mu8.500.csv")

mu9.500=Sim_SpatialGraphNonComm(mu=rep(0.9,200),N=500)
write.csv(mu9.500,"mu9.500.csv")

mu10.500=Sim_SpatialGraphNonComm(mu=rep(1,200),N=500)
write.csv(mu10.500,"mu10.500.csv")

######## Saving all 500 nodes data to csv
# mu1=read.csv("mu1.500.csv");mu2=read.csv("mu2.500.csv");mu3=read.csv("mu3.500.csv")
# mu4=read.csv("mu4.500.csv");mu5=read.csv("mu5.500.csv");mu6=read.csv("mu6.500.csv")
# mu7=read.csv("mu7.500.csv");mu8=read.csv("mu8.500.csv");mu9=read.csv("mu9.500.csv")
# mu10=read.csv("mu10.500.csv")

df_mu_500=rbind(mu1.500,mu2.500,mu3.500,mu4.500,mu5.500,mu6.500,
                mu7.500,mu8.500,mu9.500,mu10.500)
dim(df_mu_500)
node.num=500
df_mu_500<-df_mu_500 %>%dplyr::filter(order==500)
dim(df_mu_500)
write.csv(df_mu_500,"df_mu_500.csv")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# setwd("C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi-New")
# #source("SPATIAL-PIPELINE-NEW-MODEL.R")
# library(igraph)
# library(Matrix)
# library(tidyverse)
# #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # [1] Spatial Expander Propagation Graph with weighted edges, 
# # node attributes, scale free degree distribution, small world effect and community structures
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# #library(quadtree)
# # Function to generate a spatial scale-free expander graph
# # with the given parameters
# #
# # Arguments:
# #   N: number of nodes
# #   r: cutoff distance for adjacency matrix
# #   beta: power law parameter for distance-dependent probability
# #   m: number of edges added for each new node
# #   alpha: parameter controlling strength of degree effect
# #   mu: parameter controlling preference when connecting nodes. Thus, preference to care about 
# #       spatial distance between nodes when connecting them, vs caring only about vertex degree. 0:care aboutdegree
# #1: care about distance
# #   prewire: rewiring probability for small world effect
# #   node_attrs: list of node attributes to add to graph
# #   edge_weights: logical indicating whether to add edge weights
# #
# # Returns:
# #   A graph object of class igraph
# 
# spatial_scale_free_expander_noncomm <- function(N=50, beta=2, m=2, prewire = 0.1, 
#                                                 mu=0.2,add_edge_weight = FALSE, 
#                                                 add_node_attr = FALSE,r=0.1,alpha=0.2) {
#   
#   
#   ##--Generate points on a unit square
#   points <- matrix(runif(N*2), ncol=2)
#   # calculate distance matrix between all pairs of points
#   dist_mat <- as.matrix(dist(points))
#   # initialize graph with a single node and no edges
#   graph <- list()
#   graph$adj_mat <- matrix(0, nrow = N, ncol = N)
#   graph$adj_mat[1, 2] <- graph$adj_mat[2, 1] <- 1
#   graph$degrees=rep(m,m)
#   # initialize the edge_list
#   graph$edge_list <- as.data.frame(matrix(ncol=2, nrow=0,byrow = T))
#   
#   # Grow the network with new nodes added one at a time
#   # Grow the network with new nodes added one at a time
#   for (i in 3:N) {
#     # preferential attachment effect with scale free technique
#     graph$degrees= rowSums(graph$adj_mat[1:(i-1),])
#     deg_probs<- (graph$degrees^alpha)/sum(graph$degrees^alpha)
#     #--spatial distance effect incorporating short spatial distances (cut-off distance)
#     spat_dist_effect <- ifelse(dist_mat[1:(i-1), i] <= r, 1/(1+(dist_mat[1:(i-1), i])^beta), 1/N)
#     spatial_probs= spat_dist_effect /sum(spat_dist_effect )#normalise
#     # total probability of attachment which is normalised
#     TotalprobsOfAttach=(mu*spatial_probs)+((1-mu)*deg_probs) 
#     #Add m edges to existing nodes with probability proportional to their degree and distance
#     node_to_attach <- sample(1:(i-1), size = m, replace = TRUE, prob =  TotalprobsOfAttach)
#     #update edge_list
#     graph$edge_list<- rbind(graph$edge_list, c(i, node_to_attach))
#     
#     # Attach nodes
#     graph$adj_mat[i, node_to_attach] <- 1
#     graph$adj_mat[node_to_attach, i] <- 1
#     # Update degree vectors
#     graph$degrees[node_to_attach] <- graph$degrees[node_to_attach] + 1
#     graph$degrees[i] <- m
#     graph$degrees <- colSums( graph$adj_mat)
#     #graph$degrees[i] <- sum(graph$adj_mat[i, ])
#     
#     # Small-world rewiring with probability p
#     if (runif(1) < prewire) {
#       # Select random node to rewire within cutoff distance
#       neighbors <- which(graph$adj_mat[i,] == 1)
#       if (length(neighbors) > 0) {
#         new_neighbor <- sample(neighbors, 1)
#         graph$adj_mat[i, new_neighbor] <- 0
#         graph$adj_mat[new_neighbor, i] <- 0
#         available_nodes <- setdiff(1:N, c(i, neighbors))
#         new_neighbor <- sample(available_nodes, 1)
#         graph$adj_mat[i, new_neighbor] <- 1
#         graph$adj_mat[new_neighbor, i] <- 1
#       }
#     }
#   } 
#   
#   ### Graph object
#   GraphModel <- graph.adjacency(as.matrix(graph$adj_mat), mode="undirected")
#   GraphModel=igraph::simplify(GraphModel,remove.loops = TRUE,remove.multiple = TRUE)
#   graph$GraphObject <- GraphModel
#   
#   # Compute graph Laplacian for the expansion property
#   D <- diag(rowSums(graph$adj_mat))
#   L <- D - graph$adj_mat
#   # Compute eigenvalues and eigenvectors of Laplacian matrix
#   eig <- eigen(Matrix(L))
#   eig_values <- Re(eig$values)
#   eig_vectors <- eig$vectors
#   
#   # Compute the expansion ratio
#   lambda_2 <- eig_values[2]
#   lambda_n <- eig_values[nrow(L)]
#   expansion_ratio <- lambda_2 / lambda_n
#   # Check for strong expansion properties
#   #spectral_gap=eigen(Matrix::t(L) %*% L, symmetric=TRUE, only.values=TRUE)$values[n-1]
#   if(!is.connected(GraphModel)||expansion_ratio <= 0){  
#     warning("Graph may not exhibit strong expansion properties or Graph is disconnected")
#   }
#   
#   graph$ExpansionRatio=expansion_ratio
#   
#   if (add_edge_weight) {
#     graph$adj_mat[graph$adj_mat == 1] <- runif(sum(graph$adj_mat == 1))
#     graph$Weights=graph$adj_mat
#   }
#   
#   # Create node attributes if requested
#   if (add_node_attr) {
#     node_attrs <- data.frame(x = points[,1], y = points[,2])
#     graph$NodeAttributes=node_attrs
#   } 
#   return(graph)
# }      
# 
# 
# g=spatial_scale_free_expander_noncomm(N=100, beta=0, m=10, prewire = 0, 
#                                       add_edge_weight = FALSE,mu=0, 
#                                       add_node_attr = T,r=0,alpha=2)
# 
# plot(g$GraphObject,vertex.label=NA,vertex.size=2)
# degree(g$GraphObject)
# 
# #hist(igraph::degree(g$GraphObject), breaks=20, main="")
# 
# ### Testing
# SF=sample_pa(100,power = 2, m=10,directed = F,algorithm = c("psumtree"))
# plot(SF,vertex.label=NA,vertex.size=2)
# degree(SF)
# 
# #hist(igraph::degree(SF), breaks=20, main="")
# 
# lo=layout.norm(as.matrix(cbind(g$NodeAttributes$x,g$NodeAttributes$y)))
# plot(g$GraphObject,vertex.label=NA,vertex.size=2,layout=lo)
# 
# #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # [2] Spatial Expander Propagation Graph with weighted edges, 
# # node attributes, scale free degree distribution, small world effect and community structures
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# #++ Code in R to grow a spatial scale-free expander graph with the following properties. 
# #+This code should be built from scratch and not contain any igraph packages etc.
# 
# # 1) Nodes should be generated on a unit square where the initial set of nodes are distributed spatially with an underlying spatial structure, e.g. by not permitting nodes to land too near each other, or not too far from each other using an in built 'quadtree' algorithm built from scratch
# # 2) Use the distance matrix to create an adjacency matrix for the graph.o create an adjacency matrix based on the distance matrix, you can define a cutoff distance 'r' and set the adjacency matrix element A_{ij} to 1 if dist_{ij} <= 'r' and 0 otherwise . This parameter 'r' controls the strength of the spatial distance effect
# # 3) This network should be grown and nodes should be organized into 
# #non-overlapping communities ( or blocks) with nodes connecting with each 
# #other within the same community or between communities with an attachment probability 
# #favouring short spatial distances and higher degree nodes that follows scale-free sdegree distribution. 
# #Include 'alpha': the parameter that controls the strength of the degree effect
# # 4)Add a rewiring probability 'prewire' for small world effect to the resulting graph to enhance the connectivity to long range nodes. Thus you can randomly rewire each edge with probability 'prewire' to a random node within a certain distance
# # 5) Add arguments to create edge weights, node attributes etc
# # 6)Check the expansion properties of the graph by computing the Laplacian matrix. The Laplacian matrix represents the graph's connectivity and its ability to spread information. A good spatial graph should have a small eigenvalue gap, indicating strong expansion properties.
# # 7)This graph should be a single function where we can tun', N, 'alpha' , 'prewire', 'm', to generate different graph models
# 
# # count nodes for each clusters
# #NumOfNodesPerCluster <- table(graph$clusters)
# # num_nodes_per_comm <- rep(1, numofComm)
# 
# 
# spatial_scale_free_expander_comm <- function(N=50, beta=2, m=1, prewire = 0.1, numofComm=4, 
#                                              mu=0.2,add_edge_weight = FALSE,
#                                              Prob_Betw_comm=0.1,add_node_attr = FALSE,r=0.1,alpha=0.2) {
#   
#   
#   #--Step 1:Generate initial node positions
#   #--Note: Profs Quad tree implementation here
#   node_pos <- matrix(runif(N*2), ncol = 2)
#   
#   # Calculate distance matrix
#   dist_mat <- as.matrix(dist(node_pos))
#   
#   # Initialize adjacency matrix
#   graph <- list()
#   graph$adj_mat <- matrix(0, nrow = N, ncol = N)
#   
#   # Add first node to network
#   graph$adj_mat[1, 1] <- 1
#   
#   # initialize the edge_list and community membership
#   graph$edge_list <- as.data.frame(matrix(ncol=2, nrow=0))
#   
#   # create initial community for first node. Assin first node to the first community
#   graph$clusters <- rep(1, N)
#   
#   # Step 2: Create initial community/block probability matrix
#   P_ij <- matrix(0, nrow = numofComm, ncol = numofComm)
#   
#   Prob_Within_comm=1-Prob_Betw_comm
#   
#   if(is.numeric(Prob_Betw_comm)){
#     diag(P_ij) <- Prob_Within_comm
#     P_ij[lower.tri(P_ij)] <- Prob_Betw_comm
#     P_ij[upper.tri(P_ij)] <- t(P_ij)[upper.tri(P_ij)]
#     
#   }else {
#     P_ij <- matrix(runif(numofComm^2), numofComm, numofComm)
#   }
#   
#   # Step 3: Add nodes to the network one by one
#   # Grow the network by adding new nodes and edges
#   for (i in 3:N) {
#     # Step 3i: Create arbitrary number of communities/blocks
#     comm_probs<-rep(1/numofComm, numofComm)#sample.int(numofComm, (i-1), replace = TRUE)
#     graph$clusters[i] <- sample(numofComm, 1, prob = comm_probs)
#     
#     # Step 3ii: Define community probability matrix
#     P_comm <- P_ij[1:numofComm, 1:numofComm]
#     
#     ##--Step 3iii: Define probability of attachment function for within and between-community connections
#     # P_unique_within <- rep(0, (i-1))
#     # P_unique_between <- rep(0, (i-1))
#     
#     spatial_probs <- ifelse(dist_mat[1:(i-1), i] <= r, 1/(1+(dist_mat[1:(i-1), i])^beta), 1/N)
#     graph$degrees<- rowSums(graph$adj_mat[1:(i-1), ])
#     deg_probs <- (graph$degrees^alpha) / sum( graph$degrees^alpha)
#     
#     # Step 3iv: Add edges within and between communities based on community and attachment probabilty matrix 
#     P_noncomm <- (mu*spatial_probs)+((1-mu)*deg_probs) 
#     P_noncomm<-P_noncomm/sum(P_noncomm)
#     if (numofComm==1){
#       P_join <-P_comm[graph$clusters[(i-1)]] * P_noncomm
#     }else{
#       P_join <-P_comm[graph$clusters[1:(i-1)], graph$clusters[i]] * P_noncomm
#     }
#     
#     #Connect nodes based on attachment probability
#     neighbors_within <- sample(1:(i-1), m, replace = TRUE, prob = P_join)
#     graph$adj_mat[i, neighbors_within] <- 1
#     graph$adj_mat[neighbors_within, i] <- 1
#     graph$degrees[i] <- graph$degrees[i] + 1
#     graph$degrees[neighbors_within] <- graph$degrees[neighbors_within] + 1
#     graph$degrees[i] <- sum(graph$adj_mat[i, ])
#     # Bind to the empty edge_list for within communities
#     graph$edge_list <- rbind(graph$edge_list, c(i, neighbors_within))
#     
#     # Step 4: Small-world rewiring with probability p
#     if (runif(1) < prewire) {
#       # Select random node to rewire within cutoff distance
#       neighbors <- which(graph$adj_mat[i,] == 1)
#       if (length(neighbors) > 0) {
#         new_neighbor <- sample(neighbors, 1)
#         graph$adj_mat[i, new_neighbor] <- 0
#         graph$adj_mat[new_neighbor, i] <- 0
#         available_nodes <- setdiff(1:N, c(i, neighbors))
#         new_neighbor <- sample(available_nodes, 1)
#         graph$adj_mat[i, new_neighbor] <- 1
#         graph$adj_mat[new_neighbor, i] <- 1
#       }
#     }
#   }
#   
#   ### Graph object
#   GraphModel <- graph.adjacency(as.matrix(graph$adj_mat), mode="undirected")
#   GraphModel=simplify(GraphModel,remove.loops = TRUE,remove.multiple = TRUE)
#   graph$GraphObject <- GraphModel
#   
#   # Compute graph Laplacian for the expansion property
#   D <- diag(rowSums(graph$adj_mat))
#   L <- D - graph$adj_mat
#   # Compute eigenvalues and eigenvectors of Laplacian matrix
#   eig <- eigen(Matrix(L))
#   eig_values <- Re(eig$values)
#   eig_vectors <- eig$vectors
#   
#   # Compute the expansion ratio
#   lambda_2 <- eig_values[2]
#   lambda_n <- eig_values[nrow(L)]
#   expansion_ratio <- lambda_2 / lambda_n
#   # Check for strong expansion properties
#   #spectral_gap=eigen(Matrix::t(L) %*% L, symmetric=TRUE, only.values=TRUE)$values[n-1]
#   if(!is.connected(GraphModel)||expansion_ratio <= 0){  
#     warning("Graph may not exhibit strong expansion properties or Graph is disconnected")
#   }
#   
#   graph$ExpansionRatio=expansion_ratio
#   
#   if (add_edge_weight) {
#     graph$adj_mat[graph$adj_mat == 1] <- runif(sum(graph$adj_mat == 1))
#     graph$Weights=graph$adj_mat
#   }
#   
#   # Create node attributes if requested
#   if (add_node_attr) {
#     node_attrs <- data.frame(x = node_pos[,1], y = node_pos[,2])
#     graph$NodeAttributes=node_attrs
#   } 
#   return(graph)
# }      
# 
# g=spatial_scale_free_expander_comm(N=200, beta=2, m=10, prewire = 0.01, numofComm=2, 
#                                    mu=1,add_edge_weight = FALSE,
#                                    Prob_Betw_comm=0.0001,add_node_attr = T,r=0.1,alpha=0.2) 
# 
# lo=layout.norm(as.matrix(cbind(g$NodeAttributes$x,g$NodeAttributes$y)))
# 
# 
# plot(g$GraphObject,vertex.label=NA,vertex.size=2)
# 
# plot(g$GraphObject,vertex.label=NA,vertex.size=2,layout=lo)
# 
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# #+ Training a machine learning regression model with random forest
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# #Regression models
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
# ##--Data importation and processing
# library(dplyr)
# df50=read.csv("df_mu_50.csv",sep = ",", header = T)
# df100=read.csv("df_mu_100.csv",sep = ",", header = T)
# df150=read.csv("df_mu_150.csv",sep = ",", header = T)
# 
# Data=df150#rbind(df50,df100,df150)
# 
# 
# df=as_tibble(Data)
# df%>%count(GraphName,sort=T)
# df%>%View()
# 
# Data= df%>%dplyr::select(-c(connected,max_component,minDegree,maxDegree,
#                             threshold,GraphReplicate,GraphID,X,GraphName))%>%
#   mutate_if(is.character,factor)
# ##--Shuffle--data
# Data<- Data[sample(1:nrow(Data)), ]##shuffle row indices and randomly re order 
# head(Data)
# 
# 
# #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # Data prepropecessing recipe
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# library(visdat)
# 
# #Visualizing all data structure
# vis_dat(Data)
# 
# #Checking levels of factor or character nd creating a table
# library(gt)
# Data %>% 
#   count(GraphName,
#         sort = TRUE)%>%
#   gt()
# 
# set.seed(4595)
# data_split <- rsample::initial_split(Data, strata = "mu")#, prop = 0.75)
# 
# train.data <- training(data_split)
# test.data  <- testing(data_split)
# 
# dim(train.data)
# 
# ##---Creating--Data--Recipe---and--Prep--
# df_recipe=recipe(mu~., data=train.data)%>%
#   update_role(order,new_role = "numberOfNodes")%>%
#   step_zv(all_predictors())%>%
#   step_normalize(all_predictors())%>%
#   step_corr(all_predictors(), threshold = 0.7, method = "spearman")
# 
# summary(df_recipe)
# 
# df_prep <- 
#   df_recipe %>% # use the recipe object
#   prep() %>% # perform the recipe on training data
#   juice() # extract only the preprocessed dataframe 
# 
# #Take a look at the data structure:
# glimpse(df_prep)
# 
# ####------Setting--up--Cross--validation--set--to--tunes--the--models--we--are--going--to--compare
# set.seed(100)
# mu_folds <-
#   bootstraps(train.data,strata = mu)
# 
# mu.folds <- 
#   recipes::bake(
#     xgb.recipe, 
#     new_data = training(data_split)
#   ) %>%  
#   rsample::vfold_cv(v = 5)
# # mu_folds <-
# #   vfold_cv(train.data, 
# #            v = 10, 
# #            strata = mu)
# 
# library(usemodels)
# #This package contains information for options on setting up common types of models
# #eg
# use_ranger(mu~.,data=train.data)
# ##########------------SECIFYING--MODEL:SETTING--MODEL--ENGINE--AND--WORKFLOW----------#####################
# 
# # The process of specifying our models is always as follows:
# ###Pick a model type
# ###set the engine
# ###Set the mode: regression or classification
# 
# #######################--Setting---up---model--and--engine--###############################
# 
# ##--[A]--Random--Forest---model-engine and workflow---
# library(ranger)
# 
# rf.recipe=recipe(formula=mu~., data=train.data)%>%
#   update_role(order,new_role = "numberOfNodes")
# # step_zv(all_predictors())%>%
# #  step_normalize(all_predictors())%>%
# # step_corr(all_predictors(), threshold = 0.7, method = "spearman")
# 
# summary(rf.recipe)
# 
# rf.prep <- 
#   rf.recipe %>% # use the recipe object
#   prep() %>% # perform the recipe on training data
#   juice() # extract only the preprocessed dataframe 
# 
# #Take a look at the data structure:
# glimpse(rf.prep)
# 
# ### Parallel running
# #cores <- parallel::detectCores()
# #cores
# 
# rf.model <- 
#   rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>% 
#   set_engine("ranger") %>% 
#   set_mode("regression")
# 
# ##---Create--workflows
# # To combine the data preparation recipe with the model building, 
# # we use the package workflows. 
# #A workflow is an object that can bundle together your pre-processing recipe,
# #modeling, and even post-processing requests (like calculating the RMSE).   
# 
# rf.wkflow <-
#   workflow() %>%
#   add_recipe(rf.recipe) %>% 
#   add_model(rf.model) 
# 
# ## shows num of params to be tuned
# extract_parameter_set_dials(rf.wkflow)
# 
# set.seed(90223)
# doParallel::registerDoParallel()
# 
# rf.tune.params <-
#   tune_grid(rf.wkflow, 
#             resamples = mu_folds, 
#             grid = 11)
# 
# 
# ###---Explore results from tuned model
# show_best(rf.tune.params, metric="rmse")
# show_best(rf.tune.params, metric="rsq")
# 
# autoplot(rf.tune.params )
# 
# rf.final.wkflow<-rf.wkflow%>%
#   finalize_workflow(select_best(rf.tune.params ))
# 
# ##--Final fitted model
# rf.model.fit=last_fit(rf.final.wkflow,data_split )
# ##--collect metrics
# collect_metrics(rf.model.fit)
# 
# collect_predictions(rf.model.fit)%>%
#   ggplot(aes(mu,.pred))+
#   geom_abline(lty=2,color="blue")+
#   geom_point(alpha=0.5,color="red")+
#   coord_fixed()
# 
# ##---Making Predictions
# predicted.vals=predict(rf.model.fit$.workflow[[1]],test.data[,])
# 
# ##--Variable--importance
# library(vip)
# varimp.model=rf.model%>%
#   finalize_model(select_best(rf.tune.params ))%>%
#   set_engine("ranger",importance="permutation")
# 
# workflow() %>%
#   add_recipe(rf.recipe) %>% 
#   add_model(varimp.model)%>% 
#   fit(train.data)%>% 
#   pull_workflow_fit()%>% 
#   vip(aesthetics=list(alpha=0.8,fill="midnightblue"))
# 
# # rf_form_fit=
# #   rand_forest(mode = "regression", mtry = 24, trees = 500) %>%
# #   set_engine("ranger") %>%
# #   fit(mu~.,
# #       data = train.data
# #   )
# # 
# # rf_form_fit
# # 
# # ##--Predicting with formula model
# # test_results_form <- 
# #   test.data %>%
# #   select(mu) %>%
# #   bind_cols(
# #     predict(rf_form_fit, new_data = test.data[,])
# #   )
# # 
# # test_results_form
# 
# #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ##--[B]--Boosted--tree--(XGBoost)
# #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# library(xgboost)
# xgb.recipe <- 
#   recipe(formula = mu ~ ., data = train.data) %>% 
#   step_zv(all_predictors()) 
# 
# 
# xgb.model <- 
#   boost_tree(rees = tune(), min_n = tune(), tree_depth = tune(), learn_rate = tune(),
#              loss_reduction = tune(), sample_size = tune()) %>% 
#   set_engine("xgboost") %>% 
#   set_mode("regression")
# 
# 
# xgb.wkflow <-
#   workflow() %>%
#   add_recipe(df_recipe) %>% 
#   add_model(xgb.model)
# 
# set.seed(90223)
# doParallel::registerDoParallel()
# 
# xgb.tune.params <-
#   tune_grid(xgb.wkflow, 
#             resamples = mu_folds, 
#             grid = 11)
# 
# ##------General--example----###
# xgb.recipe <- 
#   recipes::recipe(mu ~ ., data = training(data_split)) %>%
#   # convert categorical variables to factors
#   recipes::step_string2factor(all_nominal()) %>%
#   # combine low frequency factor levels
#   recipes::step_other(all_nominal(), threshold = 0.01) %>%
#   # remove no variance predictors which provide no predictive information 
#   recipes::step_nzv(all_nominal()) %>%
#   recipes::step_zv(all_predictors()) %>%
#   prep()
# 
# 
# # XGBoost model specification
# xgb.model <- 
#   parsnip::boost_tree(
#     mode = "regression",
#     trees = 1000,
#     min_n = tune(),
#     tree_depth = tune(),
#     learn_rate = tune(),
#     loss_reduction = tune()
#   ) %>%
#   set_engine("xgboost", objective = "reg:squarederror")
# 
# 
# # grid specification
# xgb.params <- 
#   dials::parameters(
#     min_n(),
#     tree_depth(),
#     learn_rate(),
#     loss_reduction()
#   )
# 
# # grid space
# xgb.grid <- 
#   dials::grid_max_entropy(
#     xgb.params, 
#     size = 60
#   )
# knitr::kable(head(xgb.grid))
# 
# # Workflow
# xgb.wkflow <- 
#   workflows::workflow() %>%
#   add_model(xgb.model) %>% 
#   add_formula(mu ~ .)
# 
# # hyperparameter tuning
# xgb.tuned <- tune::tune_grid(
#   object = xgb.wkflow,
#   resamples = mu.folds,
#   grid = xgb.grid,
#   metrics = yardstick::metric_set(rmse, rsq, mae),
#   control = tune::control_grid(verbose = TRUE)
# )
# 
# # Best hyper parameter values which minimises rmse
# xgb.tuned %>%
#   tune::show_best(metric = "rmse") %>%
#   knitr::kable()
# 
# 
# #isolate the best performing hyperparameter values.
# xgb.best.params <- xgb.tuned %>%
#   tune::select_best("rmse")
# knitr::kable(xgb.best.params)
# 
# #Finalize the XGBoost model to use the best tuning parameters.
# xgb.final.model<- xgb.model %>% 
#   finalize_model(xgb.best.params)
# 
# ##--Evauating--Performance---of---model--on--train--and--test--data
# train_processed <- bake(xgb.recipe,  new_data = training(data_split))
# 
# train_prediction <- xgb.final.model %>%
#   # fit the model on all the training data
#   fit(
#     formula = mu ~ ., 
#     data    = train_processed
#   ) %>%
#   # predict mu for the training data
#   predict(new_data = train_processed) %>%
#   bind_cols(training(data_split))
# 
# xgb.score.train <- 
#   train_prediction %>%
#   yardstick::metrics(mu, .pred) %>%
#   mutate(.estimate = format(round(.estimate, 2), big.mark = ","))
# knitr::kable(xgb.score.train)
# 
# #cbind(head(train_prediction$mu),head(train.data$mu))
# 
# ##------Predicting on test set
# test_processed  <- bake(xgb.recipe, new_data = testing(data_split))
# test_prediction <- xgb.final.model %>%
#   # fit the model on all the training data
#   fit(
#     formula = mu ~ ., 
#     data    = train_processed
#   ) %>%
#   # use the training model fit to predict the test data
#   predict(new_data = test_processed) %>%
#   bind_cols(testing(data_split))
# 
# # measure the accuracy of our model using `yardstick`
# xgb.score.test <- 
#   test_prediction %>%
#   yardstick::metrics(mu, .pred) %>%
#   mutate(.estimate = format(round(.estimate, 2), big.mark = ","))
# 
# knitr::kable(xgb.score.test)
# 
# 
# 
# #To quickly check that there is not an obvious issue with our model's predictions, 
# #let's plot the test data residuals.
# 
# mu_prediction_residual <- test_prediction %>%
#   arrange(.pred) %>%
#   mutate(residual_pct = (mu - .pred) / .pred) %>%
#   select(.pred, residual_pct)
# ggplot(mu_prediction_residual, aes(x = .pred, y = residual_pct)) +
#   geom_point() +
#   xlab("Predicted mu") +
#   ylab("Residual (%)") +
#   scale_x_continuous(labels = scales::dollar_format()) +
#   scale_y_continuous(labels = scales::percent)
# 
# 
# 
# 
# 
# 
# ##--[C]--GLM--network
# glm.model <-
#   linear_reg(penalty = tune(),mixture = 1) %>% 
#   set_engine("glmnet") %>%
#   set_mode("regression")
# 
# glm.wkflow <-
#   workflow() %>%
#   add_recipe(df_recipe) %>% 
#   add_model(glm.model)
# lr_reg_grid <- tibble(penalty = 10^seq(-4, -1, length.out = 30))
# 
# 
# 
# ##############-----------------Evaluate------models--------------------############
# # Now we can use our validation 
# # set (cv_folds) to estimate the performance of our models using the fit_resamples() 
# # function to fit the models on each of the folds and store the results.
# # 
# # Note that fit_resamples() will fit our model to each resample and 
# # evaluate on the heldout set from each resample. 
# # The function is usually only used for computing performance metrics across some 
# # set of resamples to evaluate our models (like accuracy) - the models are
# # not even stored. However, in our example we save the predictions in order to
# # visualize the model fit and residuals with control_resamples(save_pred = TRUE).
# #Finally, we collect the performance metrics with collect_metrics() and 
# #pick the model that does best on the validation set.    
# 
# rf.fit <-
#   rf.wkflow %>% 
#   fit(mu~.,data = train.data)
# 
# 
# rf.res %>%  collect_metrics(summarize = TRUE)
# 
# 
# xgb.res <- 
#   xgb.wkflow %>% 
#   fit_resamples(
#     resamples = cv_folds, 
#     metrics = metric_set(
#       recall, precision, f_meas, 
#       accuracy, kap,
#       roc_auc, sens, spec),
#     control = control_resamples(save_pred = TRUE)
#   ) 
# 
# xgb.res %>% collect_metrics(summarize = TRUE)




# # # Function to generate Barabasi-Albert scale-free network with preferential attachment
# # # n: number of nodes in the network
# # # m: number of edges to attach from each new node
# # # power: power-law exponent for the degree distribution
# # generate_ba_network <- function(n, m, power) {
# #   # Create initial network with m nodes and complete graph
# #   g <- graph.empty(n = m, directed = FALSE)
# #   g <- add_edges(g, 1:m)
# #   graph=list()
# #   # Add new nodes with preferential attachment
# #   for (i in seq(m+1, n)) {
# #     # Calculate the degree of each node in the current network
# #     graph$degrees <- degree(g)
# #     
# #     # Calculate the probability of connecting to each existing node
# #     prob <- graph$degrees^power
# #     prob <- prob/sum(prob)
# #     
# #     # Choose m nodes to attach the new node to
# #     attach_to <- sample(1:vcount(g), size = m, replace = TRUE, prob = prob)
# #     
# #     # Add edges between the new node and the chosen nodes
# #     g <- add_vertices(g, 1)
# #     g <- add_edges(g, c(rep(i, m), attach_to))
# #   }
# #   
# #   # Return the generated network
# #   return(list(g,graph$degrees))
# # }
# 
# ba_model <- function(n, m) {
#   # n: number of nodes
#   # m: number of edges to attach from a new node to existing nodes
#   
#   # Create initial fully-connected graph with m+1 nodes
#   edges <- matrix(rep(0, nrow=m,ncol=N, byrow=TRUE)
#   edges[lower.tri(edges)] <- t(edges)[lower.tri(edges)]
#   edges <- edges[edges != 0]
#   degs=rep(m,N)
#   # Add nodes
#   for (i in (m+2):n) {
#     # Calculate probability distribution based on degree of existing nodes
#     degs <- colSums(matrix(edges))
#     prob <- degs / sum(degs)
#     # Select m nodes to attach to
#     targets <- sample(1:(i-1), m, replace=TRUE, prob=prob)
#     
#     # Add edges
#     edges <- c(edges, rep(i, m), targets)
#     edges <- c(edges, rep(targets, each=m), i)
#   }
#   
#   # Convert edge list to adjacency matrix
#   mat <- matrix(0, nrow=n, ncol=n)
#   for (i in 1:length(edges)) {
#     mat[edges[i], edges[i+length(edges)]] <- 1
#   }
#   
#   return(mat)
# }
# 
# ba_model(100, 4)
# 
# edges <- matrix(0, nrow=m,ncol=n)
# edges[lower.tri(edges)] <-1
# edges[upper.tri(edges)]= edges[lower.tri(edges)]
# edges
# 
# M <- matrix(0, nrow = m, ncol = n)
# # Set the upper triangular part of the matrix to be symmetric
# M[lower.tri(M)] <- 1
# M[t(upper.tri(t(M)))] <- 1
# 
# #edges[upper.tri(edges)]=edges[lower.tri(t(edges))]
# 
# mat <- matrix(0, m, n)
# # fill in the upper and lower triangles of the first m rows and columns with 1's
# mat[lower.tri(mat[1:m, m:1])] <- 1
# mat[upper.tri(mat[1:m, m:1])] <- 1
# 
# mat <- matrix(0, m, n)
# # set the diagonal entries and the symmetric entries above the diagonal
# mat[lower.tri(mat,diag = F)] <- 1
# mat1=mat
# mat1[upper.tri(mat1)] <- t(mat1)[lower.tri(mat1)]


