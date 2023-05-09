setwd("C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi-New")
#source("SPATIAL-PIPELINE-NEW.R")
spatial_scale_free_expander_noncomm <- function(N=50, beta=2, m=2, prewire = 0.1, 
                                                mu=0.2,add_edge_weight = FALSE, 
                                                add_node_attr = FALSE,r=0.1,alpha=0.2) {
  
  
  ##--Generate points on a unit square
  points <- matrix(runif(N*2), ncol=2)
  # calculate distance matrix between all pairs of points
  dist_mat <- as.matrix(dist(points))
  # Initialize network with m nodes
  graph <- list()
  graph$nodes <- 1:m
  graph$adj_mat <-matrix(0, ncol = N, nrow = m)
  graph$adj_mat[1:2]=graph$adj_mat[2:1]=1
  graph$degrees=rep(m,N)
  # initialize the edge_list
  graph$edge_list <- as.data.frame(matrix(ncol=2, nrow=0,byrow = T))
  
  # Grow the network with new nodes added one at a time
  # Grow the network with new nodes added one at a time
  for (i in (m+1):N) {
    # preferential degree attachment effect with scale free technique
    #    graph$degrees <- rowSums(graph$adj_mat[1:(i-1),])
    deg_probs<- (graph$degrees[1:(i-1)]^alpha)/sum(graph$degrees[1:(i-1)]^alpha)
    #--spatial distance effect incorporating short spatial distances (cut-off distance)
    spat_dist_effect <- ifelse(dist_mat[1:(i-1), i] <= r, 1/(1+(dist_mat[1:(i-1), i])^beta), 1/N)
    spatial_probs= spat_dist_effect /sum(spat_dist_effect )#normalise
    # total probability of attachment which is normalised
    TotalprobsOfAttach=(mu*spatial_probs)+((1-mu)*deg_probs) 
    #Add m edges to existing nodes with probability proportional to their degree and distance
    node_to_attach <- sample(1:(i-1), size = m, replace = TRUE, prob =  TotalprobsOfAttach)
    #update edge_list
    graph$edge_list<- rbind(graph$edge_list, c(i, node_to_attach))
    
    # Add new node and edges
    graph$nodes <- c( graph$nodes, i)
    graph$adj_mat <- rbind(graph$adj_mat, matrix(0, ncol = N, nrow = 1))
    graph$adj_mat[i, node_to_attach] <- 1
    graph$adj_mat[node_to_attach, i] <- 1
    # Update degree vectors
    graph$degrees[i] <- m 
    graph$degrees[node_to_attach] <- graph$degrees[node_to_attach]+1
   # graph$degrees[i] <- sum(graph$adj_mat[i, ])
    
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
  
 # graph$degrees <- rowSums( graph$adj_mat)
  
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

g=spatial_scale_free_expander_noncomm(N=100, beta=0, m=5, prewire = 0, 
                                      add_edge_weight = FALSE,mu=0, 
                                      add_node_attr = T,r=0,alpha=6)

plot(g$GraphObject,vertex.label=NA,vertex.size=2)
g$degrees

#hist(igraph::degree(g$GraphObject), breaks=20, main="")

### Testing
SF=sample_pa(100,power = 6, m=3,directed = F,algorithm = c("psumtree"))
plot(SF,vertex.label=NA,vertex.size=2)
degree(SF)

#hist(igraph::degree(SF), breaks=20, main="")

lo=layout.norm(as.matrix(cbind(g$NodeAttributes$x,g$NodeAttributes$y)))
plot(g$GraphObject,vertex.label=NA,vertex.size=2,layout=lo)

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

library(Matrix)
library(igraph)

spatial_scale_free_expander <- function(N, beta, m, prewire = 0.1, 
                                        mu=0.2,add_edge_weight = FALSE, 
                                        add_node_attr = FALSE,r=0.1,alpha=0.2) {
  #Notes:
  ##--(1) use of of N or lambda
  ##--(2) fix L or ignore
  ###--(3) mu > 0 for preferential attachment: mu captures what p does in the doi:10.1088/1367-2630/9/6/190
  ###--(4) use unit square instead of torus
  
  ##--Generate points on a unit square
  points <- matrix(runif(N*2), ncol=2)
  
  
  # calculate distance matrix between all pairs of points
  dist_mat <- as.matrix(dist(points))
  
  
  # initialize graph with a single node and no edges
  graph <- list()
  graph$adj_mat <- matrix(0, nrow = N, ncol = N)
  graph$adj_mat[1:2, 2:1] <- 1
  
  # initialize the edge_list
  graph$edge_list <- as.data.frame(matrix(ncol=2, nrow=0))
  
  # Grow the network with new nodes added one at a time
  for (i in 3:N) {
    ##--calculate probability of attaching to each existing node
    probattach <- rep(0, (i-1))
    # preferential attachment effect with scale free technique
    graph$degrees= rowSums(graph$adj_mat[1:(i-1), ])
    deg_probs<- (graph$degrees^alpha)/(sum(graph$degrees^alpha))
    
    #--spatial distance effect incorporating short spatial distances (cut-off distance)
    spatial_prob <- ifelse(dist_mat[1:(i-1), i] <= r, 1/(1+(dist_mat[1:(i-1), i])^beta), 0.01)
    
    # total probability of attachment 
    f_ij=(mu*spatial_prob) + ((1-mu)*deg_probs)
    # normalise probabilities
    probattach=f_ij/sum(f_ij)
    
    #Add m edges to existing nodes with probability proportional to their degree and distance
    node_to_attach <- sample(1:(i-1), size = m,replace = TRUE, prob = probattach)
    graph$adj_mat[i, node_to_attach] <- 1
    graph$adj_mat[node_to_attach, i] <- 1
    # graph$degrees[i-1] <- graph$degrees[i-1] + 1
    #  graph$degrees[node_to_attach-1] <-graph$degrees[node_to_attach-1] + 1
    
    #graph$degrees[c(i, node_to_attach)] <-  graph$degrees[c(i, node_to_attach)] + 1
    
    
    # Bind to the empty edge_list
    graph$edge_list <- rbind(graph$edge_list, c(i, node_to_attach))
    
    
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
  GraphModel <- igraph::graph.adjacency(as.matrix(graph$adj_mat), mode="undirected")
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

g=spatial_scale_free_expander(N=100, beta=2, m=4, prewire = 0, 
                              add_edge_weight = FALSE,mu=0, 
                              add_node_attr = T,r=0.1,alpha=6)

g$degrees
plot(g$GraphObject,vertex.label=NA,vertex.size=2)


lo=layout.norm(as.matrix(cbind(g$NodeAttributes$x,g$NodeAttributes$y)))

plot(g$GraphObject,vertex.label=NA,vertex.size=2,layout=lo)


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

# for (i in 1:numofComm) {
#   for (j in 1:numofComm) {
#     if (i == j) {
#       P_ij[i, j] <- probWithin
#     } else {
#       P_ij[i, j] <- probBetween
#     }
#   }
# }
library(Matrix)
library(stats)
library(spatstat)
spatial_scale_free_expander <- function(N=50, beta=2, m=1, prewire = 0.1, numofComm=4, 
                                        mu=0.2,add_edge_weight = FALSE,
                                        probBetween=0.1, probWithin=0.5,
                                        add_node_attr = FALSE,r=0.1,alpha=0.2) {
  
  
  #--Step 1:Generate initial node positions
  #--Note: Profs Quad tree implementation here
  node_pos <- matrix(runif(N*2), ncol = 2)
  
  # Calculate distance matrix
  dist_mat <- as.matrix(dist(node_pos))
  
  # Initialize adjacency matrix
  graph <- list()
  graph$adj_mat <- matrix(0, nrow = N, ncol = N)
  
  # Add first node to network
  graph$adj_mat[1, 1] <- 1
  
  # initialize the edge_list and community membership
  graph$edge_list <- as.data.frame(matrix(ncol=2, nrow=0))
  
  # create initial community for first node. Assin first node to the first community
  graph$clusters <- rep(1, N)
  
  # Step 2: Create initial community/block probability matrix
  P_ij <- matrix(0, nrow = numofComm, ncol = numofComm)
  if(is.numeric(probWithin) && is.numeric(probBetween)){
    diag(P_ij) <- probWithin
    P_ij[lower.tri(P_ij)] <- probBetween 
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
    P_ij_comm <- P_ij[1:numofComm, 1:numofComm]
    
    ##--Step 3iii: Define probability of attachment function for within and between-community connections
    P_unique_within <- rep(0, (i-1))
    P_unique_between <- rep(0, (i-1))
    
    spatial_probs <- ifelse(dist_mat[1:(i-1), i] <= r, 1/(1+(dist_mat[1:(i-1), i])^beta), 0)
    graph$degrees<- rowSums(graph$adj_mat[1:(i-1), ])
    deg_probs_within <- ( graph$degrees^alpha) / sum( graph$degrees^alpha)
    
    
    # Step 3iv: Add edges within the same community based on community and attachment probabilty matrix 
    P_unique_within <- (mu*spatial_probs)+((1-mu)*deg_probs_within) 
    P_unique_within<-P_unique_within/sum(P_unique_within)
    P_within <- P_ij[graph$clusters[1:(i-1)], graph$clusters[i]] * (1-prewire) + P_unique_within * prewire
    P_within[is.na(P_unique_within)] <- 0
    neighbors_within <- sample(1:(i-1), m, replace = TRUE, prob = P_within)
    graph$adj_mat[i, neighbors_within] <- 1
    graph$adj_mat[neighbors_within, i] <- 1
    graph$degrees[i] <- graph$degrees[i] + 1
    degrees[neighbors_within] <- degrees[neighbors_within] + 1
    graph$degrees[i] <- sum(graph$adj_mat[i, ])
    # Bind to the empty edge_list for within communities
    graph$edge_list <- rbind(graph$edge_list, c(i, neighbors_within))
    
    ##--Step 3v: Define probability of attachment function for between-community connections
    # Add edges between communities based on community and attachment probabilty matrix
    if (sum(P_ij_comm == probBetween) > 0) {
      deg_probs_between=deg_probs_within
      P_unique_between <-((1-prewire)*deg_probs_between)+(prewire*spatial_probs) #outer(deg_probs_between, rep(prewire, i-1), "*")  
      P_unique_between<-P_unique_between/sum(P_unique_between)
      P_unique_between[graph$clusters[1:(i-1)] == graph$clusters[i]] <- 0
      P_between <- P_ij[graph$clusters[1:(i-1)], graph$clusters[i]] * (1-prewire) + P_unique_between * prewire
      neighbors_between <- sample(1:(i-1), m, replace = TRUE, prob = P_between)
      graph$adj_mat[i, neighbors_between] <- 1
      graph$adj_mat[neighbors_between, i] <- 1
      graph$degrees[i] <- graph$degrees[i] + 1
      degrees[neighbors_between ] <- degrees[neighbors_between] + 1
      graph$degrees[i] <- sum(graph$adj_mat[i, ])
      # Bind to the empty edge_list for between communities
      graph$edge_list <- rbind(graph$edge_list, c(i, neighbors_between))
    }
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
  GraphModel=simplify(GraphModel,remove.loops = TRUE,remove.multiple = TRUE)
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
    node_attrs <- data.frame(x = node_pos[,1], y = node_pos[,2])
    graph$NodeAttributes=node_attrs
  } 
  return(graph)
}      

g=spatial_scale_free_expander(N=100, beta=2, m=1, prewire = 0.01, numofComm=2, 
                              mu=0.2,add_edge_weight = FALSE,
                              probBetween=0.001, probWithin=0.04,
                              add_node_attr = T,r=0.15,alpha=0.2)

lo=layout.norm(as.matrix(cbind(g$NodeAttributes$x,g$NodeAttributes$y)))


plot(g$GraphObject,vertex.label=NA,vertex.size=2)

plot(g$GraphObject,vertex.label=NA,vertex.size=2,layout=lo)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [1] Spatial Expander Propagation Graph with weighted edges, node attributes and community structures
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(Matrix)
library(igraph)

N=n=50;L=2;r=0.1;lambda=10;use_edge_weights=FALSE

spatial_expander_propagation_graph <- function(n, L, r, lambda,deg=4,
                                               use_edge_weights=FALSE,use_node_attributes=FALSE,
                                               use_comunity_structures=FALSE,
                                               community_groups=4) {
  
  # Generate 2D poisson point process on the L-dimensional torus
  coords <- matrix(runif(n*L), n, L)
  radius_sq <- r^2
  
  # Compute the euclidean (pairwise) distances between points
  num_coords <- nrow(coords)
  dist_mat <- matrix(0, nrow=num_coords, ncol=num_coords)
  for (i in 1:(num_coords-1)) {
    for (j in (i+1):num_coords) {
      dist_mat[i,j] <- sqrt(sum((coords[i,] - coords[j,])^2))
      dist_mat[j,i] <- dist_mat[i,j]
    }
  }
  
  #dist_mat <- as.matrix(dist(coords))
  #Create adjacency matrix with and without edge weights
  A <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      # Compute shortest distance between points on torus
      dist_sq <- sum(pmin(abs(coords[i,]-coords[j,]), 1-abs(coords[i,]-coords[j,]))^2)
      if (dist_sq <= radius_sq) {
        if (use_edge_weights) { #Uses edge weights for the adjacency matrix and
          # Assign edge weights to edges as inverse of distance
          A[i,j] <- 1/sqrt(dist_sq)
          A[j,i] <- A[i,j]#undirected graphs
        }}
      else { # uses nearest neighbour connection for the adjacency matrix (connects spatial points closest to each other)
        neighbors <- which(dist_mat[i,] < r)
        neighbors <- neighbors[neighbors != i]
        if (length(neighbors) > deg) {
          neighbors <- sample(neighbors, deg)
        }
        A[i, neighbors] <- 1 
        A[neighbors,i] <- A[i, neighbors]#undirected graphs
      }
    }
  }
  # Compute graph Laplacian for the expansion property
  D <- diag(rowSums(A))
  L <- D - A
  
  # Check for strong expansion properties
  spectral_gap <- eigen(Matrix::t(L) %*% L, symmetric=TRUE, only.values=TRUE)$values[n-1]
  if (spectral_gap <= 0) {
    warning("Graph may not exhibit strong expansion properties")
  }
  
  # Convert Laplacian to igraph object
  G <- graph.adjacency(as.matrix(A), mode="undirected")
  G=G%>%simplify(remove.loops = TRUE,remove.multiple = TRUE)
  
  # Add node attributes
  if (use_node_attributes){ 
    V(G)$coords <- coords
  }
  # Add community structures
  if (use_comunity_structures){ 
    V(G)$community <- cutree(cluster_walktrap(G),k=community_groups)
  }
  else{
    G <- graph.adjacency(as.matrix(A), mode="undirected")
    G=G%>%simplify(remove.loops = TRUE,remove.multiple = TRUE)
  }
  return(G)
}   

Graph <- spatial_expander_propagation_graph(n=500, L=1, r=0.01, 10,deg=2,
                                            use_edge_weights=T,use_node_attributes=FALSE,
                                            use_comunity_structures=FALSE,community_groups=4)
Graph
plot(Graph,vertex.label=NA,vertex.size=2)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.1] Spatial Expander Propagation Graph with weighted edges, 
# node attributes, scale free degree distribution, small world effect and community structures
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#library(quadtree)
# Function to generate a spatial scale-free expander graph
# with the given parameters
#
# Arguments:
#   N: number of nodes
#   lambda: intensity parameter for 2D Poisson point process
#   L: dimension of torus
#   r: cutoff distance for adjacency matrix
#   beta: power law parameter for distance-dependent probability
#   m: number of edges added for each new node
#   alpha: parameter controlling strength of degree effect
#   sigma: parameter controlling decay of degree effect with distance
#   mu: parameter controlling preference when connecting nodes. Thus, preference to care about 
#       spatial distance between nodes when connecting them, vs caring only about vertex degree. 0:care aboutdegree
#1: care about distance
#   prewire: rewiring probability for small world effect
#   node_attrs: list of node attributes to add to graph
#   edge_weights: logical indicating whether to add edge weights
#
# Returns:
#   A graph object of class igraph

spatial_scale_free_expander <- function(N, beta, m, prewire = 0.1, 
                                        mu=0.2,add_edge_weight = FALSE, 
                                        add_node_attr = FALSE,r=0.1,alpha=0.2) {
  #Notes:
  ##--(1) use of of N or lambda
  ##--(2) fix L or ignore
  ###--(3) mu > 0 for preferential attachment: mu captures what p does in the doi:10.1088/1367-2630/9/6/190
  ###--(4) use unit square instead of torus
  
  # Generate initial set of nodes using Poisson point process on a torus
  #L: dimension
  #L=2
  #volume <- (2 * pi)^(L/2) / gamma(L/2 + 1) # volume of unit sphere in d dimensions
  #intensity <- lambda / volume # intensity of the Poisson process
  #points <- matrix(runif(N*L), ncol=L) # generate uniform random points
  #points <- points - floor(points) # wrap around torus
  points <- matrix(runif(N*2), ncol=2)
  
  # Calculate distance matrix between all pairs of points
  dist_mat <- as.matrix(dist(points))
  
  
  # Set the adjacency matrix element A_{ij} to 1 if dist_{ij} <= r and 0 otherwise
  # r <- 0.2 # adjust this to control the degree of sparsity
  #r=sigma*(log(N)/sqrt(N))
  #sigma=r*sqrt(N)/log(N)
  adj_mat <- as.matrix(dist_mat <=r)
  #excludes self from degree 
  diag(adj_mat) <- 0 
  
  # Initialize the degree vector and the edge list
  adj_mat[1, 1] <- 1
  degree <- colSums(adj_mat)
  edge_list <- matrix(ncol=2, nrow=0)
  
  
  # Grow the network with nodes added one at a time
  for (i in 3:N) {
    # Compute the distance-dependent probability function
    dist_probs <- 1/(1+(dist_mat[1:(i-1), i])^beta) #exp(-sigma*dist_mat[i,])
    # Compute the degree-dependent probability function
    degree <- rowSums(adj_mat[1:(i-1), ])
    deg_probs<- (degree^alpha)/(sum(degree^alpha))
    
    f_ij=mu*dist_probs+(1-mu)*deg_probs
    Q_i= sum(f_ij)
    probattach=f_ij/Q_i
    
    # Compute the attachment probabilities for the new node
    #probattach<- mu * dist_probs[1:(i-1)]*deg_probs +
    # (1 - mu)* rep(sum(dist_probs[1:(i-1)]*deg_probs)/sum(deg_probs), (i-1))
    
    probattach[which(probattach == Inf)] <- 0 # remove infinite probabilities
    
    probattach=probattach/sum(probattach) #normalize attachment probabilities
    
    #Add m edges to existing nodes with probability proportional to their degree and distance
    for (j in 1:m){
      # Choose a node to attach to
      node_to_attach <- sample(1:(i-1), size=m, replace = TRUE, prob=probattach)
      # Add the edge if it doesn't already exist
      if(adj_mat[(i+1)]==0 && node_to_attach == 0) {
        edge_list <- rbind(edge_list, c(i+1, node_to_attach))
        degree[i+1] <- degree[i+1] + 1
        degree[node_to_attach] <- degree[node_to_attach] + 1
        adj_mat[i+1, node_to_attach] <- 1
        adj_mat[node_to_attach, i+1] <- 1
      }
    }
    # Small-world rewiring with probability p
    if (runif(1) < prewire) {
      # Select random node within cutoff distance
      neighbors <- which(adj_mat[i,] == 1)
      if (length(neighbors) > 0) {
        new_neighbor <- sample(neighbors, 1)
        adj_mat[i, new_neighbor] <- 0
        adj_mat[new_neighbor, i] <- 0
        available_nodes <- setdiff(1:N, c(i, neighbors))
        new_neighbor <- sample(available_nodes, 1)
        adj_mat[i, new_neighbor] <- 1
        adj_mat[new_neighbor, i] <- 1
      }
    }
  }
  
  ### Graph object
  GraphModel <- graph.adjacency(as.matrix(adj_mat), mode="undirected")
  GraphModel=igraph::simplify(GraphModel,remove.loops = TRUE,remove.multiple = TRUE)
  
  graph <- list()
  graph$A=adj_mat
  graph$GraphObject <- GraphModel
  
  # Compute graph Laplacian for the expansion property
  D <- diag(rowSums(adj_mat))
  L <- D - adj_mat
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
    adj_mat[adj_mat == 1] <- runif(sum(adj_mat == 1))
    graph$Weights=adj_mat
  }
  
  # Create node attributes if requested
  if (add_node_attr) {
    node_attrs <- data.frame(x = points[,1], y = points[,2])
    graph$NodeAttributes=node_attrs
    graph$degree=degree
  } 
  return(graph)
}      

g=spatial_scale_free_expander(N=100, beta=4, m=1, prewire = 0.1, 
                              add_edge_weight = FALSE,mu=0, 
                              add_node_attr = T,r=0.4,alpha=0)

lo=layout.norm(as.matrix(cbind(g$NodeAttributes$x,g$NodeAttributes$y)))

plot(g$GraphObject,vertex.label=NA,vertex.size=2,layout=lo)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.1] Spatial Expander Propagation Graph with weighted edges, 
# node attributes, scale free degree distribution, small world effect and community structures
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(quadtree)
# Function to generate a spatial scale-free expander graph
# with the given parameters
#
# Arguments:
#   N: number of nodes
#   lambda: intensity parameter for 2D Poisson point process
#   d: dimension of torus
#   r: cutoff distance for adjacency matrix
#   beta: power law parameter for distance-dependent probability
#   m: number of edges added for each new node
#   alpha: parameter controlling strength of degree effect
#   sigma: parameter controlling decay of degree effect with distance
#   mu: parameter controlling preference when connecting nodes. Thus, preference to care about 
#       spatial distance between nodes when connecting them, vs caring only about vertex degree
#   L: number of iterations for Barabasi-Albert model
#   p: rewiring probability for small world effect
#   node_attrs: list of node attributes to add to graph
#   edge_weights: logical indicating whether to add edge weights
#
# Returns:
#   A graph object of class igraph

spatial_scale_free_expander <- function(N, lambda, L, beta, m, prewire = 0.1, 
                                        mu=0.2,add_edge_weight = FALSE, 
                                        add_node_attr = FALSE,r=0.1,alpha=0.2) {
  
  # Set seed for reproducibility
  #  set.seed(123)
  # Generate 2D Poisson point process
  # points <- matrix(rpois(N*d, lambda), N, d)
  # # Create torus by wrapping coordinates around
  #points <- apply(points, 2, function(x) x %% 1)
  # 
  
  # Generate initial set of nodes using Poisson point process on a torus
  #L: dimension
  volume <- (2 * pi)^(L/2) / gamma(L/2 + 1) # volume of unit sphere in d dimensions
  intensity <- lambda / volume # intensity of the Poisson process
  points <- matrix(runif(N*L), ncol=L) # generate uniform random points
  points <- points - floor(points) # wrap around torus
  
  # Generate N points using a 2D Poisson point process with intensity parameter lambda
  #points <- matrix(runif(2*N, min=0, max=L), ncol=2)
  
  # Apply periodic boundary condition to ensure points are distributed on a torus
  #points <- points %% L
  min_distance = r/2
  
  ##---Quad Tree Algorithm to distribute nodes spatially-----------
  # A quadtree is a tree data structure in which each internal node has exactly four children. Quadtrees are the two-dimensional analog of octrees and are most often used to partition a 
  # two-dimensional space by recursively subdividing it into four quadrants or regions. The data associated with a leaf cell varies
  # by application, but the leaf cell represents a "unit of interesting spatial information".
  # The subdivided regions may be square or rectangular, or may have arbitrary shapes. 
  # This data structure was named a quadtree by Raphael Finkel and J.L. Bentley in 1974.[1] A similar partitioning is also known as a Q-tree.
  # All forms of quadtrees share some common features:
  # They decompose space into adaptable cells
  # Each cell (or bucket) has a maximum capacity. When maximum capacity is reached, the bucket splits
  # The tree directory follows the spatial decomposition of the quadtree. 
  
  # n <- dim(points)[1]
  # quad <- list()
  # if (n > 1) {
  #   mid_x <- median(points[,1])
  #   mid_y <- median(points[,2])
  #   quad[[1]] <- quadtree(points[points[,1] <= mid_x & points[,2] <= mid_y,], min_distance)
  #   quad[[2]] <- quadtree(points[points[,1] <= mid_x & points[,2] > mid_y,], min_distance)
  #   quad[[3]] <- quadtree(points[points[,1] > mid_x & points[,2] <= mid_y,], min_distance)
  #   quad[[4]] <- quadtree(points[points[,1] > mid_x & points[,2] > mid_y,], min_distance)
  #   # Check for nodes that are too close and remove them
  #   for (i in 1:4) {
  #     if (is.null(quad[[i]])) {
  #       next
  #     }
  #     if (dim(quad[[i]])[1] == 1) {
  #       continue
  #     }
  #     dist_matrix <- as.matrix(dist(quad[[i]]))
  #     close_nodes <- which(dist_matrix < min_distance & dist_matrix != 0, arr.ind = TRUE)
  #     if (length(close_nodes) > 0) {
  #       for (j in 1:nrow(close_nodes)) {
  #         quad[[i]] <- quadtree(quad[[i]][-close_nodes[j,2],], min_distance)
  #       }
  #     }
  #   }
  # } else {
  #   quad <- points
  # }
  
  # Construct quadtree and set nodes as points that are not too close
  nodes <- do.call(rbind, quad)
  
  # Calculate distance matrix between all pairs of points
  dist_mat <- as.matrix(dist(nodes))
  min_dist <- min(dist_mat[dist_mat > 0]) # Minimum non-zero distance
  
  # Set the adjacency matrix element A_{ij} to 1 if dist_{ij} <= r and 0 otherwise
  # r <- 0.2 # adjust this to control the degree of sparsity
  #r=sigma*(log(N)/sqrt(N))
  sigma=r*sqrt(N)/log(N)
  adj_mat <- as.matrix(dist_mat <=r)
  #excludes self from degree 
  diag(adj_mat) <- 0 
  
  # Initialize the degree vector and the edge list
  degree <- colSums(adjn <- 5
adj_mat <- matrix(0, n, n)

# initialize the upper triangle with 1s
adj_mat[upper.tri(adj_mat)] <- 1

# initialize the lower triangle with 1s
adj_mat[lower.tri(adj_mat)] <- 1

# print the adjacency matrix
adj_mat_mat)
  edge_list <- matrix(ncol=2, nrow=0)
  
  
  # Grow the network with nodes added one at a time
  for (i in 2:N) {
    # Compute the distance-dependent probability function
    dist_probs <- (1/(dist_mat[i,]))^beta * exp(-sigma*dist_mat[i,])
    # Compute the degree-dependent probability function
    deg_probs<-  ((degree[1:(i-1)])^alpha)/(sum(degree[1:(i-1)]^alpha)) }
    
    # Compute the attachment probabilities for the new node
    probattach<- mu * dist_probs[1:(i-1)]*deg_probs +
      (1 - mu)* rep(sum(dist_probs[1:(i-1)]*deg_probs)/sum(deg_probs), (i-1))
    
    probattach[which(probattach == Inf)] <- 0 # remove infinite probabilities
    
    probattach=probattach/sum(probattach) #normalize attachment probabilities
    
    #Add m edges to existing nodes with probability proportional to their degree and distance
    for (j in 1:m){
      # Choose a node to attach to
      node_to_attach <- sample(1:(i-1), size=m, replace = TRUE, prob=probattach)
      # Add the edge if it doesn't already exist
      if(adj_mat[(i+1)]==0 && node_to_attach == 0) {
        edge_list <- rbind(edge_list, c(i+1, node_to_attach))
        degree[i+1] <- degree[i+1] + 1
        degree[node_to_attach] <- degree[node_to_attach] + 1
        adj_mat[i+1, node_to_attach] <- 1
        adj_mat[node_to_attach, i+1] <- 1
      }
    }
    # Small-world rewiring with probability p
    if (runif(1) < prewire) {
      # Select random node within cutoff distance
      neighbors <- which(adj_mat[i,] == 1)
      if (length(neighbors) > 0) {
        new_neighbor <- sample(neighbors, 1)
        adj_mat[i, new_neighbor] <- 0
        adj_mat[new_neighbor, i] <- 0
        available_nodes <- setdiff(1:N, c(i, neighbors))
        new_neighbor <- sample(available_nodes, 1)
        adj_mat[i, new_neighbor] <- 1
        adj_mat[new_neighbor, i] <- 1
      }
    }
  }
  
  ### Graph object
  GraphModel <- graph.adjacency(as.matrix(adj_mat), mode="undirected")
  GraphModel=simplify(GraphModel,remove.loops = TRUE,remove.multiple = TRUE)
  
  graph <- list()
  graph$A=adj_mat
  graph$GraphObject <- GraphModel
  
  # Compute graph Laplacian for the expansion property
  D <- diag(rowSums(adj_mat))
  L <- D - adj_mat
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
    adj_mat[adj_mat == 1] <- runif(sum(adj_mat == 1))
    graph$Weights=adj_mat
  }
  
  # Create node attributes if requested
  if (add_node_attr) {
    node_attrs <- data.frame(x = points[,1], y = points[,2], degree = degree)
    graph$NodeAttributes=node_attrs
  } 
  return(graph)
}      

g=spatial_scale_free_expander(N=100, lambda=10, L=2, beta=4, m=1, prewire = 0.1, 
                              add_edge_weight = FALSE,mu=0.1, 
                              add_node_attr = T,r=0.4,alpha=0)

lo=layout.norm(as.matrix(cbind(g$NodeAttributes$x,g$NodeAttributes$y)))

plot(g$GraphObject,vertex.label=NA,vertex.size=2,layout=lo)





# Function to generate a spatial scale-free expander graph
# with the given parameters
#
# Arguments:
#   N: number of nodes
#   lambda: intensity parameter for 2D Poisson point process
#   d: dimension of torus
#   r: cutoff distance for adjacency matrix
#   beta: power law parameter for distance-dependent probability
#   m: number of edges added for each new node
#   alpha: parameter controlling strength of degree effect
#   sigma: parameter controlling decay of degree effect with distance
#   mu: parameter controlling number of communities
#   L: number of iterations for Barabasi-Albert model
#   p: rewiring probability for small world effect
#   node_attrs: list of node attributes to add to graph
#   edge_weights: logical indicating whether to add edge weights
#
# Returns:
#   A graph object of class igraph

spatial_scale_free_expander <- function(N, lambda, L, beta, m, prewire = 0.1, 
                                        mu=0.2,add_edge_weight = FALSE, 
                                        add_node_attr = FALSE,r=0.1,alpha=0.2) {
  
  # Set seed for reproducibility
  #  set.seed(123)
  
  # Generate N points using a 2D Poisson point process with intensity parameter lambda
  points <- matrix(runif(2*N, min=0, max=L), ncol=2)
  
  # Apply periodic boundary condition to ensure points are distributed on a torus
  points <- points %% L
  
  # Add spatial structure to the point distribution using jittering
  jitter_factor <- 0.1 # adjust this to control the degree of perturbation
  points <- points + matrix(runif(2*N, min=-jitter_factor, max=jitter_factor), ncol=2)
  
  
  # Calculate distance matrix between all pairs of points
  dist_mat <- as.matrix(dist(points))
  min_dist <- min(dist_mat[dist_mat > 0]) # Minimum non-zero distance
  
  # Set the adjacency matrix element A_{ij} to 1 if dist_{ij} <= d_cutoff and 0 otherwise
  # r <- 0.2 # adjust this to control the degree of sparsity
  #exp(-sigma*dist_mat[i,])
  #r=sigma*(log(N)/sqrt(N))
  sigma=r*sqrt(N)/log(N)
  adj_mat <- as.matrix(dist_mat <=r)
  diag(adj_mat) <- 0
  
  # Initialize the degree vector and the edge list
  degree <- colSums(adj_mat)
  #degree <- numeric(N) # Track degree of each node
  #degree[1] <- 1 # Start with a single node
  edge_list <- matrix(ncol=2, nrow=0)
  
  
  #dist_probs <- exp(-dist_matrix[i,] / (min_dist/3)) # Favor short spatial distances
  # pref_probs <- (degrees[1:(i-1)]^(-beta)) # Scale-free preferential attachment
  # Combine probabilities and include community structures
  # attach_probs <- mu*dist_probs[1:(i-1)]*pref_probs + (1-mu)*rep(sum(dist_probs[1:(i-1)]*pref_probs)/sum(pref_probs), (i-1)) 
  # print(attach_probs)
  # have parameter nu that is how much we care about neighborhood:
  # attach_probs <- nu*dist_probs[]... + (1-nu)*rep....
  
  
  # Grow the network with nodes added one at a time
  for (i in 2:N) {
    # Compute the distance-dependent probability function
    dist_probs <- (1/(dist_mat[i,]))^beta * exp(-sigma*dist_mat[i,])
    deg_probs<-  (degree[1:(i-1)]/(sum(degree[1:(i-1)])))^alpha 
    
    # Compute the attachment probabilities for the new node
    probattach<- mu * dist_probs[1:(i-1)]*deg_probs +
      (1 - mu)* rep(sum(dist_probs[1:(i-1)]*deg_probs)/sum(deg_probs), (i-1))
    
    probattach[which(probattach == Inf)] <- 0 # remove infinite probabilities
    #probattach[which(adj_mat[i,] == 1)] <- 0 # remove self-attachment
    
    probattach=probattach/sum(probattach) #normalize
    
    #probattach[which(!is.finite(probattach))] <- 0
    
    #Add m edges to existing nodes with probability proportional to their degree and distance
    for (j in 1:m){
      # Choose a node to attach to
      node_to_attach <- sample(1:(i-1), size=m, replace = TRUE, prob=probattach)
      # Add the edge if it doesn't already exist
      if(adj_mat[(i+1)]==0 && node_to_attach == 0) {
        edge_list <- rbind(edge_list, c(i+1, node_to_attach))
        degree[i+1] <- degree[i+1] + 1
        degree[node_to_attach] <- degree[node_to_attach] + 1
        adj_mat[i+1, node_to_attach] <- 1
        adj_mat[node_to_attach, i+1] <- 1
      }
    }
    # Small-world rewiring with probability p
    if (runif(1) < prewire) {
      # Select random node within cutoff distance
      neighbors <- which(adj_mat[i,] == 1)
      if (length(neighbors) > 0) {
        new_neighbor <- sample(neighbors, 1)
        adj_mat[i, new_neighbor] <- 0
        adj_mat[new_neighbor, i] <- 0
        available_nodes <- setdiff(1:N, c(i, neighbors))
        new_neighbor <- sample(available_nodes, 1)
        adj_mat[i, new_neighbor] <- 1
        adj_mat[new_neighbor, i] <- 1
      }
    }
  }
  
  ### Graph object
  GraphModel <- graph.adjacency(as.matrix(adj_mat), mode="undirected")
  GraphModel=simplify(GraphModel,remove.loops = TRUE,remove.multiple = TRUE)

  graph <- list()
  graph$A=adj_mat
  graph$GraphObject <- GraphModel
  
  # Compute graph Laplacian for the expansion property
  D <- diag(rowSums(adj_mat))
  L <- D - adj_mat
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
    adj_mat[adj_mat == 1] <- runif(sum(adj_mat == 1))
    graph$Weights=adj_mat
  }
  
  # Create node attributes if requested
  if (add_node_attr) {
    node_attrs <- data.frame(x = points[,1], y = points[,2], degree = degree)
    graph$NodeAttributes=node_attrs
    # return(list(adj_matrix = adj_matrix, node_attrs = node_attrs))
  } 
  return(graph)
}      

g=spatial_scale_free_expander(N=100, lambda=10, L=2, beta=4, m=1, prewire = 0.1, 
                              add_edge_weight = FALSE,mu=0.1, 
                              add_node_attr = T,r=0.4,alpha=0)

lo=layout.norm(as.matrix(cbind(g$NodeAttributes$x,g$NodeAttributes$y)))

plot(g$GraphObject,vertex.label=NA,vertex.size=2,layout=lo)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.2] Spatial Expander Propagation Graph with weighted edges, 
# node attributes, scale free degree distribution, small world effect and community structures
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to generate a Spatial Expander Graph
# with the given parameters
#
# Arguments:
#   N: number of nodes
#   lambda: intensity parameter for Poisson point process
#   L: length of torus
#   d: dimension of torus 
#   m: number of nodes to add at each instance
#   beta: power law exponent 
#   p: rewiring probability for small world effect
#   probWithin: probability of connecting nodes within communities
#   probBetween: probability of connecting nodes between communities

# Returns:
#   A Spatial Expander Graph as an igraph object

# growSpatialScaleFreeGraph <- function(N, lambda, L,m, beta, p = 0, probWithin = 0.5, probBetween = 0.1, edgeWeights = FALSE, nodeAttributes = FALSE) {
#   
#   # Generate N points using 2D Poisson point process with intensity lambda on a torus of size L^d
#   d <- 2 # Dimension
#   xy <- matrix(runif(N*d, 0, L), ncol = d) # Generate random coordinates in [0,L]^d
#   xy_jittered <- xy + runif(N*d, -lambda/2, lambda/2) # Jitter points with uniform noise of amplitude lambda
#   
# #Calculate distance matrix between all pairs of points and minimum distance between pairs of nodes
#   dist_matrix <- as.matrix(dist(xy_jittered, method = "euclidean", diag = TRUE, upper = TRUE))
#   min_dist <- min(dist_matrix[dist_matrix > 0]) # Minimum non-zero distance
#   
# #Create adjacency matrix with preferential attachment and small-world rewiring
#   adj_matrix <- matrix(0, nrow = N, ncol = N)
#   degrees <- numeric(N) # Track degree of each node
#   degrees[1] <- 1 # Start with a single node
# 
#   for (i in 2:N) {
#     # Calculate probability of attachment to each existing node
#     dist_probs <- exp(-dist_matrix[i,] / (min_dist/3)) # Favor short spatial distances
#     pref_probs <- (degrees[1:(i-1)]/(sum(degrees[1:(i-1)])))^beta # Scale-free preferential attachment
#     # Combine probabilities and include community structures
#     attach_probs <- mu*dist_probs[1:(i-1)]*pref_probs + (1-mu)*rep(sum(dist_probs[1:(i-1)]*pref_probs)/sum(pref_probs), (i-1)) 
#   print(attach_probs)
#     # have parameter nu that is how much we care about neighborhood:
#   # attach_probs <- nu*dist_probs[]... + (1-nu)*rep....
#   
#     # # normalize probabilities
#     # pa_prob[i] <-  attach_probs / sum(pa_prob[1:i])
#     
#     # Choose nodes to attach to based on the attachment probabilities
#     attach_nodes <- sample(1:(i-1), size=m,replace=T, prob = attach_probs)
#     adj_matrix[i, attach_nodes] <- 1
#     adj_matrix[attach_nodes, i] <- 1
#     degrees[c(i, attach_nodes)] <- degrees[c(i, attach_nodes)] + 1
#     
#     # Small-world rewiring with probability p
#     if (p > 0 && runif(1) < p) {
#       neighbors <- which(adj_matrix[i,] == 1)
#       if (length(neighbors) > 0) {
#         new_neighbor <- sample(neighbors, 1)
#         adj_matrix[i, new_neighbor] <- 0
#         adj_matrix[new_neighbor, i] <- 0
#         available_nodes <- setdiff(1:N, c(i, neighbors))
#         new_neighbor <- sample(available_nodes, 1)
#         adj_matrix[i, new_neighbor] <- 1
#         adj_matrix[new_neighbor, i] <- 1
#       }
#     }
#   }
#   
#   # Compute graph Laplacian for the expansion property
#   D <- diag(rowSums(adj_matrix))
#   L <- D - adj_matrix
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
#   if (expansion_ratio <= 0) {
#     warning("Graph may not exhibit strong expansion properties")
#   }
#   
#   ### Graph object
#   GraphModel <- graph.adjacency(as.matrix(adj_matrix), mode="undirected")
#   GraphModel=simplify(GraphModel,remove.loops = TRUE,remove.multiple = TRUE)
#   graph <- list()
#   graph$A=adj_mat
#   graph$GraphObject <- GraphModel
#   graph$ExpansionRatio=expansion_ratio
#  
#   if (edgeWeights) {
#     adj_matrix[adj_matrix == 1] <- runif(sum(adj_matrix == 1))
#     graph$Weights=adj_matrix
#   }
#   
#   # Create node attributes if requested
#   if (nodeAttributes) {
#     node_attrs <- data.frame(x = xy_jittered[,1], y = xy_jittered[,2], degree = degrees)
#     graph$NodeAttributes=node_attr
#     # return(list(adj_matrix = adj_matrix, node_attrs = node_attrs))
#   } 
#   return(graph)
# }      
# 
# g=growSpatialScaleFreeGraph(N=300, lambda=2, L=2,m=10, beta=0, p = 1, probWithin = 0.5, probBetween = 0.2, edgeWeights = FALSE, nodeAttributes = FALSE) 
# 
# plot(g$GraphObject,vertex.label=NA,vertex.size=2)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.3] Spatial Expander Propagation Graph with weighted edges, 
# node attributes, scale free degree distribution, small world effect and community structures
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to generate a Spatial Expander Graph
# with the given parameters
#
# Arguments:
#   N: number of nodes
#   lambda: intensity parameter for Poisson point process
#   L: length of torus
#   r: adjust this to control the degree of sparsity
#   m: number of nodes to add at each time
#   beta: power law exponent 
#   p: rewiring probability for small world effect
#   probWithin: probability of connecting nodes within communities
#   probBetween: probability of connecting nodes between communities

# Optional and may not be included
# alpha: the parameter that controls the number of edges per node
# (Included): gamma/r: the parameter that controls the strength of the spatial distance effect (r)
# mu: the parameter that controls the strength of the degree effect
# sigma (beta): the parameter that controls the decay of the degree effect with distance

#distances <- dist(poisson_points)
#adjacency_probabilities <- alpha * exp(-r * distances^2) + mu / (1 + distances)^sigma


# Returns:
#   A Spatial Expander Graph as an igraph object
# Function to generate a spatial scale-free expander graph
# spatial_scale_free_expander <- function(N, lambda, L, beta, m, prewire = 0.1, 
#                                         mu=0.2,add_edge_weight = FALSE, 
#                                         add_node_attr = FALSE,r) {
#   
#   # Set seed for reproducibility
# #  set.seed(123)
#   
#   # Generate N points using a 2D Poisson point process with intensity parameter lambda
#   points <- matrix(runif(2*N, min=0, max=L), ncol=2)
#   
#   # Apply periodic boundary condition to ensure points are distributed on a torus
#   points <- points %% L
#   
#   # Add spatial structure to the point distribution using jittering
#   jitter_factor <- 0.1 # adjust this to control the degree of perturbation
#   points <- points + matrix(runif(2*N, min=-jitter_factor, max=jitter_factor), ncol=2)
#   
#   
#   # Calculate distance matrix between all pairs of points
#   dist_mat <- as.matrix(dist(points))
#   min_dist <- min(dist_mat[dist_mat > 0]) # Minimum non-zero distance
#   
#   # Set the adjacency matrix element A_{ij} to 1 if dist_{ij} <= d_cutoff and 0 otherwise
#  # r <- 0.2 # adjust this to control the degree of sparsity
#   adj_mat <- as.matrix(dist_mat <=r)
#   
#   # Initialize the degree vector and the edge list
#   degree <- colSums(adj_mat)
#   edge_list <- matrix(ncol=2, nrow=0)
#   
#   
#   #dist_probs <- exp(-dist_matrix[i,] / (min_dist/3)) # Favor short spatial distances
#  # pref_probs <- (degrees[1:(i-1)]^(-beta)) # Scale-free preferential attachment
#   # Combine probabilities and include community structures
#  # attach_probs <- mu*dist_probs[1:(i-1)]*pref_probs + (1-mu)*rep(sum(dist_probs[1:(i-1)]*pref_probs)/sum(pref_probs), (i-1)) 
#  # print(attach_probs)
#   # have parameter nu that is how much we care about neighborhood:
#   # attach_probs <- nu*dist_probs[]... + (1-nu)*rep....
#   
#   
#   # Grow the network with nodes added one at a time
#   for (i in 1:(N-2)) {
#     # Compute the distance-dependent probability function
#     p_ij <- 1/(dist_mat[,(i+1)])^beta
# #    print(p_ij)
#     # Compute the attachment probabilities for the new node
#     p_i <- (degree/(2*sum(degree)))*(p_ij^prewire)#change p to prewire
#     p_i[which(!is.finite(p_i))] <- 0
#  #   print(p_i)
#    
#     #Add m edges to existing nodes with probability proportional to their degree and distance
#     for (j in 1:m){
#       # Choose a node to attach to
#       node_to_attach <- sample(1:N, size=m, replace = TRUE, prob=p_i)
#       # Add the edge if it doesn't already exist
#       if(adj_mat[(i+1)]==0 && node_to_attach == 0) {
#         edge_list <- rbind(edge_list, c(i+1, node_to_attach))
#         degree[i+1] <- degree[i+1] + 1
#         degree[node_to_attach] <- degree[node_to_attach] + 1
#         adj_mat[i+1, node_to_attach] <- 1
#         adj_mat[node_to_attach, i+1] <- 1
#       }
#     }
#     # Small-world rewiring with probability p
#     if (runif(1) < prewire) {
#       # Select random node within cutoff distance
#       neighbors <- which(adj_mat[i,] == 1)
#       if (length(neighbors) > 0) {
#         new_neighbor <- sample(neighbors, 1)
#         adj_mat[i, new_neighbor] <- 0
#         adj_mat[new_neighbor, i] <- 0
#         available_nodes <- setdiff(1:N, c(i, neighbors))
#         new_neighbor <- sample(available_nodes, 1)
#         adj_mat[i, new_neighbor] <- 1
#         adj_mat[new_neighbor, i] <- 1
#       }
#     }
#   }
#   
#   ### Graph object
#   GraphModel <- graph.adjacency(as.matrix(adj_mat), mode="undirected")
#   GraphModel=simplify(GraphModel,remove.loops = TRUE,remove.multiple = TRUE)
#   graph <- list()
#   graph$A=adj_mat
#   graph$GraphObject <- GraphModel
#   
#   # Compute graph Laplacian for the expansion property
#   D <- diag(rowSums(adj_mat))
#   L <- D - adj_mat
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
# if(is.connected(GraphModel)||expansion_ratio < 0){  
#     warning("Graph may not exhibit strong expansion properties or Graph is disconnected")
#   }
#  
#   graph$ExpansionRatio=expansion_ratio
#   
#   if (add_edge_weight) {
#     adj_mat[adj_mat == 1] <- runif(sum(adj_mat == 1))
#     graph$Weights=adj_mat
#   }
#   
#   # Create node attributes if requested
#   if (add_node_attr) {
#     node_attrs <- data.frame(x = points[,1], y = points[,2], degree = degree)
#     graph$NodeAttributes=node_attrs
#     # return(list(adj_matrix = adj_matrix, node_attrs = node_attrs))
#   } 
#   return(graph)
# }      
# 
# g=spatial_scale_free_expander(N=100, lambda=10, L=2, beta=0, m=1, prewire = 0.5, 
#                                         add_edge_weight = FALSE, 
#                                         add_node_attr = T,r=0.33)
# lo=layout.norm(cbind(g$NodeAttributes$x,g$NodeAttributes$y))
# 
# plot(g$GraphObject,vertex.label=NA,vertex.size=2,layout=lo)





#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.4] Spatial Expander Propagation Graph with weighted edges, 
# node attributes, scale free degree distribution, small world effect and community structures
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# grow_spatial_scale_free_graph <- function(beta, N,lambda, L, p, m, probWithin = 0.5, probBetween = 0.1, add_edge_weights = FALSE, add_node_attributes = FALSE) {
#   # Generate points using a 2D Poisson point process with intensity parameter lambda on a d-dimensional torus
#   coords <- matrix(runif(N*L, 0, L), ncol = L)
#   coords_jittered <- coords + runif(N*L, -0.5, 0.5)/lambda
#   coords_jittered <- coords_jittered %% L
#   
#   # Calculate the distance matrix between all pairs of points
#   dist_mat <- as.matrix(dist(coords_jittered))
#   
#   # Create an adjacency matrix for the graph using preferential attachment
#   adj_mat <- matrix(0, nrow = N, ncol = N)
#   degree <- numeric(N)
#   degree[1] <- 1
#   for (i in 2:N) {
#     # Compute probability of attachment based on spatial distance and preferential attachment
#     prob_attach <- (dist_mat[i,1:(i-1)] + 1e-6) ^ (-beta) * degree[1:(i-1)]
#    # prob_attach[i] <- 0 # exclude self-attachment
#     probsAttach <- prob_attach/sum(prob_attach)
#     x2=probsAttach
#     # Sample m nodes to attach to
#     attach_to <- sample(1:(i-1), size=m, replace = TRUE, prob = x2)
#     
#     # Connect to selected nodes with probability probWithin or probBetween
#     for (j in attach_to) {
#       if (runif(1) < probWithin) {
#         adj_mat[i,j] <- 1
#         adj_mat[j,i] <- 1
#       } else if (runif(1) < probBetween) {
#         available_nodes <- setdiff(1:N, c(i, j, which(adj_mat[i,] == 1), which(adj_mat[j,] == 1)))
#         if (length(available_nodes) > 0) {
#           to_node <- sample(available_nodes, 1)
#           adj_mat[i,to_node] <- 1
#           adj_mat[to_node,i] <- 1
#         }
#       }
#     }
#     
#     # Update degree vector
#     degree[i] <- sum(adj_mat[i,])
#     degree[attach_to] <- degree[attach_to] + 1
#   }
#   
#   # Rewire edges with probability p for small-world effect
#   for (i in 1:(N-1)) {
#     for (j in (i+1):N) {
#       if (adj_mat[i,j] == 1 && runif(1) < p) {
#         available_nodes <- setdiff(1:N, c(i, j, which(adj_mat[i,] == 1), which(adj_mat[j,] == 1)))
#         if (length(available_nodes) > 0) {
#           to_node <- sample(available_nodes, 1)
#           adj_mat[i,j] <- 0
#           adj_mat[j,i] <- 0
#           adj_mat[i,to_node] <- 1
#           adj_mat[to_node,i] <- 1
#         }
#       }
#     }
#   }
#   
#   
#   # Compute graph Laplacian for the expansion property
#   D <- diag(rowSums(adj_mat))
#   L <- D - adj_mat
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
#   if (expansion_ratio <= 0) {
#     warning("Graph may not exhibit strong expansion properties")
#   }
#   
#   ### Graph object
#   GraphModel <- graph.adjacency(as.matrix(adj_mat), mode="undirected")
#   GraphModel=simplify(GraphModel,remove.loops = TRUE,remove.multiple = TRUE)
#   graph <- list()
#   graph$A=adj_mat
#   graph$GraphObject <- GraphModel
#   graph$ExpansionRatio=expansion_ratio
#   
#   if (add_edge_weights) {
#     adj_mat[adj_mat == 1] <- runif(sum(adj_mat== 1))
#     graph$Weights=adj_mat
#   }
#   
#   # Create node attributes if requested
#   if (add_node_attributes) {
#     node_attrs <- data.frame(x =  coords_jittered[,1], y =  coords_jittered[,2], degree = degree)
#     graph$NodeAttributes=node_attr
#     # return(list(adj_matrix = adj_matrix, node_attrs = node_attrs))
#   } 
#   return(graph)
# }      
# 
# g=grow_spatial_scale_free_graph(N=200, lambda=10, L=2,m=1, beta=20, p = 0.5, probWithin = 0.1, probBetween = 0.1, 
#                                 add_edge_weights = FALSE, add_node_attributes = FALSE) 
# 
# plot(g$GraphObject,vertex.label=NA,vertex.size=2)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.5] Spatial Expander Propagation Graph with weighted edges, 
# node attributes, scale free degree distribution, small world effect and community structures
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to generate a Spatial Expander Graph
# with the given parameters
#
# Arguments:
#   N: number of nodes
#   lambda: intensity parameter for Poisson point process
#   L: dimension of torus
#   beta: power law exponent for degree distribution
#   prewire: rewiring probability for small world effect
#   add_communities: boolean indicating whether to add communities or not
#   prob_within: probability of connecting nodes within communities
#   prob_between: probability of connecting nodes between communities
# 
# Returns:
#   A Spatial Expander Graph as an igraph object
# spatial_expander_graph_1 <- function(N, lambda, L, beta, prewire, add_comm = FALSE, 
#                                    prob_within = 0.5, prob_between = 0.1,add_weights=F,
#                                    add_attr=F) {
#   
#   # Generate nodes using Poisson point process
#   coords <- matrix(runif(2 * N), ncol = 2)
#   coords <- coords * L
#   
#   # Calculate distance matrix
#   dist_mat <- as.matrix(dist(coords, diag = TRUE, upper = TRUE))
#   dist_mat <- pmin(dist_mat, L - dist_mat)
#   
#   # Initialize adjacency matrix
#   adj_mat <- matrix(0, nrow = N, ncol = N)
#   
#   # Iterate through each node and connect to other nodes
#   for (i in 1:N) {
#     
#     # Calculate probability of connecting to other nodes w.r.t short spatial distance and high degree
#     d_i <- sort(dist_mat[i, ]) #distance attachment
#     p_i <- d_i^(beta)
#     p_attach <- p_i / sum(p_i) 
#     
#     # Add edges to other nodes with probability proportional to d_i^beta (attachments prob)
#     for (j in 1:N) {
#       if (i != j && runif(1) < p_attach[j]) {
#         adj_mat[i, j] <- 1
#         adj_mat[j,i]<-adj_mat[i, j]
#       }
#     }
#   }
#   
#   # Add rewiring for small world effect
#   for (i in 1:(N-1)) {
#     for (j in (i+1):N) {
#       if (adj_mat[i, j] == 1 && runif(1) < prewire) {
#         adj_mat[i, j] <- 0
#         adj_mat[j, i] <- 0
#         
#         # Choose random long-range neighbor to connect to
#         k <- sample((1:N)[-c(i,j)], size = 1)
#         adj_mat[i, k] <- 1
#         adj_mat[k, i] <- 1
#       }
#     }
#   }
#   
#  # G <- graph.adjacency(as.matrix(adj_mat), mode="undirected")
# #  G=G%>%simplify(remove.loops = TRUE,remove.multiple = TRUE)
# #  plot(G,vertex.label=NA,vertex.size=2)
#  
#    # Add communities if requested
#   if (add_comm) {
#     # Divide nodes into communities
#     n_communities <- ceiling(sqrt(N))
#     comm_size <- ceiling(N / n_communities)
#     communities <- rep(1:n_communities, each = comm_size)[1:N]
#     
#     # Iterate through each node and connect to other nodes within and between communities
#   for (i in 1:(N-1)) {
#     for (j in (i+1):N) {
#         
#         # Check if nodes are in the same community
#         if (communities[i] == communities[j]) {
#           
#           # Connect with probability prob_within
#           if (runif(1) < prob_within) {
#             adj_mat[i, j] <- 1
#             adj_mat[j, i] <- 1
#           }
#           
#         } else {
#           
#           # Connect with probability prob_between
#           if (runif(1) < prob_between) {
#             adj_mat[i, j] <- 1
#             adj_mat[j,i] <-1
#           }
#         }
#       }
#     }
#   } 
#   
#   # Compute graph Laplacian for the expansion property
#   D <- diag(rowSums(adj_mat))
#   L <- D - adj_mat
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
#   if (expansion_ratio <= 0) {
#     warning("Graph may not exhibit strong expansion properties")
#   }
#   
#      ### Graph object
#     GraphModel <- graph.adjacency(as.matrix(adj_mat), mode="undirected")
#     GraphModel=simplify(GraphModel,remove.loops = TRUE,remove.multiple = TRUE)
#     graph <- list()
#     graph$adj_mat=adj_mat
#     graph$GraphObject <- GraphModel
#     graph$ExpansionRatio=expansion_ratio
#     
#     if (add_comm){
#       graph$communities=communities
#     }
#     
#     # Add edge weights and node attributes if specified
#     if (add_weights) {
#       weights <- runif(N*(N-1)/2)
#       adj_mat[upper.tri(adj_mat)] <- weights
#       adj_mat <- adj_mat + t(adj_mat)
#       graph$weights <- weights
#     }
#     if (add_attr) {
#       node_attr <- data.frame(x = coords[,1], y = coords[,2])
#       graph$node_attr <- node_attr
#     }
#   return(graph)
# }  
# 
# g=spatial_expander_graph_1(N=500, lambda=10, L=2,
#                        beta=1, prewire=0.4, add_comm= T, 
#                        prob_within = 0.1, prob_between = 0.2,add_weights=T,
#                        add_attr=T)
# 
# plot(g$GraphObject,vertex.label=NA,vertex.size=2)
# 



# Compute Laplacian matrix
# degree_matrix <- diag(colSums(adj_matrix))
# laplacian_matrix <- degree_matrix - adj_matrix
# 
# # Return graph object
# graph <- list()
# graph$adj_matrix <- adj_matrix
# graph$laplacian_matrix <- laplacian_matrix
# if (add_weights) {
#   graph$weights <- weights
# }
# if (add_attributes) {
#   graph$node_attr <- node_attr
# }
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.6] Spatial Expander Propagation Graph with weighted edges, 
# node attributes, scale free degree distribution, small world effect and community structures
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to generate a Spatial Expander Graph
# with the given parameters
#
# Arguments:
#   N: number of nodes
#   lambda: intensity parameter for Poisson point process
#   L: dimension of torus
#   beta: power law exponent for degree distribution
#   prewire: rewiring probability for small world effect
#   add_communities: boolean indicating whether to add communities or not
#   prob_within: probability of connecting nodes within communities
#   prob_between: probability of connecting nodes between communities
# 
spatial_expander_graph_2 <- function(N, lambda, L, beta, prewire, add_comm = FALSE, 
                                     prob_within = 0.5, prob_between = 0.1,add_weights=F,
                                     add_attr=F) {
  
  # Generate nodes using Poisson point process
  coords <- matrix(runif(2 * N), ncol = 2)
  coords <- coords * L
  
  # Calculate distance matrix
  dist_mat <- as.matrix(dist(coords, diag = TRUE, upper = TRUE))
  dist_mat <- pmin(dist_mat, L - dist_mat)
  
  # Initialize adjacency matrix
  adj_mat <- matrix(0, nrow = N, ncol = N)
  
  # adj_mat[1:2,1:2] <- 1
  # d_i <- sort(dist_mat[i, ]) #distance attachment
  # p_i <- d_i^(beta)
  # p_attach <- p_i / sum(p_i) #degree distribution attachment
  
  for (i in 3:N) {
    for(j in 1:(i-1)){
      # Calculate the probability of connecting to each existing node based on short distance and 
      # scale free degree distribution
      dist_probs <- 1/(dist_mat[i, j]^beta)
      deg_probs <- colSums(adj_mat[1:(i-1),1:(i-1)])
      probs <- dist_probs*deg_probs
      probs_attach <- probs / sum(probs)}}
  # Choose a random node to connect to with probability proportional to the calculated probabilities
  connected_to <- sample(1:N, size = 1, prob = probs_attach)
  adj_mat[i, connected_to] <- 1
  adj_mat[connected_to, i] <- 1
}
# Add rewiring
for (i in 1:(N-1)) {
  for (j in (i+1):N) {
    if (runif(1) < prewire) {
      if (runif(1) < 0.5) {
        adj_mat[i,j] <- 1
        adj_mat[j,i] <- 1
        adj_mat[i, which(adj_mat[j,] == 1)] <- 0
        adj_mat[which(adj_mat[i,] == 1), j] <- 0
      } else {
        adj_mat[j,i] <- 1
        adj_mat[i,j] <- 1
        adj_mat[j, which(adj_matrix[i,] == 1)] <- 0
        adj_mat[which(adj_matrix[j,] == 1), i] <- 0
      }
    }
  }
}

if(add_comm){  
  # Add community structure
  num_communities <- round(sqrt(N))
  community_assignment <- rep(1:num_communities, each = floor(N/num_communities))
  community_assignment <- c(community_assignment, rep(community_assignment[1:(N - length(community_assignment))], each = 1))
  comm_probs <- matrix(rep(probBetween, times = num_communities*num_communities), nrow = num_communities)
  diag(comm_probs) <- probWithin
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (community_assignment[i] == community_assignment[j]) {
        if (runif(1) < comm_probs[community_assignment[i], community_assignment[j]]) {
          adj_mat[i,j] <- 1
          adj_mat[j,i] <- 1
        }
      } else {
        if (runif(1) < comm_probs[community_assignment[i], community_assignment[j]]) {
          adj_mat[i,j] <- 1
          adj_mat[j,i] <- 1
        }
      }
    }
  }
}
# Compute graph Laplacian for the expansion property
D <- diag(rowSums(adj_mat))
L <- D - adj_mat
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
if (expansion_ratio <= 0) {
  warning("Graph may not exhibit strong expansion properties")
}

### Graph object
GraphModel <- graph.adjacency(as.matrix(adj_mat), mode="undirected")
GraphModel=simplify(GraphModel,remove.loops = TRUE,remove.multiple = TRUE)
graph <- list()
graph$adj_mat=adj_mat
graph$GraphObject <- GraphModel
graph$ExpansionRatio=expansion_ratio

if (add_comm){
  graph$communities=communities
}

# Add edge weights and node attributes if specified
if (add_weights) {
  weights <- runif(N*(N-1)/2)
  adj_mat[upper.tri(adj_mat)] <- weights
  adj_mat <- adj_mat + t(adj_mat)
  graph$weights <- weights
}
if (add_attr) {
  node_attr <- data.frame(x = coords[,1], y = coords[,2])
  graph$node_attr <- node_attr
}
return(graph)
}  

g=spatial_expander_graph_2(N=500, lambda=10, L=2,
                           beta=1, prewire=0.4, add_comm= T, 
                           prob_within = 0.1, prob_between = 0.2,add_weights=T,
                           add_attr=T)

plot(graph$GraphObject,vertex.label=NA,vertex.size=2)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.6] Spatial Expander Propagation Graph with weighted edges, 
# node attributes, scale free degree distribution, small world effect and community structures
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(Matrix)
library(igraph)

n=100;L=2;r=0.1;lambda=10;use_edge_weights=FALSE

spatial_expander <- function(L, n, r, beta, p_rewire, attribute) {
  # L: number of dimensions of torus
  # n: number of nodes
  # r: radius of the circle on the torus
  # beta: parameter for preferential attachment
  # p_rewire: probability of rewiring
  # attribute: vector of node attributes
  
  # generate n nodes randomly on the torus
  coords <- matrix(runif(n*L), ncol = L)
  coords <- (2*pi*r) * coords
  coords <- apply(coords, 2, function(x) (cos(x) + 1i*sin(x)))
  coords <- cbind(Re(coords), Im(coords))
  
  # calculate pairwise distances between nodes
  dist_mat <- as.matrix(dist(coords, diag = TRUE, upper = TRUE))
  
  
  # create an empty adjacency matrix
  adj_mat <- matrix(0, nrow = n, ncol = n)
  
  # iterate through each node and connect it to its k-nearest neighbors
  for (i in 1:n) {
    neighbors <- sort(dist_mat[i,])[2:(beta+1)]
    for (j in 1:beta) {
      if (runif(1) < p_rewire) {
        # rewiring with probability p_rewire
        new_neighbor <- sample(1:n, size = 1)
        adj_mat[i, new_neighbor] <- 1
        adj_mat[new_neighbor, i] <- 1
      } else {
        # preferential attachment
        p <- neighbors[j] / sum(neighbors)
        connected_nodes <- which(adj_mat[i,] == 1)
        for (k in connected_nodes) {
          p <- p + beta * (dist_mat[i,k] / sum(dist_mat[i, connected_nodes]))
        }
        p <- p / (beta + degree(graph_from_adjacency_matrix(adj_mat), i))
        if (runif(1) < p) {
          adj_mat[i, j] <- 1
          adj_mat[j, i] <- 1
        }
      }
    }
  }
  
  # create a community structure
  comm <- sample(1:3, size = n, replace = TRUE)
  
  # set node attributes
  attributes <- list(comm = comm, attr = attribute)
  
  # create the graph object
  graph <- list(adj_mat = adj_mat, attributes = attributes)
  
  return(graph)
}

graph <- spatial_expander(L = 2, n = 100, r = 0.7, beta = 1, p_rewire = 0.9, attribute = rnorm(100))
G <- graph.adjacency(as.matrix(graph$adj_mat), mode="undirected")
G=G%>%simplify(remove.loops = TRUE,remove.multiple = TRUE)
plot(G,vertex.label=NA,vertex.size=2)

# generate_spatial_graph <- function(N, lambda, L, beta, p, probWithin, probBetween, add_weights = FALSE, add_attributes = FALSE) {
#   
#   # Generate points using 2D Poisson point process on a torus
#   grid <- expand.grid(replicate(L, seq(0, 1, length.out = sqrt(N)), simplify = FALSE))
#   points <- grid[sample(nrow(grid), N, replace = FALSE), ]
#   
#   # Calculate distance matrix
#   dist_matrix <- as.matrix(dist(points, method = "euclidean", diag = TRUE, upper = TRUE))
#   
#   # Create adjacency matrix using scale-free preferential attachment with rewiring probability
#   adj_matrix <- matrix(0, N, N)
#   for (i in 1:N) {
#     neighbors <- which(dist_matrix[i,] < quantile(dist_matrix[i,], prob = beta))
#     neighbors <- neighbors[-which(neighbors == i)]
#     if (length(neighbors) > 0) {
#       prob <- 1 / (neighbors^beta)
#       prob_attach <- prob / sum(prob)
#       for (j in neighbors) {
#         if (runif(1) < p) {
#           adj_matrix[i, j] <- 1
#           adj_matrix[j, i] <- 1
#         } else {
#           k <- sample(1:N, 1, prob = prob_attach)
#           adj_matrix[i, k] <- 1
#           adj_matrix[k, i] <- 1
#         }
#       }
#     }
#   }
#   
#   # Add community structure if specified
#   if (probWithin > 0 && probBetween > 0) {
#     community <- cutree(cluster::fastgreedy.community(graph.adjacency(adj_matrix)), k = floor(sqrt(N)))
#     for (i in 1:(N-1)) {
#       for (j in (i+1):N) {
#         if (community[i] == community[j]) {
#           if (runif(1) < probWithin) {
#             adj_matrix[i, j] <- 1
#             adj_matrix[j, i] <- 1
#           }
#         } else {
#           if (runif(1) < probBetween) {
#             adj_matrix[i, j] <- 1
#             adj_matrix[j, i] <- 1
#           }
#         }
#       }
#     }
#   }
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+[2.3]
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# spatial_expander <- function(L, n, r, beta, p_rewire, attribute) {
#   # L: number of dimensions of torus
#   # n: number of nodes
#   # r: radius of the circle on the torus
#   # beta: parameter for preferential attachment
#   # p_rewire: probability of rewiring
#   # attribute: vector of node attributes
#   
#   # generate n nodes randomly on the torus
#   # coords <- matrix(runif(n*L), ncol = L)
#   # coords <- (2*pi*r) * coords
#   # coords <- apply(coords, 2, function(x) (cos(x) + 1i*sin(x)))
#   # coords <- cbind(Re(coords), Im(coords))
#   
#   # Generate 2D poisson point process on the L-dimensional torus
#   # coords <- matrix(runif(n*L), n, L)
#   # set the parameters
#   lambda <- 10   # intensity parameter
#   n <- 200       # number of points to generate
#   beta=1
#   side <- n/10     # side length of the square lattice
#   
#   # create a grid of points
#   x <- seq(0, side, length.out = sqrt(n))
#   y <- seq(0, side, length.out = sqrt(n))
#   grid <- expand.grid(x, y)
#   
#   # generate Poisson points
#   num_points <- rpois(1, lambda * side^2)
#   if (num_points > n) {
#     points <- matrix(runif(n * 2), ncol = 2)
#   } else {
#     points <- matrix(runif(num_points * 2), ncol = 2)
#   }
#   
#   # jitter the points
#   points[,1] <- points[,1] + 0.5
#   points[,2] <- points[,2] + 0.5
#   
#   # assign points to grid cells
#   cell_size <- side / sqrt(n)
#   grid$counts <- apply(grid, 1, function(p) {
#   sum(points[,1] >= p[1] & points[,1] < p[1] + cell_size &
#   points[,2] >= p[2] & points[,2] < p[2] + cell_size)
#   })
#   
#   
#   # calculate pairwise distances between nodes
#   dist_mat <- as.matrix(dist(points))
#   
#   # create an empty adjacency matrix
#   adj_mat <- matrix(0, nrow = n, ncol = n)
#   beta=20
#   # iterate through each node and connect it to its k-nearest neighbors
#  for (i in 1:n) {
#     neighbors <- sort(dist_mat[i,])[2:(beta+1)]
#     for (j in 1:beta) {
#        # preferential attachment
#       p <- neighbors[j] / sum(neighbors) #degee attachment
#       connected_nodes <- which(adj_mat[i,] == 1)#connected nodes
#       for (k in connected_nodes) {
#         p1 <- p + beta * (dist_mat[i,k] / sum(dist_mat[i, connected_nodes]))
#         } 
#       p <- p1 / (beta + degree(graph_from_adjacency_matrix(adj_mat), i))
#       if (runif(1) < p) {
#         adj_mat[i, j] <- 1
#         adj_mat[j, i] <- 1
#         }
#     }}
#   
#   # Rewire edges with probability p
#   G <- graph.adjacency(as.matrix(adj_mat), mode="undirected")
#   G=G%>%simplify(remove.loops = TRUE,remove.multiple = TRUE)
#   plot(G,vertex.label=NA,vertex.size=2)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Graph <- spatial_expander_propagation_graph(500, 1, 0.01, 10,deg=2,
#                                             use_edge_weights=T,use_node_attributes=FALSE,
#                                             use_comunity_structures=FALSE,community_groups=4)
# Graph
# for (i in 1:(n-1)) {
#   for (j in (i+1):n) {
#     if (adj_mat[i,j] == 1 && runif(1) < p_rewire) {
#       k <- sample((1:N)[-c(i,j)], 1) # choose a new node uniformly at random
#       adj_mat[i,j] <- 0
#       adj_mat[j,i] <- 0
#       adj_mat[i,k] <- 1
#       adj_mat[k,i] <- 1
#     }
#   }
# }

# #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # [2.1] Spatial Expander Propagation Graph with weighted edges, 
# # node attributes, scale free degree distribution, small world effect and community structures
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # N: the number of points to generate
# # lambda: the intensity parameter of the Poisson point process
# # L: the dimensionality of the torus (should be 2 for a 2D torus)
# # beta: the power law exponent for the preferential attachment degree distribution
# # r: the parameter controlling the strength of the preferential attachment effect
# # p: the rewiring probability for small world effect
# # probWithin: the probability of connecting nodes within communities
# # probBetween: the probability of connecting nodes between communities
# # use_weights: a logical argument indicating whether to use edge weights based on the pairwise distances between points
# 
# Spatial_Expander_Graph <- function(N, lambda, L, beta, r, p, probWithin, probBetween, use_weights = FALSE) {
#   # Generate N points using a 2D Poisson point process with intensity lambda on an L-dimensional torus
#   points <- matrix(runif(N*L), ncol=L)
#   points <- points * L
#   
#   # Calculate distance matrix between all pairs of points
#   dist_mat <- as.matrix(dist(points, diag=TRUE, upper=TRUE))
#   
#   # Create adjacency matrix with preferential attachment probability
#   adj_mat <- matrix(0, nrow=N, ncol=N)
#   for (i in 1:(N-1)) {
#     for (j in (i+1):N) {
#       if (runif(1) < probWithin) {
#         # preferential attachment within communities
#         if (runif(1) < r^(-beta * dist_mat[i,j])) {
#           adj_mat[i,j] <- 1
#           adj_mat[j,i] <- 1
#         }
#       } else {
#         # preferential attachment between communities
#         if (runif(1) < r^(-beta * dist_mat[i,j])/probBetween) {
#           adj_mat[i,j] <- 1
#           adj_mat[j,i] <- 1
#         }
#       }
#     }
#   }
#   
#   # Add rewiring probability p for small world effect
#   for (i in 1:(N-1)) {
#     for (j in (i+1):N) {
#       if (adj_mat[i,j] == 1) {
#         if (runif(1) < p) {
#           adj_mat[i,j] <- 0
#           adj_mat[j,i] <- 0
#           k <- sample(setdiff(1:N, c(i,j)), 1)
#           adj_mat[i,k] <- 1
#           adj_mat[k,i] <- 1
#           adj_mat[j,k] <- 1
#           adj_mat[k,j] <- 1
#         }
#       }
#     }
#   }
#   
#   # Compute Laplacian matrix
#   deg_vec <- apply(adj_mat, 1, sum)
#   deg_mat <- diag(deg_vec)
#   lap_mat <- deg_mat - adj_mat
#   
#   # Return graph object with desired properties
#   if (use_weights) {
#     graph_obj <- graph_from_adjacency_matrix(adj_mat, weighted=TRUE)
#     E(graph_obj)$weight <- 1/dist_mat[adjacent_vertices(graph_obj, mode="both")]
#   } else {
#     graph_obj <- graph_from_adjacency_matrix(adj_mat,mode = "undirected")
#   }
#   return(list(graph_obj=graph_obj, lap_mat=lap_mat, dist_mat=dist_mat))
# }
# 
# g=Spatial_Expander_Graph(N=500, lambda=10,
#                        L=1, beta=2, r=0.1, p=0.1, probWithin=0,
#                        probBetween=0.004, use_weights = FALSE)
# 
# plot(g$graph_obj,vertex.label=NA,vertex.size=2)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.1] Spatial Expander Propagation Graph with weighted edges, 
# node attributes, scale free degree distribution, small world effect and community structures
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# L: dimensions of torus
# N: number of nodes (points)
# r: distance contraintparameter
# beta: parameter for preferential attachment
# p: probability of rewiring
# probWithin: within community connections
# probBetween: between community connections
# attribute: vector of node attributes
#useWeightedEdges: flag to use edge weights or not
# library(Matrix)
# SpatialExpanderGraph <- function(N, lambda, L, beta, r, p, probWithin, probBetween, useWeightedEdges = FALSE) {
#   
#   # Generate points using a 2D Poisson point process on a torus
#   xy <- matrix(runif(N * L, 0, 1), N, L)
#   lambda_eff <- lambda / (2 * pi)^L # Effective intensity for torus
#   interpoint_distances <- as.matrix(dist(xy))
#   
#   # Create adjacency matrix based on spatial short distances and
#   # preferential attachment with small-world effect
#   adj_matrix <- matrix(0, N, N)
#   
#   for (i in 1:N) {
#     deg <- colSums(adj_matrix) # Calculate degrees of existing nodes
#     dists <- interpoint_distances[i,] # Distances to other nodes
#     
#     # Calculate attachment probabilities
#     if (useWeightedEdges) {
#       probs <- (dists ^ (-beta)) * (deg + 1) / sum((dists ^ (-beta)) * (deg + 1))
#     } else {
#       probs <- (1 / (dists ^ r)) * (deg + 1) ^ beta / sum((1 / (dists ^ r)) * (deg+ 1) ^ beta)
#     }
#     
#     # Connect to new node with preferential attachment
#     neighbors <- sample(1:N, size = 1, prob = probs)
#     adj_matrix[i, neighbors] <- 1
#     adj_matrix[neighbors, i] <- 1
#   
# 
# #    Small-world effect with rewiring probability p
#     if (runif(1) < p) {
#       non_neighbors <- setdiff(1:N, c(i, neighbors))
#       new_neighbor <- sample(non_neighbors, size = 1)
#       adj_matrix[i, new_neighbor] <- 1
#       adj_matrix[new_neighbor, i] <- 1
#       adj_matrix[i, neighbors] <- 0
#       adj_matrix[neighbors, i] <- 0
#     }
#   }
#   
#   # Add community structure
#   community_labels <- cutree(graph.adjacency(adj_matrix, mode = "undirected"), k = sqrt(N))
#   for (i in 1:(sqrt(N) - 1)) {
#     for (j in (i + 1):sqrt(N)) {
#       within_community <- community_labels == i & community_labels == j
#       between_communities <- community_labels == i & community_labels != j | community_labels != i & community_labels == j
#       if (runif(1) < probWithin) {
#         adj_matrix[within_community, within_community] <- 1
#       }
#       if (runif(1) < probBetween) {
#         adj_matrix[between_communities, within_community] <- 1
#         adj_matrix[within_community, between_communities] <- 1
#       }
#     }
#   }
#   
#   # Compute Laplacian matrix and eigenvalues
#   degree_matrix <- diag(colSums(adj_matrix))
#   laplacian_matrix <- degree_matrix - adj_matrix
#   eigenvalues <- eigen(laplacian_matrix)$values
#   
#   # Return adjacency matrix and eigenvalue gap
#   return(list(adj_matrix = adj_matrix, eigenvalue_gap = eigenvalues[2]))
# }
# 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.2] Spatial Expander Propagation Graph with weighted edges, 
# node attributes, scale free degree distribution, small world effect and community structures
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# spatial_expander <- function(L, n, r, beta, p_rewire, attribute) {
#   # L: number of dimensions of torus
#   # n: number of nodes
#   # r: radius of the circle on the torus
#   # beta: parameter for preferential attachment
#   # p_rewire: probability of rewiring
#   # attribute: vector of node attributes
#   
#   # generate n nodes randomly on the torus
#   coords <- matrix(runif(n*L), ncol = L)
#   coords <- (2*pi*r) * coords
#   coords <- apply(coords, 2, function(x) (cos(x) + 1i*sin(x)))
#   coords <- cbind(Re(coords), Im(coords))
#   
#   # calculate pairwise distances between nodes
#   dist_mat <- as.matrix(dist(coords, diag = TRUE, upper = TRUE))
#   
#   # create an empty adjacency matrix
#   adj_mat <- matrix(0, nrow = n, ncol = n)
#   
#   # iterate through each node and connect it to its k-nearest neighbors
#   for (i in 1:n) {
#     neighbors <- sort(dist_mat[i,])[2:(beta+1)]
#     for (j in 1:beta) {
#       if (runif(1) < p_rewire) {
#         # rewiring with probability p_rewire
#         new_neighbor <- sample(1:n, size = 1)
#         adj_mat[i, new_neighbor] <- 1
#         adj_mat[new_neighbor, i] <- 1
#       } else {
#         # preferential attachment
#         p <- neighbors[j] / sum(neighbors)
#         connected_nodes <- which(adj_mat[i,] == 1)
#         for (k in connected_nodes) {
#           p <- p + beta * (dist_mat[i,k] / sum(dist_mat[i, connected_nodes]))
#         }
#         p <- p / (beta + degree(graph_from_adjacency_matrix(adj_mat), i))
#         if (runif(1) < p) {
#           adj_mat[i, j] <- 1
#           adj_mat[j, i] <- 1
#         }
#       }
#     }
#   }
#   
#   # create a community structure
#   comm <- sample(1:3, size = n, replace = TRUE)
#   
#   # set node attributes
#   attributes <- list(comm = comm, attr = attribute)
#   
#   # create the graph object
#   graph <- list(adj_mat = adj_mat, attributes = attributes)
#   
#   return(graph)
# }
# 
# graph <- spatial_expander(L = 1, n = 300, r = 0.9, beta = 1, p_rewire = 1, attribute = rnorm(100))
# G=graph_from_adjacency_matrix(graph$adj_mat,mode="undirected")
# G=G%>%simplify(remove.loops = TRUE,remove.multiple = TRUE)
# plot(G,vertex.size=2,vertex.label=NA)
# 


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.2] Spatial Expander Propagation Graph with weighted edges, 
# node attributes, scale free degree distribution, small world effect and community structures
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# N: the number of nodes in the graph
# L: the dimensionality of the torus
# lambda: the mean of the Poisson point process used to generate the nodes
# alpha: the parameter that controls the number of edges per node
# beta: the parameter that controls the strength of the spatial distance effect
#** mu: the parameter that controls the strength of the degree effect
#** sigma: the parameter that controls the decay of the degree effect with distance
# rewiring_prob: the probability of rewiring an edge
# K:number of community groups
# use_edge_weights: whether to use edge weights or not
# q - proportion of edges within communities compared to between communities
# community_prob: Probability of within-community edge connection
# node_attributes - a matrix of size n x k containing node attributes
# Returns: a sparse adjacency matrix representing the graph
# 
# Required packages
# library(Matrix)
# library(mvtnorm)
# library(Matrix.utils) 
# 
# # Function to generate a spatial expander propagation graph
# spatial_expander_graph1 <- function(N, L, lambda, beta,
#                       rewiring_prob,use_edge_weights = TRUE,
#                       q = 0.8, node_attributes = NULL) {
#   
#   # Step 1: Generate Poisson point process
#   coords <- rmvnorm(N, mean = rep(0, L), sigma = diag(rep(1/lambda, L)))
#   coords <- coords %% 1 # L-dimensional torus
#   
#   # Step 2: Compute the euclidean (pairwise) distances between points
#   num_coords <- nrow(coords)
#   dist_mat <- matrix(0, nrow=num_coords, ncol=num_coords)
#   for (i in 1:(num_coords-1)) {
#     for (j in (i+1):num_coords) {
#       dist_mat[i,j] <- sqrt(sum((coords[i,] - coords[j,])^2))
#       dist_mat[j,i] <- dist_mat[i,j]
#     }
#   }
#   
#   # Calculate distances and probabilities
#   #distances <- dist_mat
#   #adjacency_probabilities <- alpha * exp(-beta * distances^2) + mu / (1 + distances)^sigma
#   
#   # Step 3: Create adjacency matrix based on preferential attachment (high degree) and short spatial distance
#   adj_mat <- matrix(0, nrow = N, ncol = N)
#   # Initialize node degrees
#   deg <- rep(0, N)
#   #Step 4: Construct edges with preferential attachment and within-community connection probability
#   # Add edges between nodes with probability proportional to short spatial distance and high degree nodes
#   for (i in 1:(N-1)) {
#     for (j in (i+1):N) {
#       dist <- sqrt(sum(( coords[i,] -  coords[j,] - L*round(( coords[i,] -  coords[j,])/L))^2))
#       prob <- exp(-dist) * (deg[i] + deg[j] + 1)
#       if (runif(1) < prob) {
#         adj_mat[i,j] <- 1
#         adj_mat[j,i] <- adj_mat[i,j]
#         deg[i] <- deg[i] + 1
#         deg[j] <- deg[j] + 1
#       }
#     }
#   }
#  

# for (i in 1:N) {
#   # Find nearest neighbors
#  nn_idx <- order(dist_mat[i,])[2:(beta+1)]# nearest/shortest distance of nodes
#  nn_degrees <- colSums(adj_mat[nn_idx,])#degree of nodes
#   # Preferentially attach to high degree neighbors and short spatial distance
#  new_edges <- sample(nn_idx, size=beta, replace=TRUE, prob=nn_degrees)
#   # Connect to nodes in same community with higher probability
#   for (j in new_edges) {
#     if (runif(1) < community_prob || all(nn_degrees == 0)) {
#       adj_mat[i,j] <- 1
#       adj_mat[j,i] <- 1
#     }
#   }
# }

# #Step 5: Rewire edges with small world probability
# for (i in 1:(N-1)) {
#   for (j in (i+1):N) {
#     if (adj_mat[i,j] == 1) {
#       if (runif(1) < rewire_prob) {
#         # Rewire edge
#         k <- sample(1:N, size=1)
#         while (k == i || k == j || adj_mat[i,k] == 1) {
#           k <- sample(1:N, size=1)
#         }
#         adj_mat[i,j] <- 0
#         adj_mat[j,i] <- 0
#         adj_mat[i,k] <- 1
#         adj_mat[k,i] <- 1
#         if (use_edge_weights) {
#           # Add edge weights proportional to inverse distance
#           dist_ij <- dist_mat[i,j]
#           dist_ik <- dist_mat[i,k]
#           w_ij <- 1 / dist_ij
#           w_ik <- 1 / dist_ik
#           adj_mat[j,i] <- w_ij
#           adj_mat[k,i] <- w_ik
#         } else {
#           adj_mat[j,k] <- 1
#           adj_mat[k,j] <- 1
#         }
#       }
#     }
#   }
# }

#graph=graph_from_adjacency_matrix(adjacency_matrix)
# ## Step 6: Create community structures
#   communities <- stats::cutree(cluster::walktrap.community(as.undirected(graph_from_adjacency_matrix(adj_mat))), k=2)
#   prob_within <- q / sum(communities == 1) + (1-q) / sum(communities == 2)
#   prob_between <- (1-q) / sum(communities == 1) + q / sum(communities == 2)
#   for (i in 1:nrow(adj_mat)) {
#     for (j in 1:ncol(adj_mat)) {
#       # Check if nodes i and j belong to the same community
#       if (communities[i] == communities[j]& adj_matrix[i, j] == 0) {
# 
#         # If the probability of a within-community edge is less than the threshold, add an edge
#         if (runif(1) < prob_within && adj_mat[i, j] == 0) {
#           adj_mat[i, j] <- 1
#           adj_mat[j, i] <- 1
#         }
#       } else {
# 
#        
#         # If the probability of a between-community edge is less than the threshold and the number of edges within each community is greater than the number of edges between communities, add an edge
#         if (runif(1) < prob_between && adj_mat[i, j] == 0 && sum(adj_mat[communities == communities[i], communities == communities[j]]) > sum(adj_mat[communities == communities[i],]) / num_communities && sum(adj_mat[community == community[i], communities == communities[j]]) > sum(adj_mat[communities==communities[j],]) / num_communities) {
#         {
#           adj_mat[i, j] <- 1
#           adj_mat[j, i] <- 1
#         }
#       }
#     }
#   }
#   
#     
#   #Step 7: Compute graph Laplacian for the expansion property
#   D <- diag(rowSums(adj_mat))
#   L <- D - adj_mat
#   
#   # Step 8: Check for strong expansion properties
#   spectral_gap <- eigen(Matrix::t(L) %*% L, symmetric=TRUE, only.values=TRUE)$values[n-1]
#   if (spectral_gap <= 0) {
#     warning("Graph may not exhibit strong expansion properties")
#   }
#   
#   # Add node attributes
#   if (!is.null(attributes)) {
#     node_attributes <- attributes
#   }
#   
#   # Step 9: Convert Adjacency to igraph object
#   G <- graph.adjacency(as.matrix(adj_mat), mode="undirected")
#   G=G%>%simplify(remove.loops = TRUE,remove.multiple = TRUE)
#  
#   # Return the graph as a list
#   return(list(adjacency_matrix = adj_mat,
#               node_attributes = node_attributes,
#               expansion_gap = expansion_gap))
# } 
# 
# plot(G,vertex.label=NA,vertex.size=2)



# Add node attributes
# for (i in 1:N) {
#  neighbors <- sample(1:N, size = round(alpha * N), prob = adjacency_probabilities[i, ], replace = TRUE)
#   if (use_edge_weights) {
#     edge_weights <- rexp(length(neighbors), 1)
#     adjacency_matrix[i, neighbors] <- edge_weights
#   } else {
#     adjacency_matrix[i, neighbors] <- 1
#   }
# }
#   
#   # Step 5: Rewiring
#   for (i in 1:nrow(adjacency_matrix)) {
#     for (j in 1:ncol(adjacency_matrix)) {
#       if (runif(1) < rewiring_prob) {
#         k <- sample(1:N, size = 1)
#         if (use_weights) {
#           adjacency_matrix[i, j] <- 0
#           adjacency_matrix[i, k] <- rexp(1, 1)
#         } else {
#           adjacency_matrix[i, j] <- 0
#           adjacency_matrix[i, k] <- 1
#         }
#       }
#     }
#   }
# 
#   
#   #graph=graph_from_adjacency_matrix(adjacency_matrix)
#   ## Step 6: Create community structures
#   communities <- cutree(cluster::walktrap.community(as.undirected(graph_from_adjacency_matrix(adjacency_matrix))), k = 2)
#   prob_within <- q / sum(communities == 1) + (1-q) / sum(communities == 2)
#   prob_between <- (1-q) / sum(communities == 1) + q / sum(communities == 2)
#   for (i in 1:nrow(adjacency_matrix)) {
#     for (j in 1:ncol(adjacency_matrix)) {
#       # Check if nodes i and j belong to the same community
#       if (community[i] == community[j]) {
#         
#         # If the probability of a within-community edge is less than the threshold, add an edge
#         if (runif(1) < prob_within && adjacency_matrix[i, j] == 0) {
#           adjacency_matrix[i, j] <- 1
#           adjacency_matrix[j, i] <- 1
#         }
#       } else {
#         
#         # If the probability of a between-community edge is less than the threshold, add an edge
#         if (runif(1) < prob_between && adjacency_matrix[i, j] == 0) {
#           adjacency_matrix[i, j] <- 1
#           adjacency_matrix[j, i] <- 1
#         }
#       }
#     }
#   }



# spatial_expander_propagation_graph <- function(n, L, r, lambda,deg=4,
#                                                use_edge_weights=FALSE,use_node_attributes=FALSE,
#                                                use_comunity_structures=FALSE,
#                                                community_groups=4) {
#   # Generate 2D poisson point process on the L-dimensional torus
#   coords <- matrix(runif(n*L), n, L)
#   radius_sq <- r^2
#   
#   # Compute the euclidean (pairwise) distances between points
#   num_coords <- nrow(coords)
#   dist_mat <- matrix(0, nrow=num_coords, ncol=num_coords)
#   for (i in 1:(num_coords-1)) {
#     for (j in (i+1):num_coords) {
#       dist_mat[i,j] <- sqrt(sum((coords[i,] - coords[j,])^2))
#       dist_mat[j,i] <- dist_mat[i,j]
#     }
#   }
#   
#   # Compute weights based on pairwise distances (optional)
#   if (use_weights) {
#     weights <- exp(-dist_matrix^2 / r^2)
#     adj_matrix <- as.matrix(weights)
#   } else {
#     adj_matrix <- as.matrix(dist_matrix <= r)
#   }
# 
# 
# 
# Graph <- spatial_expander_propagation_graph(500, 1, 0.01, 10,deg=2,
#                                             use_edge_weights=T,use_node_attributes=FALSE,
#                                             use_comunity_structures=FALSE,community_groups=4)
# Graph
# plot(Graph,vertex.label=NA,vertex.size=2)



# for (i in 1:(n-1)) {
#   for (j in (i+1):n) {
#     neighbors <- which(dist_mat[i,] < r)
#     neighbors <- neighbors[neighbors != i]
#     if (length(neighbors) > deg) {
#       neighbors <- sample(neighbors, deg)
#     }
#     adj_mat[i, neighbors] <- 1
#     adj_mat[neighbors,i] <- adj_mat[i, neighbors]
#     
#     # Compute shortest distance between points on torus
#     dist_sq <- sum(pmin(abs(coords[i,]-coords[j,]), 1-abs(coords[i,]-coords[j,]))^2)
#     #dist_sq=dist_mat
#     if (dist_sq <= radius_sq) {
#       # Assign edge weight as inverse of distance
#       W[i,j] <- 1/sqrt(dist_sq)
#       W[j,i] <- W[i,j]
#     }
#     
#     #check if the adjacency matrix has edge weights or not and select which matrix to use for graph
#     if (sum(adj_mat) == 0) {
#       print("No edge weights in the adjacency matrix")
#       # Convert adjacency matrix to igraph object
#       G <- graph_from_adjacency_matrix(as.matrix(W),mode="undirected")
#       
#     } 
#     else if (sum(adj_mat) != 0) {
#       print("The adjacency matrix has edge weights")
#       # Convert adjacency matrix to igraph object
#       G <- graph_from_adjacency_matrix(as.matrix(W),mode="undirected")
#     }
#   }
#   }
# return(adj_mat)
#} 

# Generate a spatial expander propagation graph with 500 nodes on a 4D torus

#  # Create adjacency matrix with edge weights
#    W <- Matrix(0, n, n)
#    for (i in 1:(n-1)) {
#      for (j in (i+1):n) {
#        # Compute shortest distance between points on torus
#        dist_sq <- sum(pmin(abs(coords[i,]-coords[j,]), 1-abs(coords[i,]-coords[j,]))^2)
#        if (dist_sq <= radius_sq) {
#          # Assign edge weight as inverse of distance
#          W[i,j] <- 1/sqrt(dist_sq)
#          W[j,i] <- W[i,j]
#          adj_mat=W
#        }
#      }
#    }
#  




# Convert adjacency matrix to igraph object
#G <- graph_from_adjacency_matrix(as.matrix(W),mode="undirected")

# Add node attributes
#if (node.attribute==True){
# V(G)$coords <- coords
#  V(G)$community <- cutree(cluster_walktrap(G))
#}
#} 
# Check for strong expansion properties
#spectral_gap <- eigen(Matrix::t(W) %*% W, symmetric=TRUE, only.values=TRUE)$values[n-1]
#if (spectral_gap <= 0) {
# warning("Graph may not exhibit strong expansion properties")
#}

#   return(adj_mat)
# }




# Plot the graph with node colors based on community structure
#plot(G, vertex.color=V(G)$community)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sample Code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#library(Matrix)
#library(spatstat)

# Function to construct a spatial expander propagation graph on an L-dimensional torus
# Arguments:
# - n: number of nodes
# - L: torus dimension
# - r: radius of connection (shortest spatial distance)
# - lambda: intensity parameter for Poisson point process
# Returns:
# - an adjacency matrix representing the spatial expander propagation graph




#Generate Poisson point process on the torus
#torus <- owin(rep(0, L), rep(1, L))
# pp <- rpoispp(lambda, win=torus)
# points <- pp$marks

# Compute distances between points
#dist_mat <- as.matrix(nndist(points), torus)

#  adjacency <- matrix(0,n,n)
# for (i in 1:n) {
#   neighbors <- which(distances[i,] < r)
#   neighbors <- neighbors[neighbors != i]
#   if (length(neighbors) > deg) {
#     neighbors <- sample(neighbors, deg)
#   }
#   adjacency[i,neighbors] <- 1
# }


# spatial_expander_propagation_graph <- function(n, L, r,deg,lambda) {  
#   # Generate the 2D Poisson point process
#   set.seed(123)
#   lambda <- n^(1/L) #intensity parameter
#   points <- matrix(runif(n*L), ncol=L)
#   points <- points - floor(points) #wrap around torus edges
#   points <- points[runif(n) < lambda, ]
#   
#   # Compute the euclidean (pairwise) distances between points
#   dist_mat <- as.matrix(dist(points))
#   
#   # Construct adjacency matrix
#   adj_mat <- matrix(0, n, n)
#   for (i in 1:n) {
#     neighbors <- which(dist_mat[i,] < r)
#     neighbors <- neighbors[neighbors != i]
#     if (length(neighbors) > deg) {
#       neighbors <- sample(neighbors, deg)
#     }
#     adj_mat[i, neighbors] <- 1
#    # adj_mat[neighbors,i] <- 1
#   }
#   diag(adj_mat) <- 0
#   adj_mat <- as(adj_mat, "dgCMatrix")
#   expander_r <- (2*L)^2*log(n)/(r^L)
#   while(min(eigen(adj_mat, symmetric = TRUE)$values) <= expander_r) {
#     # Remove edges randomly until we get an expander graph
#     edges <- which(adj_mat != 0, arr.ind = TRUE)
#     num_edges <- length(edges[,1])
#     if (num_edges == 0) {
#       return(adj_mat)
#     }
#     to_remove <- sample(num_edges, size = 1)
#     adj_mat[edges[to_remove,]] <- 0
#   }
#   return(adj_mat)
# } 
# 
# 
# 
# adj_mat <- spatial_expander_propagation_graph(n = 5, L = 2,deg=4, r = 0.1, lambda = 10)
# 
#   # Ensure graph is connected and an expander graph
#   
#   



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Spatial Expandergraph
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Explanation of the code:
# We define a function torus_ppp that generates points on the L-dimensional torus using a 
# 2D Poisson point process with intensity lambda. We use the runif function to generate 
# uniform random numbers in the unit cube, and then map these points to the torus using 
# the modulo operator (%%).
# We define a function dist_torus that computes the pairwise distances between points on
# the torus. We use a nested loop to compute the distance between every pair of points, 
# taking into account the periodic boundary conditions of the torus.
# We define a function spatial_expander_propagation_torus that constructs a spatial 
# expander propagation graph on the L-dimensional torus. This function takes three 
# arguments: lambda (the intensity of the Poisson point process), 
# L (the dimension of the torus), and r (the distance threshold for connecting two points).
# The function first generates points on the torus using torus_ppp, then computes 
# the pairwise distances between points using dist_torus, and constructs an 
# adjacency matrix based on the distance threshold r. The function then 
# checks that the resulting graph has strong expansion properties 
# (i.e., every subset of vertices of size less than lambda/2 has at least lambda/2 
# neighbors), and returns the graph as a list with three elements: 
# the adjacency matrix (adj), the number of vertices


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [1.1] Spatial Expander graph
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#library(Matrix)

#sp.exp.graph=function(n=100,d=2,lambda=10){
# Define the number of nodes and the dimension of the torus

# n : number of nodes
# d : dimension of torus
# lambda : intensity parameter of ppp

# Generate the Poisson point process on the torus
# 
# X <- matrix(runif(n*d, 0, 1), ncol=d)
# for (i in 1:d) {
#   X[,i] <- X[,i] * (2*pi)
# }
# Y <- matrix(rpois(n, lambda), ncol=1)
# nodes <- matrix(0, nrow=sum(Y), ncol=d)
# index <- 1
# for (i in 1:n) {
#   for (j in 1:Y[i]) {
#     nodes[index,] <- X[i,]
#     index <- index + 1
#   }
# }
# 
# # Compute the distances between nodes on the torus
# distances <- matrix(0, nrow=n, ncol=n)
# for (i in 1:(n-1)) {
#   for (j in (i+1):n) {
#     d1 <- abs(nodes[i,] - nodes[j,])
#     d2 <- 2*pi - d1
#     distances[i,j] <- sqrt(sum(pmin(d1, d2)^2))
#     distances[j,i] <- distances[i,j]
#   }
# }
# 
# # Compute the Cheeger constant
# vol <- sum(Y)
# vol_left <- 0
# vol_right <- vol
# perm <- order(Y, decreasing=TRUE)
# h <- Inf
# for (i in 1:(n-1)) {
#   vol_left <- vol_left + Y[perm[i]]
#   vol_right <- vol_right - Y[perm[i]]
#   cut_size <- sum(distances[perm[1:i], perm[(i+1):n]])
#   conductance <- cut_size / min(vol_left, vol_right)
#   if (conductance < h) {
#     h <- conductance
#     perm_star <- perm[1:i]
#   }
# }
# 
# # Construct the adjacency matrix of the expander graph
# A <- matrix(0, nrow=n, ncol=n)
# for (i in 1:n) {
#   for (j in perm_star) {
#     if (i != j) {
#       d1 <- abs(nodes[i,] - nodes[j,])
#       d2 <- 2*pi - d1
#       d <- sqrt(sum(pmin(d1, d2)^2))
#       if (d <= h/2) {
#         A[i,j] <- 1
#         A[j,i] <- 1
#       }
#     }
#   }
# }
# return(A)
# }
# 
# adjMatrix=sp.exp.graph(n=500,d=1,lambda=10)
# A=adjMatrix
# # Construct the combinatorial Laplacian matrix of the expander graph
# D <- diag(rowSums(A))
# L <- D - A
# 
# # Check the expansion properties of the expander graph
# lambda_2 <- eigen(L, symmetric=TRUE, only.values=TRUE)$values[2]
# if (lambda_2 > 0.1*h^2) {
#   print("The expander graph has strong expansion properties.")
# } else {
#   print("The expander graph does not have strong expansion properties.")
# }
# 
# G=graph_from_adjacency_matrix(A, mode = "undirected")
# plot(G,vertex.label=NA, vertex.size=2,edge.arrow.size=0.01,
#      edge.color = "darkgrey", edge.width = 1,
#      edge.labels = NA)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [1.2] Spatial Expander graph
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#1.Generate random points on the torus using a 2D Poisson point process. 
#We can use the spatstat package in R to generate these points. Specifically, 
#we can use the rpoispp() function to generate a Poisson point pattern on a 
#rectangular grid, and then transform the points to lie on the torus using the torusify()
#function.

#2.Calculate the pairwise distances between the points using the Euclidean distance metric.

#3.Construct the adjacency matrix of the graph by connecting each point to its
#k nearest neighbors, where k is a parameter that determines the degree of the graph.
#To do this, we can sort the distances for each point and take the first 
#k distances as the distances to its k nearest neighbors. We then set the 
#corresponding entries in the adjacency matrix to 1.

#4.To ensure that the graph has strong expansion properties, 
#we can use a spectral partitioning algorithm to partition the graph into two 
#approximately equal-sized sets. To do this, we can compute the eigenvectors 
#corresponding to the k smallest eigenvalues of the Laplacian matrix of the graph,
#and then use these eigenvectors to perform k-means clustering on the points. 
#The two resulting clusters will be our partition of the graph

# library(Matrix)
# 
# SpatModel1=function(d=2,n=300,r=0.5,deg=4,epsilon=0.2){
# 
# # Define the parameters
# # d:  dimension of torus
# # n:  number of vertices
# # r:  radius for neighbour search and connection with neighbour (threshold that controls the connection) 
# # deg: degree for each vertex
# # epsilon: mixing parameter
# 
# # Generate the 2D Poisson point process
# set.seed(123)
# lambda <- n^(1/d) #intensity parameter
# points <- matrix(runif(n*d), ncol=d)
# points <- points - floor(points) #wrap around torus edges
# points <- points[runif(n) < lambda, ]
# 
# # Compute the euclidean (pairwise) distances between points
# distances <- as.matrix(dist(points))
# # coords <- coords %% 1 # wrap around torus
# # coords <- coords * (1-2*radius) + radius # shift to avoid boundary effects
# 
# # Generate the adjacency matrix and connects points to k-nearest neighbours
# #connect two points if their Euclidean distance is less than or equal to r
# adjacency <- matrix(0,n,n)
# for (i in 1:n) {
#   neighbors <- which(distances[i,] < r)
#   neighbors <- neighbors[neighbors != i]
#   if (length(neighbors) > deg) {
#     neighbors <- sample(neighbors, deg)
#   }
#   adjacency[i,neighbors] <- 1
# }
# 
# # Ensure the graph is connected whiles repeating the step above for every node
# # and connecting them to their k nearest neighbours
# #while (!is_connected(graph_from_adjacency_matrix(adjacency))) {
#   for (i in 1:n) {
#     neighbors <- which(adjacency[i,] == 1)
#     if (length(neighbors) < deg) {
#       candidates <- setdiff(1:n, c(i,neighbors))
#       if (length(candidates) > deg - length(neighbors)) {
#         candidates <- sample(candidates, deg - length(neighbors))
#       }
#       adjacency[i,candidates] <- 1
#     }
#   }
# #}
# return(adjacency)
# }
# 
# spgraph=SpatModel1(d=2,n=300,r=0.5,deg=3,epsilon=0.2)
# 
# g=graph_from_adjacency_matrix(spgraph, mode = "undirected")
# plot(g,vertex.label=NA, vertex.size=2,edge.arrow.size=0.01,
#      edge.color = "darkgrey", edge.width = 1,
#      edge.labels = NA)
# 
# 
# # # Compute the degree matrix
#  degMatrix <- diag(colSums(adjacency))
# # 
# # # Compute the Laplacian matrix
# laplacian <- degMatrix- adjacency
# # 
# # # Compute the eigenvalues and eigenvectors of the Laplacian matrix
#  eig <- eigen(laplacian)
# eigenvalues <- eig$values
#  eigenvectors <- eig$vectors
# # 
# # # Compute the expander graph
#  graph <- eigenvectors[,2:nrow(degMatrix+1)] %*% diag(sqrt(eigenvalues[2:nrow(degMatrix+1)]))
#  graph.1 <- graph/sqrt(graph)
# # 
# # # Add noise for strong expansion properties
#  noise <- matrix(rnorm(n*deg), nrow=n, ncol=deg)
#  graph.2 <-  graph.1 + 1:nrow(noise * epsilon)
# # 
# # # Ensure the graph is normalized
# graph.3 <- Matrix(as.numeric(graph.2), sparse=TRUE)
# graph.4 <- graph.3 / norm(graph.3, "F")
# z=as.matrix(graph.4)
# # Visualize the graph
# g = graph_from_adjacency_matrix(z, mode = "undirected")



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Spatial Expander graphs with scale free degree distribution, small world
# and community structures with expansion properties (simple example)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(spatstat)
library(igraph)

# Parameters
N <- 1000 # Number of points
m <- 12 # Number of edges to attach to each new vertex in Barabasi-Albert model
p <- 0.3 # Probability of adding a random edge
ncommunities <- 13 # Number of communities
community_size <- 100 # Number of vertices in each community

# Spatial distribution of points
x <- runifpoint(N)

#add.vertices(N,g)

##---Make a Barabasi-Albert random graph and rewire with rewiring prop p for small world effect
# g <- make_empty_graph(N,directed = FALSE)
g1=barabasi.game(N,power=2,m, directed = FALSE,algorithm="psumtree")%>%
  rewire(each_edge(p = .1, loops = FALSE)) 
#degree <- degree(g1)

# Divide vertices into communities
communities <- sample(rep(1:ncommunities, each = community_size), N, replace = TRUE)

##Connect vertices within each community
within_community_prob <- 0.2
for (i in 1:ncommunities) {
  community_vertices <- which(communities == 3)
  within_community_edges <- sample(community_vertices, size = length(community_vertices) * (length(community_vertices) - 1) * within_community_prob / 2, replace = T)
  g2 <- add.edges(g1, within_community_edges, verbose = FALSE)
}

# Compute the first k eigenvectors of the graph Laplacian
k <- 15
# Compute the adjacency matrix
A <- get.adjacency(g2)
# Compute the laplacian matrix
L <- laplacian_matrix(g2)
eigenvectors <- eigen(L)$vectors[, 1:k]

# Perform spectral partitioning (with the number of communities/clusters) on the eigenvectors
spectral_partition <- kmeans(eigenvectors, ncommunities)
between_community_edges <- spectral_partition$cluster
# Add edges between different clusters (comunities) to create an expander graph
g3 <- add.edges(g2, between_community_edges, verbose = FALSE)%>%
  simplify(remove.multiple = T,remove.loops = T)
# Visualize the graph
plot.igraph(g3, vertex.col = communities,
            vertex.label=NA, vertex.size=2,edge.arrow.size=0.01,
            edge.color = "darkgrey", edge.width = 1,
            edge.labels = NA)


# # Add edges between different clusters (comunities) to create an expander graph
# for (i in 1:nrow(A)) {
#   for (j in (i+1):ncol(A)) {
#     if (A[i,j] == 0 && clusters$cluster[i] != clusters$cluster[j]) {
#       add_edges(g2, c(i,j))
#     }
#   }
# }

# cutsize <- 0.1 * N * (ncommunities - 1)
# adjacency_matrix <- as.matrix(as_adjacency_matrix(g2))
# spectral_partition <- spantree.partition(adjacency_matrix, cutsize)
# between_community_edges <- spectral_partition$cut
# g3 <- add.edges(g2, between_community_edges, verbose = FALSE)

###--Add edges between vertices in different communities to create an expander graph
###--First, calculate the normalized Laplacian matrix of the graph
# L <- laplacian_matrix(g2)
# D <- diag(vcount(g2))
# Dinv <- solve(D)
# Lnorm <- Dinv %*% L #normalize laplacian matrix
# Then, use spectral partitioning to find a cut that minimizes edge expansion
# spectral_partition <- spectral_partition(Lnorm, K)
# cut_edges <- get.cut_edges(g2, spectral_partition$membership)
# # Finally, add the cut edges to the graph
# g3 <- add.edges(g2, cut_edges,verbose = FALSE)




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Spatial Expander Propagation graph Model
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Goal:To Construct a Spatial Expander Propagation Graph that has the following properties.

###--[1]--The graph should have expansion properties of an expander graph on an L-dimesional torus.

###--[2]--The graph is grown on the torus and  each node added to the network is given 
# a spatial position that follows a 2D poisson point process and are connected to existing 
# vertices with a probability favouring short spatial distances and high degrees. 

###--[3]--This graph should have properties of scale free degree distribution, small world and
# community structures and should be constructed without igraph and be built from scratch.

###--[4]--This graph should be a single function where we can tune few parameters to 
#generate different types of models 
#(eg graphs with larger network sizes, sparse and dense connections,random graphs etc)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Creating Helper functions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# library(Matrix)
# library(RANN)
# 
# # Set parameters for the spatial point process
# n_points <- 100 # number of nodes to add
# lambda <- 1 # density of the Poisson point process
# 
# # Set parameters for the graph growth
# p_degree <- 0.1 # probability of connecting to existing vertices based on degree
# p_distance <- 0.9 # probability of connecting to existing vertices based on spatial distance
# r_distance <- 0.1 # radius of spatial neighborhood for connecting based on distance
# 
# ### Generate spatial positions of nodes using Poisson point process with intensity parameter lambda
# poisonn_pp=function(lambda,xmax,ymax){
#   n=rpois(1,lambda*xmax*ymax)
#   x=runif(n,0,xmax)
#   y=runif(n,0,ymax)
#   return(data.frame(x,y))
# }
# 
# poisonn_pp(2,3,3)
# 
# # # Generate spatial positions of nodes using Poisson point process on torus
# # x <- runif(n_points)
# # y <- runif(n_points)
# # xy <- cbind(x, y)
# # theta <- 2*pi*(x+y)
# # xy_torus <- cbind(cos(theta), sin(theta))
# 
# # Initialize adjacency matrix and degrees vector
# adj_mat <- Matrix(0, n_points, n_points)
# deg <- rep(0, n_points)
# 
# # Add nodes to the graph one at a time
# for (i in 1:n_points) {
#   # Connect to existing vertices based on degree and distance
#   neighbors <- which(deg > 0)
#   if (length(neighbors) > 0) {
#     distances <- RANN::nn2(xy[i,], xy[neighbors,], k = 1)$nn.dist
#     p_distance_neigh <- exp(-distances/r_distance)
#     p_degree_neigh <- deg[neighbors]/sum(deg[neighbors])
#     p_neighbors <- p_degree*p_degree_neigh + p_distance*p_distance_neigh
#     neighbors <- sample(neighbors, size = 5, prob = p_neighbors)
#     adj_mat[i,neighbors] <- 1
#     adj_mat[neighbors,i] <- 1
#   }
#   
#   # Add node to graph and update degrees
#   deg[i] <- sum(adj_mat[i,])
# }



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Algorithm: Steps
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  

# 1.Initialize the graph with a small number of vertices randomly distributed on the torus.

# 2.Add new vertices to the graph one by one, following a 2D Poisson point process.
# 
# 3.Connect the new vertex to k existing vertices, where k is determined by the
# degree distribution of the graph.
# 
# 4.Choose the existing vertices to connect to by giving preference to vertices that are 
# close to the new vertex in terms of spatial distance and have high degrees.
# 
# 5.Repeat steps 2-4 until the desired number of vertices is reached.
# 
# To ensure that the resulting graph has properties of a scale-free degree distribution,
# small world, and community structures, the following modifications can be
# made to the algorithm:
#   
# i)Use a preferential attachment mechanism to assign degrees to new vertices.
# Specifically, new vertices should be connected to existing vertices with a
# probability proportional to their degree.
# 
# ii)Introduce random rewiring of edges to promote small-worldness. 
# Specifically, with a small probability, existing edges can be 
# rewired to connect to vertices that are farther away in terms of 
# spatial distance but have higher degrees.
# 
# iii)Use community detection algorithms to identify and extract communities within the graph.
# 
# To implement the algorithm, you can use the R programming language and the packages spatstat,
# igraph, and community.

# Function to generate a scale-free network
generate_scale_free_network <- function(N, alpha, m) {
  # Create an initial connected network with m+1 nodes
  G <- matrix(0, nrow=N, ncol=N)
  G[1:(m+1), 1:(m+1)] <- 1
  diag(G) <- 0
  
  # Create the degree sequence of the initial network
  deg_seq <- rep(m+1, m+1)
  
  # Add nodes to the network
  for (i in (m+2):N) {
    # Calculate the degree distribution
    p <- deg_seq / sum(deg_seq)
    
    # Select m nodes to connect to
    targets <- sample(1:(i-1), m, prob=p, replace=TRUE)
    
    # Connect the new node to the selected nodes
    G[i,targets] <- 1
    G[targets,i] <- 1
    
    # Update the degree sequence
    deg_seq[i] <- m
    deg_seq[targets] <- deg_seq[targets] + 1
  }
  
  # Calculate the degree distribution of the final network
  deg_dist <- degree(G)
  
  # Check if the degree distribution follows a power law
  fit <- power.law.fit(deg_dist)
  if (fit$alpha < alpha) {
    # If the power law parameter is less than alpha, try again
    G <- generate_scale_free_network(N, alpha, m)
  }
  
  return(G)
}


