setwd("C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi-New")
#source("SPATIAL-PIPELINE-NEW-MODEL.R")
library(igraph)
library(Matrix)
library(tidyverse)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [1] Spatial Expander Propagation Graph with weighted edges, 
# node attributes, scale free degree distribution, small world effect and community structures
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to generate a spatial scale-free graph
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
  #+++++++++++++++++++++++++++++++++++++++++++
  #+quad tree part goes here. This is where I 
  #will call the c++ function
  #+++++++++++++++++++++++++++++++++++++++++++
  points <- matrix(runif(N*2), ncol=2)
  
  
  ##--calculate distance matrix between all pairs of points
  dist_mat <- as.matrix(dist(points))
  ##--Set up graph object
  graph <- list()
  graph$nodes <- 1:m
  ##--initialize the edge_list and adjacency matrix
  graph$edge_list <- data.frame(from=NULL, to=NULL)
  graph$adj_mat <- matrix(0, nrow = N, ncol = N)
  ##--initialize graph with initial set of nodes fully connected to each other
  graph$adj_mat[1:m,m:1]<-1
  diag(graph$adj_mat) <- 0 #exclude self connections for initialialized nodes
  ##--initialize the degree vector
  graph$degrees=rep(m,N)#degree sequence have degree of at least m
  #graph$degrees[1:m]=m-1# Initial set of nodes have degree m-1
  
  ##--Grow the network with new nodes added one at a time using scale-free technique
  for (i in (m+1):N) {
    #preferential attachment probability
    #graph$degrees= rowSums(graph$adj_mat[1:(i-1), ])
    deg_probs<- graph$degrees[1:(i-1)]^alpha/sum(graph$degrees[1:(i-1)]^alpha)
    #spatial attachment probability
    spat_dist_effect <- ifelse(dist_mat[1:(i-1), i] <= r, 1/(1+(dist_mat[1:(i-1), i])^beta), 
                               1+((dist_mat[1:(i-1), i])-N)^beta)
    
    spatial_probs <- spat_dist_effect /sum(spat_dist_effect )#normalise
    print(spatial_probs)
    # total probability of attachment which is normalised
    TotalprobsOfAttach <-(mu*spatial_probs)+((1-mu)*deg_probs)
    #Add m edges between new node to existing nodes with probability proportional to their degree and distance
    targetnodes<- sample(1:(i-1), size = m, replace = TRUE, prob =  TotalprobsOfAttach)
    
    #Attach nodes
    graph$nodes <- c(graph$nodes, i)
   # graph$adj_mat <- rbind(graph$adj_mat, matrix(0, ncol = N, nrow = 1))
    graph$adj_mat[i,  targetnodes] <- 1
    graph$adj_mat[targetnodes, i] <- 1
    
    #Update degree vectors
    graph$degrees[targetnodes] <- graph$degrees[targetnodes]+1
    graph$degrees[i] <- m
    
    #update edge_list
    new_edges <- data.frame(from=i, to=targetnodes)
    graph$edge_list<- rbind(graph$edge_list, new_edges)
    
    #Small-world rewiring with probability prewire
    if (runif(1) < prewire) {
      #Select random node to rewire within cutoff distance
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

  ##--Graph object
  GraphModel <- graph.adjacency(as.matrix(graph$adj_mat), mode="undirected")
  GraphModel=igraph::simplify(GraphModel,remove.loops = TRUE,remove.multiple = TRUE)
  graph$GraphObject <- GraphModel

  ##--Create edge weights if requested
  if (add_edge_weight) {
    graph$adj_mat[graph$adj_mat == 1] <- runif(sum(graph$adj_mat == 1))
    graph$Weights=graph$adj_mat
  }

  ##--Create node attributes if requested
  if (add_node_attr) {
    node_attrs <- data.frame(x = points[,1], y = points[,2])
    graph$NodeAttributes=node_attrs
  }
  return(graph)
}

g=spatial_scale_free_expander_noncomm (N=100, beta=0, m=3, prewire = 0, 
                                           add_edge_weight = FALSE,mu=1, 
                                           add_node_attr = T,r=0.9,alpha=3)
plot(g$GraphObject,vertex.label=NA,vertex.size=2)
g$degrees


### Testing with Scale Free barabassi model
SFF=sample_pa(100,power = 1.5, m=4,directed = F,algorithm = c("psumtree"))
SF=igraph::simplify(SFF,remove.loops = TRUE,remove.multiple = TRUE)
degree(SF)
plot(SF,vertex.label=NA,vertex.size=2)
#hist(igraph::degree(SF), breaks=20, main="")

# lo=layout.norm(as.matrix(cbind(g$NodeAttributes$x,g$NodeAttributes$y)))
# plot(g$GraphObject,vertex.label=NA,vertex.size=2,layout=lo)

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


# spatial_scale_free_expander_comm <- function(N=50, beta=2, m=1, prewire = 0.1, numofComm=4, 
#                                              mu=0.2,add_edge_weight = FALSE,
#                                              Prob_Betw_comm=0.1,add_node_attr = FALSE,r=0.1,alpha=0.2) {
#   
#   
# 
#   ##--Generate spatial points on a unit square
#   points <- matrix(runif(N*2), ncol=2)
#   # calculate distance matrix between all pairs of points
#   dist_mat <- as.matrix(dist(points))
#   # Set up graph ogject
#   graph <- list()
#   graph$nodes <- 1:m
#   # initialize the edge_list and adjacency matrix
#   graph$edge_list <- data.frame(from=NULL, to=NULL)
#   graph$adj_mat <- matrix(0, nrow = N, ncol = N)
#   # initialize graph with an initial node and no edges
#   graph$adj_mat[1]<-1
#   # initialize the degree vector
#   #graph$degrees <- numeric(N)
#   graph$degrees=rep(m,N)
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
#    # graph$degrees<- rowSums(graph$adj_mat[1:(i-1), ])
#     deg_probs <- (graph$degrees[1:(i-1)]^alpha) / sum( graph$degrees[1:(i-1)]^alpha)
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
#     targetnodes <- sample(1:(i-1), m, replace = TRUE, prob = P_join)
#     graph$adj_mat[i,targetnodes] <- 1
#     graph$adj_mat[targetnodes, i] <- 1
#     graph$degrees[targetnodes] <- graph$degrees[targetnodes] + 1
#     graph$degrees[i] <- m
#     
#     #update edge_list
#     new_edges <- data.frame(from=i, to=targetnodes)
#     graph$edge_list<- rbind(graph$edge_list, new_edges)
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
#   GraphModel<-igraph::simplify(GraphModel,remove.loops = TRUE,remove.multiple = TRUE)
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
#     node_attrs <- data.frame(x =   points[,1], y =   points[,2])
#     graph$NodeAttributes=node_attrs
#   } 
#   return(graph)
# }      
# 
# g=spatial_scale_free_expander_comm(N=100, beta=0, m=2, prewire = 0, numofComm=1, 
#                                    mu=0,add_edge_weight = FALSE,
#                                    Prob_Betw_comm=0.0001,add_node_attr = T,r=0,alpha=2) 
# 
# lo=layout.norm(as.matrix(cbind(g$NodeAttributes$x,g$NodeAttributes$y)))
# 
# 
# plot(g$GraphObject,vertex.label=NA,vertex.size=2)
# 
# plot(g$GraphObject,vertex.label=NA,vertex.size=2,layout=lo)
# g$degrees
# 



###++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The probability attachment function based on the given properties can be defined as follows:
#   
#   Let p(u,v) be the probability of an edge between nodes u and v, then
# 
# p(u,v) = f(u,v)/sum(w(x,v)) where x belongs to N(v) and N(v) is the set of all nodes adjacent to v, and w(x,v) = f(x,v) if d(x,v) <= r, else w(x,v) = alpha * f(x,v) where alpha is a constant such that sum(w(x,v)) <= 1.
# 
# Here, we are using a normalization factor to ensure that the sum of probabilities for all adjacent nodes of v is equal to 1.
# 
# The function f(u,v) is given by:
#   
#   If d(u,v) <= r, then f(u,v) = 1/(1+(d(u,v)^beta))
# 
# If d(u,v) > r, then f(u,v) = c/(1+(d(u,v)^beta))
# 
# where c is a constant such that f(u,v) < 1/(1+(d(u,v)^beta)) when d(u,v) > r.
# 
# The step function for f(u,v) is as follows:
#   
# For d(u,v) <= r, f(u,v) is a monotonically decreasing function of d(u,v) with a maximum value of 1.
# 
# For d(u,v) > r, f(u,v) is a monotonically decreasing function of d(u,v) with a maximum value of c/(1+(r^beta)).
# 
# The value of beta determines the rate at which f(u,v) decreases with distance. A larger value of beta means that the edge probability decreases more rapidly with distance.
# 
# The cutoff distance parameter 'r' determines the distance beyond which the edge probability decreases more rapidly with distance, and the constant 'c' is chosen such that the maximum value of f(u,v) is less than 1/(1+(r^beta)).
