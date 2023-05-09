# Function to simulate a scale-free network using preferential attachment
# Arguments:
#   n: number of nodes in the network
#   m: number of edges to attach from a new node to existing nodes
# Returns:
#   adjacency matrix of the network
scale_free_pa <- function(n, m) {
  # Initialize the network with a single node
  network <- matrix(0, nrow = n, ncol = n)
  network[1:2, 2:1] <- 1
  
  # Loop through each new node and attach edges
  for (i in 3:n) {
    # Calculate the degree distribution of the existing nodes
    degrees <- rowSums(network[1:(i-1), ])
    probs <- degrees / sum(degrees)
    
    # Select m nodes to attach to with probability proportional to degree
    new_edges <- sample(1:(i-1), size = m, replace = TRUE, prob = probs)
    
    # Add the new edges to the network
    network[i, new_edges] <- 1
    network[new_edges, i] <- 1
  }
  return(network)
}

# Example usage: simulate a scale-free network with 100 nodes and 3 edges attached from each new node
network <- scale_free_pa(100, 2)
g=graph_from_adjacency_matrix(network,mode="undirected")
g=igraph::simplify(g,remove.multiple = T,remove.loops = T)
plot(g, vertex.label=NA,vertex.size=2)
degree(g)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Initialize an empty graph with one node
N <- 100  # Number of nodes
m <- 3    # Number of edges to add at each instance
alpha <- 1 # Power law or preferential attachment parameter


# Add edges to the graph with preferential attachment
sf=function(N=100,m=3,alpha=1){
  graph <- matrix(0, nrow = N, ncol = N)
  graph[1,1]=1
  
for (i in 3:N) {
  # Choose m target nodes proportional to their degree
  deg=rowSums(graph[1:(i-1),])
  deg=deg^alpha/sum(deg^alpha)
  
  targets <- sample(1:(i-1), size = m, replace = TRUE, prob = deg)
  
  # Add edges from the new node to the target nodes
  graph[i,targets] <- 1
  graph[targets,i] <- 1
}
  return(graph)
}
# Calculate the degree distribution and plot it
grph=sf(N=100,m=2,alpha=6)
grph=igraph::graph_from_adjacency_matrix(grph,mode="undirected")
g.graph=igraph::simplify(grph,remove.multiple = T,remove.loops = T)
plot(g.graph, vertex.label=NA,vertex.size=2)

degree(g.graph)
#degrees <- colSums(graph)
hist(degrees, breaks = 30, freq = FALSE, main = "Degree distribution", xlab = "Degree")

# Fit a power law distribution to the degree distribution
library("poweRlaw")
fit <- displ$new(degrees)
fit$setXmin(1)
fit$setPars(list(alpha = alpha))
fit$setVerbose(FALSE)
fit$estimate()
fit$plot()


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+# Define the number of nodes and edges to add at each instance
N <- 100 # number of nodes
m <- 2 # number of edges to add at each instance
alpha <- 2 # preferential attachment parameter (power law exponent)

# Initialize the graph with m nodes and m edges
sfg = function(N=100,m=3,alpha=1){
  
edges <- matrix(0, nrow = m, ncol = 2)
for (i in 1:m) {
  edges[i,] <- c(1,i+1)
}

adj_mat <- matrix(0, nrow = m+1, ncol = m+1)
for (i in 1:m) {
  adj_mat[edges[i,1], edges[i,2]] <- 1
  adj_mat[edges[i,2], edges[i,1]] <- 1
}

# Iterate through remaining nodes and add edges preferentially
for (i in (m+1):N) {
  # Calculate the degree distribution of existing nodes
  degree_dist <- rowSums(adj_mat[1:(i-1), ])
  
  # Calculate the preferential attachment probabilities for each node
  pa_probs <- (degree_dist/ sum(degree_dist))^alpha
  
  # Choose m nodes to connect to
  new_edges <- sample(1:(i-1), m, replace = T, prob = pa_probs)}
  
  # Add the new edges to the adjacency matrix
  adj_mat[i, new_edges] <- 1}
  adj_mat[new_edges, i] <- 1
  }
}
return(adj_mat)
}
# Plot the degree distribution of the resulting graph
# Calculate the degree distribution and plot it
grph=sfg(N=100,m=2,alpha=6)

grph=igraph::graph_from_adjacency_matrix(grph,mode="undirected")
g.graph=igraph::simplify(grph,remove.multiple = T,remove.loops = T)
plot(g.graph, vertex.label=NA,vertex.size=2)


degree_dist <- rowSums(adj_mat)
hist(degree_dist, breaks = seq(0,max(degree_dist), by = 1), main = "Degree distribution")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to generate scale-free network
preferential_attachment <- function(N, alpha, m) {
  # Initialize the graph with m nodes
  nodes <- 1:m
  edges <- matrix(0, nrow = m, ncol = N)
  edges[lower.tri(edges)] <-1
  #edges[upper.tri(edges)]<-1
  edges <- edges + t(edges)# inititialise adj matrix
  degree <- rowSums(edges)
  
  # Grow the network
  for (i in (m+1):N) {
    # Calculate the probability of attaching to each existing node
    prob <- (degree^alpha) / sum(degree^alpha)
    # Choose m nodes to attach to based on their probability
    attach_to <- sample(nodes, m, replace = TRUE, prob = prob)
    # Add the new node and edges to the graph
    #nodes <- c(nodes, i)
    #edges <- rbind(edges, rep(0, (i-1)))
    #edges <- cbind(edges, c(rep(0, i), 1))
    edges[i, attach_to] <- 1
    edges[attach_to, i] <- 1
    # Update the degree of each node
    # degree <- rowSums(edges)
  }
  
  # Return the graph as an adjacency matrix
  return(edges)
}

gg=preferential_attachment(N=50, alpha=1, m=3)

# This function takes three parameters: N (the number of nodes in the final network), 
# alpha (the power law or preferential attachment parameter), 
# and m (the number of edges to add at each step). It initializes
# the graph with m nodes and m(m-1)/2 edges,
# and then grows the network by adding one node at a 
# time and m edges to attach to existing nodes.
# 
# The algorithm calculates the probability of attaching to each existing node
# based on its degree raised to the power of alpha. It then chooses m nodes to 
# attach to based on their probabilities, and adds the new node and edges to the graph. 
# Finally, it updates the degree of each node and repeats the process until 
# the graph has N nodes.
#The function returns the final graph as an adjacency matrix. 
#You can visualize the graph using the igraph package or any other graph visualization tool.



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Define the parameters
N <- 100# number of nodes
alpha <- 6 # power law or preferential attachment parameter
m <- 4 # number of edges to add at each instance

# Create a vector to store the degrees of the nodes
degrees <- numeric(N)

# Create the initial graph with two connected nodes
degrees[1:2] <- m
edges <- matrix(c(1, 2), ncol = 2, byrow = TRUE)

# Loop over the remaining nodes
for (i in 3:N) {
  # Compute the probability of attaching to each existing node
  probabilities <- degrees[1:(i-1)]^alpha / sum(degrees[1:(i-1)]^alpha)
  
  # Choose m nodes to attach to with replacement, according to the probabilities
  targets <- sample(1:(i-1), size = m, replace = TRUE, prob = probabilities)
  
  # Update the degrees of the chosen nodes
  degrees[targets] <- degrees[targets] + 1
  degrees[i]=m
  
  # Add the new edges to the edge matrix
  new_edges <- matrix(c(i, targets), ncol = 2, byrow = TRUE)
  edges <- rbind(edges, new_edges)
}

# Plot the degree distribution
hist(degrees, breaks = seq(0.5, max(degrees) + 0.5, by = 1),
     col = "steelblue", xlab = "Degree", ylab = "Frequency", main = "Degree Distribution")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Define function to generate the network
preferential_attachment <- function(N, alpha, m) {
  # Create initial fully-connected graph with m nodes
  nodes <- 1:m
  edges <- t(combn(nodes, 2))
  # Create vector to hold the degrees of each node
  degrees <- rep(m, m)
  
  # Iterate through remaining nodes and add edges
  for (i in (m+1):N) {
    # Compute probability distribution for attachment
    prob <- (degrees[edges[,1]]^alpha) * (degrees[edges[,2]]^alpha)
    prob <- prob / sum(prob)
    # Sample m edges to attach to new node
    new_edges <- sample(edges, m, prob = prob, replace = TRUE)
    # Add new node and edges to network
    nodes <- c(nodes, i)
    edges <- rbind(edges, cbind(rep(i, m), new_edges[,2]))
    degrees <- c(degrees, rep(m, m))
    degrees[new_edges] <- degrees[new_edges] + 1
  }
  
  # Return adjacency matrix of network
  adj_matrix <- matrix(0, N, N)
  adj_matrix[edges] <- 1
  return(adj_matrix)
}

preferential_attachment(100,6, 4)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to generate scale-free preferential attachment network
# with parameters N, alpha, and m
scale_free_network <- function(N, alpha, m) {
  
  # Create an initial connected network of m+1 nodes
  nodes <- 1:(m)
  edges <- matrix(nodes, ncol=2)
  edges <- edges[lower.tri(edges)]
  graph <- igraph::graph_from_edgelist(as.matrix(edges,ncol=2))
  
  # Calculate initial degree sequence
  degree_seq <- igraph::degree(graph)
  
  # Loop to add remaining N-m-1 nodes
  for (i in (m+2):N) {
    
    # Calculate probability of connecting to each existing node
    prob <- (degree_seq^alpha) / sum(degree_seq^alpha)
    
    # Choose m nodes to connect to with preferential attachment
    new_edges <- numeric()
    for (j in 1:m) {
      chosen_node <- sample(1:length(prob), 1, prob=prob)
      new_edges <- c(new_edges, chosen_node)
      degree_seq[chosen_node] <- degree_seq[chosen_node] + 1
    }
    
    # Add new node and its edges to graph
    nodes <- c(nodes, i)
    edges <- rbind(edges, c(rep(i, m), new_edges))
    degree_seq <- c(degree_seq, m)
    
    #graph <- igraph::graph_from_edgelist(edges)
  }
  
  return(edges)
}

# Generate a scale-free network with 100 nodes, alpha=2, and m=2
graph <- scale_free_network(100, 2, 2)

# Plot the degree distribution
hist(igraph::degree(graph), breaks=20, main="")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Define parameters
N <- 100 # Number of nodes
alpha <- 2 # Power law parameter
m <- 2 # Number of edges to add at each step

# Initialize the network
edges <- matrix(data = NA, nrow = 0, ncol = 2) # Empty adjacency matrix
nodes <- c(1) # Start with one node

# Define preferential attachment function
preferential_attachment <- function(nodes, edges, m, alpha) {
  # Calculate degrees of all nodes
  degrees <- table(c(edges[,1], edges[,2]))
  
  # Add m edges at each step
  for (i in 1:m) {
    # Calculate probability of each node being selected
    probs <- (degrees ^ alpha) / sum(degrees ^ alpha)
    
    # Select target node based on probabilities
    target_node <- sample(nodes, size = 1, prob = probs)
    
    # Add edge to target node
    edges <- rbind(edges, c(target_node, max(nodes) + 1))
    
    # Add new node to network
    nodes <- c(nodes, max(nodes) + 1)
  }
  
  # Return updated nodes and edges
  return(list(nodes = nodes, edges = edges))
}

# Add edges to network iteratively
for (i in 2:N) {
  network <- preferential_attachment(nodes, edges, m, alpha)
  nodes <- network$nodes
  edges <- network$edges
}

# Plot degree distribution
degrees <- table(c(edges[,1], edges[,2]))
plot(degree_distribution, log = "xy", type = "h", main = "Degree Distribution")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Initialize the network with one node
network <- data.frame(from = 1, to = 1)
num_nodes <- 1

# Define the parameters
N <- 1000  # Number of nodes to add
alpha <- 1  # Power law or preferential attachment parameter
m <- 2  # Number of edges to add at each instance

# Function to sample nodes with probability proportional to their degree
sample_nodes <- function(degree) {
  probs <- degree / sum(degree)
  sample(num_nodes, m, replace = TRUE, prob = probs)
}

# Add nodes and edges
for (i in 2:(N+1)) {
  # Add the new node to the network
  num_nodes <- i
  network <- rbind(network, data.frame(from = i, to = i))
  
  # Calculate the degree of each node
  degree <- table(c(network$from, network$to))
  
  # Sample nodes to connect to
  to_connect <- sample_nodes(degree)
  
  # Add edges to the selected nodes
  for (j in to_connect) {
    network <- rbind(network, data.frame(from = i, to = j))
  }
}

# Calculate the degree distribution
degree_dist <- table(c(network$from, network$to))
degree_dist <- degree_dist / sum(degree_dist)

# Plot the degree distribution
plot(degree_dist, log = "xy", xlab = "Degree", ylab = "Frequency")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Define the number of nodes and the number of edges to add at each instance
N <- 100
m <- 3
alpha <- 2.5

# Initialize the network with one node
edges <- matrix(0, ncol = 2, nrow = 1)
degree <- c(1)

# Loop over the remaining nodes
for (i in 2:N) {
  # Add m edges to the network
  for (j in 1:m) {
    # Compute the probability of connecting to each existing node
    prob <- degree^alpha / sum(degree^alpha)
    # Sample an existing node to connect to
    target <- sample(1:(i-1), size = m, replace = T, prob = prob)
    # Add an edge between the new node and the selected target node
    edges <- rbind(edges, c(i, target))
    # Increase the degree of the new node and the target node
    degree[c(i, target)] <- degree[c(i, target)] + 1
  }
}

# Plot the degree distribution of the resulting network
library(ggplot2)
degree_dist <- table(degree) / length(degree)
ggplot(data = data.frame(degree = as.numeric(names(degree_dist)), 
                         freq = as.numeric(degree_dist)),
       aes(x = degree, y = freq)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  xlab("Degree") +
  ylab("Frequency")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to generate scale-free network
scale_free_network <- function(N, alpha, m) {
  # Initialize network with m nodes
  nodes <- 1:m
  edges <- matrix(0, ncol = N, nrow = m)
  degree <- rep(m, N)
  
  # Loop through the remaining nodes
  for (i in (m+1):N) {
    # Choose m edges to attach to existing nodes, based on preferential attachment
    new_edges <- sample(1:(i-1), m, replace = TRUE, prob = degree[1:(i-1)]^alpha/sum(degree[1:(i-1)]^alpha))
  
    # Add new node and edges
    nodes <- c(nodes, i)
    edges <- rbind(edges, matrix(0, ncol = N, nrow = 1))
    edges[i, new_edges] <- 1
    edges[new_edges, i] <- 1
    
    # Update degree of existing nodes
    degree[new_edges] <- degree[new_edges] + 1
    degree[i] <- m
  }
  
  # Return adjacency matrix and node degrees
  degree_seq <- degree#colSums(edges)
  return(list(adjacency = edges, degree_sequence = degree_seq))
}

# Generate a scale-free network with N=1000, alpha=2.0, and m=5
network <- scale_free_network(100, 6, 2)
net=graph.adjacency(as.matrix(network$adjacency), mode="undirected")%>%
igraph::simplify(remove.loops = TRUE,remove.multiple = TRUE)

plot(net,vertex.label=NA,vertex.size=2)
network$degree_sequence


# Plot degree distribution
library(ggplot2)
degree_df <- data.frame(degree = network$degree_sequence)
ggplot(degree_df, aes(x = degree)) + 
  geom_histogram(binwidth = 1, fill = "blue", color = "black") + 
  scale_x_log10() + scale_y_log10() + 
  xlab("Degree") + ylab("Count")
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to generate a scale-free network using the Barabasi-Albert model
# n: the number of nodes
# m: the number of edges to attach from a new node to existing nodes
# Set up the initial graph
# Initialize the graph with two nodes and an edge between them
adj_matrix <- matrix(0, nrow = n, ncol = n)
adj_matrix[1, 2] <- 1
adj_matrix[2, 1] <- 1

# Initialize the degree sequence and degree exponent
degree_seq <- c(1, 1)
alpha <- 6

# Add nodes to the graph
for (i in 3:n) {
  # Calculate the probability of connecting to each existing node based on degree
  prob <- (degree_seq^alpha) / sum(degree_seq^alpha)
  
  # Choose m nodes to connect to
  neighbors <- sample(1:(i-1), size = m, replace = TRUE, prob = prob)
  
  # Add edges to the graph
  adj_matrix[i, neighbors] <- 1
  adj_matrix[neighbors, i] <- 1
  
  # Update the degree sequence
  degree_seq[i] <- m
  degree_seq[neighbors] <- degree_seq[neighbors] + 1
}

# Convert the adjacency matrix to an igraph object for visualization and analysis
library(igraph)
g <- graph_from_adjacency_matrix(adj_matrix)

net11=graph.adjacency(as.matrix(graph), mode="undirected")%>%
  igraph::simplify(remove.loops = TRUE,remove.multiple = TRUE)
net11

# Function to generate a Barabasi-Albert scale-free network with preferential attachment parameter alpha
# N: number of nodes
# m: number of nodes to add at each instance
# alpha: preferential attachment parameter

barabasi_albert <- function(N, m, alpha) {
  # Initialize the network with m nodes and complete graph
  edges <- matrix(1:m, ncol=2)
  degrees <- rep(m, m)
  
  # Loop over remaining nodes
  for (i in (m+1):N) {
    # Compute probability distribution for preferential attachment
    prob <- degrees / sum(degrees) 
    prob_higher <- prob^alpha
    prob_higher <- prob_higher / sum(prob_higher) # renormalize
    prob_final <- (1-alpha)*prob + alpha*prob_higher # mix with uniform distribution
    
    # Choose m nodes to connect to
    neighbors <- sample(1:(i-1), m, replace=FALSE, prob=prob_final)
    
    # Add edges to new node
    edges <- rbind(edges, cbind(rep(i,m), neighbors))
    
    # Update degree distribution
    degrees[neighbors] <- degrees[neighbors] + 1
    degrees <- c(degrees, m)
  }
  
  # Create igraph object from edge list
  library(igraph)
  graph <- graph_from_edgelist(edges, directed=FALSE)
  return(graph)
}

barabasi_albert(100, 4, 6)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set parameters
N <- 1000    # Number of nodes
alpha <- 2   # Power-law exponent
m <- 2       # Number of edges to add at each step

# Initialize the network with m0 nodes and m0(m0-1)/2 edges
m0 <- m + 1
edges <- matrix(0, nrow = m0, ncol = m0)
for (i in 1:m0) {
  for (j in (i+1):m0) {
    edges[i, j] <- 1
    edges[j, i] <- 1
  }
}

# Loop over remaining nodes
for (i in (m0+1):N) {
  # Compute node degrees and probability distribution
  degrees <- colSums(edges)
  prob <- degrees^alpha / sum(degrees^alpha)
  
  # Add m edges using preferential attachment
  new_edges <- sample(1:(i-1), m, replace = TRUE, prob = prob)}
  edges[i, new_edges] <- 1
  edges[new_edges, i] <- 1
}

# Compute degree distribution
degrees <- colSums(edges)
hist(degrees, breaks = 50, col = "blue", xlab = "Degree", ylab = "Frequency", main = "Degree Distribution")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Plot degree distribution
degree_dist <- degree[degree >= m]
hist(degree_dist, breaks = 50, col = "blue", xlab = "Degree", ylab = "Frequency", main = "Degree Distribution")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Required packages
library(Matrix)
library(pracma)
library(spatstat)

# Function to generate Spatial Expander Graph
Spatial_Expander_Graph <- function(N, lambda, L, beta, p, probWithin, probBetween, edge_weight=FALSE, node_attribute=FALSE){
  
  # Generate nodes using Poisson point process
  coords <- matrix(runif(2 * N), ncol = 2)
  coords <- coords * L
  
  # Calculate distance matrix
  dist_mat <- as.matrix(dist(coords, diag = TRUE, upper = TRUE))
  dist_mat <- pmin(dist_mat, L - dist_mat)
  
  # Create adjacency matrix using preferential attachment model
  adj_mat <- matrix(0, nrow = N, ncol = N)
  for(i in 1:N){
    if(i == 1){
      adj_mat[1, 2] <- 1
      adj_mat[2, 1] <- 1
    }else{
      deg <- colSums(adj_mat[1:i-1, 1:i-1])
      for(j in 1:i-1){
        prob_ij <- ((dist_mat[i, j]^(-beta))*(deg[j]^p))/sum((dist_mat[i, 1:i-1]^(-beta))*(deg[1:i-1]^p))
        if(runif(1) < prob_ij){
          adj_mat[i, j] <- 1
          adj_mat[j, i] <- 1
        }
      }
    }
  }
  
  # Add rewiring probability for small world effect
  for(i in 1:N){
    for(j in 1:i-1){
      if(runif(1) < p){
        adj_mat[i, j] <- 0
        adj_mat[j, i] <- 0
        avail_nodes <- which(colSums(adj_mat)==0)
        rand_node <- sample(avail_nodes, 1)
        adj_mat[i, rand_node] <- 1
        adj_mat[rand_node, i] <- 1
      }
    }
  }
  
  # Add community structures
  community_size <- floor(N/4)
  communities <- rep(1:4, each = community_size)
  adj_mat_within <- matrix(0, nrow = N, ncol = N)
  for(i in 1:4){
    community_nodes <- which(communities == i)
    for(j in community_nodes){
      for(k in community_nodes){
        if(j < k){
          if(runif(1) < probWithin){
            adj_mat_within[j, k] <- 1
            adj_mat_within[k, j] <- 1
          }
        }
      }
    }
  }
  
  adj_mat_between <- matrix(0, nrow = N, ncol = N)
  for(i in 1:4){
    for(j in 1:4){
      if(i < j){
        community_i_nodes <- which(communities == i)
        community_j_nodes <- which(communities == j)
        for(k in community_i_nodes){
          for(l in community_j_nodes){
            if(runif(1) < probBetween){
              adj_mat_between[k, l] <- 1
              adj_mat_between[l, k] <- 1
            }
          }
        }
      }
    }
  }
  
  adj_mat <- adj_mat + adj_mat_within + adj_mat_between
  
  # Create weighted adjacency matrix and node attributes
  if(edge_weight){
    adj_mat[adj_mat == 1] <- runif(sum(adj_mat == 1))
  }
  if(node_attribute){
    
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to create a Spatial Expander Graph
Spatial_Expander_Graph2 <- function(N, lambda, L, beta, prewire,
                          prob_within, prob_between, add_weights=FALSE, add_attr=FALSE,
                          add_comm){
  
  # Generate nodes using a 2D Poisson point process on an L-dimensional torus
  coords <- matrix(runif(N*L), ncol=L)
  coords <- coords - floor(coords) # wrap the points on the torus
  
  # Calculate the distance matrix between all pairs of points
  dist_mat <- as.matrix(dist(coords))
  
  # Create the adjacency matrix
  adj_mat <- matrix(0, nrow=N, ncol=N)
  
  # Add edges to the graph with a probability that favors short spatial distances and a high degree that follows a scale-free preferential attachment with power law parameter beta
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      if(runif(1) < (lambda^beta)/((1 + dist_mat[i,j])^beta)){
        adj_mat[i,j] <- 1
        adj_mat[j,i] <- 1
      }
    }
  }
  
  # Rewire the graph with probability p for small world effect
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      if(adj_mat[i,j] == 1 && runif(1) < prewire){
        adj_mat[i,j] <- 0
        adj_mat[j,i] <- 0
        
        # choose a random node and connect to it
        k <- sample(1:N, 1)
        while(k == i || adj_mat[i,k] == 1){
          k <- sample(1:N, 1)
        }
        adj_mat[i,k] <- 1
        adj_mat[k,i] <- 1
      }
    }
  }
  
  # Add community structure to the graph by connecting nodes within communities with probability probWithin more than between communities with probability probBetween
  n_communities <- ceiling(sqrt(N))
  comm_size <- ceiling(N / n_communities)
  communities <- rep(1:n_communities, each = comm_size)[1:N]
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      if(communities[i] == communities[j]){
        if(runif(1) < probWithin){
          adj_mat[i,j] <- 1
          adj_mat[j,i] <- 1
        }
      } else {
        if(runif(1) < probBetween){
          adj_mat[i,j] <- 1
          adj_mat[j,i] <- 1
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

g=Spatial_Expander_Graph2(N=10, lambda=10, L=2,
                           beta=1, prewire=0.4, add_comm=F, 
                           prob_within = 0.1, prob_between = 0.2,add_weights=T,
                           add_attr=T)

plot(g$GraphObject,vertex.label=NA,vertex.size=2)

  
                
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SpatialExpanderGraph <- function(N, lambda, L, beta, r, p, probWithin, probBetween, useWeightedEdges = FALSE) {
  library(Matrix)
  
  # Generate points using a 2D Poisson point process on a torus
  xy <- matrix(runif(N * L, 0, 1), N, L)
  lambda_eff <- lambda / (2 * pi)^L # Effective intensity for torus
  interpoint_distances <- as.matrix(dist(xy))
  
  # Create adjacency matrix based on preferential attachment with small-world effect
  adj_matrix <- matrix(0, N, N)
  for (i in 1:N) {
    degrees <- colSums(adj_matrix) # Calculate degrees of existing nodes
    dists <- interpoint_distances[i,] # Distances to other nodes
    
    # Calculate attachment probabilities
    if (useWeightedEdges) {
      probs <- (dists ^ (-beta)) * (degrees + 1) / sum((dists ^ (-beta)) * (degrees + 1))
    } else {
      probs <- (1 / (dists ^ r)) * (degrees + 1) ^ beta / sum((1 / (dists ^ r)) * (degrees + 1) ^ beta)
    }
    
    # Connect to new node with preferential attachment
    neighbors <- sample(1:N, size = 1, prob = probs)
    adj_matrix[i, neighbors] <- 1
    adj_matrix[neighbors, i] <- 1
    
    # Small-world effect with rewiring probability p
    if (runif(1) < p) {
      non_neighbors <- setdiff(1:N, c(i, neighbors))
      new_neighbor <- sample(non_neighbors, size = 1)
      adj_matrix[i, new_neighbor] <- 1
      adj_matrix[new_neighbor, i] <- 1
      adj_matrix[i, neighbors] <- 0
      adj_matrix[neighbors, i] <- 0
    }
  }
  
  # Add community structure
  community_labels <- cutree(graph.adjacency(adj_matrix, mode = "undirected"), k = sqrt(N))
  for (i in 1:(sqrt(N) - 1)) {
    for (j in (i + 1):sqrt(N)) {
      within_community <- community_labels == i & community_labels == j
      between_communities <- community_labels == i & community_labels != j | community_labels != i & community_labels == j
      if (runif(1) < probWithin) {
        adj_matrix[within_community, within_community] <- 1
      }
      if (runif(1) < probBetween) {
        adj_matrix[between_communities, within_community] <- 1
        adj_matrix[within_community, between_communities] <- 1
      }
    }
  }
  
  # Compute Laplacian matrix and eigenvalues
  degree_matrix <- diag(colSums(adj_matrix))
  laplacian_matrix <- degree_matrix - adj_matrix
  eigenvalues <- eigen(laplacian_matrix)$values
  
  # Return adjacency matrix and eigenvalue gap
  return(list(adj_matrix = adj_matrix, eigenvalue_gap = eigenvalues[2]))
}

###+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###+sample code
###++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to generate Spatial Expander Graph
spatial_expander_graph <- function(N, lambda, L, beta, r, p, probWithin, probBetween, use_edge_weights = TRUE) {
  library(Matrix)
  
  # Generate points using 2D Poisson point process on a torus
  set.seed(123) # for reproducibility
  x <- runif(N, 0, L)
  y <- runif(N, 0, L)
  dx <- outer(x, x, function(x1, x2) min(abs(x1 - x2), L - abs(x1 - x2)))
  dy <- outer(y, y, function(y1, y2) min(abs(y1 - y2), L - abs(y1 - y2)))
  d <- sqrt(dx^2 + dy^2)
  
  # Calculate the distance matrix
  D <- as.matrix(d)
  
  # Create adjacency matrix using preferential attachment
  A <- Matrix(0, N, N, sparse = TRUE)
  for (i in 1:N) {
    for (j in 1:(i-1)) {
      if (runif(1) < r*(d[j,i]^(-beta))) {
        A[i,j] <- 1
        A[j,i] <- 1
      }
    }
  }
  
  # Add small world effect to the graph
  for (i in 1:N) {
    for (j in 1:(i-1)) {
      if (A[i,j] == 1 && runif(1) < p) {
        k <- sample(1:N, 1)
        if (k != i && k != j && A[i,k] == 0) {
          A[i,j] <- 0
          A[j,i] <- 0
          A[i,k] <- 1
          A[k,i] <- 1
        }
      }
    }
  }
  
  # Add community structures to the graph
  comms <- rep(1:(N/10), each = 10)
  C <- outer(comms, comms, function(x, y) x == y)
  for (i in 1:N) {
    for (j in 1:(i-1)) {
      if (C[i,j] && runif(1) < probWithin) {
        A[i,j] <- 1
        A[j,i] <- 1
      } else if (!C[i,j] && runif(1) < probBetween) {
        A[i,j] <- 1
        A[j,i] <- 1
      }
    }
  }
  
  # Create weighted adjacency matrix if specified
  if (use_edge_weights) {
    w <- exp(-D^2/(2*lambda^2))
    W <- Matrix(w, sparse = TRUE)
    A <- A * W
  }
  
  # Compute Laplacian matrix
  L <- Diagonal(N) - A
  eig <- eigen(L, symmetric = TRUE, only.values = TRUE)
  eig_vals <- sort(eig$values)
  
  # Plot graph
  plot(graph_from_adjacency_matrix(A))
  
  # Return adjacency matrix, Laplacian matrix, and eigenvalue gap
  return(list(A = A, L = L, eigenvalue_gap = eig_vals[2] - eig_vals[1]))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

spatial_expander_graph <- function(N, lambda, L, beta, r, p, probWithin, probBetween, useEdgeWeights = TRUE) {
  library(Matrix)
  
  # Step 1: Generate N points using a 2D Poisson point process on an L-dimensional torus
  points <- matrix(runif(N*L), ncol=L)
  points <- points - floor(points) # wrap around torus
  
  # Step 2: Calculate the distance matrix
  dist_matrix <- as.matrix(dist(points))
  
  # Step 3: Create an adjacency matrix using the distance matrix
  adj_matrix <- as.matrix(dist_matrix <= r)
  if (useEdgeWeights) {
    adj_matrix[adj_matrix == 1] <- exp(-dist_matrix[adj_matrix == 1])
  }
  
  # Step 4: Preferential attachment with probability favoring short spatial distances and high degree
  for (i in 1:(N-2)) {
    degrees <- colSums(adj_matrix)
    p_new <- (degrees^beta) * exp(-dist_matrix[, i]/r)
    p_new[i] <- 0
    p_new[1:(i-1)] <- 0
    p_new <- p_new/sum(p_new)
    j <- sample(1:N, 1, prob=p_new)
    adj_matrix[i, j] <- 1
    adj_matrix[j, i] <- 1
  }
  
  # Step 5: Small world effect with rewiring probability p
  for (i in 1:(N-2)) {
    for (j in (i+1):N) {
      if (adj_matrix[i,j] == 0) {
        if (runif(1) < p) {
          adj_matrix[i,j] <- 1
          adj_matrix[j,i] <- 1
        }
      }
    }
  }
  
  # Step 6: Community structure
  community_size <- floor(N/2)
  communities <- rep(1:2, each=community_size)
  within_prob <- probWithin
  between_prob <- probBetween
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (communities[i] == communities[j]) {
        if (runif(1) < within_prob) {
          adj_matrix[i,j] <- 1
          adj_matrix[j,i] <- 1
        }
      } else {
        if (runif(1) < between_prob) {
          adj_matrix[i,j] <- 1
          adj_matrix[j,i] <- 1
        }
      }
    }
  }
  
  # Step 7: Compute Laplacian matrix and check expansion properties
  degrees <- colSums(adj_matrix)
  lap_matrix <- diag(degrees) - adj_matrix
  eig_vals <- eigen(lap_matrix)$values
  gap_size <- eig_vals[2] - eig_vals[1]
  
  # Return adjacency matrix and Laplacian matrix
  return(list(adj_matrix = adj_matrix, lap_matrix = lap_matrix, gap_size = gap_size))
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
spatial_expander_graph <- function(N, lambda, L, beta, r, p, probWithin, probBetween, use_edge_weights = FALSE) {
  
  library(Matrix)
  
  # Step 1: Generate N points using a 2D Poisson point process with intensity parameter lambda on an L-dimensional torus
  coords <- matrix(runif(N * L), ncol = L)
  coords <- coords %% 1  # wrap around to torus
  
  # Step 2: Calculate the distance matrix between all pairs of points
  dists <- as.matrix(dist(coords))
  
  # Step 3: Create adjacency matrix based on distance matrix with probability favoring short spatial distances and high degree
  A <- matrix(0, N, N)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (runif(1) < r * exp(-beta * dists[i,j]) || i == j) {
        if (use_edge_weights) {
          A[i,j] <- exp(-dists[i,j])
          A[j,i] <- exp(-dists[i,j])
        } else {
          A[i,j] <- 1
          A[j,i] <- 1
        }
      }
    }
  }
  
  # Step 4: Rewire edges with probability p for small world effect
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (A[i,j] == 1 && runif(1) < p) {
        # find a random node that is not a neighbor of i and rewire edge to it
        possible_neighbors <- which(A[i,] == 0 & seq_along(A[i,]) != i)
        if (length(possible_neighbors) > 0) {
          new_neighbor <- sample(possible_neighbors, size = 1)
          A[i,j] <- 0
          A[j,i] <- 0
          A[i,new_neighbor] <- 1
          A[new_neighbor,i] <- 1
        }
      }
    }
  }
  
  # Step 5: Add community structure to the graph
  num_communities <- ceiling(sqrt(N))
  community_sizes <- rep(N/num_communities, num_communities)
  community_sizes[1:(N %% num_communities)] <- community_sizes[1:(N %% num_communities)] + 1
  community_assignments <- rep(1:num_communities, community_sizes)
  within_probs <- diag(rep(probWithin, num_communities))
  between_probs <- (1-probBetween) / (num_communities-1) * (matrix(1, num_communities, num_communities) - diag(rep(1, num_communities)))
  community_probs <- within_probs + between_probs
  
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (community_assignments[i] == community_assignments[j] && runif(1) < probWithin) {
        A[i,j] <- 1
        A[j,i] <- 1
      } else if (community_assignments[i] != community_assignments[j] && runif(1) < probBetween) {
        A[i,j] <- 1
        A[j,i] <- 1
      }
    }
  }
  
  # Step 6: Compute Laplacian matrix and eigenvalues
  D <- Matrix(0, N, N
              
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
spatial_expander_graph <- function(N, lambda, L, beta, r, p, probWithin, probBetween, use_edge_weights = FALSE) {
  # Step 1: Generate N points using a 2D Poisson point process on an L-dimensional torus
  points <- matrix(runif(N*L), ncol = L)
  for (i in 1:L) {
    points[,i] <- points[,i] * L
  }
  
  # Step 2: Calculate the distance matrix between all pairs of points
  dist_mat <- as.matrix(dist(points, method = "euclidean"))
  
  # Step 3: Create an adjacency matrix for the graph using the distance matrix
  adj_mat <- matrix(0, nrow = N, ncol = N)
  if (use_edge_weights) {
    for (i in 1:(N-1)) {
      for (j in (i+1):N) {
        w <- exp(-r*dist_mat[i,j]) + rgamma(1, shape = beta, rate = 1/beta)
        if (runif(1) < p) {
          adj_mat[i,j] <- w
          adj_mat[j,i] <- w
        }
      }
    }
  } else {
    for (i in 1:(N-1)) {
      for (j in (i+1):N) {
        if (runif(1) < exp(-r*dist_mat[i,j]) + rgamma(1, shape = beta, rate = 1/beta)) {
          if (runif(1) < p) {
            adj_mat[i,j] <- 1
            adj_mat[j,i] <- 1
          }
        }
      }
    }
  }
  
  # Step 4: Add community structures to the graph
  comm_sizes <- c(rep(ceiling(N/2), 2), rep(floor(N/2), 2))
  comm_sizes <- comm_sizes[1:4]
  comm_assignments <- rep(1:4, times = comm_sizes)
  within_probs <- matrix(probWithin, nrow = 4, ncol = 4)
  within_probs[1,2] <- within_probs[2,1] <- probBetween
  within_probs[3,4] <- within_probs[4,3] <- probBetween
  within_probs[1,3] <- within_probs[3,1] <- 0
  within_probs[1,4] <- within_probs[4,1] <- 0
  within_probs[2,4] <- within_probs[4,2] <- 0
  within_probs[2,3] <- within_probs[3,2] <- 0
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (comm_assignments[i] == comm_assignments[j]) {
        if (runif(1) < within_probs[comm_assignments[i], comm_assignments[j]]) {
          adj_mat[i,j] <- 1
          adj_mat[j,i] <- 1
        }
      } else {
        if (runif(1) < probBetween) {
          adj_mat[i,j] <- 1
          adj_mat[j,i] <- 1
        }
      }
    }
  }
  

  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

spatial_expander_graph <- function(N, lambda, L, beta, p, probWithin, probBetween, use_weights = TRUE) {
  
  # Generate points using a 2D Poisson point process on an L-dimensional torus
  points <- matrix(runif(N*L), ncol = L)
  points <- points - floor(points)
  
  # Calculate the distance matrix between all pairs of points
  dist_matrix <- as.matrix(dist(points, method = "euclidean"))
  
  # Create adjacency matrix with probability favoring short spatial distances and high degree
  adj_matrix <- matrix(0, nrow = N, ncol = N)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (runif(1) < exp(-beta * dist_matrix[i,j])) {
        if (use_weights) {
          adj_matrix[i,j] <- adj_matrix[j,i] <- dist_matrix[i,j]
        } else {
          adj_matrix[i,j] <- adj_matrix[j,i] <- 1
        }
      }
    }
  }
  
  # Add rewiring probability for small world effect
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (runif(1) < p) {
        if (use_weights) {
          adj_matrix[i,j] <- adj_matrix[j,i] <- dist_matrix[i,j]
        } else {
          adj_matrix[i,j] <- adj_matrix[j,i] <- 1
        }
      }
    }
  }
  
  # Add community structures to the graph
  community_sizes <- c(rep(ceiling(N/2), 2), rep(floor(N/2), 2))
  num_communities <- length(community_sizes)
  communities <- rep(1:num_communities, times = community_sizes)
  
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (communities[i] == communities[j]) {
        if (runif(1) < probWithin) {
          adj_matrix[i,j] <- adj_matrix[j,i] <- 1
        }
      } else {
        if (runif(1) < probBetween) {
          adj_matrix[i,j] <- adj_matrix[j,i] <- 1
        }
      }
    }
  }
  
  # Compute Laplacian matrix and eigenvalues
  degree_matrix <- diag(rowSums(adj_matrix))
  laplacian_matrix <- degree_matrix - adj_matrix
  eigenvalues <- eigen(laplacian_matrix)$values
  
  # Return the adjacency matrix and eigenvalues
  return(list(adj_matrix = adj_matrix, eigenvalues = eigenvalues))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

spatialExpanderGraph <- function(N, lambda, L, beta, p, probWithin, probBetween, useEdgeWeights = FALSE) {
  
  # Generate points using 2D Poisson point process on a torus
  points <- matrix(runif(N * L, 0, 1), ncol = L) # generate uniform points
  points <- points - floor(points) # wrap around torus
  lambdaEff <- lambda / (2^L) # effective intensity for torus
  
  # Calculate distance matrix between all pairs of points
  distMat <- as.matrix(dist(points, method = "euclidean"))
  
  # Create adjacency matrix using preferential attachment with small world effect
  A <- matrix(0, nrow = N, ncol = N)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (runif(1) < exp(-beta * distMat[i,j])) {
        A[i,j] <- 1
        A[j,i] <- 1
      }
    }
  }
  if (p > 0) {
    for (i in 1:(N-1)) {
      for (j in (i+1):N) {
        if (A[i,j] == 1 && runif(1) < p) {
          # rewiring edge
          k <- sample(setdiff(1:N, c(i,j)), 1)
          A[i,j] <- 0
          A[j,i] <- 0
          A[i,k] <- 1
          A[k,i] <- 1
        }
      }
    }
  }
  
  # Add community structure to the graph
  communitySizes <- table(sample(1:(N/probWithin), N, replace = TRUE)) * probWithin
  communityLabels <- rep(1:length(communitySizes), communitySizes)
  withinCommunityPairs <- t(combn(1:N, 2)[communityLabels[combn(1:N, 2)[,1]] == communityLabels[combn(1:N, 2)[,2]],])
  betweenCommunityPairs <- t(combn(1:N, 2)[communityLabels[combn(1:N, 2)[,1]] != communityLabels[combn(1:N, 2)[,2]],])
  for (i in 1:nrow(withinCommunityPairs)) {
    if (runif(1) < probWithin) {
      A[withinCommunityPairs[i,1], withinCommunityPairs[i,2]] <- 1
      A[withinCommunityPairs[i,2], withinCommunityPairs[i,1]] <- 1
    }
  }
  for (i in 1:nrow(betweenCommunityPairs)) {
    if (runif(1) < probBetween) {
      A[betweenCommunityPairs[i,1], betweenCommunityPairs[i,2]] <- 1
      A[betweenCommunityPairs[i,2], betweenCommunityPairs[i,1]] <- 1
    }
  }
  
  # Create weighted adjacency matrix if specified
  if (useEdgeWeights) {
    W <- A * exp(-beta * distMat)
    diag(W) <- 0
    return(W)
  }
  
  # Create Laplacian matrix and return
  D <- diag(rowSums(A))
  L <- D - A
  return(L)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to generate Spatial Expander Graph
SpatialExpanderGraph <- function(N, lambda, L, beta, r, p, probWithin, probBetween, useEdgeWeight = FALSE) {
  
  # Generate N points using 2D Poisson point process on an L dimensional torus
  points <- matrix(runif(N*L), ncol = L) * 2 * pi
  points <- t(t(points) %% (2 * pi))
  
  # Calculate distance matrix between all pairs of points
  distMat <- as.matrix(dist(points, diag = TRUE, upper = TRUE))
  
  # Create adjacency matrix for the graph with arguments to use edge weights or not
  adjMat <- ifelse(distMat < r, 1, 0)
  if (useEdgeWeight) {
    adjMat <- ifelse(adjMat == 1, r - distMat, 0)
  }
  
  # Apply preferential attachment with power law parameter 'beta'
  degrees <- apply(adjMat, 1, sum)
  for (i in 1:(N - 1)) {
    prob <- (degrees[i] + beta) / sum(degrees + beta)
    for (j in (i+1):N) {
      if (runif(1) < prob && adjMat[i,j] == 0) {
        adjMat[i,j] <- 1
        adjMat[j,i] <- 1
        degrees[i] <- degrees[i] + 1
        degrees[j] <- degrees[j] + 1
      }
    }
  }
  
  # Apply rewiring probability 'p' for small world effect
  for (i in 1:(N - 1)) {
    for (j in (i+1):N) {
      if (adjMat[i,j] == 0 && runif(1) < p) {
        adjMat[i,j] <- 1
        adjMat[j,i] <- 1
      }
    }
  }
  
  # Apply community structure to the graph
  communityLabels <- cutree(cluster::walktrap.community(graph.adjacency(adjMat)))
  communityMat <- ifelse(outer(communityLabels, communityLabels, "=="), 1, 0)
  for (i in 1:(N - 1)) {
    for (j in (i+1):N) {
      if (communityMat[i,j] == 1 && runif(1) < probWithin) {
        adjMat[i,j] <- 1
        adjMat[j,i] <- 1
        degrees[i] <- degrees[i] + 1
        degrees[j] <- degrees[j] + 1
      } else if (communityMat[i,j] == 0 && runif(1) < probBetween) {
        adjMat[i,j] <- 1
        adjMat[j,i] <- 1
        degrees[i] <- degrees[i] + 1
        degrees[j] <- degrees[j] + 1
      }
    }
  }
  
  # Compute Laplacian matrix
  degreeMat <- diag(degrees)
  laplacianMat <- degreeMat - adjMat
  
  # Return adjacency matrix and Laplacian matrix
  return(list(adjMat = adjMat, laplacianMat = laplacianMat))
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SpatialExpanderGraph <- function(beta, r, N, lambda, L, p, probWithin, probBetween) {
  
  # Generate points using a 2D Poisson point process with intensity parameter lambda
  X <- matrix(runif(2*N), ncol = 2)
  lambdaArea <- lambda*L^2 # Adjust intensity for torus
  R <- sqrt(runif(N, 0, 1/lambdaArea))
  Theta <- runif(N, 0, 2*pi)
  X[,1] <- L*sqrt(R)*cos(Theta)
  X[,2] <- L*sqrt(R)*sin(Theta)
  
  # Calculate distance matrix between all pairs of points
  distMat <- as.matrix(dist(X, method = "euclidean"))
  
  # Create adjacency matrix
  adjMat <- matrix(0, nrow = N, ncol = N)
  
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (runif(1) < exp(-distMat[i,j]^beta/r)) { # preferential attachment with distance favoring
        adjMat[i,j] <- 1
        adjMat[j,i] <- 1
      }
    }
  }
  
  # Rewire edges with probability p
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (adjMat[i,j] == 1 && runif(1) < p) {
        k <- sample((1:N)[-c(i,j)], 1) # choose a new node uniformly at random
        adjMat[i,j] <- 0
        adjMat[j,i] <- 0
        adjMat[i,k] <- 1
        adjMat[k,i] <- 1
      }
    }
  }
  
  # Add community structures
  # communityLabels <- cutree(walktrap.community(graph.adjacency(adjMat)))
  # 
  # for (i in 1:(N-1)) {
  #   for (j in (i+1):N) {
  #     if (communityLabels[i] == communityLabels[j]) { # within community
  #       if (runif(1) < probWithin) {
  #         adjMat[i,j] <- 1
  #         adjMat[j,i] <- 1
  #       }
  #     } else { # between communities
  #       if (runif(1) < probBetween) {
  #         adjMat[i,j] <- 1
  #         adjMat[j,i] <- 1
  #       }
  #     }
  #   }
  # }
  # 
  # # Compute Laplacian matrix
  # degreeVec <- rowSums(adjMat)
  # lapMat <- diag(degreeVec) - adjMat
  # 
  # Return adjacency matrix, Laplacian matrix, and community labels
  #return(list(adjMat = adjMat, lapMat = lapMat, communityLabels = communityLabels))
  return(adjMat)
}

# Generate a Spatial Expander Graph with N = 100, beta = 2, r = 0.1, lambda = 0.1, L = 10, p = 0.01,
# probWithin = 0.5, probBetween = 0.01
graph <- SpatialExpanderGraph(beta = 1, r = 0.1, N = 100, lambda = 0.5, L = 2, p = 0.01,
                              probWithin = 0.2,probBetween=0.1)
                              
g <- graph.adjacency(as.matrix(graph), mode="undirected")
g=g%>%simplify(remove.loops = TRUE,remove.multiple = TRUE)
plot(g,vertex.label=NA,vertex.size=2)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+samplecode
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Spatial_Expander_Graph <- function(beta, r, N, lambda, L, p, probWithin, probBetween, use_edge_weights = TRUE) {
  
  # Generate points using 2D Poisson point process on a torus
  points <- matrix(runif(2*N), ncol = 2)
  points <- points * L
  distance_matrix <- matrix(0, nrow = N, ncol = N)
  
  # Calculate distance matrix
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        dx <- abs(points[i,1] - points[j,1])
        dx <- min(dx, L - dx)
        dy <- abs(points[i,2] - points[j,2])
        dy <- min(dy, L - dy)
        distance_matrix[i,j] <- sqrt(dx^2 + dy^2)
      }
    }
  }
  
  # Create adjacency matrix with preferential attachment
  adj_matrix <- matrix(0, nrow = N, ncol = N)
  for (i in 1:N) {
    for (j in 1:i) {
      if (distance_matrix[i,j] < r) {
        weight <- 1/distance_matrix[i,j]^beta
        adj_matrix[i,j] <- weight
        adj_matrix[j,i] <- weight
      }
    }
  }
  
  # Add small world effect
  for (i in 1:N) {
    for (j in (i+1):N) {
      if (runif(1) < p) {
        adj_matrix[i,j] <- 0
        adj_matrix[j,i] <- 0
        neighbors <- which(adj_matrix[i,] > 0)
        new_neighbor <- sample(setdiff(1:N, c(i,neighbors)), 1)
        adj_matrix[i,new_neighbor] <- 1
        adj_matrix[new_neighbor,i] <- 1
      }
    }
  }
  
  # Add community structure
  community_size <- ceiling(N/2)
  community_labels <- rep(1:2, each = community_size)
  prob_matrix <- matrix(probBetween, nrow = 2, ncol = 2)
  prob_matrix[1,1] <- probWithin
  prob_matrix[2,2] <- probWithin
  for (i in 1:N) {
    for (j in (i+1):N) {
      if (community_labels[i] == community_labels[j]) {
        if (runif(1) < prob_matrix[community_labels[i], community_labels[j]]) {
          adj_matrix[i,j] <- 1
          adj_matrix[j,i] <- 1
        }
      } else {
        if (runif(1) < prob_matrix[community_labels[i], community_labels[j]]) {
          adj_matrix[i,j] <- 1
          adj_matrix[j,i] <- 1
        }
      }
    }
  }
  
  # Create Laplacian matrix
  if (use_edge_weights) {
    degree_vector <- rowSums(adj_matrix)
    diagonal_matrix <- diag(degree_vector)
    laplacian_matrix <- diagonal_matrix - adj_matrix
  } else {
    adj_matrix[adj_matrix > 0] <- 1
    degree_vector <- rowSums(adj_matrix)
    diagonal_matrix <- diag(degree_vector)
    laplacian_matrix <- diagonal_matrix - adj_matrix
  }
  
  # Return graph as a list
  graph_list <- list()
  graph_list$adj
  

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  spatial_expander_graph <- function(N, lambda, L, beta, r, p, probWithin, probBetween, use_weights = FALSE) {
    
    # Step 1: Generate points using a 2D Poisson point process on an L-dimensional torus
    points <- matrix(runif(N*L, 0, 1), ncol = L)
    lambda_vol <- lambda^L
    num_points <- rpois(1, lambda_vol)
    if (num_points < N) {
      points <- matrix(runif(num_points*L, 0, 1), ncol = L)
    }
    dist_matrix <- as.matrix(dist(points, method = "euclidean", diag = TRUE, upper = TRUE))
    
    # Step 2: Calculate the distance matrix between all pairs of points
    # Step 3: Create an adjacency matrix for the graph using preferential attachment
    adjacency_matrix <- matrix(0, nrow = N, ncol = N)
    degree_vector <- rep(0, N)
    for (i in 1:N) {
      for (j in (i+1):N) {
        if (runif(1) < exp(-r * dist_matrix[i,j])) {
          adjacency_matrix[i,j] <- 1
          adjacency_matrix[j,i] <- 1
          if (use_weights) {
            adjacency_matrix[i,j] <- exp(-r * dist_matrix[i,j])
            adjacency_matrix[j,i] <- exp(-r * dist_matrix[i,j])
          }
          degree_vector[i] <- degree_vector[i] + 1
          degree_vector[j] <- degree_vector[j] + 1
        }
      }
    }
    
    # Calculate the degree distribution and its cumulative distribution
    degree_dist <- table(degree_vector)/N
    cum_degree_dist <- cumsum(degree_dist)
    
    # Use the degree distribution to add edges with preferential attachment
    for (i in 1:N) {
      for (j in (i+1):N) {
        if (runif(1) < p) {
          if (runif(1) < exp(-r * dist_matrix[i,j])) {
            adjacency_matrix[i,j] <- 1
            adjacency_matrix[j,i] <- 1
            if (use_weights) {
              adjacency_matrix[i,j] <- exp(-r * dist_matrix[i,j])
              adjacency_matrix[j,i] <- exp(-r * dist_matrix[i,j])
            }
          }
          else {
            # Select a neighbor based on degree distribution
            prob <- cum_degree_dist - degree_dist[j]
            prob[prob < 0] <- 0
            prob <- prob/sum(prob)
            neighbor <- sample(1:N, 1, prob = prob)
            adjacency_matrix[i,neighbor] <- 1
            adjacency_matrix[neighbor,i] <- 1
            if (use_weights) {
              adjacency_matrix[i,neighbor] <- exp(-r * dist_matrix[i,neighbor])
              adjacency_matrix[neighbor,i] <- exp(-r * dist_matrix[i,neighbor])
            }
            degree_vector[i] <- degree_vector[i] + 1
            degree_vector[neighbor] <- degree_vector[neighbor] + 1
          }
        }
      }
    }
    
    # Step 5: Add community structure to the graph
    num_communities <- round(sqrt(N))
    community_sizes <- rep(ceiling(N/num_communities), num_communities)
    community_sizes[1:(N %% num_communities)] <- community_sizes[1:(N %% num_comm
                                                                    
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Load required packages
library(Matrix)
library(Matrix.utils)

# Function to generate a 2D Poisson point process with intensity parameter lambda on an L dimensional torus
generate_points <- function(N, lambda, L) {
  points <- matrix(runif(N*L), ncol=L)
  for (i in 1:L) {
    points[,i] <- points[,i]*lambda^(1/L)
  }
  return(points %% 1)
}

# Function to calculate the distance matrix between all pairs of points
calculate_distance_matrix <- function(points, L) {
  N <- nrow(points)
  distance_matrix <- matrix(0, nrow=N, ncol=N)
  for (i in 1:N) {
    for (j in 1:N) {
      d <- min(abs(points[i,] - points[j,]), 1 - abs(points[i,] - points[j,]))
      distance_matrix[i,j] <- sqrt(sum(d^2))
    }
  }
  return(distance_matrix)
}

# Function to create an adjacency matrix for the graph
create_adjacency_matrix <- function(distance_matrix, r, p, beta) {
  N <- nrow(distance_matrix)
  adjacency_matrix <- matrix(0, nrow=N, ncol=N)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (distance_matrix[i,j] <= r && runif(1) <= exp(-distance_matrix[i,j]^beta)) {
        if (runif(1) <= p) {
          # Small-world rewiring
          k <- sample((1:N)[-c(i,j)][distance_matrix[i,-c(i,j)] <= r], 1)
          if (k < i) {
            adjacency_matrix[i,k] <- 1
          } else if (k > i && k < j) {
            adjacency_matrix[k,i] <- 1
          } else if (k > j) {
            adjacency_matrix[j,k] <- 1
          }
        } else {
          adjacency_matrix[i,j] <- 1
          adjacency_matrix[j,i] <- 1
        }
      }
    }
  }
  return(adjacency_matrix)
}

# Function to create community structures in the graph
create_communities <- function(adjacency_matrix, probWithin, probBetween) {
  N <- nrow(adjacency_matrix)
  communities <- sample(x=1:10, size=N, replace=TRUE)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (communities[i] == communities[j]) {
        if (runif(1) <= probWithin) {
          adjacency_matrix[i,j] <- 1
          adjacency_matrix[j,i] <- 1
        }
      } else {
        if (runif(1) <= probBetween) {
          adjacency_matrix[i,j] <- 1
          adjacency_matrix[j,i] <- 1
        }
      }
    }
  }
  return(adjacency_matrix)
}

# Function to generate the spatial expander graph
generate_spatial_expander_graph <- function(N, lambda, L, r, p, probWithin, probBetween, beta) {
  # Generate points
  points <- generate_points(N, lambda, L)
  
  # Calculate distance matrix
  distance_matrix <- calculate_distance_matrix(points, L)
  
  # Create adjacency matrix
  adjacency_matrix <- create_adjacency_matrix(distance_matrix, r, p, beta)
  
  # Add community structures to the graph
  adjacency_matrix
  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  library(Matrix)
  
  generate_spatial_graph <- function(N, lambda, L, beta, r, p = 0, probWithin = 0.9, probBetween = 0.1) {
    
    # Generate N points using a 2D Poisson point process on an L-dimensional torus
    coords <- matrix(runif(N*L), ncol=L)
    coords <- coords * (2*pi*r)
    coords <- cbind(cos(coords), sin(coords))
    
    # Calculate distance matrix
    dist_mat <- as.matrix(dist(coords, diag = TRUE, upper = TRUE))
    
    # Create adjacency matrix
    adj_mat <- matrix(0, nrow = N, ncol = N)
    for (i in 1:N) {
      for (j in 1:N) {
        if (i == j) {
          next
        }
        else if (runif(1) < lambda * exp(-dist_mat[i,j])) {
          adj_mat[i,j] <- 1
          adj_mat[j,i] <- 1
        }
      }
    }
    
    # Add small world effect with probability p
    if (p > 0) {
      for (i in 1:N) {
        for (j in (i+1):N) {
          if (adj_mat[i,j] == 1 && runif(1) < p) {
            # Rewire edge i-j to a random node k
            k <- sample(setdiff(1:N, c(i,j)), 1)
            adj_mat[i,j] <- 0
            adj_mat[j,i] <- 0
            adj_mat[i,k] <- 1
            adj_mat[k,i] <- 1
          }
        }
      }
    }
    
    # Add community structure
    comm_sizes <- sample(1:(N/4), 4, replace = TRUE)
    comm_sizes <- comm_sizes / sum(comm_sizes) * N
    comm_sizes <- round(comm_sizes)
    communities <- rep(1:4, comm_sizes)
    within_comm_prob <- probWithin / (mean(comm_sizes) - 1)
    between_comm_prob <- probBetween / (N - mean(comm_sizes))
    for (i in 1:N) {
      for (j in (i+1):N) {
        if (communities[i] == communities[j] && runif(1) < within_comm_prob * exp(-dist_mat[i,j])) {
          adj_mat[i,j] <- 1
          adj_mat[j,i] <- 1
        }
        else if (communities[i] != communities[j] && runif(1) < between_comm_prob * exp(-dist_mat[i,j])) {
          adj_mat[i,j] <- 1
          adj_mat[j,i] <- 1
        }
      }
    }
    
    # Add scale-free structure with power law parameter beta
    degree_seq <- rowSums(adj_mat)
    for (i in (N+1):(2*N)) {
      prob_attach <- degree_seq / sum(degree_seq)
      attach_node <- sample(1:N, 1, prob = prob_attach)
      adj_mat[attach_node, i] <- 1
      adj_mat[i, attach_node] <- 1
      degree_seq[attach_node] <- degree_seq[attach_node] + 1
      degree_seq[i] <- 1
      prob_attach[i] <- 0
      prob_attach[attach_node] <- 0
      prob_attach <- prob_attach^beta
      prob_attach <- prob_attach / sum(prob
                                       
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  
# Function to generate a spatial expander graph
# Arguments:
#   N: number of nodes
#   lambda: intensity parameter of Poisson point process
#   L: number of dimensions of torus
#   r: radius parameter for connecting nodes
#   p: rewiring probability
#   probWithin: within-node connection probability for community structure
#   probBetween: between-node connection probability for community structure
#   beta: preferential attachment power law parameter
#   use_edge_weights: logical indicating whether to use edge weights in adjacency matrix
# Returns:
#   adjacency matrix of graph

generate_spatial_expander_graph <- function(N, lambda, L, r, p, probWithin, probBetween, beta, use_edge_weights) {
  # Generate N points using Poisson point process on L-dimensional torus
  points <- matrix(runif(N*L), ncol=L)
  points <- points * (2*pi) # scale to torus
  
  # Calculate distance matrix
  dist_mat <- as.matrix(dist(points, method="euclidean", diag=TRUE, upper=TRUE))
  dist_mat <- pmin(dist_mat, dist_mat[t(1:N), 1:N] + matrix(rep(2*pi, N*N), ncol=N))
  
  # Create adjacency matrix using spatial distance and rewiring probability
  adj_mat <- matrix(0, nrow=N, ncol=N)
  for (i in 1:N) {
    for (j in (i+1):N) {
      if (dist_mat[i,j] < r) {
        if (runif(1) < (1-p)) {
          adj_mat[i,j] <- 1
          adj_mat[j,i] <- 1
          if (use_edge_weights) {
            adj_mat[i,j] <- 1/dist_mat[i,j]
            adj_mat[j,i] <- 1/dist_mat[i,j]
          }
        } else {
          # Rewire edge
          rewired <- FALSE
          while (!rewired) {
            k <- sample(N, 1)
            if (k != i && k != j && adj_mat[i,k] == 0 && adj_mat[j,k] == 0 && runif(1) < p) {
              adj_mat[i,j] <- 0
              adj_mat[j,i] <- 0
              adj_mat[i,k] <- 1
              adj_mat[k,i] <- 1
              adj_mat[j,k] <- 1
              adj_mat[k,j] <- 1
              if (use_edge_weights) {
                adj_mat[i,k] <- 1/dist_mat[i,k]
                adj_mat[k,i] <- 1/dist_mat[i,k]
                adj_mat[j,k] <- 1/dist_mat[j,k]
                adj_mat[k,j] <- 1/dist_mat[j,k]
              }
              rewired <- TRUE
            }
          }
        }
      }
    }
  }
  
  # Add community structure
  num_communities <- 4 # can be tuned
  community_size <- ceiling(N / num_communities)
  community_assignments <- rep(1:num_communities, each=community_size)[1:N]
  for (i in 1:N) {
    for (j in (i+1):N) {
      if (community_assignments[i] == community_assignments[j]) {
        if (runif(1) < probWithin) {
          adj_mat[i,j] <- 1
          adj_mat[j,i] <- 1
          if (use_edge_weights) {
            
            
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
spatial_expander_graph <- function(N, lambda, L, r, p, probWithin, probBetween, beta){
  # Step 1: Generate a set of N points using a 2D Poisson point process with intensity parameter 'lambda' on an 'L' dimensional torus
  points <- matrix(runif(N*L), ncol=L)
  points <- points - floor(points)
  points <- matrix(qnorm(points, lambda), ncol=L)
  
  # Step 2: Calculate the distance matrix between all pairs of points
  dist_matrix <- as.matrix(dist(points, diag=TRUE, upper=TRUE))
  
  # Step 3: Use the distance matrix to create an adjacency matrix for the graph
  adj_matrix <- matrix(0, nrow=N, ncol=N)
  for (i in 1:N) {
    for (j in (i+1):N) {
      if (dist_matrix[i, j] <= r && runif(1) > p) {
        adj_matrix[i, j] <- 1
        adj_matrix[j, i] <- 1
      }
    }
  }
  
  # Step 4: Add a rewiring probability 'p' for small world effect
  for (i in 1:N) {
    for (j in (i+1):N) {
      if (adj_matrix[i, j] == 0 && runif(1) < p) {
        adj_matrix[i, j] <- 1
        adj_matrix[j, i] <- 1
      }
    }
  }
  
  # Step 5: Add community structures to the graph
  community_size <- floor(N / 2)
  community_assignment <- rep(1:2, each=community_size)[1:N]
  for (i in 1:N) {
    for (j in (i+1):N) {
      if (community_assignment[i] == community_assignment[j]) {
        if (runif(1) < probWithin) {
          adj_matrix[i, j] <- 1
          adj_matrix[j, i] <- 1
        }
      } else {
        if (runif(1) < probBetween) {
          adj_matrix[i, j] <- 1
          adj_matrix[j, i] <- 1
        }
      }
    }
  }
  
  # Step 6: Add a preferential attachment power law parameter 'beta'
  degree_seq <- colSums(adj_matrix)
  for (i in (N+1):(2*N)) {
    prob_vec <- degree_seq / sum(degree_seq)
    new_node_degree <- rpower(1, beta) * (N-1)
    neighbors <- sample(1:N, size=new_node_degree, replace=TRUE, prob=prob_vec)
    adj_matrix[i, neighbors] <- 1
    adj_matrix[neighbors, i] <- 1
    degree_seq <- colSums(adj_matrix)
  }
  
  # Step 7: Check the expansion properties of the graph by computing the Laplacian matrix
  degree_matrix <- diag(degree_seq)
  laplacian_matrix <- degree_matrix - adj_matrix
  eigenvalues <- eigen(laplacian_matrix, symmetric=TRUE)$values
  eigenvalue_gap <- abs(max(eigenvalues)) - abs(min(eigenvalues))
  
  # Return the adjacency matrix and Laplacian matrix
  return(list(adj_matrix=adj_matrix, laplacian_matrix=laplacian_matrix, eigenvalue_gap=eigenvalue_gap))
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(Matrix)
library(RANN)

build_spatial_graph <- function(N, lambda, L, r, beta, p, probWithin, probBetween, use_weight=TRUE){
  
  # Generate set of points using Poisson point process
  X <- matrix(runif(N*L), ncol = L)
  X <- X - floor(X) # map to torus
  
  # Calculate distance matrix between all pairs of points
  D <- RANN::nn2(X, X, k = N, searchtype = "default", algorithm = "kd_tree", 
                 delta = 1e-04, eps = 0, verbose = FALSE)$nn.dist
  D <- D + diag(N)*100 # set diagonal to high value to prevent self-loops
  
  # Create adjacency matrix with edge weights favoring short spatial distance and high degree nodes
  A <- matrix(0, nrow=N, ncol=N)
  for(i in 1:N){
    for(j in i:N){
      if(i != j){
        if(D[i,j] <= r){
          if(runif(1) <= p){
            A[i,j] <- 1
            A[j,i] <- 1
          }
        }else{
          prob <- 1/(D[i,j]^beta)
          if(runif(1) <= prob){
            A[i,j] <- 1
            A[j,i] <- 1
          }
        }
      }
    }
  }
  
  # Add community structures to the graph
  communities <- cutree(igraph::graph_from_adjacency_matrix(A), no = round(sqrt(N)))
  for(i in 1:N){
    for(j in i:N){
      if(communities[i] == communities[j]){
        if(runif(1) <= probWithin){
          A[i,j] <- 1
          A[j,i] <- 1
        }
      }else{
        if(runif(1) <= probBetween){
          A[i,j] <- 1
          A[j,i] <- 1
        }
      }
    }
  }
  
  # Add rewiring probability for small world effect
  for(i in 1:N){
    for(j in (i+1):N){
      if(A[i,j] == 1){
        if(runif(1) <= p){
          A[i,j] <- 0
          A[j,i] <- 0
          candidates <- which(A[i,] == 0 & D[i,] <= r)
          if(length(candidates) > 0){
            A[i,sample(candidates, 1)] <- 1
            A[sample(candidates, 1),i] <- 1
          }
        }
      }
    }
  }
  
  # Create sparse matrix representation of the adjacency matrix
  if(use_weight){
    A[A>0] <- 1/D[A>0]
  }
  A <- Matrix(A, sparse = TRUE)
  
  # Compute Laplacian matrix to check expansion properties
  Lap <- Matrix::laplacianMatrix(A, normalized = TRUE)
  lambda <- eigen(Lap, symmetric = TRUE, only.values = TRUE)$values
  
  return(list(A = A, lambda = lambda))
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

generate_spatial_expander_graph <- function(N, lambda, L, r, p, probWithin, probBetween, beta) {
  
  # Generate N points using 2D Poisson point process with intensity parameter lambda on an L-dimensional torus
  points <- matrix(runif(N*L, min=0, max=1), ncol=L)
  distances <- dist(points, method="euclidean", diag=TRUE, upper=TRUE)
  
  # Create the adjacency matrix using a function that connects nodes with a probability favoring short spatial distance
  adjacency_matrix <- ifelse(distances <= r, 1, 0)
  degree_vector <- rowSums(adjacency_matrix)
  for (i in 1:N) {
    for (j in i+1:N) {
      if (runif(1) < p) {
        if (runif(1) < probWithin && degree_vector[i] > 1 && degree_vector[j] > 1) {
          adjacency_matrix[i,j] <- 1
          adjacency_matrix[j,i] <- 1
        } else if (runif(1) < probBetween && degree_vector[i] > 1 && degree_vector[j] > 1) {
          adjacency_matrix[i,j] <- 1
          adjacency_matrix[j,i] <- 1
        }
      }
    }
  }
  
  # Add a preferential attachment power law parameter to the graph to create a scale-free structure that follows the power law distribution
  for (i in seq(3, N)) {
    weights <- degree_vector[1:(i-1)]^beta
    probs <- weights / sum(weights)
    j <- sample(1:(i-1), size=1, prob=probs)
    adjacency_matrix[i,j] <- 1
    adjacency_matrix[j,i] <- 1
  }
  
  # Compute Laplacian matrix to check expansion properties of the graph
  degree_matrix <- diag(degree_vector)
  laplacian_matrix <- degree_matrix - adjacency_matrix
  
  # Return the adjacency matrix and Laplacian matrix
  return(list(adjacency_matrix=adjacency_matrix, laplacian_matrix=laplacian_matrix))
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to generate a spatial expander graph
spatial_expander_graph <- function(N, lambda, L, r, p, probWithin, probBetween, beta, use_edge_weights = FALSE) {
  
  # Generate N points using a 2D Poisson point process on an L-dimensional torus
  points <- matrix(runif(N*L), ncol = L)
  points <- points - floor(points) # Wrap around the torus
  
  # Calculate the distance matrix between all pairs of points
  distances <- as.matrix(dist(points, diag = TRUE, upper = TRUE))
  
  # Use the distance matrix to create an adjacency matrix for the graph
  adjacency <- matrix(0, nrow = N, ncol = N)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (distances[i,j] <= r && runif(1) > p) { # Connect nodes with a probability favoring short spatial distance
        if (use_edge_weights) { # Use edge weights based on distances
          adjacency[i,j] <- distances[i,j]
          adjacency[j,i] <- distances[i,j]
        } else { # Use binary edge weights
          adjacency[i,j] <- 1
          adjacency[j,i] <- 1
        }
      }
    }
  }
  
  # Add community structures to the graph
  community_sizes <- rpois(ceiling(sqrt(N)), N/2) # Poisson distribution for number of nodes in each community
  cum_sizes <- c(0, cumsum(community_sizes)) # Cumulative sizes of communities
  community_labels <- rep(1:length(community_sizes), times = community_sizes) # Labels for each node's community
  within_probs <- matrix(0, nrow = N, ncol = N)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (community_labels[i] == community_labels[j]) { # Connect nodes within communities
        within_probs[i,j] <- probWithin
        within_probs[j,i] <- probWithin
      } else { # Connect nodes between communities
        within_probs[i,j] <- probBetween
        within_probs[j,i] <- probBetween
      }
    }
  }
  within_probs[which(within_probs == 0)] <- NA # Set non-connections to NA
  within_probs[which(adjacency == 1)] <- NA # Set existing connections to NA
  within_probs[which(is.na(within_probs))] <- 0 # Set NA values to 0
  for (i in 1:N) {
    within_probs[i,i] <- 0 # Set diagonal to 0
  }
  adjacency <- adjacency + (runif(N*N) < within_probs) # Add connections based on within node probability
  
  # Add preferential attachment power law parameter
  degrees <- apply(adjacency, 1, sum)
  for (i in (N+1):N+50) { # Add 50 new nodes
    degree_probs <- degrees^beta / sum(degrees^beta)
    new_neighbors <- sample(1:N, size = 5, prob = degree_probs, replace = TRUE) # Choose 5 existing nodes to connect to
    adjacency <- rbind(adjacency, rep(0, N)) # Add new row
    adjacency <- cbind(adjacency, rep(0, N+1)) # Add new column
    adjacency[i,new_neighbors] <-
      
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+smple code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Generate N points using a 2D Poisson point process with intensity parameter lambda on an L-dimensional torus. You can use the spatstat package in R to generate a 2D Poisson point process with a toroidal boundary. The code snippet for generating the points is as follows:
      R
    Copy code
    library(spatstat)
    set.seed(123)
    n <- 100
    lambda <- 0.1
    L <- c(10, 10) # Dimensions of the torus
    pts <- rToroidalPoisson(n, lambda, L)
    Calculate the distance matrix between all pairs of points. You can use the dist function in R to calculate the Euclidean distance between all pairs of points.
    R
    Copy code
    dist_mat <- dist(pts)
    Use the distance matrix to create an adjacency matrix for the graph with arguments to use edge weights if specified. You can use the following function to create the adjacency matrix based on the distance matrix.
    R
    Copy code
    adj_matrix <- function(dist_mat, r, p = 0) {
      n <- nrow(dist_mat)
      adj_mat <- matrix(0, n, n)
      for (i in 1:(n-1)) {
        for (j in (i+1):n) {
          if (dist_mat[i,j] < r && runif(1) > p) {
            adj_mat[i,j] <- adj_mat[j,i] <- 1
          }
        }
      }
      return(adj_mat)
    }
    Add a rewiring probability p for small world effect to the graph to enhance the connectivity. You can use the following function to rewire the edges in the graph with a probability p.
    R
    Copy code
    rewire_edges <- function(adj_mat, p) {
      n <- nrow(adj_mat)
      for (i in 1:(n-1)) {
        for (j in (i+1):n) {
          if (adj_mat[i,j] == 1 && runif(1) < p) {
            # Find a random neighbor of node i that is not already connected to i
            neighbor_i <- which(adj_mat[i,] == 0 & seq_len(n) != i)
            if (length(neighbor_i) > 0) {
              idx <- sample(neighbor_i, 1)
              # Swap the edge (i,j) with (i,idx)
              adj_mat[i,j] <- adj_mat[j,i] <- 0
              adj_mat[i,idx] <- adj_mat[idx,i] <- 1
            }
          }
        }
      }
      return(adj_mat)
    }
    Add community structures to the graph by connecting nodes within communities with within node probability probWithin more than between communities with probability probBetween. You can use the following function to create a community structure in the graph.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Generate a set of N points using a 2D Poisson point process on an L-dimensional torus
generate_points <- function(N, lambda, L) {
  # Generate uniformly distributed points on the torus
  x <- matrix(runif(N*L), ncol=L)
  # Apply the torus periodicity
  x <- x %% 1
  # Generate a Poisson point process on the torus
  ppp <- spatstat::rpoispp(lambda, win=spatstat::owin(rep(0, L), rep(1, L)), nsim=N)
  # Combine the uniform and Poisson points
  points <- rbind(x, spatstat::as.data.frame.ppp(ppp))
  return(points)
}

# Calculate the distance matrix between all pairs of points
distance_matrix <- function(points, L) {
  N <- nrow(points)
  # Compute the distance matrix between all pairs of points
  d <- as.matrix(dist(points, diag=TRUE, upper=TRUE))
  # Apply the torus periodicity
  for (i in 1:N) {
    for (j in (i+1):N) {
      for (k in 1:L) {
        d[i,j] <- min(d[i,j], abs(points[i,k]-points[j,k]), abs(1-points[i,k]+points[j,k]), abs(points[i,k]-1-points[j,k]))
      }
      d[j,i] <- d[i,j]
    }
  }
  return(d)
}

# Create the adjacency matrix for the graph
create_adjacency_matrix <- function(d, r, p) {
  N <- nrow(d)
  # Create an empty adjacency matrix
  A <- matrix(0, nrow=N, ncol=N)
  # Connect nodes with a probability favoring short spatial distance with radius parameter 'r'
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (d[i,j] < r && runif(1) < exp(-d[i,j]/r)) {
        A[i,j] <- 1
        A[j,i] <- 1
      }
    }
  }
  # Add a rewiring probability 'p' for small world effect
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (A[i,j] == 0 && runif(1) < p) {
        A[i,j] <- 1
        A[j,i] <- 1
      }
    }
  }
  return(A)
}

# Add community structures to the graph
add_community_structure <- function(A, probWithin, probBetween, numCommunities) {
  N <- nrow(A)
  # Compute the community membership for each node
  membership <- rep(1:numCommunities, each=ceiling(N/numCommunities))[1:N]
  # Connect nodes within communities with within node probability 'probWithin'
  for (c in 1:numCommunities) {
    members <- which(membership == c)
    for (i in members) {
      for (j in members) {
        if (i < j && A[i,j] == 0 && runif(1) < probWithin) {
          A[i,j] <- 1
          A[j,i] <- 1
        }
      }
    }
  }
  # Connect nodes between communities with between node probability 'prob
  



code in R to construct spatial expander propagation graph on an L-dimensional torus. This code should be a single function built from scratch with the following properties                                                       
(1)  Nodes distributed spatially with underlying spatial structure on the torus using 2D poisson point process. 
(2) Graph should have arguments to suggest if to use edge weights or no edge weights for the creation of adjacency matrix of points
(3)
(3) nodes in graph should be connected to each other with a probability favouring short spatial distance and higher degree and also have some rewiring probability for small world effect and have community structures (where we have more connections within communities than between communities using a rewiring probability).
(4) graph should have preferential attachment parameter for scale free ,  distance parameter for short spatial distances, rewiring probability parameter for small world, and community structures parameters that we can tune to generate different types of models.
(5) over all, graph should  be an expander graph with these properties and demonstrate strong expansion properties with check with laplacian matrix
(6)  This code should not use igraph. and must be built from scratch


#########################################################################################
# Python
#########################################################################################
import numpy as np
import networkx as nx

def generate_graph(N, L, r, beta, p, probWithin, probBetween, edge_weight=None):
  # Step 1: Generate points with a 2D Poisson point process
  lambd = N / (L ** 2)  # intensity parameter
points = np.random.uniform(low=0, high=L, size=(N, 2))  # 2D torus

# Step 2: Calculate the distance matrix
dist_matrix = np.zeros((N, N))
for i in range(N):
  for j in range(i+1, N):
  dx = np.abs(points[i][0] - points[j][0])
dy = np.abs(points[i][1] - points[j][1])
dx = np.minimum(dx, L-dx)  # torus distance
dy = np.minimum(dy, L-dy)
dist_matrix[i][j] = dist_matrix[j][i] = np.sqrt(dx**2 + dy**2)

# Step 3: Create the adjacency matrix with preferential attachment
G = nx.Graph()
degree = np.zeros(N)
for i in range(N):
  G.add_node(i)
for j in range(i):
  if np.random.rand() < r * np.exp(-beta * dist_matrix[i][j]):
  if edge_weight:
  G.add_edge(i, j, weight=edge_weight(dist_matrix[i][j]))
else:
  G.add_edge(i, j)
degree[i] += 1
degree[j] += 1

# Step 4: Rewire edges with probability p
for i in range(N):
  for j in G.neighbors(i):
  if np.random.rand() < p:
  k = np.random.randint(N)
while k == i or G.has_edge(i, k):
  k = np.random.randint(N)
G.remove_edge(i, j)
if edge_weight:
  G.add_edge(i, k, weight=edge_weight(dist_matrix[i][k]))
else:
  G.add_edge(i, k)
degree[i] += 1
degree[j] -= 1

# Step 5: Add community structure
num_communities = int(np.sqrt(N))
community_size = int(N / num_communities)
for i in range(num_communities):
  community = np.arange(i*community_size, (i+1)*community_size)
for u in community:
  for v in community:
  if u < v and np.random.rand() < probWithin:
  if edge_weight:
  G.add_edge(u, v, weight=edge_weight(dist_matrix[u][v]))
else:
  G.add_edge(u, v)
degree[u] += 1
degree[v] += 1
for j in range(i):
  other_community = np.arange(j*community_size, (j+1)*community_size)
for u in community:
  for v in other_community:
  if np.random.rand() < probBetween:
  if edge_weight:
  G.add_edge(u, v, weight=edge_weight(dist_matrix[u][v]))
else:
  G.add_edge(u, v)
degree[u] += 1
degree[v] += 1

# Step 6: Compute the Laplacian matrix and check expansion properties
laplacian = nx.laplacian_matrix(G).todense()
eigenvalues = np.linalg.eigvals(laplacian)
eigenvalues.sort()


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(Matrix)
library(igraph)

set.seed(123) # For reproducibility

# Function to generate N points using a 2D Poisson point process on an L-dimensional torus
generate_points <- function(N, lambda, L) {
  d <- L/2
  unlist(lapply(1:L, function(x) seq(-d, d, length.out = sqrt(N))))
  coords <- matrix(sample(runif(N*2*d), N*2), ncol = 2)
  coords <- coords - 0.5
  coords <- coords * L
  return(coords)
}

# Function to create an adjacency matrix for the graph with scale-free preferential attachment
# and rewiring probability p
create_adjacency_matrix <- function(N, beta, distance_matrix, p) {
  A <- matrix(0, nrow = N, ncol = N)
  for (i in 1:N) {
    for (j in 1:N) {
      if (i == j) next
      d_ij <- distance_matrix[i,j]
      p_ij <- 1/(d_ij^beta)
      if (runif(1) < p_ij) {
        A[i,j] <- 1
        A[j,i] <- 1
      }
    }
  }
  A <- igraph::graph_from_adjacency_matrix(A, mode = "undirected")
  A <- igraph::simplify(A)
  A <- igraph::rewire(A, each_edge(prob = p))
  A <- as_adjacency_matrix(A)
  return(A)
}

# Function to create a community structure in the graph
create_communities <- function(A, probWithin, probBetween, num_communities) {
  V <- igraph::V(igraph::graph_from_adjacency_matrix(A))
  comm <- sample(x = 1:num_communities, size = igraph::vcount(V), replace = TRUE)
  for (i in 1:num_communities) {
    V_in <- igraph::V(A)[comm == i]
    V_out <- igraph::V(A)[comm != i]
    E <- get.edgelist(subgraph.edges(A, which(V_in %--% V_out)))
    for (j in 1:nrow(E)) {
      if (runif(1) < probBetween) {
        A[E[j,1], E[j,2]] <- 0
        A[E[j,2], E[j,1]] <- 0
      }
    }
    E <- get.edgelist(subgraph.edges(A, which(V_in %--% V_in)))
    for (j in 1:nrow(E)) {
      if (runif(1) < probWithin) {
        A[E[j,1], E[j,2]] <- 1
        A[E[j,2], E[j,1]] <- 1
      }
    }
  }
  return(A)
}

# Function to compute the Laplacian matrix
compute_laplacian <- function(A) {
  D <- diag(x = colSums(A))
  L <- D - A
  return(L)
}

# Function to generate the graph
generate_graph <- function(N, lambda, L, beta, p, probWithin, probBetween, num_communities) {
  coords <- generate_points(N, lambda, L)
  distance_matrix <- as.matrix(dist(coords, method = "euclidean", diag = TRUE, upper = TRUE))
  A <- create_adjacency_matrix(N, beta, distance_matrix, p)
  A <- create_communities(A, probWithin,
                          
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(spatstat)
library(Matrix)

# Function to generate a Spatial Expander Graph
spatial_expander_graph <- function(N, lambda, L, beta, r, p, probWithin, probBetween, use_edge_weights = FALSE) {
  
  # Generate a 2D Poisson point process on an L-dimensional torus
  points <- rToroidalPoissDisc(N, lambda = lambda, dim = L)
  
  # Calculate the distance matrix between all pairs of points
  dist_mat <- as.matrix(dist(points))
  
  # Create an adjacency matrix for the graph
  adj_mat <- Matrix(0, nrow = N, ncol = N, sparse = TRUE)
  
  # Connect nodes with a probability favoring short spatial distance and a scale free preferential attachment degree distribution
  for (i in 1:N) {
    for (j in (i+1):N) {
      if (runif(1) < r / (dist_mat[i, j]^beta)) {
        if (use_edge_weights) {
          adj_mat[i, j] <- adj_mat[j, i] <- dist_mat[i, j]
        } else {
          adj_mat[i, j] <- adj_mat[j, i] <- 1
        }
      }
    }
  }
  
  # Add rewiring probability p for small world effect
  for (i in 1:N) {
    for (j in (i+1):N) {
      if (runif(1) < p) {
        if (use_edge_weights) {
          adj_mat[i, j] <- adj_mat[j, i] <- dist_mat[i, j]
        } else {
          adj_mat[i, j] <- adj_mat[j, i] <- 1
        }
      }
    }
  }
  
  # Add community structures to the graph
  community_labels <- cutree(cluster_walktrap(graph_from_adjacency_matrix(adj_mat)), k = N/10) # use walktrap algorithm to detect communities
  community_members <- split(1:N, community_labels)
  for (i in 1:length(community_members)) {
    for (j in (i+1):length(community_members)) {
      if (runif(1) < probBetween) {
        # connect between communities
        for (k in community_members[[i]]) {
          for (l in community_members[[j]]) {
            if (use_edge_weights) {
              adj_mat[k, l] <- adj_mat[l, k] <- dist_mat[k, l]
            } else {
              adj_mat[k, l] <- adj_mat[l, k] <- 1
            }
          }
        }
      } else {
        # connect within communities
        for (k in community_members[[i]]) {
          for (l in community_members[[i]]) {
            if (k != l && runif(1) < probWithin) {
              if (use_edge_weights) {
                adj_mat[k, l] <- adj_mat[l, k] <- dist_mat[k, l]
              } else {
                adj_mat[k, l] <- adj_mat[l, k] <- 1
              }
            }
          }
        }
        for (k in community_members[[j]]) {
          for (l in community_members[[j]]) {
            if (k != l && runif(1) < probWithin) {
              if (use_edge_weights) {
                adj_mat[k,
                        
                        
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(Matrix)
library(RANN)
library(Matrix.utils)

# Function to generate a 2D Poisson point process on an L-dimensional torus
generate_points <- function(N, lambda, L) {
  # Calculate the size of the torus
  torus_size <- rep(1, L)
  # Generate points using the Poisson point process
  points <- matrix(rpois(N * L, lambda), ncol = L, byrow = TRUE) / lambda
  # Map points to the torus
  points <- points %% torus_size
  # Return points as a matrix
  return(points)
}

# Function to calculate the distance matrix
calculate_distance_matrix <- function(points, L) {
  # Calculate the pairwise distances between points
  distances <- RANN::nn2(points, k = N)$nn.dist
  # Convert distances to a distance matrix
  distance_matrix <- matrix(distances, nrow = N, ncol = N)
  # Adjust for the torus topology
  for (i in 1:N) {
    for (j in 1:N) {
      for (k in 1:L) {
        diff <- abs(points[i, k] - points[j, k])
        if (diff > 0.5) {
          diff <- 1 - diff
        }
        distance_matrix[i, j] <- distance_matrix[i, j] + diff ^ 2
      }
    }
  }
  distance_matrix <- sqrt(distance_matrix)
  # Return the distance matrix
  return(distance_matrix)
}

# Function to create an adjacency matrix
create_adjacency_matrix <- function(distance_matrix, beta, r, p) {
  # Create an empty adjacency matrix
  adjacency_matrix <- Matrix(0, nrow = N, ncol = N, sparse = TRUE)
  # Calculate the degree of each node
  degree <- colSums(distance_matrix < r)
  # Calculate the probabilities of connecting to each node
  probabilities <- (degree + 1) ^ (-beta)
  probabilities <- probabilities / sum(probabilities)
  # Connect nodes with a probability favoring short spatial distance
  for (i in 1:N) {
    neighbors <- which(distance_matrix[i, ] < r)
    probabilities_i <- probabilities[neighbors]
    probabilities_i <- probabilities_i / sum(probabilities_i)
    connect <- sample(neighbors, size = 1, prob = probabilities_i)
    if (runif(1) < p) {
      neighbors <- setdiff(neighbors, connect)
      probabilities_i <- probabilities_i[-match(connect, neighbors)]
      probabilities_i <- probabilities_i / sum(probabilities_i)
      connect <- sample(neighbors, size = 1, prob = probabilities_i)
    }
    if (connect > 0) {
      adjacency_matrix[i, connect] <- 1
      adjacency_matrix[connect, i] <- 1
    }
  }
  # Return the adjacency matrix
  return(adjacency_matrix)
}

# Function to add community structures
add_community_structure <- function(adjacency_matrix, probWithin, probBetween, num_communities) {
  # Create community labels
  community_labels <- rep(1:num_communities, each = floor(N / num_communities))
  community_labels <- c(community_labels, rep(num_communities, N - length(community_labels)))
  # Calculate the community matrix
  community_matrix <- outer(community_labels, community_labels, "==")
  # Add edges within communities
  adjacency_matrix[community_matrix & adjacency_matrix == 0] <- runif(sum(community_matrix & adjacency_matrix == 0)) < probWithin
  # Add
  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  library(Matrix)
  library(spatstat)
  
  # Set parameters
  beta <- 1.5
  r <- 0.1
  N <- 500
  lambda <- 10
  L <- c(1,1)
  p <- 0.1
  probWithin <- 0.8
  probBetween <- 0.2
  
  # Generate points using Poisson point process
  pts <- rpoispp(lambda*N, win = owin(rep(0,2),L))
  pts <- as.data.frame(coordinates(pts))
  
  # Calculate distance matrix
  dist_mat <- as.matrix(dist(pts))
  
  # Create adjacency matrix using preferential attachment and small world rewiring
  adj_mat <- matrix(0, nrow=N, ncol=N)
  adj_mat[1,2] <- 1
  for (i in 3:N){
    degree <- rowSums(adj_mat[1:(i-1),1:(i-1)])
    prob <- (dist_mat[1:(i-1),i]^(-beta))*degree/(sum(degree)+r*i^(1+p))
    prob_within <- prob*probWithin
    prob_between <- prob*(1-probWithin)/(N*(N-1)/2-(i-1)*degree*(1-probWithin)/(2*(1-probWithin)))
    prob_mat <- matrix(c(rep(prob_within, i-2), prob_between), ncol=i-1)
    for (j in 1:(i-1)){
      if (runif(1) < p) prob_mat[,j] <- prob_between
      prob_mat[j,j] <- 0
    }
    chosen <- sample.int(i-1, size=sum(rbinom(i-1, 1, prob_within)+rbinom((i-1)*(i-2)/2, 1, prob_between)))
    adj_mat[chosen,i] <- 1
    adj_mat[i,chosen] <- 1
  }
  
  # Add community structure to the graph
  num_communities <- 5
  community_size <- floor(N/num_communities)
  communities <- rep(1:num_communities, each=community_size)
  if (N > num_communities*community_size) communities <- c(communities, rep(num_communities, N-num_communities*community_size))
  for (i in 1:N){
    for (j in (i+1):N){
      if (communities[i] == communities[j]){
        if (runif(1) < probWithin) {
          adj_mat[i,j] <- 1
          adj_mat[j,i] <- 1
        }
      } else {
        if (runif(1) < probBetween) {
          adj_mat[i,j] <- 1
          adj_mat[j,i] <- 1
        }
      }
    }
  }
  
  # Create Laplacian matrix and check eigenvalue gap
  diag_deg <- diag(rowSums(adj_mat))
  lap_mat <- diag_deg - adj_mat
  eigenvals <- eigen(lap_mat)$values
  eigen_gap <- min(diff(sort(eigenvals)))
  
  # Convert adjacency matrix to sparse matrix format
  sparse_adj_mat <- Matrix(adj_mat, sparse=TRUE)
  
  # Print results
  print(paste0("Number of nodes: ", N))
  print(paste0("Number of edges: ", sum(adj_mat)/2))
  print(paste0("Minimum eigenvalue gap: ", eigen_gap))
  