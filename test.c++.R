setwd("C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi-New")
library(Rcpp)
library(devtools)
#library(RcppArmadillo)
#include <Rcpp.h>
#Rcpp::sourceCpp('c++functions.cpp')
Rcpp::sourceCpp('manet.cpp')



mat <- matrix(0, m, 3)

# generate the upper triangle and lower triangle
upper_tri <- upper.tri(mat)
lower_tri <- lower.tri(mat)

# randomly select m nodes to fill in
nodes <- sample(sum(upper_tri), m, replace = FALSE)

# fill in the selected nodes with random values
mat[upper_tri][nodes] <- runif(m)
mat[lower_tri][nodes] <- mat[upper_tri][nodes]

# make the matrix symmetric
mat <- t(mat)


m <- 5  # number of nodes
n <- 7  # number of columns

# create empty matrix with n columns and m rows
mat <- matrix(0, nrow = m, ncol = n)
mat[lower.tri(mat, diag = FALSE)] <- 1
mat[upper.tri(mat)]=mat[t(lower.tri(mat))]

mat[upper.tri(mat)] <- t(mat[upper.tri(t(mat))])
mat


m <- 4 # number of rows (nodes)
n <- 7 # number of columns

# create an empty m x n matrix
mat <- matrix(0, m, n)

# assign 1 to nodes 1:m and m:1
mat[1:m, m:1] <- 1

# make the matrix symmetric
mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]

# print the resulting matrix
mat

# 
# m <- 3
# n <- 6
# # create an empty adjacency matrix of size m x n
# adj_matrix <- matrix(0, m, n)
# # assign nodes 1:m and m:1 the number 1
# adj_matrix[1:m, (n-(m+1)):n] <- 1
# adj_matrix[(n-(m+1)):n, 1:m] <- 1
# 
# # set the diagonal to zero
# diag(adj_matrix) <- 0
# # print the resulting matrix
# adj_matrix
# 
# 
# M <- 3 # number of rows
# N <- 5 # number of columns
# 
# # create an empty M x N adjacency matrix
# adj_mat <- matrix(0, M, N)
# adj_mat[1:m,m:1]=1
# diag(adj_mat) <- 0
# # assign edge values to the upper triangular part
# # adj_mat[lower.tri(adj_mat, diag = FALSE)] <- 1
# # # assign edge values to the lower triangular part (excluding the diagonal)
# # adj_mat[upper.tri(adj_mat)] <- upper.tri(t(adj_mat))[lower.tri(adj_mat)]
# # # assign diagonal values
# # diag(adj_mat) <- 0
# # adj_mat
# 
# adj_mat <- matrix(0, m, n)
# # initialize the upper triangle with 1s
# adj_mat[upper.tri(adj_mat)] <- t(adj_mat)[upper.tri(adj_mat)]
# 
# adj_mat[upper.tri(adj_mat)] <- 1
# # initialize the lower triangle with 1s
# adj_mat[lower.tri(adj_mat)] <- 1
# 
# edges <- matrix(0, nrow = m, ncol = n)
# edges[lower.tri(edges)] <- 1
# edges <- edges[lower.tri(edges)] #inititialise adj matrix
# 
# # res <- matrix(0, m, m)
# # res.upper <- 1
# # res[lower.tri(res)] <- 1
# # rm(res.upper)
# # diag(res) <- 0
# # res[upper.tri(res)]  <- t(res)[upper.tri(res)]
# # 
# 
# 
# # Define the number of nodes in the network
# N <- 100
# 
# # Define the power law parameter (typically between 2 and 3)
# alpha <- 2.5
# 
# # Define the number of nodes to add at each instance
# m <- 5
# 
# # Create an initial connected network of m nodes
# adj_matrix <- matrix(0, nrow = m, ncol = N)
# for (i in 1:m) {
#   for (j in 1:m) {
#     if (i != j) {
#       adj_matrix[i, j] <- 1
#     }
#   }
# }
# degree_seq <- rep(m-1, m-1)
# 
# # Add nodes to the network
# for (i in (m+1):N) {
#   # Calculate the probability of attaching to each existing node
#   prob <- degree_seq^alpha / sum(degree_seq^alpha)
#   # Choose m nodes to attach to
#   attach_to <- sample(1:(i-1), m, replace = TRUE, prob = prob)
#   # Add the new node to the adjacency matrix and degree sequence
#   adj_matrix <- rbind(adj_matrix, matrix(0, ncol = N, nrow = 1))}
#   adj_matrix[i, attach_to] <- 1
#   adj_matrix[attach_to, i] <- 1
# }
#   degree_seq[i] <- m
#   degree_seq[attach_to] <- degree_seq[attach_to] + 1
# }
# 
# # Check the degree sequence of the final network
# hist(degree_seq, breaks = max(degree_seq))

# Function to generate a scale-free network
# with N nodes, power law parameter 'alpha', and
# the number of nodes to add at each instance 'm'
scale_free_network <- function(N, alpha, m) {
  
  # Initialize the network with m nodes
  nodes <- 1:m
  edges <- matrix(0, ncol = N, nrow = N)
  for (i in 1:m-1) {
    for (j in (i+1):m) {
      edges[i,j] <- edges[j,i] <- 1
    }
  }
  
  # Generate the rest of the network using the BA model
  for (i in (m+1):N) {
    # Calculate the degree distribution of existing nodes
    deg <- apply(edges[,1:i-1], 1, sum)
    # Calculate the degree probabilities
    prob <- deg^(-alpha)
    prob <- prob / sum(prob)
    # Choose m nodes to connect to
    neighbors <- sample(1:(i-1), m, replace = TRUE, prob = prob)
    # Add edges between the new node and its neighbors
    edges[i,neighbors] <- edges[neighbors,i] <- 1
  }
  
  
  
  
  # Return the adjacency matrix
  return(edges)
}

# Example usage
N <- 50
alpha <- 2.5
m <- 2
network <- scale_free_network(N, alpha, m)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set the number of nodes and the initial number of nodes to connect to
N <- 100 # number of nodes
m <- 5 # initial number of nodes to connect to

# Create a vector to store the degree of each node
degrees <- rep(0, N)

# Create an adjacency matrix for the first m nodes
adj_mat <- matrix(0, nrow = m, ncol = m)
for (i in 1:m-1) {
  for (j in (i+1):m) {
    adj_mat[i, j] <- 1
    adj_mat[j, i] <- 1
  }
  degrees[i] <- m - 1
}

# Add nodes to the network one at a time
for (i in (m+1):N) {
  # Calculate the probability of connecting to each existing node
  prob <- degrees[1:(i-1)] / sum(degrees[1:(i-1)])
  
  # Choose m nodes to connect to, with probability proportional to their degree
  neighbors <- sample(1:(i-1), m, replace = FALSE, prob = prob)
  
  # Add the new node and its connections to the adjacency matrix
  adj_mat[i, neighbors] <- 1
  adj_mat[neighbors, i] <- 1
  
  # Update the degree of each node
  degrees[i] <- m
  degrees[neighbors] <- degrees[neighbors] + 1
}

# Convert the adjacency matrix to an igraph object
library(igraph)
g <- graph_from_adjacency_matrix(adj_mat)

# Plot the graph
plot(g)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set the values for N, alpha, and m
N <- 100 # number of nodes
alpha <- 2 # power law parameter
m <- 5 # number of nodes to add at each instance

# Create a graph with m nodes
adj_mat <- matrix(0, nrow = m, ncol = m)
for (i in 1:m) {
  for (j in (i+1):m) {
    adj_mat[i, j] <- 1
    adj_mat[j, i] <- 1
  }
}

# Generate the degree sequence using the power law distribution
degrees <- numeric(N)
degrees[1:m] <- m - 1
for (i in (m+1):N) {
  degree_sum <- sum(degrees[1:i-1]^(-alpha))
  prob <- (m-1) / degree_sum
  for (j in 1:i-1) {
    if (runif(1) < prob) {
      adj_mat[i, j] <- 1
      adj_mat[j, i] <- 1
      degrees[i] <- degrees[i] + 1
      degrees[j] <- degrees[j] + 1
    }
  }
}

# Create the adjacency matrix for the full graph
adj_mat_full <- matrix(0, nrow = N, ncol = N)
adj_mat_full[1:m, 1:m] <- adj_mat

# Add remaining nodes to the graph
for (i in (m+1):N) {
  degrees_subset <- degrees[1:i-1]
  available_nodes <- which(degrees_subset >= m)
  if (length(available_nodes) < m) {
    m_new <- m - length(available_nodes)
    degrees_subset_sorted <- sort(degrees_subset, decreasing = TRUE)
    high_degree_nodes <- which(degrees_subset == degrees_subset_sorted[1:m_new])
    available_nodes <- c(available_nodes, high_degree_nodes)
  }
  selected_nodes <- sample(available_nodes, m, replace = FALSE)
  adj_mat_full[i, selected_nodes] <- 1
  adj_mat_full[selected_nodes, i] <- 1
  degrees[i] <- length(selected_nodes)
  degrees[selected_nodes] <- degrees[selected_nodes] + 1
}

# Convert the adjacency matrix to an igraph object
library(igraph)
g <- graph_from_adjacency_matrix(adj_mat_full)

# Plot the graph
plot(g)


