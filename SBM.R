set.seed(123) # set seed for reproducibility

# function to generate stochastic block matrix
stochastic_block_matrix <- function(n, k, p) {
  # n: number of nodes
  # k: number of clusters
  # p: matrix of probabilities
  
  # create empty adjacency matrix
  adj_matrix <- matrix(0, n, n)
  
  # assign nodes to clusters
  nodes_per_cluster <- rep(n %/% k, k)
  nodes_left <- n %% k
  for (i in 1:nodes_left) {
    nodes_per_cluster[i] <- nodes_per_cluster[i] + 1
  }
  clusters <- rep(1:k, nodes_per_cluster)
  
  # generate edges
  for (i in 1:n) {
    for (j in i:n) {
      if (i != j) {
        # probability of edge between i and j
        prob <- p[clusters[i], clusters[j]]}}
        # generate edge with probability prob
        if (runif(1) < prob) {
          adj_matrix[i, j] <- 1
          adj_matrix[j, i] <- 1
        }
      }
    }
  }
  
  return(adj_matrix)
}

# example usage
n <- 100 # number of nodes
k <- 4 # number of clusters
p <- matrix(c(0.8, 0.2, 0.2, 0.1,
              0.2, 0.8, 0.1, 0.2,
              0.2, 0.1, 0.8, 0.2,
              0.1, 0.2, 0.2, 0.8), nrow=k, ncol=k) # matrix of probabilities
sbm <- stochastic_block_matrix(n, k, p)
print(sbm)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set the number of nodes and the number of communities
n <- 100
K <- 2

# Set the probability matrix P for each community
P <- matrix(c(0.2, 0.1, 0.1, 0.3), nrow = K, ncol = K)

# Generate the block matrix
Z <- matrix(sample(1:K, n, replace = TRUE), nrow = n, ncol = 1)
B <- matrix(0, nrow = n, ncol = n)
for (i in 1:n) {
  for (j in i:n) {
    if (Z[i] == Z[j]) {
      B[i,j] <- rbinom(1, 1, P[Z[i], Z[j]])
      B[j,i] <- B[i,j]
    }
  }
}

# Visualize the block matrix
image(B, col = gray((0:10)/10), axes = FALSE, xlab = "", ylab = "")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Define the parameters
n <- 100  # Number of nodes
k <- 4    # Number of blocks
p <- matrix(runif(k^2), k, k)  # Probability matrix

# Generate the block matrix
block_mat <- matrix(0, n, n)
block_sizes <- rep(n/k, k)
cumulative_sizes <- c(0, cumsum(block_sizes))
for (i in 1:k) {
  for (j in 1:k) {
    block_mat[cumulative_sizes[i]+1:cumulative_sizes[i+1],
              cumulative_sizes[j]+1:cumulative_sizes[j+1]] <- rbinom(block_sizes[i]*block_sizes[j], 1, p[i,j])
  }
}
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+    library(igraph)

# Define the number of nodes in the graph
n <- 100

# Define the block sizes
block_sizes <- c(30, 30, 40)

# Define the probability matrix
P <- matrix(c(0.1, 0.02, 0.01,
              0.02, 0.15, 0.03,
              0.01, 0.03, 0.2), nrow = length(block_sizes), byrow = TRUE)

# Generate the community memberships
community_memberships <- rep(1:length(block_sizes), times = block_sizes)

# Generate the adjacency matrix
A <- matrix(0, nrow = n, ncol = n)
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    if (community_memberships[i] == community_memberships[j]) {
      A[i,j] <- rbinom(1, 1, P[community_memberships[i], community_memberships[j]])
      A[j,i] <- A[i,j]
    }
  }
}

# Convert the adjacency matrix to a graph object
G <- graph.adjacency(A, mode = "undirected")

# Plot the graph
plot(G)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
generate_spatial_scale_free_graph <- function(N, numofComm, probWithin, probBetween, beta, m, r, prewire) {
  # Step 1: Initialize network with a single node
  adj_mat <- matrix(0, nrow = N, ncol = N)
  adj_mat[1, 1] <- 1
  degrees <- rep(0, N)
  community_assignments <- rep(1, N)
  num_communities <- 1
  
  # Step 2: Create initial community/block probability matrix
  P_ij <- matrix(0, nrow = numofComm, ncol = numofComm)
  for (i in 1:numofComm) {
    for (j in 1:numofComm) {
      if (i == j) {
        P_ij[i, j] <- probWithin
      } else {
        P_ij[i, j] <- probBetween
      }
    }
  }
  
  # Step 3-6: Add nodes to the network one by one
  for (i in 2:N) {
    # Step 3i: Create arbitrary number of communities/blocks
    if (i > numofComm) {
      num_communities <- sample(1:numofComm, 1)
    }
    
    # Step 3ii: Define community probability matrix
    P_ij_i <- P_ij[1:num_communities, 1:num_communities]}
    
    # Step 3iii: Define probability of attachment function for within-community connections
    spat_distprobs <- 1 / (1 + (dist_mat[1:(i-1), i])^beta)
    spat_distprobs[dist_mat[1:(i-1), i] > r] <- 0
    degrees_i <- degrees[1:(i-1)]
    deg_probs_within <- (degrees_i^alpha) / sum(degrees_i^alpha)
    P_unique_within <- outer(deg_probs_within, spat_distprobs, "*")
    
    # Step 3iv: Define probability of attachment function for between-community connections
    deg_probs_between <- deg_probs_within
    P_unique_between <- outer(deg_probs_between, rep(prewire, i-1), "*")
    
    # Step 3v: Assign new node to a community using 2D Poisson point process
    if (num_communities == 1) {
      community_assignments[i] <- 1
    } else {
      coords <- matrix(runif(num_communities * 2), ncol = 2)
      distances <- sqrt((coords[,1] - coords[i-1,1])^2 + (coords[,2] - coords[i-1,2])^2)
      lambda <- num_communities * exp(-distances^2 / (2 * r^2))
      community_assignments[i] <- sample(1:num_communities, 1, prob = lambda / sum(lambda))
    }
    
    # Step 3vi: Add edges to the new node based on community/block assignments and probability matrix
    for (j in 1:(i-1)) {
      if (community_assignments[i] == community_assignments[j]) {
        if (runif(1) < P_unique_within[j, i-1] * P_ij_i[community_assignments[i], community_assignments[j]]) {
          adj_mat[i, j] <- 1
          adj_mat[j, i] <- 1
          degrees[i] <- degrees[i] + 1
          degrees[j] <- degrees[j] + 1
            
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          generateSpatialScaleFreeGraph <- function(N, numofComm, beta, m, r, probWithin, probBetween, prewire) {
            # initialize adjacency matrix and distance matrix
            adj_mat <- matrix(0, nrow = N, ncol = N)
            dist_mat <- matrix(0, nrow = N, ncol = N)
            
            # randomly distribute the first node
            x <- runif(1)
            y <- runif(1)
            xy <- c(x, y)
            
            # calculate distances between nodes
            for (i in 1:N) {
              for (j in 1:N) {
                dx <- xy[i, 1] - xy[j, 1]
                dy <- xy[i, 2] - xy[j, 2]
                dist_mat[i, j] <- sqrt(dx^2 + dy^2)
              }
            }
            
            # initialize first node in its own community
            community <- rep(1, N)
            numofComm <- 1
            
            # add new nodes to the graph
            for (i in 2:N) {
              # create random partition of nodes into communities
              community <- kmeans(xy, centers = numofComm)$cluster
              
              # initialize probability matrix for each community
              P_ij <- matrix(0, nrow = numofComm, ncol = numofComm)
              
              # populate probability matrix based on parameters
              for (j in 1:numofComm) {
                for (k in 1:numofComm) {
                  if (j == k) {
                    P_ij[j, k] <- probWithin
                  } else {
                    P_ij[j, k] <- probBetween
                  }
                }
              }
              
              # calculate degree and spatial probabilities for within-community connections
              degrees <- rowSums(adj_mat[1:(i-1), ])
              spat_distprobs <- 1 / (1 + (dist_mat[1:(i-1), i])^beta) * (dist_mat[1:(i-1), i] <= r)
              deg_probs <- (degrees^m) / sum(degrees^m)
              P_unique_within <- P_ij[community[i], community[1:(i-1)]] * deg_probs[1:(i-1)] * spat_distprobs
              
              # calculate degree probability for between-community connections
              deg_probs <- (degrees^m) / sum(degrees^m)
              P_unique_between <- P_ij[community[i],] * deg_probs * prewire
              
              # add edges for within-community connections
              for (j in 1:m) {
                if (sum(P_unique
                        
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
generate_spatial_scalefree_graph <- function(N, numofComm, probWithin, probBetween, beta, m, prewire, r, alpha) {
  # N: number of nodes in the network
  # numofComm: number of communities/blocks
  # probWithin: probability of an edge between nodes in the same community/block
  # probBetween: probability of an edge between nodes in different communities/blocks
  # beta: power law parameter for preferential attachment
  # m: number of edges to add at each step
  # prewire: probability of rewiring an edge
  # r: maximum distance for within-community attachment
  # alpha: scaling parameter for degree probabilities
  
  # initialize the adjacency matrix and the distance matrix
  adj_mat <- matrix(0, nrow = N, ncol = N)
  dist_mat <- matrix(0, nrow = N, ncol = N)
  
  # generate the spatial coordinates of the nodes and compute the distance matrix
  coords <- matrix(runif(N * 2), nrow = N, ncol = 2)
  for (i in 1:N) {
    for (j in (i+1):N) {
      dist_mat[i,j] <- sqrt(sum((coords[i,] - coords[j,])^2))
      dist_mat[j,i] <- dist_mat[i,j]
    }
  }
  
  # assign the first node to a community
  comm_assign <- rep(1, N)
  
  # add nodes to the network
  for (i in 2:N) {
    # create numofComm communities
    comm_probs <- matrix(probBetween, nrow = numofComm, ncol = numofComm)
    diag(comm_probs) <- probWithin
    comm_assign[i] <- sample(numofComm, 1)
    
    # compute the within-community and between-community attachment probabilities
    comm_degrees <- rowSums(adj_mat * (comm_assign == comm_assign[i]))
    comm_spat_dists <- dist_mat[1:(i-1),i] * (comm_assign[1:(i-1)] == comm_assign[i])
    comm_spat_probs <- 1 / (1 + comm_spat_dists^beta)
    comm_spat_probs[comm_spat_dists > r] <- 0
    comm_spat_probs[is.na(comm_spat_probs)] <- 0
    comm_spat_probs <- comm_spat_probs / sum(comm_spat_probs)
    comm_deg_probs <- (comm_degrees^alpha) / sum(comm_degrees^alpha)
    
    within_comm_probs <- comm_spat_probs * comm_deg_probs
    within_comm_probs <- within_comm_probs / sum(within_comm_probs)
    
    between_comm_probs <- matrix(0, nrow = numofComm, ncol = 1)
    for (j in 1:num
         s
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         spatial_scale_free_graph <- function(N, numofComm, probWithin, probBetween, m, beta, alpha, prewire, r) {
           # Initialize the adjacency matrix and assign the first node to a community
           adj_mat <- matrix(0, nrow = N, ncol = N)
           comm_assignments <- rep(1, N)
           degrees <- numeric(N)
           degrees[1] <- 1
           
           # Generate the initial positions of the nodes on a unit square
           positions <- matrix(runif(N*2), ncol = 2)
           dist_mat <- as.matrix(dist(positions))
           
           # Generate the community probability matrix
           P_ij <- matrix(0, nrow = numofComm, ncol = numofComm)
           diag(P_ij) <- probWithin
           P_ij[lower.tri(P_ij)] <- probBetween
           P_ij[upper.tri(P_ij)] <- t(P_ij)[upper.tri(P_ij)]
           
           # For each new node, add edges to existing nodes based on community assignments
           for (i in 3:N) {
             # Step 1: Create arbitrary number of communities
             if (i == 2) {
               comm_assignments[i] <- 1
             } else {
               comm_assignments[i] <- sample(1:numofComm, 1)
             }
             
             # Step 2: Define the community probability matrix
             P_ij_comm <- P_ij[comm_assignments[1:(i-1)], comm_assignments[i]]
             
             # Step 3: Add edges within the same community
             degrees_comm <- degrees[comm_assignments == comm_assignments[i]]}
             spat_distprobs <- 1 / (1 + (dist_mat[1:(i-1), i]^beta))
             spat_distprobs[dist_mat[1:(i-1), i] > r] <- 0
             deg_probs <- (degrees_comm^alpha) / sum(degrees_comm^alpha)
             P_unique_within <- deg_probs * spat_distprobs
             P_unique_within[is.na(P_unique_within)] <- 0}
             neighbors_within <- sample(1:(i-1), m, replace = TRUE, prob = P_unique_within)
             adj_mat[i, neighbors_within] <- 1
             adj_mat[neighbors_within, i] <- 1
             degrees[i] <- sum(adj_mat[i, ])
             }
             
             # Step 4: Add edges between different communities
             if (sum(P_ij_comm == probBetween) > 0) {
               deg_probs <- (degrees^alpha) / sum(degrees^alpha)
               P_unique_between <- deg_probs * prewire
               P_unique_between[comm_assignments[1:(i-1)] == comm_assignments[i]] <- 0
               neighbors_between <- sample(1:(i-1), m, replace = TRUE, prob = P_unique_between)
               adj_mat[i, neighbors_between] <- 1
               adj_mat[neighbors_between, i] <- 1
               degrees[i] <- sum(adj_mat[i, ])
             }
             
             # Step 5: Generate community assignments for new node
             if (sum(P_ij_comm == probBetween) > 0) {
               prob_comm <- P_ij_comm / sum(P_ij_comm)
               comm_assignments[i] <- sample(numofComm, 1, prob = prob_comm)
             } else {
               comm_assignments[i] <- comm_assignments[sample(1:(i-1), 1)]
             }
             
             # Step 6:
             
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  spatial_scale_free_graph <- function(N, numofComm, beta, m, prewire, probWithin, probBetween) {
    # initialize adjacency matrix and community assignment vector
    adj_mat <- matrix(0, nrow = N, ncol = N)
    comm_assign <- rep(1, N)
    
    # initialize first node and set its community assignment
    adj_mat[1, 1] <- 1
    comm_assign[1] <- 1
    
    # initialize distance matrix
    coords <- matrix(runif(N * 2), ncol = 2)
    dist_mat <- as.matrix(dist(coords))
    
    # add new nodes to the graph
    for (i in 2:N) {
      # Step 1: create numofComm communities
      comm_assign[i] <- sample(numofComm, 1)
      
      # Step 2: create community probability matrix
      P_ij <- matrix(probBetween, nrow = numofComm, ncol = numofComm)
      diag(P_ij) <- probWithin
      
      # Step 3: within community attachment
      comm_nodes <- which(comm_assign[1:(i-1)] == comm_assign[i])
      spat_distprobs <- ifelse(dist_mat[comm_nodes, i] <= r, 1 / (1 + dist_mat[comm_nodes, i]^beta), 0)
      degrees <- rowSums(adj_mat[comm_nodes, 1:(i-1)])
      deg_probs <- degrees^alpha / sum(degrees^alpha)
      P_unique_within <- deg_probs * spat_distprobs
      P_unique_within[is.na(P_unique_within)] <- 0
      if (length(comm_nodes) > m) {
        selected_nodes <- sample(comm_nodes, m, prob = P_unique_within, replace = TRUE)
        adj_mat[i, selected_nodes] <- 1
        adj_mat[selected_nodes, i] <- 1
      } else {
        adj_mat[i, comm_nodes] <- 1
        adj_mat[comm_nodes, i] <- 1
        
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        spatial_scale_free_graph <- function(N, numofComm, probWithin, probBetween, m, beta, prewire, r, alpha) {
          
          # Create an adjacency matrix
          adj_mat <- matrix(0, nrow = N, ncol = N)
          
          # Create a distance matrix
          coords <- matrix(runif(N*2), ncol = 2)
          dist_mat <- as.matrix(dist(coords))
          
          # Create a vector to store node degrees
          degrees <- rep(0, N)
          
          # Assign the first node to a community
          communities <- rep(1, N)
          
          # Add new nodes to the network
          for (i in 2:N) {
            
            # Create numofComm communities
            community_sizes <- rep(1, numofComm)
            community_sizes[1:(numofComm-1)] <- sample(1:(N-i+1), numofComm-1)
            community_sizes <- community_sizes / sum(community_sizes) * (N-i+1)
            community_sizes <- cumsum(community_sizes)
            communities[i] <- findInterval(i, community_sizes)
            
            # Create the community probability matrix
            P_ij <- matrix(0, nrow = numofComm, ncol = numofComm)
            diag(P_ij) <- probWithin
            P_ij[lower.tri(P_ij)] <- probBetween
            P_ij[upper.tri(P_ij)] <- t(P_ij)[upper.tri(P_ij)]
            
            # Create the probability matrix for within-community attachment
            spat_distprobs <- 1 / (1 + (dist_mat[1:(i-1), i])^beta)
            spat_distprobs[dist_mat[1:(i-1), i] > r] <- 0
            deg_probs <- (degrees[1:(i-1)]^alpha) / sum(degrees[1:(i-1)]^alpha)
            P_unique_within <- deg_probs * spat_distprobs
            P_within <- P_ij[communities[1:(i-1)], communities[i]] * (1-prewire) + P_unique_within * prewire
            
            # Create the probability matrix for between-community attachment
            deg_probs <- (degrees[1:(i-1)]^alpha) / sum(degrees[1:(i-1)]^alpha)
            P_unique_between <- deg_probs * prewire
            P_between <- P_ij[communities[1:(i-1)], communities[i]] * (1-prewire) + P_unique_between * prewire
            
            # Connect the new node to existing nodes based on community assignments
            for (j in 1:(m/2)) {
              within_comm <- sample(1:(i-1), 1, prob = P_within)
              adj_mat[i, within_comm] <- 1
              adj_mat[within_comm, i] <- 1
              degrees[i] <- degrees[i] + 1
              degrees[within_comm] <- degrees[within_comm] + 1
              
              between_comm <- sample(1:(i-1), 1, prob = P_between)
              adj_mat[i, between_comm] <- 1
              adj_mat[between_comm, i] <- 1
              degrees[i] <- degrees[i] + 1
              degrees[between_comm] <- degrees[between_comm] + 1
            }
          }
          
          # Rewire edges with probability prewire
          for (i in 1:N) {
            for (j in (i+1):
                 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
generate_community_network <- function(N, numofComm, probWithin, probBetween, beta, m, alpha, prewire, r) {
  # Step 1: Initialize the network with a single node
  adj_mat <- matrix(0, nrow = N, ncol = N)
  comm_assignments <- rep(1, N)
  degrees <- rep(0, N)
  adj_mat[1, 1] <- 1
  
  # Step 2: Assign the first node to a community
  comm_assignments[1] <- sample(1:numofComm, 1)
  
  # Loop through the remaining nodes
  for (i in 2:N) {
    # Step 3: Create probability matrices for within and between community attachment
    P_within <- matrix(0, nrow = numofComm, ncol = numofComm)
    P_within[lower.tri(P_within)] <- probWithin
    P_within[upper.tri(P_within)] <- probWithin
    diag(P_within) <- probWithin
    
    P_between <- matrix(probBetween, nrow = numofComm, ncol = numofComm)
    diag(P_between) <- probWithin
    
    # Calculate distance probability for within community attachment
    dist_probs <- ifelse(dist_mat[1:(i-1), i] <= r, 1 / (1 + dist_mat[1:(i-1), i]^beta), 0)
    
    # Calculate scale free degree probability for within community attachment
    deg_probs <- degrees[1:(i-1)]^alpha / sum(degrees[1:(i-1)]^alpha)
    
    # Calculate unique attachment probability for within community attachment
    P_unique_within <- outer(deg_probs, dist_probs, `*`)
    P_unique_within <- P_unique_within / sum(P_unique_within)
    
    # Step 4: Calculate unique attachment probability for between community attachment
    deg_probs <- degrees[1:(i-1)]^alpha / sum(degrees[1:(i-1)]^alpha)
    P_unique_between <- outer(deg_probs, rep(prewire, numofComm), `*`)
    P_unique_between <- P_unique_between / sum(P_unique_between)
    
    # Step 5: Generate community assignment for new node
    comm_probs <- ifelse(1:numofComm == comm_assignments[i-1], probWithin, probBetween)
    comm_probs <- comm_probs / sum(comm_probs)
    comm_assignments[i] <- sample(1:numofComm, 1, prob = comm_probs)
    
    
    
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#     To construct a network with an underlying spatial structure, you can follow these steps:
#       
#       Define a spatial domain: Start by defining the spatial domain where you want to construct your network. This could be a bounded region like a city or a park or an unbounded region like a road network.
#     
#     Generate points: Generate a set of points in the spatial domain that will serve as nodes for the network. You can generate these points using a spatial point process like a Poisson process or a clustered point process that has an underlying spatial structure.
#     
#     Define the network topology: Once you have your set of nodes, you need to define the topology of the network, which is the set of rules that govern how the nodes are connected. This can be done using a distance-based criterion or a more complex spatial algorithm like a Delaunay triangulation.
#     
#     Create edges: Using the topology defined in step 3, create edges that connect the nodes to form the network. The edges can be weighted or unweighted, directed or undirected, depending on the nature of the underlying spatial structure.
#     
#     Add new nodes: To add new nodes to the network, you can use a point process simulation that takes into account the existing network structure. This could be a spatial birth-death process or a preferential attachment model that considers the existing nodes' degree and location.
# 
# Validate the network: Finally, validate the network by checking if it satisfies the desired properties like connectivity, resilience, and degree distribution. You can use various network analysis tools like centrality measures, clustering coefficients, and path lengths to evaluate the network's performance.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++