#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# generate_spatial_scale_free_expander_graph <- function(N, beta, alpha, r, prewire, m) {
#   # Generate initial set of nodes with spatial structure
#   nodes <- data.frame(x = runif(N), y = runif(N))
#   
#   # Initialize adjacency matrix
#   adj_mat <- matrix(0, N, N)
#   
#   # Initialize node degrees
#   degrees <- rep(0, N)
#   
#   # Add initial node to the graph
#   adj_mat[1, 1] <- 1
#   
#   # Grow the rest of the network
#   for (i in 2:N) {
#     # Assign nodes to communities using stochastic block matrix techniques
#     community_assignments <- sample(1:ceiling(sqrt(N)), 1, prob = rep(1:ceiling(sqrt(N)), each = ceiling(N/ceiling(sqrt(N)))))
#     
#     # Calculate attachment probabilities
#     attachment_probs <- rep(0, i-1)
#     for (j in 1:(i-1)) {
#       d <- sqrt((nodes$x[j] - nodes$x[i])^2 + (nodes$y[j] - nodes$y[i])^2)
#       attachment_probs[j] <- (degrees[j]^alpha) / (1 + beta*d^r)
#     }
#     attachment_probs <- attachment_probs / sum(attachment_probs)
#     
#     # Choose m neighbors with preferential attachment
#     neighbors <- sample(1:(i-1), m, replace = TRUE, prob = attachment_probs)
#     
#     # Connect node to its neighbors
#     adj_mat[i, neighbors] <- 1
#     adj_mat[neighbors, i] <- 1
#     
#     # Update node degrees
#     degrees[i] <- m
#     degrees[neighbors] <- degrees[neighbors] + 1
#     
#     # Rewire edges with probability prewire
#     for (j in 1:m) {
#       if (runif(1) < prewire) {
#         k <- sample(1:N, 1)
#         d <- sqrt((nodes$x[k] - nodes$x[i])^2 + (nodes$y[k] - nodes$y[i])^2)
#         if (k != i && adj_mat[i, k] == 0 && d < r) {
#           adj_mat[i, neighbors[j]] <- 0
#           adj_mat[neighbors[j], i] <- 0
#           adj_mat[i, k] <- 1
#           adj_mat[k, i] <- 1
#           neighbors[j] <- k
#         }
#       }
#     }
#   }
#   
#   # Return adjacency matrix and node coordinates
#   return(list(adj_mat = adj_mat, nodes = nodes))
# }

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# library(Matrix)
# library(Rcpp)
# cppFunction("
# NumericMatrix StochasticBlockModel(NumericVector s, NumericMatrix P) {
#   int n = s.size();
#   NumericMatrix adj(n, n);
#   for(int i=0; i<n; i++) {
#     for(int j=0; j<n; j++) {
#       double u = R::runif(0.0, 1.0);
#       if(u < P(i,j)) {
#         adj(i,j) = 1;
#         adj(j,i) = 1;
#       }
#     }
#   }
#   return adj;
# }")
# 
# generate_scale_free_expander_graph <- function(N, beta, alpha, prewire, r, m) {
#   # Generate initial coordinates
#   coords <- matrix(runif(N*2), ncol=2)
#   
#   # Generate adjacency matrix
#   adj <- Matrix(0, nrow=N, ncol=N, sparse=TRUE)
#   
#   # Add first node
#   adj[1,1] <- 1
#   
#   # Stochastic block matrix
#   s <- rep(1, N)
#   sbm_p <- matrix(1, nrow=N, ncol=N)
#   sbm_p[lower.tri(sbm_p)] <- runif(N*(N-1)/2)
#   sbm_p <- sbm_p + t(sbm_p)
#   
#   # Degree distribution
#   deg <- rep(1, N)
#   
#   # Add new nodes one by one
#   for (i in 2:N) {
#     # Categorize node into a community
#     comm <- sample(1:N, 1, prob=s)
#     
#     # Compute attachment probabilities
#     dist <- sqrt(rowSums((coords - coords[comm,])^2))
#     prob <- (deg^alpha / sum(deg^alpha)) * (1 + beta*deg[comm]^2 / (1 + beta*sum(deg^2)))
#     prob[dist > r] <- 0
#     prob[comm] <- 0
#     prob <- prob / sum(prob)
#     
#     # Choose m neighbors with preferential attachment
#     if (sum(prob) > 0) {
#       neighbor_probs <- rep(0, N)
#       neighbor_probs[sample(1:N, m, prob=prob, replace=FALSE)] <- 1
#       
#       # Add edges to the adjacency matrix
#       adj[i,neighbor_probs > 0] <- 1
#       adj[neighbor_probs > 0,i] <- 1
#       
#       # Update degree distribution
#       deg[i] <- sum(adj[i,])
#       deg[neighbor_probs > 0] <- deg[neighbor_probs > 0] + 1
#     }
#     
#     # Rewire edges with probability prewire
#     for (j in 1:N) {
#       if (j != i & adj[i,j] > 0) {
#         if (runif(1) < prewire) {
#           dist <- sqrt((coords[i,1] - coords[j,1])^2 + (coords[i,2] - coords[j,2])^2)
#           if (dist <= r) {
#             choices <- which(dist <= r & adj[j,] == 0)
#             if (length(choices) > 0) {
#               rewired <- sample(choices, 1)
#               adj[i,j] <- 0
#               adj[j,i] <- 0
#               adj[i,rewired] <- 1
#               adj[rewired,i] <- 1
#             }
#           }
#         }
#       }
#     }
#     
#     #
#     
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#     library(igraph)
#     
#     generate_spatial_scalefree <- function(N, beta, alpha, r, m, prewire) {
#       # Create initial node and block
#       nodes <- data.frame(x = runif(1), y = runif(1), block = 1)
#       edges <- data.frame(from = integer(), to = integer())
#       
#       # Initialize block/community matrix
#       block_probs <- matrix(rep(alpha/m, m^2), ncol = m)
#       diag(block_probs) <- alpha/(m-1)
#       
#       # Create graph iteratively
#       for (i in 2:N) {
#         # Step 1: Create arbitrary number of blocks
#         num_blocks <- sample(1:m, 1, prob = alpha/(1:m))
#         blocks <- sample(1:num_blocks, i, replace = TRUE)
#         
#         # Step 2: Define block probability attachment matrix
#         block_probs <- matrix(rep(alpha/m, m^2), ncol = m)
#         diag(block_probs) <- alpha/(m-1)
#         dist_matrix <- as.matrix(dist(nodes[,1:2]))
#         for (b in 1:num_blocks) {
#           in_block <- blocks == b
#           block_size <- sum(in_block)
#           block_dist <- dist_matrix[in_block, in_block]
#           block_probs[b,] <- (1-alpha)*rowSums(in_block)/sum(in_block) + 
#             alpha*rowSums(block_dist <= r)/block_size^(beta-1)
#         }
#         
# #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# #+sample code
# #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#         grow_spatial_sf_graph <- function(beta, N, prewire, m, r) {
#           # create initial node
#           nodes <- data.frame(x = runif(1), y = runif(1), community = 1)
#           edges <- data.frame(from = integer(), to = integer())
#           
#           # initialize community/block probability attachment matrix
#           prob_mat <- matrix(0, ncol = 1, nrow = 1)
#           
#           for (i in 2:N) {
#             # create arbitrary number of communities/blocks
#             num_communities <- sample(2:5, 1)
#             prob_mat <- matrix(runif(num_communities^2), ncol = num_communities, nrow = num_communities)
#             
#             # generate community/block assignments for new node
#             new_community <- rpois(1, num_communities)
#             new_node <- data.frame(x = runif(1), y = runif(1), community = new_community)
#             
#             # compute distances between new node and existing nodes
#             distances <- sqrt((nodes$x - new_node$x)^2 + (nodes$y - new_node$y)^2)
#             
#             # compute probability of edge based on distance and community/block assignment
#             probabilities <- rep(0, nrow(nodes))
#             for (j in 1:nrow(nodes)) {
#               if (nodes$community[j] == new_community) {
#                 probabilities[j] <- prob_mat[new_community, new_community]
#               } else {
#                 probabilities[j] <- prob_mat[nodes$community[j], new_community]
#               }
#               probabilities[j] <- probabilities[j] * (distances[j] <= r)
#             }
#             
#             # apply scale-free degree preferential attachment
#             if (i <= m + 1) {
#               to_connect <- sample(1:(i-1), 1)
#             } else {
#               degrees <- table(edges$from)
#               degree_probs <- (degrees+1)^beta / sum((degrees+1)^beta)
#               to_connect <- sample(names(degree_probs), size = m, replace = TRUE, prob = degree_probs)
#             }
#             
#             # add edges to existing nodes based on probability and degree preferential attachment
#             for (j in 1:length(to_connect)) {
#               if (runif(1) < probabilities[to_connect[j]]) {
#                 edges <- rbind(edges, data.frame(from = i, to = to_connect[j]))
#               }
#             }
#             
#             # randomly rewire edges for small world effect
#             for (j in 1:nrow(edges)) {
#               if (runif(1) < prewire) {
#                 distances <- sqrt((nodes$x - new_node$x[edges$to[j]])^2 + (nodes$y - new_node$y[edges$to[j]])^2)
#                 possible_rewires <- which(distances <= r & nodes$community == new_node$community[edges$to[j]])
#                 if (length(possible_rewires) > 1) {
#                   new_to <- sample(possible_rewires[possible_rewires != edges$from[j]], 1)
#                   edges$to[j] <- new_to
#                 }
#               }
#             }
#             
#             # add new node to node list
#             nodes <- rbind(nodes, new_node)
#           }
#           
#           # create igraph object from edge list
#           g <- igraph::graph_from_data_frame(edges, directed = FALSE)
#           layout <- igraph::layout_with_fr(as.matrix(nodes[, c("x", "y")]), grid = FALSE)
#           igraph
#           
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          # generate_spatial_scale_free_graph <- function(N, r, beta, m, prewire) {
          #   # Initialize the graph with a single node in a community/block
          #   nodes <- data.frame(x = runif(1), y = runif(1), community = 1)
          #   edges <- data.frame(from = numeric(), to = numeric())
          #   
          #   # Create the remaining N-1 nodes
          #   for (i in 2:N) {
          #     # Step 1: Create an arbitrary number of communities/blocks
          #     num_communities <- sample(1:5, 1)
          #     communities <- sample(1:num_communities, i, replace = TRUE)
          #     
          #     # Step 2: Define the community/block probability attachment matrix
          #     # First compute the spatial distances
          #     distances <- as.matrix(dist(nodes[, c("x", "y")]))
          #     distances[upper.tri(distances)] <- Inf
          #     distances <- distances < r
          #     # Then compute the degree distribution
          #     degrees <- degree(graph_from_data_frame(edges, directed = FALSE), loops = FALSE)
          #     degrees <- degrees[order(-degrees)]
          #     # Finally compute the probability attachment matrix
          #     attachment_matrix <- matrix(0, num_communities, num_communities)
          #     for (j in 1:num_communities) {
          #       for (k in 1:num_communities) {
          #         if (j == k) {
          #           attachment_matrix[j, k] <- sum(degrees[communities == j]^(-beta))
          #         } else {
          #           attachment_matrix[j, k] <- sum(degrees[communities == j]^(-beta)) *
          #             sum(degrees[communities == k]^(-beta)) / sum(degrees^(-beta))
          #         }
          #       }
          #     }
          #     # Normalize the attachment matrix
          #     attachment_matrix <- attachment_matrix / sum(attachment_matrix)
          #     
          #     # Step 3: Generate community/block assignments for the new node
          #     community_probabilities <- rep(1/num_communities, num_communities)
          #     community_assignment <- sample(1:num_communities, 1, prob = community_probabilities)
          #     
          #     # Step 4: (spatial distance effect interplay) Add edges between the new node and existing nodes based on community/block assignments
          #     node_distances <- distances[communities == community_assignment, ]
          #     node_probabilities <- rowSums(node_distances * degrees^(-beta))
          #     node_probabilities <- node_probabilities / sum(node_probabilities)
          #     new_node_index <- nrow(nodes) + 1
          #     new_edges <- sample(which(node_distances), m, replace = TRUE, prob = node_probabilities)
          #     new_edges <- cbind(rep(new_node_index, m), which(node_distances)[new_edges])
          #     edges <- rbind(edges, new_edges)
          #     
          #     # Small world effect: Rewire edges with probability prewire
          #     for (j in 1:nrow(edges)) {
          #       if (runif(1) < prewire) {
          #         from <- edges[j, 1]
          #         to <- edges[j, 2]
          #         available_nodes <- setdiff(1:N, c(from, to))
          #         node_distances <- sqrt((nodes$x[available_nodes] - nodes$x[from])^2 +
          #                                  (nodes$y[available_nodes] - nodes$y[from])^2)
          #         rewire_probabilities <- node_distances < r
          #         rewire_probabilities[to] <- FALSE
          #         if (any(rewire_probabilities)) {
          #           to <- sample(available_nodes[rewire_probabilities], 1)
          #           edges[j, 2] <- to
          #         }
          #       }
          #     }
          #     
          #     # Add the new node
          #     
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
generateSpatialScaleFreeGraph <- function(N=50, beta=2, numofComm=10, r=0.4, m=2, prewire=0.3) {
                
                # create initial node
                coords <- runif(2, 0, 1)
                adjMat <- matrix(0, 1, 1)
                comm <- 1
                
                # grow graph
                for (i in 2:N) {
                  
                  # create communities
                  commProb <- matrix(runif(numofComm^2), numofComm, numofComm)
                  commProb <- commProb/sum(commProb)
                  
                  # generate community assignment for new node
                  commAssignment <- rpois(1, numofComm)
                  
                  # compute distances to existing nodes
                  distToExisting <- as.matrix(dist(coords))
                  
                  # compute within community attachment probability
                  withinCommProb <- (distToExisting <= r)*(1 + m*degree(adjMat[1:(i-1),1:(i-1)]))^beta
                  withinCommProb[is.na(withinCommProb)] <- 0
                  withinCommProb <- withinCommProb/sum(withinCommProb)
                  
                  # compute between community attachment probability
                  betweenCommProb <- rep(0, i-1)
                  for (j in 1:(i-1)) {
                    if (commAssignment == comm[j]) {
                      betweenCommProb[j] <- withinCommProb[j]
                    } else {
                      betweenCommProb[j] <- (distToExisting[j] <= r)*(1 + m*degree(adjMat[1:(i-1),j]))^beta
                    }
                  }
                  betweenCommProb[is.na(betweenCommProb)] <- 0
                  betweenCommProb <- betweenCommProb/sum(betweenCommProb)
                  
                  # add edges to new node
                  newEdges <- c()
                  for (j in 1:m) {
                    if (runif(1) < prewire) {
                      # rewire to random node within certain distance
                      newDist <- runif(1, 0, r)
                      eligibleNodes <- which(distToExisting <= newDist)
                      eligibleNodes <- eligibleNodes[eligibleNodes != (i-1)]
                      if (length(eligibleNodes) > 0) {
                        newEdge <- sample(eligibleNodes, 1)
                        newEdges <- c(newEdges, newEdge)
                      }
                    } else {
                      # attach to within community node or between community node
                      if (runif(1) < commProb[commAssignment, comm[i-1]]) {
                        newEdge <- sample(1:(i-1), 1, prob = withinCommProb)
                        newEdges <- c(newEdges, newEdge)
                      } else {
                        newEdge <- sample(1:(i-1), 1, prob = betweenCommProb)
                        newEdges <- c(newEdges, newEdge)
                      }
                    }
                  }
                  
                  # add new node to graph
                  adjMat <- rbind(adjMat, c(rep(0, i-1), newEdges))
                  adjMat[newEdges, i] <- 1
                  coords <- rbind(coords, runif(2, 0, 1))
                  comm <- c(comm, commAssignment)
                }
                
                # return adjacency matrix and coordinates
                return(list(adjMat = adjMat))
                }
              generateSpatialScaleFreeGraph(N=50, beta=2, numofComm=10, r=0.4, m=2, prewire=0.3)                            
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(Matrix)
library(stats)

generate_spatial_scalefree_graph <- function(N, beta, numofComm, r, m, prewire) {
  # Step 1: Generate random coordinates for nodes
  coordinates <- matrix(runif(N * 2), ncol = 2)
  
  # Step 2: Initialize the network with a single node and place it in an initial community / block
  adj_matrix <- spMatrix(1, 1, 0)
  community_assignment <- c(1)
  
  # Step 3: For each new node now added to the network, do the following:
  for (i in 2:N) {
    # Step 3(i): Create an arbitrary number of communities/blocks to partition or divide the network into using the parameter "numofComm"
    comm_prob <- matrix(runif(numofComm^2), nrow = numofComm)
    
    # Step 3(ii): Define the community probability attachment matrix
    comm_prob[lower.tri(comm_prob)] <- 0
    comm_prob <- comm_prob + t(comm_prob) - diag(diag(comm_prob))
    
    # Step 3(iii): Generate within-community attachment probability matrix
    dist_mat <- as.matrix(dist(coordinates[1:i-1, ]))
    within_prob <- (dist_mat <= r) * (1 + dist_mat / r) ^ (-beta)
    within_prob <- within_prob * (1 - diag(numofComm))
    within_prob <- within_prob / rowSums(within_prob)
    
    # Step 3(iv): Generate between-community attachment probability matrix
    between_prob <- (dist_mat <= r) * (1 + dist_mat / r) ^ (-beta)
    between_prob <- between_prob * comm_prob
    between_prob <- between_prob / rowSums(between_prob)
    
    # Step 3(v): Generate community assignment for the new node
    comm_assignment_prob <- rep(1 / numofComm, numofComm)
    comm_assignment <- sample.int(numofComm, 1, prob = comm_assignment_prob)
    
    # Step 3(vi): Add edges between the new node and existing nodes based on the community/block assignments and community/block probability matrix for within and between communities
    new_row <- rep(0, i - 1)
    comm_assignment_indices <- which(community_assignment == comm_assignment)
    for (j in 1:m) {
      # within-community attachment
      within_indices <- sample(comm_assignment_indices, 2, replace = TRUE, prob = within_prob[comm_assignment, ])
      new_row[within_indices] <- new_row[within_indices] + 1
      
      # between-community attachment
      between_indices <- sample(setdiff(1:(i-1), comm_assignment_indices), 2, replace = TRUE, prob = between_prob[comm_assignment, ])
      new_row[between_indices] <- new_row[between_indices] + 1
    }
    
    # Add rewiring probability
    for (j in 1:(i-1)) {
      if (runif(1) < prewire && new_row[j] == 0) {
        new_indices <- sample(setdiff(1:(i-1), c(j)), 1)
        new_row[c(j, new_indices)] <- c(1, 1)
      }
    }
    
    # Add the new row and column to the adjacency matrix
    adj_matrix <- bdiag(adj_matrix, spMatrix(rep(0, i-1), new_row))
    community_assignment <- c(community_assignment
                              
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to generate a spatial scale-free graph
# spatial_scale_free_graph <- function(N, beta, r, numofComm, m, prewire) {
#   
#   # Initialize the network with a single node
#   nodes <- data.frame(x = runif(1), y = runif(1), community = 1)
#   edges <- data.frame(from = numeric(0), to = numeric(0))
#   
#   # Generate nodes and edges
#   for (i in 2:N) {
#     # Step 1: Create an arbitrary number of communities/blocks
#     communities <- sample.int(numofComm, i-1, replace = TRUE)
#     
#     # Step 2: Define community probability attachment matrix
#     P_within <- matrix(0, nrow = numofComm, ncol = numofComm)
#     P_between <- matrix(0, nrow = numofComm, ncol = numofComm)
#     for (c1 in 1:numofComm) {
#       for (c2 in 1:numofComm) {
#         if (c1 == c2) {
#           P_within[c1,c2] <- 1
#         } else {
#           P_within[c1,c2] <- 1/(sum(communities == c1)*sum(communities == c2))
#           P_between[c1,c2] <- prewire*r*exp(-beta*r^2)
#         }
#       }
#     }
#     P_between <- P_between/sum(P_between)
#     
#     # Step 3: Generate degree sequence and within-community edges
#     degree_sequence <- numeric(i-1)
#     for (j in 1:(i-1)) {
#       dist <- sqrt((nodes$x[j] - nodes$x[i])^2 + (nodes$y[j] - nodes$y[i])^2)
#       if (dist < r) {
#         degree_sequence[j] <- sum(edges$from == j) + sum(edges$to == j)
#       }
#     }
#     degree_sequence <- degree_sequence^beta
#     degree_sequence <- degree_sequence/sum(degree_sequence)
#     for (j in 1:m) {
#       k <- sample.int(i-1, size = 1, prob = degree_sequence)
#       while (k %in% edges$from[edges$to == i] | k %in% edges$to[edges$from == i]) {
#         k <- sample.int(i-1, size = 1, prob = degree_sequence)
#       }
#       edges <- rbind(edges, data.frame(from = i, to = k))
#     }
#     
#     # Step 4: Generate between-community edges
#     for (j in 1:m) {
#       c1 <- sample.int(numofComm, size = 1, prob = P_between[,communities[i-1]])
#       c2 <- sample.int(numofComm, size = 1, prob = P_within[c1,])
#       k <- sample(which(communities == c2), size = 1)
#       while ((i %in% edges$from & k == edges$to[i]) | (i %in% edges$to & k == edges$from[i])) {
#         k <- sample(which(communities == c2), size = 1)
#       }
#       edges <- rbind(edges, data.frame(from = i, to = k))
#     }
#     
#     # Step 5: Assign node to a community
#     lambd <- matrix(0, nrow = 10, ncol = 10)
#     for (j in 1:(i-1)) {
#       lambd[ceiling(nodes$x
                    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(Matrix)
library(stats)
library(spatstat)

spatial_scale_free_graph <- function(N, beta, numofComm, prewire, m, r) {
  # Generate initial node randomly distributed in the unit square
  nodes <- matrix(runif(N*2), ncol=2)
  
  # Initialize the community assignments for each node
  community_assignments <- matrix(1, ncol=numofComm, nrow=N)
  
  # Initialize the adjacency matrix
  adj_matrix <- Matrix(0, ncol=N, nrow=N)
  
  # Create initial node
  adj_matrix[1,1] <- 1
  community_assignments=NULL
  # Create communities
  for (i in 2:N) {
    # Create numofComm communities
    communities <- sample(numofComm, 1)
    community_assignments[i,] <- sample(communities, N, replace=TRUE)}
    
    # Calculate within community attachment probability matrix
    prob_within_comm <- matrix(0, ncol=numofComm, nrow=numofComm)
    for (c1 in 1:numofComm) {
      for (c2 in 1:numofComm) {
        if (c1 == c2) {
          prob_within_comm[c1,c2] <- 1
        } else {
          # Calculate the number of nodes in each community
          nodes_in_c1 <- sum(community_assignments[,c1] == c1)
          nodes_in_c2 <- sum(community_assignments[,c2] == c2)
          # Calculate the distance matrix between nodes in each community
          dist_within_c1 <- as.matrix(dist(nodes[community_assignments[,c1] == c1,]))
          dist_within_c2 <- as.matrix(dist(nodes[community_assignments[,c2] == c2,]))
          # Calculate the probability of attachment
          prob_within_comm[c1,c2] <- sum((1 / (1 + dist_within_c1^r))^(beta) * nodes_in_c1 * nodes_in_c2) / (nodes_in_c1 * nodes_in_c2)
        }
      }
    }
    
    # Calculate between community attachment probability matrix
    prob_between_comm <- matrix(0, ncol=numofComm, nrow=numofComm)
    for (c1 in 1:numofComm) {
      for (c2 in 1:numofComm) {
        if (c1 == c2) {
          prob_between_comm[c1,c2] <- 0
        } else {
          # Calculate the number of nodes in each community
          nodes_in_c1 <- sum(community_assignments[,c1] == c1)
          nodes_in_c2 <- sum(community_assignments[,c2] == c2)
          # Calculate the distance matrix between nodes in each community
          dist_between_c1c2 <- as.matrix(dist(nodes[community_assignments[,c1] == c1,], nodes[community_assignments[,c2] == c2,]))
          # Calculate the probability of attachment
          prob_between_comm[c1,c2] <- sum((1 / (1 + dist_between_c1c2^r))^(beta) * nodes_in_c1 * nodes_in_c2) / (nodes_in_c1 * nodes_in_c2)
        }
      }
    }
    
    # Add edges for new node based on community assignments and probability matrices
    for (j in 1:(i-1)) {
      if (community_assignments[i,j] == community_assignments[i,i]) {
        if (runif(1) < prob_within_comm[
          
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          # Function to generate a spatial scale-free graph
          # Inputs:
          #   beta: power law parameter for degree distribution
          #   N: number of nodes to generate
          #   numofComm: number of communities/blocks to partition the network into
          #   prewire: rewiring probability for small world effect
          #   m: number of edges added at each instance
          #   r: parameter that controls the degree of sparsity
          # Output:
          #   A spatial scale-free graph object
          spatialScaleFreeGraph <- function(beta, N, numofComm, prewire, m, r) {
            # Create initial node and assign it to a block/community
            blocks <- sample(1:numofComm, 1)
            x <- runif(1)
            y <- runif(1)
            nodes <- data.frame(id = 1, block = blocks, x = x, y = y, degree = 0)
            
            # Generate remaining nodes
            for (i in 2:N) {
              # Step 1: Create communities
              blocks <- sample(1:numofComm, i, replace = TRUE)
              
              # Step 2: Define probability attachment matrix
              P <- matrix(0, nrow = numofComm, ncol = numofComm)
              for (j in 1:numofComm) {
                for (k in 1:numofComm) {
                  if (j == k) {
                    P[j, k] <- 1 - prewire
                  } else {
                    P[j, k] <- prewire / (numofComm - 1)
                  }
                }
              }
              
              # Step 3: Generate within community attachment
              withinBlocks <- blocks
              for (j in 1:numofComm) {
                # Select nodes within the same block
                blockNodes <- nodes[nodes$block == j,]
                if (nrow(blockNodes) > 0) {
                  # Calculate distances
                  distances <- sqrt((blockNodes$x - x[i-1])^2 + (blockNodes$y - y[i-1])^2)
                  
                  # Calculate attachment probabilities based on distance and degree
                  probs <- (distances <= r) * (blockNodes$degree + 1)^beta
                  probs <- probs / sum(probs)
                  
                  # Add edges to selected nodes
                  selected <- sample(blockNodes$id, m, replace = TRUE, prob = probs)
                  edges <- data.frame(from = i, to = selected)
                  nodes$degree[nodes$id %in% selected] <- nodes$degree[nodes$id %in% selected] + 1
                  nodes$degree[i] <- nodes$degree[i] + length(selected)
                }
              }
              
              # Step 4: Generate between community attachment
              for (j in 1:numofComm) {
                # Select nodes in other blocks
                blockNodes <- nodes[nodes$block != j,]
                if (nrow(blockNodes) > 0) {
                  # Calculate distances
                  distances <- sqrt((blockNodes$x - x[i-1])^2 + (blockNodes$y - y[i-1])^2)
                  
                  # Calculate attachment probabilities based on distance and degree
                  probs <- (distances <= r) * (blockNodes$degree + 1)^beta
                  probs <- probs / sum(probs)
                  
                  # Add edges to selected nodes
                  selected <- sample(blockNodes$id, m, replace = TRUE, prob = probs)
                  edges <- rbind(edges, data.frame(from = i, to
                                                   
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# grow_spatial_scale_free_graph <- function(N, numofComm, beta, m, r, prewire) {
#   # Initialize the graph with a single node and place it in a random position in the unit square
#   nodes <- data.frame(x = runif(1), y = runif(1), community = 1)
#   edges <- data.frame(from = integer(), to = integer())
#   
#   # Define the within-community attachment probability matrix
#   P_within <- matrix(rep(1 / numofComm, numofComm^2), nrow = numofComm)
#   diag(P_within) <- 1
#   
#   # Define the between-community attachment probability matrix
#   P_between <- matrix(rep((1 - P_within) / (numofComm - 1), numofComm^2), nrow = numofComm)
#   diag(P_between) <- 0
#   
#   # Define a function to compute the distance between two nodes
#   dist_nodes <- function(node1, node2) {
#     sqrt((node1$x - node2$x)^2 + (node1$y - node2$y)^2)
#   }
#   
#   # Add nodes iteratively
#   for (i in 2:N) {
#     # Step 1: Create numofComm communities/blocks
#     communities <- sample.int(numofComm, i, replace = TRUE)
#     
#     # Step 2: Define the community probability attachment matrix
#     P <- P_within
#     P[communities != communities[i]] <- P_between[communities != communities[i]]
#     
#     # Step 3: Compute the degree of existing nodes and the probability of attachment
#     degrees <- degree(nodes, mode = "out")
#     probs_within <- degrees^beta
#     probs_within[communities != communities[i]] <- 0
#     probs_within <- probs_within / sum(probs_within)
#     
#     # Step 4: Compute the probability of attachment between communities
#     probs_between <- rep(1, i-1)
#     for (j in 1:(i-1)) {
#       if (communities[j] != communities[i]) {
#         d <- dist_nodes(nodes[i,], nodes[j,])
#         if (d < r) {
#           probs_between[j] <- (d/r)^beta
#         } else {
#           probs_between[j] <- 0
#         }
#       }
#     }
#     probs_between <- probs_between / sum(probs_between)
#     
#     # Step 5: Generate a community/block assignment for the new node
#     community_i <- sample.int(numofComm, 1)
#     
#     # Step 6: Add edges based on the community/block assignment and probability matrix
#     for (j in 1:(i-1)) {
#       if (communities[j] == community_i) {
#         if (runif(1) < probs_within[j]) {
#           edges <- rbind(edges, data.frame(from = j, to = i))
#         }
#       } else {
#         if (runif(1) < probs_between[j] * prewire) {
#           edges <- rbind(edges, data.frame(from = j, to = i))
#         }
#       }
#     }
#     
#     # Add the new node to the network
#     new_node <- data.frame(x = runif(1), y = runif(1), community = community_i)
#     nodes <- rbind(nodes, new_node)
#   }
#   
#   # Add small-world effect by randomly rewiring edges with probability prewire
#   for (i in 1:nrow(edges)) {
#     if (runif(1) < prewire) {
#       dists <- s
      
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # spatial_scale_free_graph <- function(N, beta, numofComm, r, m, prewire) {
      #   # Create initial node
      #   nodes <- data.frame(x = runif(1), y = runif(1))
      #   nodes$community <- 1
      #   
      #   # Create initial adjacency matrix
      #   adjacency <- matrix(0, nrow = 1, ncol = 1)
      #   
      #   # Create probability attachment matrix
      #   P_within <- matrix(0, nrow = numofComm, ncol = numofComm)
      #   P_between <- matrix(0, nrow = numofComm, ncol = numofComm)
      #   
      #   # Create within-community probability attachment matrix
      #   for (i in 1:numofComm) {
      #     for (j in 1:numofComm) {
      #       if (i == j) {
      #         P_within[i, j] <- 1
      #       } else {
      #         P_within[i, j] <- 1 / (1 + r^beta)
      #       }
      #     }
      #   }
      #   
      #   # Create between-community probability attachment matrix
      #   for (i in 1:numofComm) {
      #     for (j in 1:numofComm) {
      #       if (i == j) {
      #         P_between[i, j] <- 0
      #       } else {
      #         P_between[i, j] <- prewire / (1 + r^beta)
      #       }
      #     }
      #   }
      #   community=NULL
      #   # Generate additional nodes
      #   for (i in 2:N) {
      #     # Create communities for node i
      #     community_probs <- rep(1/numofComm, numofComm)
      #     community[i] <- sample(1:numofComm, 1, prob = community_probs)
      #     
      #     # Create within-community edges for node i
      #     within_probs <- rep(1/m, i-1)
      #     within_probs[nodes$community[-i] != nodes$community[i]] <- 0
      #     within_probs[which(pmin(dist(nodes[c("x","y")]), r) == 0)] <- 0
      #     within_probs <- within_probs / sum(within_probs)
      #     within_edges <- sample(1:(i-1), m, replace = TRUE, prob = within_probs)
      #     
      #     # Create between-community edges for node i
      #     between_probs <- rep(1/m, i-1)
      #     between_probs[nodes$community[-i] == nodes$community[i]] <- 0
      #     between_probs[which(pmin(dist(nodes[c("x","y")]), r) == 0)] <- 0
      #     between_probs <- between_probs / sum(between_probs)
      #     between_edges <- sample(1:(i-1), m, replace = TRUE, prob = between_probs)
      #     
      #     # Add edges to adjacency matrix
      #     new_row <- c(rep(0, i-1), within_edges, between_edges)
      #     adjacency <- rbind(adjacency, new_row)
      #     new_col <- c(adjacency[-1, i-1], within_edges, between_edges)
      #     adjacency <- cbind(adjacency, new_col)
      #     
      #     # Add new node to nodes data frame
      #     new_node <- data.frame(x = runif(1), y = runif(1), community = nodes$community[i])
      #     nodes <- rbind(nodes, new_node)
      #     
      #     # Rewire edges with small world probability
      #     for (j in 1:(i-1)) {
      #       if (runif(1) < prewire) {
      #         dist_ij <- sqrt((nodes$x[i] - nodes$x[j])^2 + (nodes$y[i
                                                                     
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to generate a spatial scale-free graph
# with arbitrary number of communities
# Inputs:
# - N: number of nodes
# - numofComm: number of communities
# - beta: power law parameter for degree distribution
# - m: number of edges added at each instance
# - r: parameter that controls the degree of sparsity
# - alpha: parameter that controls the degree preference
# - prewire: rewiring probability for small world effect
# Output: adjacency matrix of the resulting graph

generate_spatial_scalefree_graph <- function(N, numofComm, beta, m, r, alpha, prewire) {
  
  # Initialize adjacency matrix with first node
  adj_mat <- matrix(0, nrow = N, ncol = N)
  adj_mat[1, 1] <- 1
  
  # Define distance matrix between nodes
  coords <- matrix(runif(N*2), ncol = 2)
  dist_mat <- as.matrix(dist(coords))
  
  # Loop over the remaining nodes
  for (i in 2:N) {
    # Step 1: Create arbitrary number of communities
    comm_assign <- sample(numofComm, 1)
    # Step 2: Define community probability attachment matrix
    P_ij <- matrix(0, nrow = numofComm, ncol = numofComm)
    for (c1 in 1:numofComm) {
      for (c2 in 1:numofComm) {
        if (c1 == c2) {
          P_ij[c1, c2] <- 1
        } else {
          P_ij[c1, c2] <- 1/numofComm
        }
      }
    }
    # Step 3: Define within community attachment probability matrix
    spat_distprobs <- 1/(1+(dist_mat[1:(i-1), i])^beta)
    spat_distprobs[dist_mat[1:(i-1), i] > r] <- 0
    degrees <- rowSums(adj_mat[1:(i-1), ])
    deg_probs <- (degrees^alpha)/(sum(degrees^alpha))
    deg_probs[is.na(deg_probs)] <- 0
    within_probs <- P_ij[comm_assign, comm_assign] * spat_distprobs * deg_probs
    within_probs[which(within_probs == 0)] <- NA
    # Step 4: Define between community attachment probability matrix
    between_probs <- matrix(0, nrow = i-1, ncol = numofComm)
    for (j in 1:(i-1)) {
      if (dist_mat[j, i] <= r) {
        between_probs[j, ] <- P_ij[, comm_assign] * (1-spat_distprobs[j]) * deg_probs[j]
      }
    }
    between_probs[which(between_probs == 0)] <- NA
    # Step 5: Generate community/block assignment for new node
    block_assign <- sample(numofComm, 1)
    # Step 6: Add edges based on community/block assignment and probability matrix
    for (j in 1:m) {
      within_probs_norm <- within_probs/sum(within_probs, na.rm = TRUE)
      to_node <- sample(1:(i-1), 1, prob = within_probs_norm)
      to_block <- block_assign
      if (runif(1) <
          
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
generate_spatial_sf_graph <- function(N, beta, numofComm, r, m, prewire, alpha) {
  # Initialize the adjacency matrix
  adj_mat <- matrix(0, nrow = N, ncol = N)
  
  # Initialize the community assignment matrix
  comm_mat <- matrix(0, nrow = N, ncol = 1)
  
  # Create an initial node and assign it to a community
  comm_mat[1] <- sample(1:numofComm, 1)
  
  # Define a vector to keep track of the degrees of nodes
  degrees <- rep(0, N)
  
  # Add new nodes to the network
  for (i in 2:N) {
    # Step 1: Create communities
    comm_prob <- matrix(0, nrow = numofComm, ncol = numofComm)
    diag(comm_prob) <- 1
    
    # Step 2: Define within community attachment probability matrix
    spat_distprobs <- ifelse(dist_mat[1:(i-1), i] <= r, 1/(1+(dist_mat[1:(i-1), i])^beta), 0)
    deg_probs <- (degrees[1:(i-1)]^alpha) / sum(degrees[1:(i-1)]^alpha)
    
    for (k in 1:numofComm) {
      for (l in 1:numofComm) {
        comm_prob[k, l] <- sum(spat_distprobs[comm_mat[1:(i-1)] == k] * deg_probs[comm_mat[1:(i-1)] == l])
      }
    }
    
    # Step 3: Define between community attachment probability matrix
    spat_distprobs <- ifelse(dist_mat[1:(i-1), i] <= r, 1/(1+(dist_mat[1:(i-1), i])^beta), 0)
    deg_probs <- (degrees[1:(i-1)]^alpha) / sum(degrees[1:(i-1)]^alpha)
    
    comm_prob <- comm_prob * (1 - prewire)
    comm_prob <- comm_prob + (prewire / numofComm)
    
    # Step 4: Generate community assignments for the new node
    comm_mat[i] <- sample(1:numofComm, 1, prob = colSums(comm_prob))
    
    # Step 5: Add edges based on community assignments and community probability matrix
    for (j in 1:(i-1)) {
      if (comm_mat[i] == comm_mat[j]) {
        # within community attachment
        attach_prob <- spat_distprobs[j] * deg_probs[j]
        if (runif(1) < attach_prob) {
          adj_mat[i, j] <- 1
          adj_mat[j, i] <- 1
          degrees[i] <- degrees[i] + 1
          degrees[j] <- degrees[j] + 1
        }
      } else {
        # between community attachment
        attach_prob <- spat_distprobs[j] * deg_probs[j] * comm_prob[comm_mat[j], comm_mat[i]]
        if (runif(1) < attach_prob) {
          adj_mat[i, j] <- 1
          adj_mat[j, i] <- 1
          degrees[i] <- degrees[i] + 1
          degrees[j] <- degrees[j] + 1

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
          generate_graph <- function(N, beta, numofComm, prewire, m, r, alpha) {
            # Generate initial coordinates for nodes on a unit square
            coords <- matrix(runif(N*2), ncol = 2)
            
            # Calculate the distance matrix between nodes
            dist_mat <- as.matrix(dist(coords))
            
            # Initialize the adjacency matrix with a single node
            adj_mat <- matrix(0, nrow = N, ncol = N)
            adj_mat[1, 1] <- 1
            
            # Assign the first node to the first community
            comm_assign <- rep(1, N)
            
            # Loop over the remaining nodes
            for (i in 2:N) {
              # Step 1: Create an arbitrary number of communities/blocks
              comm_probs <- rep(1/numofComm, numofComm)
              comm_assign[i] <- sample(numofComm, 1, prob = comm_probs)}
              
              # Step 2: Define the community probability attachment matrix
              P <- matrix(0, nrow = numofComm, ncol = numofComm)
              for (c1 in 1:numofComm) {
                for (c2 in 1:numofComm) {
                  if (c1 == c2) {
                    # Within community attachment probability
                    dist_probs <- ifelse(dist_mat[1:(i-1), i] <= r,
                                         1/(1+(dist_mat[1:(i-1), i])^beta), 0)
                    degrees <- rowSums(adj_mat[1:(i-1), ])
                    deg_probs <- (degrees^alpha)/(sum(degrees^alpha))
                    within_probs <- deg_probs * dist_probs
                    P[c1, c2] <- sum(within_probs)
                  } else {
                    # Between community attachment probability
                    dist_probs <- ifelse(dist_mat[1:(i-1), i] <= r,
                                         1/(1+(dist_mat[1:(i-1), i])^beta), 0)
                    degrees <- rowSums(adj_mat[1:(i-1), ])
                    deg_probs <- (degrees^alpha)/(sum(degrees^alpha))
                    between_probs <- deg_probs * dist_probs * prewire
                    P[c1, c2] <- sum(between_probs)
                  }
                }
              }
              
              # Step 3: Generate a community assignment for the new node
              comm_probs <- P[
                
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
                generate_spatial_scalefree_graph <- function(N, numofComm, beta, m, prewire, r, alpha) {
                  # N: number of nodes
                  # numofComm: number of communities
                  # beta: power law exponent for degree distribution
                  # m: number of edges added for each new node
                  # prewire: rewiring probability for small world effect
                  # r: spatial cutoff distance for within-community edges
                  # alpha: parameter controlling degree distribution preference
                  
                  # initialize adjacency matrix and degree vector
                  adj_mat <- matrix(0, nrow = N, ncol = N)
                  degrees <- rep(0, N)
                  
                  # create initial community for first node
                  comm_assignment <- rep(1, N)
                  num_nodes_per_comm <- rep(1, numofComm)
                  
                  # loop over remaining nodes
                  for (i in 2:N) {
                    # create arbitrary number of communities
                    num_nodes_per_comm <- c(num_nodes_per_comm, 1)
                    comm_assignment[i] <- sample(1:numofComm, 1, replace = TRUE, prob = num_nodes_per_comm)
                    
                    # calculate probability matrix for within- and between-community edges
                    within_prob_mat <- matrix(0, nrow = numofComm, ncol = numofComm)
                    between_prob_mat <- matrix(0, nrow = numofComm, ncol = numofComm)
                    for (comm1 in 1:numofComm) {
                      for (comm2 in 1:numofComm) {
                        if (comm1 == comm2) {
                          # within-community probability
                          comm_nodes <- which(comm_assignment[1:(i-1)] == comm1)
                          if (length(comm_nodes) > 0) {
                            # compute distance probability
                            spat_distprobs <- 1 / (1 + (dist_mat[comm_nodes, i])^beta)
                            spat_distprobs[dist_mat[comm_nodes, i] > r] <- 0
                            # compute degree probability
                            comm_degrees <- degrees[comm_nodes]
                            deg_probs <- (comm_degrees^alpha) / sum(comm_degrees^alpha)
                            # compute overall probability
                            within_prob_mat[comm1, comm2] <- sum(spat_distprobs * deg_probs)
                          }
                        } else {
                          # between-community probability
                          # compute distance probability
                          spat_distprobs <- 1 / (1 + (dist_mat[1:(i-1), i])^beta)
                          spat_distprobs[dist_mat[1:(i-1), i] > r] <- 0
                          # compute degree probability
                          comm_nodes <- which(comm_assignment[1:(i-1)] == comm2)
                          comm_degrees <- degrees[comm_nodes]
                          deg_probs <- (comm_degrees^alpha) / sum(comm_degrees^alpha)
                          # compute overall probability
                          between_prob_mat[comm1, comm2] <- sum(spat_distprobs * deg_probs)
                        }
                      }
                    }
                    
                    # add m edges to existing nodes based on probability matrix
                    for (j in 1:m) {
                      # choose destination community
                      dest_comm <- sample(1:numofComm, 1, prob = within_prob_mat[comm_assignment[i], ])
                      # choose destination node within community
                      dest_nodes <- which(comm_assignment[1:(i-1)] == dest_comm)
                      if (length(dest_nodes) > 0) {
                        dest_node <- sample(dest
                                            
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
generate_spatial_scale_free_graph <- function(N, beta, numofComm, r, m, prewire) {
  # Generate initial node positions
  node_pos <- matrix(runif(N*2), ncol = 2)
  
  # Calculate distance matrix
  dist_mat <- as.matrix(dist(node_pos))
  
  # Initialize adjacency matrix
  adj_mat <- matrix(0, nrow = N, ncol = N)
  
  # Add first node to network
  adj_mat[1, 1] <- 1
  block_assignments <- rep(1, N)
  
  # Loop over remaining nodes
  for (i in 3:N) {
    # Create communities/blocks
    blocks <- sample(1:numofComm, i, replace = TRUE)
    block_probs <- matrix(0, nrow = numofComm, ncol = numofComm)
    
    # Define within community attachment probability matrix
    within_probs <- matrix(0, nrow = i-1, ncol = i-1)
    for (j in 1:(i-1)) {
      for (k in 1:(i-1)) {
        if (j == k) {
          within_probs[j, k] <- 0
        } else if (dist_mat[j, i] <= r && dist_mat[k, i] <= r) {
          within_probs[j, k] <- 1 / (1 + (dist_mat[j, i]^beta))
        } else {
          within_probs[j, k] <- 0
        }
      }
    }
    
    # Calculate degree probabilities for within community attachment
    degrees <- rowSums(adj_mat[1:(i-1), ])
    deg_probs <- (degrees^beta) / sum(degrees^beta)
    
    # Normalize probabilities
    within_probs_norm <- within_probs * deg_probs
    within_probs_norm <- within_probs_norm / sum(within_probs_norm)}
    
    # Assign probabilities to blocks
    for (j in 1:numofComm) {
      for (k in 1:numofComm) {
        if (j == k) {
          block_probs[j, k] <- sum(within_probs_norm[blocks == j, blocks == k])
        } else {
          block_probs[j, k] <- sum(within_probs_norm[blocks == j, blocks != k])
        }
      }
    }
    
    # Define between community attachment probability matrix
    between_probs <- matrix(0, nrow = i-1, ncol = numofComm)
    for (j in 1:(i-1)) {
      for (k in 1:numofComm) {
        if (dist_mat[j, i] <= r) {
          between_probs[j, k] <- 1 / (1 + (dist_mat[j, i]^beta))
        } else {
          between_probs[j, k] <- 0
        }
      }
    }
    
    # Calculate degree probabilities for between community attachment
    degrees <- rowSums(adj_mat[1:(i-1), ])
    deg_probs <- (degrees^beta) / sum(degrees^beta)
    
    # Normalize probabilities
    between_probs_norm <- between_probs * deg_probs
    between_probs_norm <- between_probs_norm / sum(between_probs_norm)
    
    # Generate block assignments for new node
    block_assignments[i] <- sample(1:numofComm, 1, prob = colSums(block_probs))
    
    # Add edges to existing nodes based on block assignments and probabilities
    for (j in 1:(i-1)) {
      
      
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
      # Function to generate a spatial scale-free graph
      # Arguments:
      #   N: number of nodes
      #   numofComm: number of communities
      #   beta: power law parameter for degree distribution
      #   m: number of edges added at each instance
      #   r: distance threshold for within-community attachment
      #   alpha: exponent for degree attachment probability
      #   prewire: rewiring probability for small world effect
      
      spatial_scale_free_graph <- function(N, numofComm, beta, m, r, alpha, prewire) {
        
        # Initialize adjacency matrix with a single node
        adj_mat <- matrix(0, nrow = N, ncol = N)
        adj_mat[1, 1] <- 1
        
        # Initialize block assignment for first node
        block_assignment <- rep(1, N)
        
        # Initialize degree sequence with single node
        degrees <- rep(0, N)
        
        # Initialize distance matrix between nodes
        coords <- matrix(runif(N*2), ncol = 2)
        dist_mat <- as.matrix(dist(coords))
        
        # Loop over remaining nodes to add
        for (i in 2:N) {
          
          # Step 1: Create numofComm communities/blocks
          block_probs <- rep(1/numofComm, numofComm)
          block_assignment[i] <- sample(1:numofComm, 1, prob = block_probs)
          
          # Step 2: Define community probability matrix P_ij
          P_ij <- matrix(0, nrow = numofComm, ncol = numofComm)
          for (b1 in 1:numofComm) {
            for (b2 in 1:numofComm) {
              if (b1 == b2) {
                # within-community attachment
                spat_distprobs <- 1/(1 + (dist_mat[1:(i-1), i]/r)^beta)
                deg_probs <- (degrees[1:(i-1)]^alpha) / sum(degrees[1:(i-1)]^alpha)
                P_ij[b1, b2] <- sum(spat_distprobs * deg_probs)
              } else {
                # between-community attachment
                spat_distprobs <- 1/(1 + (dist_mat[1:(i-1), i]/r)^beta)
                deg_probs <- (degrees[1:(i-1)]^alpha) / sum(degrees[1:(i-1)]^alpha)
                rewiring_probs <- as.numeric(dist_mat[1:(i-1), i] <= r)
                rewiring_probs[which(block_assignment[1:(i-1)] == b1)] <- 0
                P_ij[b1, b2] <- sum(spat_distprobs * deg_probs * rewiring_probs)
              }
            }
          }
          
          # Step 3: Assign node to a community
          block_assignment[i] <- sample(1:numofComm, 1, prob = block_probs)
          
          # Step 4: Add edges to existing nodes
          for (j in 1:(i-1)) {
            if (runif(1) < P_ij[block_assignment[i], block_assignment[j]]) {
              adj_mat[i, j] <- 1
              adj_mat[j, i] <- 1
              degrees[i] <- degrees[i] + 1
              degrees[j] <- degrees[j] + 1
            }
          }
          
          # Step 5: Add rewiring for small world effect
          for (j in 1:N
               
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
grow_spatial_scalefree <- function(N, beta, numofComm, r, m, prewire) {
  # Create distance matrix
  xy <- matrix(runif(N*2), ncol=2)
  dist_mat <- as.matrix(dist(xy))
  
  # Initialize adjacency matrix and community assignments
  adj_mat <- matrix(0, nrow = N, ncol = N)
  comm_assign <- rep(1, N)
  
  # Add first node to the network
  adj_mat[1, 1] <- 1
  
  for (i in 2:N) {
    # Step 1: Create communities/blocks
    comm_probs <- rep(1/numofComm, numofComm)
    comm_assign[i] <- sample(1:numofComm, 1, prob=comm_probs)
    
    # Step 2: Define community probability attachment matrix
    P_ij <- matrix(0, nrow = numofComm, ncol = numofComm)
    for (c1 in 1:numofComm) {
      for (c2 in 1:numofComm) {
        if (c1 == c2) {
          P_ij[c1, c2] <- 1 / (1 + (dist(comm_assign == c1, method="euclidean")[1:(i-1)] / r)^beta)
        } else {
          P_ij[c1, c2] <- 1
        }
      }
    }
    
    # Step 3: Compute degree probabilities for within-community attachment
    degrees <- rowSums(adj_mat[1:(i-1), ])
    deg_probs <- (degrees^beta) / sum(degrees^beta)
    
    # Add edges within community
    for (j in 1:(i-1)) {
      if (comm_assign[j] == comm_assign[i] && dist_mat[j, i] <= r) {
        attach_prob <- deg_probs[j] * P_ij[comm_assign[j], comm_assign[i]]
        if (runif(1) < attach_prob) {
          adj_mat[i, j] <- 1
          adj_mat[j, i] <- 1
        }
      }
    }
    
    # Step 4: Compute degree probabilities for between-community attachment
    between_comm_probs <- rep(0, numofComm)
    for (c in 1:numofComm) {
      degree_sum <- sum(degrees[comm_assign == c])
      between_comm_probs[c] <- (degree_sum^beta) / sum(degrees^beta)
    }
    
    # Add edges between communities
    for (j in 1:(i-1)) {
      if (comm_assign[j] != comm_assign[i] && dist_mat[j, i] <= r) {
        attach_prob <- between_comm_probs[comm_assign[j]] * P_ij[comm_assign[j], comm_assign[i]]
        if (runif(1) < attach_prob) {
          adj_mat[i, j] <- 1
          adj_mat[j, i] <- 1
        }
      }
    }
    
    # Step 5: Generate community assignments for new node
    new_comm_probs <- between_comm_probs * P_ij[, comm_assign[i]]
    new_comm_probs[comm_assign[i]] <- 0
    new_comm_probs <- new_comm_probs / sum(new_comm_probs)
    comm_assign[i] <- sample(1:numofComm, 1, prob=new_comm_probs)
    
    # Step 6: Add edges between new node and existing nodes based on community assignments and probability matrix
    for
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
    generate_spatial_sf_graph <- function(N, beta, numofComm, r, m, prewire, alpha) {
      
      # Step 1: Generate nodes on a unit square with spatial structure
      xy <- runif(N*2, 0, 1)
      x <- xy[1:N]
      y <- xy[(N+1):(2*N)]
      dist_mat <- as.matrix(dist(cbind(x, y)))
      
      # Initialize adjacency matrix and assign first node to a community/block
      adj_mat <- matrix(0, nrow = N, ncol = N)
      block_assignments <- rep(1, N)
      
      # Step 2: Create arbitrary number of communities/blocks
      community_prob_mat <- matrix(1/numofComm, nrow = numofComm, ncol = numofComm)
      
      # Add edges for each new node added to the network
      for (i in 2:N) {
        # Step 3: Calculate within-community probabilities
        degrees <- rowSums(adj_mat[1:(i-1), ])
        spat_distprobs <- ifelse(dist_mat[1:(i-1), i] <= r, 1/(1+(dist_mat[1:(i-1), i])^beta), 0)
        deg_probs <- (degrees^alpha) / sum(degrees^alpha)
        within_probs <- spat_distprobs * deg_probs
        
        # Step 4: Calculate between-community probabilities
        degrees <- rowSums(adj_mat[1:(i-1), ])
        spat_distprobs <- ifelse(dist_mat[1:(i-1), i] <= r, 1/(1+(dist_mat[1:(i-1), i])^beta), 0)
        deg_probs <- (degrees^alpha) / sum(degrees^alpha)
        between_probs <- rep(0, numofComm)
        for (j in 1:numofComm) {
          nodes_in_block <- which(block_assignments[1:(i-1)] == j)
          if (length(nodes_in_block) > 0) {
            within_block_degrees <- rowSums(adj_mat[nodes_in_block, nodes_in_block])
            within_block_probs <- (within_block_degrees^alpha) / sum(within_block_degrees^alpha)
            between_probs[j] <- sum(within_block_probs * within_probs[nodes_in_block])
          }
        }
        between_probs <- between_probs / sum(between_probs)
        
        # Step 5: Generate community assignments for the new node
        block_probs <- community_prob_mat[block_assignments[1:(i-1)], ]
        block_probs[cbind(1:(i-1), block_assignments[1:(i-1)])] <- 0
        block_probs <- apply(block_probs, 2, sum) / sum(block_probs)
        block_assignments[i] <- sample(1:numofComm, size = 1, prob = block_probs)
        
        # Step 6: Add edges based on community assignments and probabilities
        for (j in 1:(i-1)) {
          if (block_assignments[j] == block_assignments[i]) {
            if (runif(1) < within_probs[j]) {
              adj_mat[i,j] <- 1
              adj_mat[j,i] <- 1
            }
          } else {
            if (runif(1) < between_probs[block_assignments[j]]) {
              adj_mat[i,j] <- 1
              adj_mat[j,i] <- 1
            }
          }
          
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            

  
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
          generate_spatial_sf_graph <- function(N, numofComm, beta, m, r, alpha, prewire) {
            # initialize adjacency matrix and degree vector
            adj_mat <- matrix(0, nrow = N, ncol = N)
            degrees <- rep(0, N)
            
            # generate coordinates for nodes on a unit square
            coords <- matrix(runif(N*2), ncol = 2)
            
            # generate distance matrix between nodes
            dist_mat <- as.matrix(dist(coords))
            
            # assign first node to a random community
            comm_assignments <- sample(1:numofComm, 1)
            
            # add nodes to the network
            for (i in 2:N) {
              # generate community assignments for new node
              spat_probs <- rep(1, numofComm)
              spat_probs[comm_assignments] <- 1 + (r^beta)/(1+(dist_mat[1:(i-1), i]/r)^beta)
              comm_probs <- spat_probs * (degrees^alpha)/(sum(degrees^alpha))
              comm_probs <- comm_probs/sum(comm_probs)
              comm_assignments[i] <- sample(1:numofComm, 1, prob = comm_probs)
              
              # generate edge probabilities based on community assignments
              within_comm_probs <- rep(0, i-1)
              within_comm_probs[comm_assignments[1:(i-1)] == comm_assignments[i]] <-
                1 + (r^beta)/(1+(dist_mat[1:(i-1), i]/r)^beta)
              within_comm_probs <- within_comm_probs * (degrees[1:(i-1)]^alpha)
              within_comm_probs <- within_comm_probs/sum(within_comm_probs)
              between_comm_probs <- rep(0, i-1)
              between_comm_probs[comm_assignments[1:(i-1)] != comm_assignments[i]] <-
                1 + (r^beta)/(1+(dist_mat[1:(i-1), i]/r)^beta)
              between_comm_probs <- between_comm_probs * (degrees[1:(i-1)]^alpha)
              between_comm_probs <- between_comm_probs/sum(between_comm_probs)
              
              # add edges based on probabilities
              for (j in 1:m) {
                if (runif(1) < prewire) {
                  # rewire edge to random node within certain distance
                  dist_to_node <- dist_mat[1:(i-1), i]
                  nodes_within_range <- which(dist_to_node <= r)
                  new_neighbor <- sample(nodes_within_range, 1)
                  adj_mat[i, new_neighbor] <- 1
                  adj_mat[new_neighbor, i] <- 1
                  degrees[i] <- degrees[i] + 1
                  degrees[new_neighbor] <- degrees[new_neighbor] + 1
                } else if (runif(1) < sum(within_comm_probs)) {
                  # attach to node in same community based on degree preferential attachment
                  new_neighbor <- sample(1:(i-1), 1, prob = within_comm_probs)
                  adj_mat[i, new_neighbor] <- 1
                  adj_mat[new_neighbor, i] <- 1
                  degrees[i] <- degrees[i] + 1
                  degrees[new_neighbor] <- degrees[new_neighbor] + 1
                } else {
                  # attach to node in different community based on degree preferential attachment
                  new_neighbor <- sample(1:(i-1), 1, prob = between_comm_probs)
                  adj
                  
  
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
                  grow_spatial_scale_free_graph <- function(N, beta, numofComm, m, r, prewire) {
                    
                    # Create distance matrix
                    coords <- matrix(runif(N*2), ncol=2)
                    dist_mat <- as.matrix(dist(coords))
                    
                    # Initialize adjacency matrix and community assignments
                    adj_mat <- matrix(0, nrow = N, ncol = N)
                    comm_assignments <- rep(1, N)
                    
                    # Add the first node and assign it to the first community
                    adj_mat[1, 1] <- 1
                    
                    # Loop over the remaining nodes
                    for (i in 2:N) {
                      
                      # Create arbitrary number of communities/blocks
                      comm_prob_mat <- matrix(runif(numofComm*numofComm), ncol=numofComm)
                      comm_prob_mat <- comm_prob_mat / rowSums(comm_prob_mat)
                      
                      # Compute probabilities for within-community attachment
                      spat_distprobs <- ifelse(dist_mat[1:(i-1), i] <= r, 1/(1+(dist_mat[1:(i-1), i])^beta), 0)
                      degrees <- rowSums(adj_mat[1:(i-1), ])
                      deg_probs <- (degrees^beta)/(sum(degrees^beta))
                      
                      # Compute probabilities for between-community attachment
                      between_probs <- matrix(0, nrow = i-1, ncol = numofComm)
                      for (j in 1:numofComm) {
                        indices <- which(comm_assignments[1:(i-1)] == j)
                        if (length(indices) > 0) {
                          spat_distprobs_j <- ifelse(dist_mat[indices, i] <= r, 1/(1+dist_mat[indices, i]^beta), 0)
                          degrees_j <- rowSums(adj_mat[indices, 1:(i-1)])
                          deg_probs_j <- (degrees_j^beta)/(sum(degrees_j^beta))
                          between_probs[indices, j] <- spat_distprobs_j %*% t(matrix(deg_probs_j, nrow=1)) * comm_prob_mat[j,]
                        }
                      }
                      between_probs <- between_probs / rowSums(between_probs)
                      
                      # Assign new node to a community
                      comm_probs <- rowSums(between_probs)
                      comm_probs <- comm_probs / sum(comm_probs)
                      comm_assignments[i] <- sample(1:numofComm, size=1, prob=comm_probs)
                      
                      # Add edges to existing nodes
                      within_probs <- spat_distprobs * deg_probs
                      within_probs <- within_probs / sum(within_probs)
                      for (j in 1:(i-1)) {
                        if (comm_assignments[j] == comm_assignments[i]) {
                          if (runif(1) < within_probs[j]) {
                            adj_mat[i, j] <- 1
                            adj_mat[j, i] <- 1
                          }
                        } else {
                          between_prob <- between_probs[j, comm_assignments[i]]
                          if (runif(1) < between_prob) {
                            adj_mat[i, j] <- 1
                            adj_mat[j, i] <- 1
                          }
                        }
                      }
                    }
                    
                    # Randomly rewire edges with probability prewire
                    for (i in 1:(N-1)) {
                      for (j in (i+1):N) {
                        if (adj_mat[i, j] == 1) {
                          if (runif(1) < prewire)
                            
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
spatial_scale_free_graph <- function(N, beta, numofComm, r, m, prewire) {
  # Initialize the adjacency matrix with a single node
  adj_mat <- matrix(0, nrow = N, ncol = N)
  adj_mat[1, 1] <- 1
  
  # Initialize the spatial coordinates of the nodes
  coords <- matrix(runif(2), nrow = N, ncol = 2)
  
  # Initialize the degree and community assignment of the nodes
  degrees <- rep(0, N)
  comm_assignment <- rep(1, N)
  
  # Initialize the distance matrix between nodes
  dist_mat <- as.matrix(dist(coords))
  
  # Loop through the remaining nodes
  for (i in 2:N) {
    # Step 1: Create arbitrary number of communities/blocks
    comm_vec <- sample(1:numofComm, i, replace = TRUE)
    
    # Step 2: Define the community probability attachment matrix
    P_ij <- matrix(0, nrow = numofComm, ncol = numofComm)
    for (c1 in 1:numofComm) {
      for (c2 in 1:numofComm) {
        if (c1 == c2) {
          P_ij[c1, c2] <- 1
        } else {
          P_ij[c1, c2] <- runif(1)
        }
      }
    }
    
    # Step 3: Calculate the within community attachment probability matrix
    spat_distprobs <- 1 / (1 + (dist_mat[1:(i-1), i])^beta)
    spat_distprobs[dist_mat[1:(i-1), i] > r] <- 0
    degrees <- rowSums(adj_mat[1:(i-1), ])
    deg_probs <- (degrees^m) / sum(degrees^m)
    within_probs <- outer(deg_probs, deg_probs) * spat_distprobs
    
    # Step 4: Calculate the between community attachment probability matrix
    between_probs <- matrix(0, nrow = i-1, ncol = numofComm)
    for (j in 1:(i-1)) {
      if (comm_vec[j] == comm_vec[i]) {
        between_probs[j, ] <- 0
      } else {
        spat_distprob <- 1 / (1 + (dist_mat[j, i])^beta)
        if (dist_mat[j, i] <= r) {
          degrees <- rowSums(adj_mat[1:(i-1), ])
          deg_probs <- (degrees^m) / sum(degrees^m)
          between_probs[j, ] <- deg_probs * spat_distprob
        } else {
          between_probs[j, ] <- 0
        }
      }
    }
    between_probs <- matrix(colSums(between_probs), nrow = 1)
    between_probs <- matrix(rep(between_probs, i-1), nrow = i-1, ncol = numofComm)
    
    # Step 5: Assign the node to a community/block using a 2D poison point process
    comm_probs <- colSums(P_ij[comm_vec[1:(i-1)], ])
    comm_probs <- comm_probs / sum(comm_probs)
    comm_assignment[i] <- sample(1:numofComm, 1, prob = comm_probs)
    
    # Step 6: Add edges based on the community/block assignments and probability matrices
    for (j in 1:(i-
                 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
spatial_scale_free_graph <- function(N, numofComm, beta, m, r, prewire) {
  # Generate initial set of nodes distributed spatially
  coords <- matrix(runif(N * 2), ncol = 2)
  dist_mat <- as.matrix(dist(coords))
  
  # Initialize adjacency matrix and community assignment
  adj_mat <- matrix(0, nrow = N, ncol = N)
  comm_assign <- rep(1, N)
  
  # Add first node to initial community
  adj_mat[1, 1] <- 1
  
  # Loop over remaining nodes to add to the graph
  for (i in 2:N) {
    # Step 1: Create numofComm communities/blocks
    comm_prob_mat <- matrix(0, nrow = numofComm, ncol = numofComm)
    for (j in 1:numofComm) {
      for (k in 1:numofComm) {
        if (j == k) {
          comm_prob_mat[j, k] <- runif(1, 0.8, 1)
        } else {
          comm_prob_mat[j, k] <- runif(1, 0, 0.2)
        }
      }
    }
    comm_assign[i] <- sample(numofComm, 1, prob = rep(1/numofComm, numofComm))
    
    # Step 2: Define community probability attachment matrix
    spat_dist_probs <- 1 / (1 + (dist_mat[1:(i-1), i])^beta)
    spat_dist_probs[dist_mat[1:(i-1), i] > r] <- 0
    degrees <- rowSums(adj_mat[1:(i-1), ])
    deg_probs <- (degrees^beta) / sum(degrees^beta)
    within_comm_probs <- outer(deg_probs, deg_probs, '*') * spat_dist_probs
    
    # Step 3: Add within-community edges
    for (j in 1:(i-1)) {
      if (comm_assign[i] == comm_assign[j] && runif(1) < sum(within_comm_probs[j, (j+1):(i-1)])) {
        adj_mat[i, j] <- 1
        adj_mat[j, i] <- 1
      }
    }
    
    # Step 4: Add between-community edges
    between_comm_probs <- rep(0, numofComm)
    for (j in 1:numofComm) {
      if (j == comm_assign[i]) {
        between_comm_probs[j] <- 0
      } else {
        spat_dist_probs <- 1 / (1 + (dist_mat[1:(i-1), i])^beta)
        spat_dist_probs[dist_mat[1:(i-1), i] > r] <- 0
        degrees <- rowSums(adj_mat[1:(i-1), ])
        deg_probs <- (degrees^beta) / sum(degrees^beta)
        between_comm_probs[j] <- sum(deg_probs[comm_assign == j] * deg_probs[comm_assign == comm_assign[i]] * spat_dist_probs[comm_assign == j])
      }
    }
    between_comm_probs <- between_comm_probs / sum(between_comm_probs)
    comm_sample <- sample(numofComm, 1, prob = between_comm_probs)
    within_comm_probs <- outer(deg_probs, (rowSums(adj_mat[comm_assign == comm_sample, ])^beta) / sum(rowSums(adj_mat[comm_assign == comm
                                                                                                                      
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to generate a spatial scale-free graph
spatial_scalefree <- function(N, beta, numofComm, prewire, m, r, alpha) {
  
  # Generate initial node
  nodes <- matrix(runif(2), nrow = 1)
  adj_mat <- matrix(0, nrow = 1, ncol = 1)
  comm_assignments <- rep(1, 1)
  degrees <- rep(0, 1)
  
  # Generate remaining nodes
  for (i in 2:N) {
    
    # Step 1: Create communities/blocks
    comm_assignments <- c(comm_assignments, sample(1:numofComm, 1))
    
    # Step 2: Define community probability attachment matrix
    P_ij <- matrix(0, nrow = numofComm, ncol = numofComm)
    for (j in 1:numofComm) {
      for (k in 1:numofComm) {
        if (j == k) {
          P_ij[j, k] <- sum(degrees[comm_assignments == j]^alpha)/sum(degrees^alpha)
        } else {
          P_ij[j, k] <- 1/numofComm
        }
      }
    }
    
    # Step 3: Define within-community probability attachment matrix
    spat_distprobs <- rep(0, i-1)
    deg_probs <- rep(0, i-1)
    for (j in 1:(i-1)) {
      if (dist_mat[j, i] <= r) {
        spat_distprobs[j] <- 1/(1 + (dist_mat[j, i]^beta))
        deg_probs[j] <- degrees[j]^alpha
      }
    }
    if (sum(spat_distprobs) == 0) {
      within_probs <- rep(1/i, i-1)
    } else {
      within_probs <- (spat_distprobs * (deg_probs/sum(deg_probs))) + ((1 - spat_distprobs) * (1/i))
    }
    
    # Step 4: Define between-community probability attachment matrix
    spat_distprobs <- rep(0, i-1)
    deg_probs <- rep(0, i-1)
    for (j in 1:(i-1)) {
      if (dist_mat[j, i] <= r) {
        spat_distprobs[j] <- 1/(1 + (dist_mat[j, i]^beta))
        deg_probs[j] <- degrees[j]^alpha
      }
    }
    if (sum(spat_distprobs) == 0) {
      between_probs <- rep(1/i, i-1)
    }
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    generate_spatial_scalefree_graph <- function(N, numofComm, beta, r, alpha, m, prewire) {
      # Step 1: Generate initial node positions
      node_positions <- matrix(runif(N*2), ncol=2)
      
      # Calculate distance matrix
      dist_mat <- as.matrix(dist(node_positions))
      
      # Step 2: Initialize adjacency matrix and community assignments
      adj_mat <- matrix(0, nrow=N, ncol=N)
      comm_assignments <- rep(1, N)
      
      # Add initial node to first community
      comm_assignments[1] <- 1
      
      # Step 3 - 6: Loop over new nodes and add edges
      for (i in 2:N) {
        # Step 3i: Create random community assignments
        comm_probs <- rep(1/numofComm, numofComm)
        comm_assignments[i] <- sample(1:numofComm, size=1, prob=comm_probs)
        
        # Step 3ii: Define community probability attachment matrix
        P_ij <- matrix(0, nrow=numofComm, ncol=numofComm)
        for (j in 1:numofComm) {
          for (k in 1:numofComm) {
            if (j == k) {
              # Within-community attachment probability
              spat_distprobs <- ifelse(dist_mat[1:(i-1),i] <= r, 1/(1+(dist_mat[1:(i-1),i])^beta), 0)
              degrees <- rowSums(adj_mat[1:(i-1),][comm_assignments[1:(i-1)]==j,])
              deg_probs <- ifelse(sum(degrees)==0, rep(1/sum(degrees), length(degrees)), (degrees^alpha)/(sum(degrees^alpha)))
              P_ij[j,k] <- sum(spat_distprobs * deg_probs)
            } else {
              # Between-community attachment probability
              spat_distprobs <- ifelse(dist_mat[1:(i-1),i] <= r, 1/(1+(dist_mat[1:(i-1),i])^beta), 0)
              degrees <- rowSums(adj_mat[1:(i-1),][comm_assignments[1:(i-1)]==j,])
              deg_probs <- ifelse(sum(degrees)==0, rep(1/sum(degrees), length(degrees)), (degrees^alpha)/(sum(degrees^alpha)))
              between_probs <- sum(spat_distprobs * deg_probs)
              P_ij[j,k] <- between_probs * prewire
            }
          }
        }
        
        # Step 3v: Sample community assignment for new node using 2D Poisson process
        new_node_position <- node_positions[i,]
        comm_probs <- apply(P_ij, 2, sum)
        new_comm_assignment <- sample(1:numofComm, size=1, prob=comm_probs)
        
        # Step 3vi: Add edges between new node and existing nodes
        for (j in 1:(i-1)) {
          if (comm_assignments[j] == new_comm_assignment) {
            # Within-community attachment
            if (runif(1) < P_ij[new_comm_assignment,new_comm_assignment]) {
              adj_mat[i,j] <- 1
              adj_mat[j,i] <- 1
            }
          } else {
            # Between-community attachment
            if (runif(1) < P_ij[comm_assignments[j],new_comm
                                
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

grow_scale_free_expander <- function(N, beta, alpha, r, prewire, m) {
  # initialize network with a single node at the center of the unit square
  nodes <- data.frame(id = 1, x = 0.5, y = 0.5, community = 1, degree = 0)
  
  # initialize the attachment probability matrix
  P <- matrix(0, nrow = 1, ncol = 1)
  
  for (i in 2:N) {
    # determine the community for the new node
    communities <- unique(nodes$community)
    num_communities <- length(communities)
    community_probs <- rep(1/num_communities, num_communities)
    new_community <- sample(communities, size = 1, prob = community_probs)
    # determine the attachment probabilities for each existing node
    distances <- sqrt((nodes$x - nodes[(i-1), "x"])^2 + (nodes$y - nodes[(i-1), "y"])^2)}
    spatial_probs <- as.numeric(distances <= r)
    degree_probs <- (nodes$degree + alpha)^(-beta)}
    if (sum(degree_probs) == 0) {
      degree_probs <- rep(1/nrow(P), nrow(P))
    } else {
      degree_probs <- degree_probs / sum(degree_probs)
    }
    P_new <- outer(degree_probs, spatial_probs, "*")
    P_new <- P_new / sum(P_new)
    P <- rbind(P, cbind(P_new, rep(0, nrow(P_new))))
    P <- cbind(P, c(rep(0, ncol(P_new)), 0))
    
    # choose m nodes to connect to
    connected_nodes <- sample(1:(i-1), size = m, prob = P[i, 1:(i-1)])
    
    # add edges to the network
    for (
      
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
grow_spatial_sf_expander <- function(N, beta, alpha, prewire, r, m) {
  # Initialize the graph with a single node
  x <- runif(1)
  y <- runif(1)
  nodes <- data.frame(id = 1, x = x, y = y, comm = 1, deg = 0)
  adj_mat <- matrix(0, nrow = 1, ncol = 1)
  
  # Add nodes to the graph one at a time
  for (i in 2:N) {
    # Determine the community for the new node
    p_comm <- numeric(i-1)
    for (j in 1:(i-1)) {
      p_comm[j] <- sum(adj_mat[j,]) / (m*(i-1))
    }
    comm <- sample.int(i-1, 1, prob = p_comm) + 1
    
    # Determine the attachment probabilities for the new node
    dist <- sqrt((nodes$x - nodes[comm, "x"])^2 + (nodes$y - nodes[comm, "y"])^2)
    dist[which(dist > r)] <- Inf
    P <- (dist^(-beta)) * (nodes$deg^alpha)
    P[which(is.infinite(P))] <- 0
    P[comm] <- 0
    P <- P / sum(P)
    
    # Determine the edges for the new node
    attach <- sample.int(i-1, m, replace = TRUE, prob = P) + 1
    edges <- cbind(rep(i, m), attach)
    edges <- edges[order(edges[,2]),]
    edges <- unique(edges, by = "V2")
    edges <- edges[order(edges[,1]),]
    
    # Add the new node and edges to the graph
    x <- runif(1)
    y <- runif(1)
    nodes <- rbind(nodes, data.frame(id = i, x = x, y = y, comm = comm, deg = 0))
    for (j in 1:nrow(edges)) {
      adj_mat[edges[j,1], edges[j,2]] <- 1
      adj_mat[edges[j,2], edges[j,1]] <- 1
      nodes[edges[j,1], "deg"] <- nodes[edges[j,1], "deg"] + 1
      nodes[edges[j,2], "deg"] <- nodes[edges[j,2], "deg"] + 1
    }
    
    # Rewire edges for small-world effect
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Function to create a spatial scale-free expander graph
    create_spatial_scale_free_expander_graph <- function(N, beta, alpha, prewire, r, m) {
      # Generate initial set of nodes
      coords <- matrix(runif(N*2, 0, 1), ncol = 2)
      while (min(dist(coords)) < r) {
        coords <- matrix(runif(N*2, 0, 1), ncol = 2)
      }
      
      # Create adjacency matrix based on distance matrix
      dist_mat <- as.matrix(dist(coords))
      adj_mat <- as.matrix(dist_mat <= r)
      
      # Create community structure using stochastic block matrix
      num_communities <- 5
      community_sizes <- rep(N/num_communities, num_communities)
      community_sizes <- community_sizes + (N - sum(community_sizes))
      block_matrix <- matrix(0, nrow = N, ncol = N)
      for (i in 1:num_communities) {
        for (j in 1:num_communities) {
          if (i == j) {
            block_matrix[sum(community_sizes[1:(i-1)])+1:sum(community_sizes[1:i]),
                         sum(community_sizes[1:(j-1)])+1:sum(community_sizes[1:j])] <- 1
          } else {
            block_matrix[sum(community_sizes[1:(i-1)])+1:sum(community_sizes[1:i]),
                         sum(community_sizes[1:(j-1)])+1:sum(community_sizes[1:j])] <- runif(1, 0, 1)
          }
        }
      }
      P_ij <- block_matrix
      
      # Initialize graph
      g <- graph.adjacency(adj_mat, mode = "undirected")
      V(g)$coords <- coords
      
      # Grow graph
      for (i in 1:(N-m)) {
        # Calculate attachment probabilities based on degree and distance
        degree_probs <- degree(g)**alpha
        distance_probs <- exp(-beta*dist_mat[i,])
        attachment_probs <- degree_probs*distance_probs
        attachment_probs[i] <- 0
        attachment_probs <- attachment_probs/sum(attachment_probs)
        
        # Choose m nodes to attach to
        attach_indices <- sample(N, m, prob = attachment_probs, replace = TRUE)
        
        # Add edges to graph
        for (j in attach_indices) {
          if (runif(1, 0, 1) < prewire) {
            # Rewire edge
            candidates <- which(dist_mat[i,] > r & dist_mat[i,] <= 2*r)
            rewire_index <- sample(candidates, 1)
            g <- delete.edges(g, which(get.edgelist(g) == c(i,j) | get.edgelist(g) == c(j,i)))
            g <- add.edges(g, c(i, rewire_index))
          } else {
            # Add edge
            if (runif(1, 0, 1) < P_ij[i,j]) {
              g <- add.edges(g, c(i, j))
            }
          }
        }
      }
      
      return(g)
    }
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Required packages
    library(Matrix)
    library(stats)
    
    # Function to generate spatial scale-free expander graph
    generateGraph <- function(N, beta, alpha, prewire, r, m) {
      
      # Generate initial set of nodes with spatial structure
      coords <- data.frame(x = runif(N), y = runif(N))
      dists <- as.matrix(dist(coords))
      
      # Set adjacency matrix based on distance matrix
      A <- ifelse(dists <= r, 1, 0)
      
      # Generate block matrix for community structure
      nCommunities <- 4
      communitySizes <- rep(floor(N/nCommunities), nCommunities)
      communitySizes[1:N%%nCommunities] <- communitySizes[1:N%%nCommunities] + 1
      blockProb <- matrix(0.5, nCommunities, nCommunities)
      diag(blockProb) <- 0.9
      P <- kronecker(blockProb, matrix(1, nrow = m, ncol = m)) # preferential attachment matrix
      
      # Initialize degree vector and edges
      degree <- rowSums(A)
      edges <- which(A == 1, arr.ind = TRUE)
      
      # Grow the graph
      for (i in 1:(N - m)) {
        
        # Select community for new node
        community <- sample(1:nCommunities, 1)
        
        # Calculate attachment probability based on degree and spatial distance
        dist_i <- dists[i,]
        P_i <- (degree^alpha) * exp(-beta*dist_i) * (1 - A[i,]) # attachment probability for node i
        P_i[which(P_i < 0)] <- 0 # replace negative probabilities with 0
        P_i <- P_i / sum(P_i) # normalize probabilities
        
        # Select m nodes to attach to
        attachNodes <- sample(1:(N-m), m, replace = TRUE, prob = P_i)
        A[i, attachNodes] <- 1
        A[attachNodes, i] <- 1
        edges <- rbind(edges, cbind(rep(i, m), attachNodes))
        
        # Update degree vector
        degree[attachNodes] <- degree[attachNodes] + 1
        degree[i] <- m
        
        # Small-world effect: rewire edges with probability prewire
        for (j in 1:m) {
          if (runif(1) < prewire) {
            # Find random node within certain distance
            candidateNodes <- which(dists[attachNodes[j],] <= r & degree < N)
            candidateNodes <- setdiff(candidateNodes, c(i, attachNodes))
            if (length(candidateNodes) > 0) {
              newAttachNode <- sample(candidateNodes, 1)
              A[i, attachNodes[j]] <- 0
              A[attachNodes[j], i] <- 0
              A[i, newAttachNode] <- 1
              A[newAttachNode, i] <- 1
              edges <- rbind(edges, c(i, newAttachNode))
              degree[attachNodes[j]] <- degree[attachNodes[j]] - 1
              degree[newAttachNode] <- degree[newAttachNode] + 1
            }
          }
        }
        
      }
      
      # Create sparse adjacency matrix
      A <- Matrix(A, sparse = TRUE)
      
      # Return graph object
      return(graph_from_adjacency_matrix(A, mode = "undirected"))
      

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        
    
      # Define parameters
      beta <- 1.5   # power law parameter for degree distribution
      m <- 4        # number of edges to add for each new node
      alpha <- 0.5  # parameter for degree effect in attachment probability matrix
      prewire <- 0.1 # probability of rewiring an edge
      
      N <- 100      # number of nodes in the network
      r <- 0.2      # cutoff distance for spatial distance effect
      nc <- 5       # number of communities
      
      # Initialize the network
      x <- runif(N)  # generate random x-coordinates for nodes
      y <- runif(N)  # generate random y-coordinates for nodes
      d <- as.matrix(dist(cbind(x, y))) # compute distance matrix
      A <- matrix(0, nrow=N, ncol=N)   # adjacency matrix
      W <- matrix(0, nrow=N, ncol=N)   # attachment probability matrix
      
      # Assign nodes to communities
      communities <- rep(1:nc, each=N/nc)
      set.seed(123)  # for reproducibility
      communities <- sample(communities)
      
      # Compute attachment probability matrix
      for (i in 1:N) {
        for (j in 1:N) {
          if (i != j) {
            if (communities[i] == communities[j]) {
              W[i,j] <- exp(-alpha*d[i,j])  # within-community attachment
            } else {
              W[i,j] <- exp(-r*alpha*d[i,j])  # between-community attachment
            }
          }
        }
      }
      
      # Add nodes to the network
      for (i in 2:N) {
        # preferential attachment
        degrees <- colSums(A)
        probs <- (degrees + m*W[i,]) / (sum(degrees) + m)
        new_edges <- sample(1:(i-1), m, replace=TRUE, prob=probs)
        A[i,new_edges] <- 1
        A[new_edges,i] <- 1
        
        # rewiring
        for (j in 1:m) {
          if (runif(1) < prewire) {
            rewired <- FALSE
            while (!rewired) {
              k <- sample(1:N, 1)
              if (d[i,k] < r && A[i,k] == 0) {
                A[i,new_edges[j]] <- 0
                A[new_edges[j],i] <- 0
                A[i,k] <- 1
                A[k,i] <- 1
                rewired <- TRUE
              }
            }
          }
        }
      }
      
      
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
      # Function to generate a spatial scale-free expander graph
      # beta: power law parameter for degree distribution
      # N: number of nodes to add to the graph
      # alpha: parameter for degree preferential attachment
      # prewire: probability of rewiring an edge to a random node
      # r: cutoff distance for spatial distance effect
      # m: number of edges to add for each new node
      # num_comm: number of communities
      # P_ij: attachment probability matrix for nodes in different communities
      spatial_scalefree_graph <- function(beta, N, alpha, prewire, r, m, num_comm, P_ij) {
        # Generate initial set of nodes randomly on unit square
        nodes <- runif(N * 2, 0, 1)
        nodes <- matrix(nodes, ncol=2)
        
        # Apply spatial constraints to initial set of nodes
        for (i in 1:N) {
          while (any(sqrt(rowSums((nodes[i,] - nodes)^2))) < 0.05) {
            nodes[i,] <- runif(2, 0, 1)
          }
        }
        
        # Generate initial adjacency matrix based on distance matrix and cutoff distance r
        dist_mat <- as.matrix(dist(nodes))
        adj_mat <- ifelse(dist_mat <= r, 1, 0)
        
        # Initialize stochastic block matrix for community structure
        block_mat <- matrix(rep(0, N * N), ncol=N)
        
        # Add new nodes to the graph one by one
        for (i in 1:(N-1)) {
          # Initialize new row and column for adjacency matrix
          adj_mat <- rbind(adj_mat, rep(0, N))
          adj_mat <- cbind(adj_mat, rep(0, N+1))
          
          # Calculate degree distribution of existing nodes
          degrees <- rowSums(adj_mat[1:i,])
          prob_vec <- (degrees / sum(degrees))^alpha
          prob_vec <- prob_vec / sum(prob_vec)
          
          # Generate m edges for new node using preferential attachment
          for (j in 1:m) {
            # Select a node to attach to based on degree preferential attachment
            attach_node <- sample(1:i, 1, prob=prob_vec)
            
            # Calculate distance between new node and selected node
            distance <- sqrt(sum((nodes[i+1,] - nodes[attach_node,])^2))
            
            # Determine whether to connect based on attachment probability matrix P_ij
            if (sample(1:2, 1, prob=c(P_ij[block_mat[i+1,attach_node]+1,]), replace=TRUE) == 1) {
              # Add edge to adjacency matrix
              adj_mat[i+1, attach_node] <- 1
              adj_mat[attach_node, i+1] <- 1
              
              # Update block matrix for new node and selected node
              block_mat[i+1, attach_node] <- sample(1:num_comm, 1)
              block_mat[attach_node, i+1] <- block_mat[i+1, attach_node]
              
              # Rewire edge with probability prewire
              if (runif(1) < prewire) {
                # Find all nodes within distance r
                dist_vec <- dist_mat[i+1,]
                potential_nodes <- which(dist_vec <= r & dist_vec > 0)
                
                # Remove nodes that are already connected to new node
                potential_nodes <- setdiff(potential_nodes, which(adj_mat[i+
                                                                            
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        
library(Matrix)
library(pracma)
library(mvtnorm)

grow_spatial_scalefree_expander <- function(N, beta, alpha, prewire, r, m) {
  # Generate initial nodes with spatial structure
  coords <- matrix(0, nrow = N, ncol = 2)
  coords[1, ] <- runif(2)
  for (i in 2:N) {
    repeat {
      proposed_coords <- runif(2)
      if (min(dist(rbind(coords[1:(i-1),], proposed_coords))) > 0.2) {
        coords[i,] <- proposed_coords
        break
      }
    }
  }
  
  # Create adjacency matrix
  dist_mat <- pdist2(coords, coords)
  adj_mat <- ifelse(dist_mat <= r, 1, 0)
  
  # Generate stochastic block model for community structure
  block_sizes <- sample(1:ceiling(sqrt(N)), floor(N/2), replace = TRUE)
  block_assignments <- rep(1:length(block_sizes), times = block_sizes)
  num_blocks <- length(block_sizes)
  P <- matrix(0, nrow = num_blocks, ncol = num_blocks)
  for (i in 1:num_blocks) {
    for (j in 1:num_blocks) {
      if (i == j) {
        P[i, j] <- rbinom(1, 1, 0.9)
      } else {
        P[i, j] <- rbeta(1, 1, 5)
      }
    }
  }
  P <- P + t(P)
  for (i in 1:num_blocks) {
    P[i, i] <- 0
    P[i, ] <- P[i, ] / sum(P[i, ])
  }
  block_prob_mat <- P[block_assignments, block_assignments]
  
  # Initialize degree sequence
  degree_seq <- rep(m, N)
  
  # Loop over nodes to add
  for (i in (m+1):N) {
    # Generate attachment probabilities
    node_probs <- rep(0, N)
    for (j in 1:N) {
      dist_prob <- exp(-alpha*dist_mat[i,j])
      comm_prob <- block_prob_mat[block_assignments[i], block_assignments[j]]
      degree_prob <- degree_seq[j]^beta
      node_probs[j] <- dist_prob * comm_prob * degree_prob
    }
    node_probs[i] <- 0
    node_probs <- node_probs / sum(node_probs)
    
    # Attach to m nodes with preferential attachment
    new_edges <- sample(1:N, m, replace = FALSE, prob = node_probs)
    adj_mat[i, new_edges] <- 1
    adj_mat[new_edges, i] <- 1
    degree_seq[new_edges] <- degree_seq[new_edges] + 1
    degree_seq[i] <- m
    
    # Rewire with probability prewire
    for (j in new_edges) {
      if (runif(1) < prewire) {
        dist_vec <- dist_mat[i,] - dist_mat[j,]
        candidate_edges <- which(abs(dist_vec) <= r & dist_vec != 0)
        if (length(candidate_edges) > 0) {
          rewire_to <- sample(candidate_edges, 1)
          adj_mat[i, j] <- 0
          adj_mat[j, i] <- 0
          adj_mat[i, rewire_to] <- 1
          adj_mat[rewire_to, i] <- 1
          
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
          # function to create a spatial scale-free expander graph
          # with strong community structure and small-world effect
          spatial_scale_free_expander_graph <- function(N, beta, alpha, r, m, prewire){
            
            # generate N nodes on a unit square
            nodes <- matrix(runif(N*2), ncol=2)
            
            # create a distance matrix based on the spatial locations of the nodes
            dist_mat <- as.matrix(dist(nodes))
            
            # create an adjacency matrix based on the distance matrix and cutoff distance r
            adj_mat <- ifelse(dist_mat <= r, 1, 0)
            
            # initialize the community membership of each node
            comm_mem <- rep(1, N)
            
            # initialize the degree sequence with m nodes, where each node is connected to all other nodes in its community
            deg_seq <- rep(0, N)
            for(i in 1:m){
              deg_seq[i] <- m - 1
              for(j in 1:m){
                if(i != j && comm_mem[i] == comm_mem[j]){
                  adj_mat[i, j] <- 1
                  adj_mat[j, i] <- 1
                  deg_seq[i] <- deg_seq[i] + 1
                  deg_seq[j] <- deg_seq[j] + 1
                }
              }
            }
            
            # loop through the remaining nodes and add them to the network
            for(i in (m+1):N){
              # calculate the probability of attaching to each existing node
              attach_probs <- rep(0, N)
              for(j in 1:N){
                if(comm_mem[i] == comm_mem[j]){
                  attach_probs[j] <- (deg_seq[j] + alpha) / (2*m + i - 1 + alpha*N)
                } else {
                  attach_probs[j] <- (deg_seq[j] + 1) / (2*m + i - 1)
                }
              }
              # sample m nodes to attach to with probability proportional to their degree
              attach_nodes <- sample(1:N, m, replace=FALSE, prob=attach_probs)
              # add edges between the new node and the selected nodes
              for(j in attach_nodes){
                adj_mat[i, j] <- 1
                adj_mat[j, i] <- 1
                deg_seq[i] <- deg_seq[i] + 1
                deg_seq[j] <- deg_seq[j] + 1
              }
              # randomly assign the new node to a community
              comm_mem[i] <- sample(1:max(comm_mem), 1)
              # rewire edges with probability prewire
              for(j in 1:N){
                if(runif(1) < prewire && adj_mat[i, j] == 1){
                  # find all nodes within distance r that are not in the same community
                  eligible_nodes <- which(dist_mat[i,] <= r & comm_mem != comm_mem[i] & seq_len(N) != i)
                  if(length(eligible_nodes) > 0){
                    # randomly select one of the eligible nodes to rewire to
                    new_node <- sample(eligible_nodes, 1)
                    adj_mat[i, j] <- 0
                    adj_mat[j, i] <- 0
                    adj_mat[i, new_node] <- 1
                    adj_mat[new_node, i] <- 1
                  }
                }
              }
            }
            
            # convert adjacency matrix to an igraph object for visualization and analysis
            g <- igraph::graph_from_adjacency_matrix(adj
                                                     
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to generate a spatial scale-free expander graph
# with the given parameters:
# beta: power-law exponent for degree distribution
# N: number of nodes to generate
# alpha: parameter controlling strength of degree effect
# r: spatial cutoff distance
# prewire: rewiring probability
# m: number of edges to add per node
spatial_sf_expander <- function(beta, N, alpha, r, prewire, m) {
  # Generate initial set of nodes with spatial structure
  coords <- matrix(runif(2*N), ncol=2)
  for (i in 2:N) {
    # Add new node
    new_coord <- runif(2)
    # Check distance to existing nodes
    while (min(sqrt(rowSums((coords - new_coord)^2))) < 0.1) {
      new_coord <- runif(2)
    }
    coords <- rbind(coords, new_coord)
  }
  # Compute distance matrix
  dist_mat <- as.matrix(dist(coords))
  # Create adjacency matrix based on distance matrix
  adj_mat <- (dist_mat <= r)
  # Initialize list of communities
  comm_list <- list(1:N)
  # Add nodes and edges to graph
  for (i in (N+1):(2*N)) {
    # Compute attachment probabilities
    attach_prob <- rep(0, i-1)
    for (j in 1:(i-1)) {
      # Compute distance-based attachment probability
      if (dist_mat[i,j] <= r) {
        dist_term <- 1 / (dist_mat[i,j]^alpha)
      } else {
        dist_term <- 0
      }
      # Compute degree-based attachment probability
      if (length(comm_list[[j]]) > 0) {
        deg_term <- (length(comm_list[[j]]))^-beta
      } else {
        deg_term <- 0
      }
      # Compute overall attachment probability
      attach_prob[j] <- dist_term * deg_term
    }
    # Normalize probabilities and sample m nodes to connect to
    attach_prob <- attach_prob / sum(attach_prob)
    conn_nodes <- sample((1:(i-1)), size=m, replace=FALSE, prob=attach_prob)
    # Add new node to its own community
    comm_list[[i]] <- i
    # Connect new node to selected nodes
    for (j in conn_nodes) {
      # Add edge between nodes
      adj_mat[i,j] <- 1
      adj_mat[j,i] <- 1
      # Update communities
      comm_i <- comm_list[[i]]
      comm_j <- comm_list[[j]]
      if (!identical(comm_i, comm_j)) {
        comm_list[[i]] <- union(comm_i, comm_j)
        comm_list[[j]] <- NULL
      }
    }
    # Rewire edges with probability prewire
    for (j in 1:m) {
      if (runif(1) < prewire) {
        # Choose a random node within a certain distance to rewire to
        dist_i <- dist_mat[i,]
        dist_i[comm_list[[i]]] <- Inf
        max_dist <- min(dist_i)
        rewire_nodes <- which(dist_i <= max_dist*prewire)
        if (length(rewire_nodes) > 0) {
          rewire_node <- sample(rewire_nodes, size=1)
          # Remove edge between nodes
          adj_mat[i,conn_nodes[j]] <- 0
          adj_mat[conn_nodes[j],
                  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
grow_spatial_scale_free_expander_graph <- function(N, beta, alpha, r, prewire, m) {
  # Generate initial set of nodes on a unit square with random spatial coordinates
  coords <- matrix(runif(N * 2), ncol = 2)
  # Compute pairwise distances between nodes
  dists <- as.matrix(dist(coords))
  # Compute attachment probability matrix based on spatial distance and degree preference
  P <- matrix(0, nrow = N, ncol = N)
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        d <- dists[i, j]
        P[i, j] <- (d <= r) * (1 + alpha) * (m + 1)^(beta - alpha) * degree[i]^alpha / sum(degree^beta)
      }
    }
  }
  # Grow network by iteratively adding nodes and connecting them to existing nodes
  edges <- matrix(0, nrow = 0, ncol = 2)
  degree <- rep(0, N)
  for (i in 1:N) {
    # Categorize node i into a community based on its spatial coordinate
    x <- coords[i, 1]
    y <- coords[i, 2]
    community <- floor(x * sqrt(m)) + floor(y * sqrt(m)) * sqrt(m)
    # Select m existing nodes to connect to based on attachment probability matrix
    for (j in 1:m) {
      while (TRUE) {
        k <- sample(N, 1, prob = P[i, ])
        if (k != i & degree[k] < m) break
      }
      # Add edge between i and k
      edges <- rbind(edges, c(i, k))
      degree[i] <- degree[i] + 1
      degree[k] <- degree[k] + 1
      # Update attachment probability matrix for k based on new degree
      for (l in 1:N) {
        if (l != i & l != k) {
          d <- dists[k, l]
          if (d <= r) {
            P[k, l] <- P[k, l] + (1 + alpha) * (m + 1)^(beta - alpha) * degree[k]^alpha / sum(degree^beta)
            P[l, k] <- (1 - prewire) * P[l, k] + prewire * (1 + alpha) * (m + 1)^(beta - alpha) * degree[l]^alpha / sum(degree^beta)
          }
        }
      }
    }
  }
  # Rewire edges with probability prewire
  for (i in 1:nrow(edges)) {
    if (runif(1) <= prewire) {
      k <- sample(N, 1)
      if (k != edges[i, 1] & k != edges[i, 2]) {
        edges[i, 2] <- k
      }
    }
  }
  # Create igraph object from edge list
  g <- igraph::graph_from_edgelist(edges)
  # Add spatial coordinates as vertex attributes
  igraph::V(g)$x <- coords[, 1]
  igraph::V(g)$y <- coords[, 2]
  # Return resulting graph
  return(g)
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to generate a spatial scale-free expander graph
# Parameters: beta - power law parameter for degree distribution
#             alpha - parameter for the strength of degree effect
#             prewire - rewiring probability for small-world effect
#             r - parameter controlling the degree of sparsity
#             m - number of edges added per new node
#             N - total number of nodes to generate
spatial_scale_free_expander <- function(beta, alpha, prewire, r, m, N) {
  
  # Generate initial set of nodes distributed spatially on a unit square
  nodes <- data.frame(x = runif(N), y = runif(N))
  
  # Initialize adjacency matrix
  adjacency_matrix <- matrix(0, nrow = N, ncol = N)
  
  # Initialize community labels for each node
  community_labels <- 1:N
  
  # For each new node
  for (i in 1:(N-1)) {
    
    # Determine the community of the new node
    community <- sample(unique(community_labels), size = 1)
    
    # Calculate the attachment probabilities to all existing nodes based on spatial distance and degree
    attachment_probs <- rep(0, N)
    for (j in 1:N) {
      if (community_labels[j] == community & j != i) {
        # Calculate spatial distance
        dist_ij <- sqrt((nodes$x[i] - nodes$x[j])^2 + (nodes$y[i] - nodes$y[j])^2)
        # Calculate degree effect
        degree_effect <- (degree(adjacency_matrix, j)^alpha) / sum(degree(adjacency_matrix)^alpha)
        # Calculate attachment probability
        attachment_probs[j] <- (dist_ij^(-beta)) * degree_effect
      }
    }
    attachment_probs <- attachment_probs / sum(attachment_probs)
    
    # Choose m nodes to connect to
    nodes_to_connect <- sample(1:N, size = m, prob = attachment_probs, replace = TRUE)
    
    # Connect the new node to the chosen nodes
    for (j in nodes_to
         
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to generate a spatial scale-free expander graph
spatial_scale_free_expander_graph <- function(N, beta, alpha, r, prewire, m) {
  
  # Generate initial set of nodes with spatial structure
  x <- runif(N)
  y <- runif(N)
  
  # Calculate pairwise distances between nodes
  dist_mat <- as.matrix(dist(cbind(x, y)))
  
  # Initialize adjacency matrix and community labels
  A <- matrix(0, N, N)
  C <- rep(1, N)
  
  # Grow network
  for (i in 1:N) {
    # Calculate attachment probabilities
    d <- dist_mat[i, ]
    P <- (d <= r) * (1 + alpha) * (seq(N) ^ (-beta))
    P[i] <- 0
    P <- P / sum(P)
    
    # Add edges to nodes in same community
    comm_nodes <- which(C == C[i])
    comm_nodes <- comm_nodes[-which(comm_nodes == i)]
    comm_probs <- P[comm_nodes]
    comm_edges <- sample(comm_nodes, m, replace = TRUE, prob = comm_probs)
    A[i, comm_edges] <- 1
    A[comm_edges, i] <- 1
    
    # Add edges to nodes in different communities
    other_nodes <- which(C != C[i])
    other_probs <- P[other_nodes]
    other_edges <- sample(other_nodes, m, replace = TRUE, prob = other_probs)
    A[i, other_edges] <- 1
    
    # Update community labels
    if (length(other_edges) > 0) {
      C[other_edges] <- max(C) + 1
    }
    
    # Rewire edges with small world effect
    rewire_edges <- sample(which(A[i, ] == 1), sum(A[i, ] == 1) * m, replace = TRUE)
    for (j in rewire_edges) {
      if (runif(1) < prewire) {
        # Find candidate nodes to rewire to
        candidates <- which(dist_mat[j, ] > r & C != C[i])
        if (length(candidates) > 0) {
          # Rewire to random candidate node
          rewired_node <- sample(candidates, 1)
          A[i, j] <- 0
          A[j, i] <- 0
          A[i, rew
            
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# function to generate a spatial scale-free expander graph
# input: N - number of nodes to add, beta - power law parameter, alpha - degree effect parameter,
#        prewire - rewiring probability, m - number of edges to attach, r - spatial distance parameter
# output: adjacency matrix of the graph
generate_spatial_scale_free_expander <- function(N, beta, alpha, prewire, m, r){
  # initialize coordinates and adjacency matrix
  coordinates <- matrix(runif(N*2), ncol = 2)
  adjacency <- matrix(0, nrow = N, ncol = N)
  
  # add nodes one by one
  for (i in 1:N){
    # compute attachment probabilities to other nodes
    distances <- sqrt(rowSums((coordinates - coordinates[i,])^2))
    degrees <- rowSums(adjacency)
    attachment_probs <- (distances <= r) * (degrees + alpha)^(beta)
    attachment_probs[i] <- 0 # self-loops not allowed
    
    # select m nodes to connect to
    if (sum(attachment_probs) == 0){
      # if no nodes are eligible, connect to a random node
      target_nodes <- sample(1:N, m, replace = FALSE)
    } else {
      # select target nodes with probability proportional to attachment probabilities
      attachment_probs <- attachment_probs / sum(attachment_probs)
      target_nodes <- sample(1:N, m, replace = FALSE, prob = attachment_probs)
    }
    
    # add edges to target nodes
    adjacency[i, target_nodes] <- 1
    adjacency[target_nodes, i] <- 1
    
    # rewire edges with probability prewire
    rewire_probs <- matrix(runif(m^2), ncol = m)
    rewire_mask <- (rewire_probs <= prewire) & (distances[target_nodes] > r)
    rewire_from <- which(rewire_mask, arr.ind = TRUE)[,1]
    rewire_to <- sample(setdiff(1:N, c(i, target
                                       
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to generate a spatial scale-free expander graph
# with the specified properties
generate_graph <- function(N, beta, alpha, r, m, prewire) {
  
  # Generate initial set of nodes with spatial structure
  coords <- matrix(runif(N*2), ncol = 2)
  
  # Compute distance matrix
  dist_mat <- as.matrix(dist(coords))
  
  # Generate stochastic block matrix for community structure
  p1 <- matrix(runif(m*m), nrow = m)
  p2 <- matrix(runif(m*m), nrow = m)
  P <- rbind(cbind(p1, p2), cbind(p2, p1))
  P <- P / sum(P)
  
  # Initialize adjacency matrix
  A <- matrix(0, nrow = N, ncol = N)
  
  # Add nodes to the network
  for (i in 1:N) {
    # Determine community membership for new node
    comm_prob <- apply(A, 1, function(x) sum(x)^(alpha*beta) / sum(x^beta))
    comm_prob[i] <- 0
    comm_prob <- comm_prob / sum(comm_prob)
    comm <- sample(1:m, size = 1, prob = comm_prob)
    
    # Determine attachment probabilities for new node
    attach_prob <- rep(0, N)
    for (j in 1:N) {
      if (j == i) {
        next
      }
      if (dist_mat[i,j] > r) {
        next
      }
      if (A[i,j] == 1) {
        next
      }
      attach_prob[j] <- sum(A[which(A[,comm] == 1), j]^beta)^(alpha*beta) / sum(A[which(A[,comm] == 1), ]^beta)
    }
    attach_prob <- attach_prob / sum(attach_prob)
    
    # Choose neighbors for new node
    neighbors <- sample(1:N, size = m, replace = TRUE, prob = attach_prob)
    
    # Add edges to adjacency matrix
    A[i,neighbors] <- 1
    A[neighbors,i] <- 1
    
    # Rewire edges with probability prewire
    for (j in neighbors) {
      if (runif(1) < prewire) {
        dist_j <- dist_mat[j,]
        
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # function to generate spatial scale-free expander graph
        generate_spatial_scalefree_expander_graph <- function(N, r, m, beta, alpha, prewire) {
          # generate initial node coordinates on a unit square
          coords <- matrix(runif(N*2), ncol = 2)
          
          # initialize adjacency matrix and degree vector
          adj_mat <- matrix(0, nrow = N, ncol = N)
          deg_vec <- rep(0, N)
          
          # add first node to the network
          nodes <- 1
          
          # add remaining nodes
          for (i in 2:N) {
            # calculate attachment probability matrix
            dist_mat <- as.matrix(dist(coords[1:i-1,]))
            P_ij <- (dist_mat < r) * (1 + alpha * deg_vec[1:i-1])^(-beta)
            
            # generate community assignment using stochastic block matrix techniques
            communities <- sample(1:ceiling(sqrt(i)), i, replace = TRUE)
            
            # add m edges to the network
            for (j in 1:m) {
              # choose target node based on attachment probability
              probs <- P_ij[i,] * (communities == communities[i])
              probs[i] <- 0
              probs <- probs / sum(probs)
              target <- sample(1:i-1, 1, prob = probs)
              
              # add edge to target node
              adj_mat[i,target] <- 1
              adj_mat[target,i] <- 1
              
              
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              generate_spatial_scale_free_graph <- function(N, beta, alpha, r, prewire, m) {
                
                # Generate initial set of nodes with spatial structure
                coords <- matrix(runif(N*2), ncol=2)
                
                # Calculate distances between nodes
                dists <- as.matrix(dist(coords))
                
                # Initialize adjacency matrix
                adj_matrix <- matrix(0, nrow=N, ncol=N)
                
                # Add first node to the network
                node_list <- 1
                
                # Add remaining nodes to the network
                for (i in 2:N) {
                  
                  # Initialize community membership for new node
                  community_memberships <- sample(1:max(node_list), size=m, replace=TRUE)
                  
                  # Compute attachment probabilities to existing nodes
                  dists_to_new_node <- dists[i, node_list]
                  degrees_of_existing_nodes <- colSums(adj_matrix[node_list,])
                  degree_power_law_effect <- degrees_of_existing_nodes^(-beta*alpha)
                  degree_power_law_effect <- degree_power_law_effect / sum(degree_power_law_effect)
                  spatial_effect <- exp(-(dists_to_new_node/r)^2)
                  spatial_effect <- spatial_effect / sum(spatial_effect)
                  attachment_probabilities <- degree_power_law_effect * spatial_effect
                  
                  # Add new node to the network
                  node_list <- c(node_list, i)
                  
                  # Connect new node to existing nodes based on attachment probabilities
                  for (j in 1:m) {
                    neighbor <- sample(node_list, size=1, prob=attachment_probabilities)
                    adj_matrix[i, neighbor] <- 1
                    adj_matrix[neighbor, i] <- 1
                  }
                  
                  # Randomly rewire edges with probability prewire
                  for (j in 1:i-1) {
                    if (adj_matrix[i,j] == 1 && runif(1) < prewire) {
                      candidate_neighbors <- which(dists[i,] > r & adj_matrix[j,] == 0)
                      if (length(candidate_neighbors) > 0) {
                        new_neighbor <- sample(candidate_neighbors, size=1)
                        adj_matrix[i,j] <- 0
                        adj_matrix[j,i] <- 0
                        adj_matrix[i,new_neighbor] <- 1
                        adj_matrix[new_neighbor,i] <- 1
                      }
                    }
                  }
                  
                }
                
                # Return adjacency matrix
                return(adj_matrix)
              }
              
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              # Function to generate a spatial scale-free expander graph
              # with specified parameters
              # Input:
              #   - beta: power law parameter for degree distribution
              #   - alpha: parameter that controls the strength of the degree effect
              #   - prewire: rewiring probability for small world effect
              #   - m: number of edges to attach from a new node to existing nodes
              #   - N: number of nodes to generate
              #   - r: spatial distance parameter
              # Output:
              #   - adjacency matrix of the generated graph
              
              spatial_scale_free_expander_graph <- function(beta, alpha, prewire, m, N, r) {
                # Generate initial set of nodes with spatial structure
                x <- runif(N)
                y <- runif(N)
                
                # Initialize adjacency matrix
                A <- matrix(0, nrow=N, ncol=N)
                
                # Add first node
                A[1,1] <- 1
                
                # Initialize degree sequence
                d <- rep(1, N)
                
                # Add nodes one at a time
                for (i in 2:N) {
                  # Categorize nodes into communities using stochastic block matrix
                  if (i <= N/2) {
                    comm <- rep(1, i)
                  } else {
                    comm <- rep(2, i)
                  }
                  
                  # Compute attachment probabilities to existing nodes
                  dist <- sqrt((x[i]-x)^2 + (y[i]-y)^2)
                  dist[dist > r] <- Inf
                  deg <- d^(alpha*beta)
                  deg[i] <- 0
                  prob <- deg*exp(-beta*log(dist))
                  prob[prob < 0] <- 0
                  prob[prob > 1] <- 1
                  
                  # Preferential attachment to existing nodes
                  for (j in 1:m) {
                    # Find node to attach to
                    node <- sample(1:N, 1, prob=prob)
                    A[i,node] <- 1
                    A[node,i] <- 1
                    d[i] <- d[i] + 1
                    d[node] <- d[node] + 1
                  }
                  
                  # Rewire edges with small world probability
                  for (j in 1:(i-1)) {
                    if (runif(1) < prewire) {
                      if (A[i,j] == 1) {
                        # Rewire outgoing edge
                        dist <- sqrt((x[j]-x)^2 + (y[j]-y)^2)
                        dist[dist > r] <- Inf
                        deg <- d^(alpha*beta)
                        deg[i] <- 0
                        prob <- deg*exp(-beta*log(dist))
                        prob[prob < 0]
                        
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        # Function to generate a spatial scale-free expander graph
                        # with the given parameters
                        #
                        # Parameters:
                        # N - number of nodes
                        # beta - power law parameter for degree distribution
                        # alpha - parameter that controls strength of degree effect
                        # r - spatial distance constraint parameter
                        # prewire - rewiring probability
                        # m - number of edges to attach from a new node to existing nodes
                        #
                        # Returns:
                        # an adjacency matrix representing the graph
                        generate_graph <- function(N, beta, alpha, r, prewire, m) {
                          # Initialize the adjacency matrix
                          adj_mat <- matrix(0, nrow = N, ncol = N)
                          
                          # Generate the initial set of nodes with a spatial structure
                          coords <- matrix(runif(N * 2), nrow = N, ncol = 2)
                          
                          # Add the first node
                          nodes_added <- 1
                          adj_mat[nodes_added, nodes_added] <- 1
                          
                          # Add the rest of the nodes
                          while (nodes_added < N) {
                            # Calculate the probability matrix for attachment
                            d <- as.matrix(dist(coords[1:nodes_added, ]))
                            P <- 1 / (1 + (d / r)^alpha)
                            P[upper.tri(P)] <- 0
                            
                            # Add the new node
                            nodes_added <- nodes_added + 1
                            adj_mat[nodes_added, nodes_added] <- 1
                            
                            # Attach edges to existing nodes
                            for (i in 1:m) {
                              # Choose a community for the new node
                              comm_probs <- rep(1/nodes_added, nodes_added)
                              comm <- sample(1:nodes_added, 1, prob = comm_probs)
                              
                              # Choose a node to attach to within the community
                              node_probs <- P[comm, 1:nodes_added]
                              node_probs[adj_mat[comm, ] == 1] <- 0
                              node_probs[comm] <- 0
                              node_probs <- node_probs / sum(node_probs)
                              node <- sample(1:nodes_added, 1, prob = node_probs)
                              
                              # Attach the edge
                              adj_mat[nodes_added, node] <- 1
                              adj_mat[node, nodes_added] <- 1
                            }
                            
                            # Rewire edges with probability prewire
                            for (i in 1:nodes_added) {
                              for (j in (i+1):nodes_added) {
                                if (adj_mat[i, j] == 1 && runif(1) < prewire) {
                                  # Choose a random node to rewire to
                                  k <- sample(1:nodes_added, 1)
                                  while (k == i || k == j || adj_mat[i, k] == 1) {
                                    k <- sample(1:nodes_added, 1)
                                  }
                                  
                                  # Rewire the edge
                                  adj_mat[i, j] <- 0
                                  adj_mat[j, i] <- 0
                                  adj_mat[i, k] <- 1
                                  adj_mat[k, i] <- 1
                                }
                              }
                            }
                          }
                          
                          return(adj_mat)
                        }
                        
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        spatialScaleFreeExpanderGraph <- function(N, beta, alpha, r, prewire, m) {
                          # Create initial node in the center of the unit square
                          nodes <- data.frame(x = 0.5, y = 0.5)
                          
                          # Create adjacency matrix for the initial node
                          adjMat <- matrix(0, nrow = 1, ncol = 1)
                          
                          # Initialize stochastic block matrix for community assignment
                          commProb <- matrix(1, nrow = 1, ncol = 1)
                          
                          # Grow the rest of the network
                          for (i in 2:N) {
                            # Assign node to a community using stochastic block matrix technique
                            community <- sample(1:ncol(commProb), 1, prob = commProb[1,])
                            
                            # Compute attachment probability matrix
                            distMat <- as.matrix(dist(nodes))
                            P_ij <- 1 / (1 + (distMat/r)^alpha)
                            P_ij[lower.tri(P_ij)] <- 0
                            P_ij[which(P_ij > 0)] <- P_ij[which(P_ij > 0)]^beta
                            
                            # Select m nodes to attach to
                            neighbors <- sample(which(adjMat[community,] == 1), m, replace = TRUE)
                            nonNeighbors <- sample(which(adjMat[community,] == 0), m - length(neighbors), replace = TRUE)
                            neighbors <- c(neighbors, nonNeighbors)
                            
                            # Attach node to selected nodes with probability P_ij
                            newAdj <- rep(0, N)
                            for (j in neighbors) {
                              attachProb <- P_ij[i,j] * (1 - prewire) + prewire/N
                              if (runif(1) < attachProb) {
                                newAdj[j
                                       
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to grow spatial scale-free expander graph
# Inputs:
#  N: number of nodes
#  beta: power law parameter for degree distribution
#  alpha: parameter controlling the strength of the degree effect
#  r: parameter controlling the degree of sparsity
#  prewire: rewiring probability for small world effect
#  m: number of edges to add per new node
# Output:
#  adjacency matrix of the resulting graph

grow_graph <- function(N, beta, alpha, r, prewire, m) {
  # Generate initial set of nodes with spatial structure
  coords <- runif(N, 0, 1)
  x <- coords
  y <- coords
  # Initialize adjacency matrix
  A <- matrix(0, nrow=N, ncol=N)
  # Add first node to graph
  A[1,1] <- 1
  # Initialize degree vector
  deg <- rep(1, N)
  # Grow rest of the graph
  for (i in 2:N) {
    # Compute community assignment probabilities based on spatial distance
    d <- sqrt((x[i]-x)^2 + (y[i]-y)^2)
    w <- exp(-d^2/r^2)
    # Normalize probabilities
    w <- w/sum(w)
    # Assign node to community
    comm <- sample(1:N, 1, prob=w)
    # Compute attachment probabilities
    p <- rep(0, N)
    for (j in 1:N) {
      if (j != i) {
        # Compute spatial distance
        d_ij <- sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2)
        # Compute degree effect
        deg_effect <- (deg[j]+alpha)^(beta/2) / (deg[i]+deg[j]+alpha*N^(beta/2))
        # Compute attachment probability
        p[j] <- w[j] * deg_effect
      }
    }
    # Normalize probabilities
    p <- p/sum(p)
    # Choose m edges to add
    edges <- sample(1:N, m, prob=p, replace=FALSE)
    # Add edges
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Function to generate a spatial scale-free expander graph
    # Parameters:
    #   - N: number of nodes
    #   - beta: power law parameter for degree distribution
    #   - alpha: parameter controlling the strength of the degree effect
    #   - r: spatial constraint parameter
    #   - m: number of edges to attach from a new node to existing nodes
    #   - prewire: rewiring probability for small world effect
    # Output:
    #   - adjacency matrix of the generated graph
    spatial_scale_free_expander <- function(N, beta, alpha, r, m, prewire) {
      
      # Generate the initial set of nodes on a unit square
      x <- runif(N)
      y <- runif(N)
      
      # Initialize adjacency matrix
      A <- matrix(0, N, N)
      
      # Add first node
      A[1,1] <- 0
      
      # Initialize degree vector
      k <- rep(0, N)
      k[1] <- 1
      
      # Add remaining nodes
      for (i in 2:N) {
        
        # Categorize nodes into communities
        community <- sample(1:ceiling(sqrt(N)), 1, prob = rep(1:ceiling(sqrt(N)), each = ceiling(N/ceiling(sqrt(N)))))
        
        # Compute attachment probabilities
        d <- sqrt((x - x[i])^2 + (y - y[i])^2)
        d[i] <- Inf
        P <- (d < r) * (k^alpha) / sum((d < r) * (k^alpha))
        P[i,] <- 0
        
        # Generate m edges based on attachment probabilities
        for (j in sample(1:N, m, replace = FALSE, prob = P)) {
          A[i,j] <- 1
          A[j,i] <- 1
          k[i] <- k[i] + 1
          k[j] <- k[j] + 1
        }
        
        # Small world rewiring
        for (j in which(A[i,] == 1)) {
          if (runif(1) < prewire) {
            dist <- sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2)
            A[i,j] <- 0
            A[j,i] <- 0
            candidates <- which(d < r & A[j,] == 0 & d > dist)
            if (length(candidates) > 0) {
              new_j <- sample(candidates, 1)
              A[i,new_j] <- 1
              A[new_j,i] <- 1
            } else {
              A[i,j] <- 1
              A[j,i] <- 1
            }
          }
        }
        
        # Add new node
        x[i] <- runif(1)
        y[i] <- runif(1)
        A[i,i] <- 0
        k[i] <- 0
        
      }
      
      return(A)
      
    }
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ 
spatial_scale_free_expander_graph <- function(N, alpha, beta, prewire, r) {
  # Generate spatial distribution of nodes
  coords <- matrix(runif(2 * N), ncol = 2)
  while (min(dist(coords)) < 0.1) { # prevent nodes from being too close
    coords <- matrix(runif(2 * N), ncol = 2)
  }
  
  # Create initial adjacency matrix based on spatial constraints
  dists <- as.matrix(dist(coords))
  adj <- ifelse(dists <= r, 1, 0)
  
  # Initialize degree sequence and community membership
  degrees <- rowSums(adj)
  comm <- rep(1, N)
  
  # Grow the graph by adding new nodes one by one
  for (i in (N + 1):(2 * N)) {
    # Determine the community to which the new node belongs
    comm_probs <- (degrees[1:(i-1)] + alpha) ^ beta
    comm_probs[comm_probs == Inf] <- 0
    comm_probs <- comm_probs * (comm == comm[i-1])
    comm_probs <- comm_probs / sum(comm_probs)
    comm[i] <- sample(1:(i-1), 1, prob = comm_probs)
    
    # Determine the neighbors of the new node
    dists_to_i <- dists[i, 1:(i-1)]
    neighbor_probs <- ((1 / (dists_to_i + 1)) ^ alpha) * (comm[1:(i-1)] == comm[i])
    neighbor_probs[i] <- 0 # cannot connect to itself
    neighbor_probs <- neighbor_probs / sum(neighbor_probs)
    
    # Attach the new node to the existing network
    adj[i,] <- rbinom(i-1, 1, neighbor_probs)
    adj[,i] <- adj[i,-i]
    degrees[i] <- sum(adj[i,])
    
    # Rewire edges for small-world effect
    rewire_probs <- rep(0, i-1)
    for (j in 1:(i-1)) {
      if (adj[i,j] == 1 && runif(1) < prewire) {
        rewirable <- (comm[j] == comm[i] | dists[i,j] <= r) & (degrees[j] > 1)
        if (any(rewirable)) {
          rewire_probs[j] <- 1 / sum(rewirable)
        }
      }
    }
    if (sum(rewire_probs) > 0) {
      rewire_to <- sample(1:(i-1), sum(rewire_probs), replace = TRUE, prob = rewire_probs)
      for (j in 1:length(rewire_to)) {
        adj[i,rewire_to[j]] <- 0
        adj[rewire_to[j],i] <- 0
        new_edge_probs <- ((1 / (dists_to_i[rewirable] + 1)) ^ alpha) * (comm[rewirable] == comm[i])
        new_edge_probs <- new_edge_probs / sum(new_edge_probs)
        new
        
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        spatial_scale_free_expander <- function(N, alpha, beta, prewire, r) {
          # Generate initial set of nodes with spatial structure
          coords <- matrix(runif(N * 2), ncol = 2)
          # while (min(dist(coords)) < 0.1) {
          #   coords <- matrix(runif(N * 2), ncol = 2)
          # }
          
          # Initialize adjacency matrix
         # A<- matrix(0, nrow = N, ncol = N)
          # Create initial adjacency matrix based on spatial constraints
          dist_mat <- as.matrix(dist(coords))
          adj <- ifelse(dist_mat <= r, 1, 0)
          
          # Initialize degree sequence and community membership
          degrees <- rowSums(adj)
          comm <- rep(1, N)
          
          # Add nodes to the graph one by one
          for (i in 2:N) {
            probs=prob = degrees^(alpha)*(1/dist_mat^beta)}
            # Assign new node to a community
            comm <- sample(1:N, 1, prob=probs)}
            
            # Connect new node to existing nodes in its community
            comm_nodes <- which(A[comm,] == 1)
            new_edges <- sample(comm_nodes, min(degree, length(comm_nodes)), replace = FALSE)
            A[i, new_edges] <- 1
            A[new_edges, i] <- 1
            
            # Connect new node to nodes in other communities
            other_nodes <- setdiff(1:(i-1), comm_nodes)
            for (j in other_nodes) {
              if (runif(1) < P) {
                if (dist(coords[i,], coords[j,]) <= r) {
                  A[i, j] <- 1
                  A[j, i] <- 1
                }
              }
            }
            
            # Rewire edges with small world probability
            for (j in 1:(i-1)) {
              if (runif(1) < prewire) {
                if (dist(coords[i,], coords[j,]) <= r) {
                  A[i, j] <- 0
                  A[j, i] <- 0
                  new_edge <- sample(setdiff(1:N, c(i, j)), 1)
                  A[i, new_edge] <- 1
                  A[new_edge, i] <- 1
                }
              }
            }
            
            # Update coordinates of new node
            coords[i,] <- runif(2)
            while (min(dist(coords[1:i,])) < 0.1) {
              coords[i,] <- runif(2)
            }
          }
          
          # Return adjacency matrix
          A
        }
        
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Function to generate a spatial scale-free expander graph
    spatial_sf_graph <- function(N, alpha, beta, r, prewire) {
      # Create initial set of nodes distributed spatially with underlying structure
      coords <- matrix(runif(N*2), ncol=2)
      dist_mat <- as.matrix(dist(coords))
      
      # Generate adjacency matrix based on distance matrix and cutoff distance r
      adj_mat <- ifelse(dist_mat <= r, 1, 0)
      
      # Add nodes one at a time, organizing them into non-overlapping communities
      for (i in 2:N) {
        # Calculate probability of attachment to existing nodes
        dist_vec <- dist_mat[1:i-1, i]
        degree_vec <- rowSums(adj_mat[1:i-1, 1:i-1])
        degree_vec[degree_vec == 0] <- 1 # avoid division by zero
        prob_vec <- (dist_vec^(-alpha)) * (degree_vec^beta)
        prob_vec <- prob_vec / sum(prob_vec) # normalize probabilities
        
        # Select community and node to attach to based on probability vector
        comm_vec <- cutree(cluster::hclust(dist(coords[1:i, ])))
        comm_prob_vec <- rep(0, max(comm_vec))
        for (j in 1:max(comm_vec)) {
          comm_prob_vec[j] <- sum(prob_vec[comm_vec == j])
        }
        comm_prob_vec <- comm_prob_vec / sum(comm_prob_vec) # normalize probabilities
        comm <- sample(1:max(comm_vec), size=1, prob=comm_prob_vec)
        nodes_in_comm <- which(comm_vec == comm)
        attach_node <- sample(nodes_in_comm, size=1, prob=prob_vec[nodes_in_comm])
        
        # Add edge between new node and selected node
        adj_mat[i, attach_node] <- 1
        adj_mat[attach_node, i] <- 1
        
        # Rewire edges with probability prewire
        for (j in 1:(i-1)) {
          if (adj_mat[i, j] == 1 && runif(1) < prewire) {
            # Find all nodes within a certain distance and rewire edge to one of them
            rewire_dist_vec <- dist_mat[, i]
            rewire_dist_vec[rewire_dist_vec > r] <- Inf
            rewire_nodes <- which(rewire_dist_vec <= r)
            rewire_nodes <- setdiff(rewire_nodes, c(i, j, which(degree_vec == max(degree_vec))))
            if (length(rewire_nodes) > 0) {
              rewire_node <- sample(rewire_nodes, size=1)
              adj_mat[i, j] <- 0
              adj_mat[j, i] <- 0
              adj_mat[i, rewire_node] <- 1
              adj_mat[rewire_node, i] <- 1
            }
          }
        }
        
        # Update distance matrix
        dist_mat[i, 1:i-1] <- sqrt((coords[i, 1] - coords[1:i-1, 1])^2 + (coords[i, 2] - coords[1:i-1, 2])^2)
        dist_mat[1:i-1, i] <- dist_mat[i, 1:i-1]
      }
      
      # Convert adjacency matrix to igraph object
      graph <- igraph::graph_from_adjacency_matrix(adj_mat, mode="undirected")
      
      
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      library(Matrix)
      
      generate_graph <- function(N, alpha, beta, prewire, r) {
        # Create initial nodes with spatial structure
        coords <- matrix(runif(2*N), ncol=2)
        while (min(dist(coords)) < 0.05) {
          coords <- matrix(runif(2*N), ncol=2)
        }
        
        # Calculate distance matrix
        dists <- as.matrix(dist(coords))
        
        # Create adjacency matrix
        adj <- as.matrix(dists <= r)
        
        # Initialize degree sequence with initial nodes
        deg_seq <- rep(1, N)
        
        # Loop over new nodes and add to graph
        for (i in 2:N) {
          # Calculate attachment probabilities
          probs <- (dists[i,] <= r) * (deg_seq^(-beta)) * (1: N)^alpha
          probs[i] <- 0
          probs <- probs / sum(probs)
          
          # Choose community based on probabilities
          community <- sample.int(N, 1, prob=probs)
          
          # Choose node within community to connect to
          community_nodes <- which(adj[community,])
          community_probs <- (deg_seq[community_nodes]^(-beta))
          community_probs <- community_probs / sum(community_probs)
          node <- sample.int(community_nodes, 1, prob=community_probs)
          
          # Connect new node to chosen node
          adj[i, node] <- 1
          adj[node, i] <- 1
          deg_seq[i] <- deg_seq[i] + 1
          deg_seq[node] <- deg_seq[node] + 1
          
          # Rewire edges with small world probability
          if (runif(1) < prewire) {
            neighbors <- which(adj[i,])
            if (length(neighbors) > 1) {
              new_neighbor <- sample(setdiff(1:N, c(i, neighbors)), 1)
              adj[i, neighbors] <- 0
              adj[neighbors, i] <- 0
              adj[i, new_neighbor] <- 1
              adj[new_neighbor, i] <- 1
            }
          }
        }
        
        # Return graph as a sparse matrix
        graph <- as(adj, "dgCMatrix")
        return(graph)
      }
      
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      generate_graph <- function(N, alpha, beta, r, prewire) {
        # Generate spatial coordinates for nodes
        coords <- matrix(runif(2*N), ncol = 2)
        while (min(dist(coords)) < r/2) {
          coords <- matrix(runif(2*N), ncol = 2)
        }
        
        # Calculate distance matrix
        dist_mat <- as.matrix(dist(coords))
        
        # Create adjacency matrix based on distance matrix and cutoff distance r
        adj_mat <- ifelse(dist_mat <= r,1,0)
        
        # Initialize graph with first node
        graph <- matrix(0, nrow = N, ncol = N)
        graph[1, 1] <- 1
        
        # Initialize degrees and community assignments
        degrees <- rep(0, N)
        communities <- rep(1, N)
        
        # Add remaining nodes to the graph
        for (i in 2:N) {
          # Calculate attachment probabilities
          probs <- rep(0, i-1)
          for (j in 1:(i-1)) {
            if (adj_mat[i, j] == 1) {
              # Calculate probability based on spatial distance and degree
              prob_spatial <- exp(-alpha*dist_mat[i, j]^2)
              prob_degree <- degrees[j]^(-beta)
              probs[j] <- prob_spatial * prob_degree
            }
          }}
          probs <- probs / sum(probs)
          
          # Choose community based on community assignment probabilities
          comm_probs <- table(communities[1:(i-1)]) / (i-1)
          comm <- sample.int(names(comm_probs), 1, prob = comm_probs)
          
          # Choose neighbor node based on attachment probabilities
          neighbor <- sample.int((i-1), 1, prob = probs)
          
          # Connect node to neighbor in community
          graph[i, neighbor] <- 1
          graph[neighbor, i] <- 1
          degrees[i] <- degrees[i] + 1
          degrees[neighbor] <- degrees[neighbor] + 1
          communities[i] <- comm
          
          # Rewire edges with probability prewire
          if (runif(1) < prewire) {
            for (j in 1:N) {
              if (graph[i, j] == 1 && runif(1) < prewire) {
                # Choose new neighbor based on distance and degree
                avail_nodes <- which(adj_mat[i, ] == 1 & degrees > 0)
                avail_probs <- exp(-alpha*dist_mat[i, avail_nodes]^2) * degrees[avail_nodes]^(-beta)
                avail_probs <- avail_probs / sum(avail_probs)
                new_neighbor <- sample.int(avail_nodes, 1, prob = avail_probs)
                
                # Rewire edge
                graph[i, j] <- 0
                graph[j, i] <- 0
                graph[i, new_neighbor] <- 1
                graph[new_neighbor, i] <- 1
                degrees[j] <- degrees[j] - 1
                degrees[new_neighbor] <- degrees[new_neighbor] + 1
              }
            }
          }
        }
        
        # Return adjacency matrix
        return(graph)
      }
      
      length(c(30,30,40))
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      library(igraph)
      
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
      
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Function to generate a spatial scale-free expander graph
      # Input:
      #   - N: number of nodes
      #   - beta: power-law exponent for degree distribution
      #   - alpha: parameter controlling the strength of the degree effect
      #   - r: cutoff distance for spatial constraint
      #   - prewire: rewiring probability for small-world effect
      #   - m: number of edges to attach from a new node
      # Output:
      #   - adjacency matrix of the resulting graph
      grow_spatial_scalefree_graph <- function(N, beta, alpha, r, prewire, m) {
        # Create initial set of nodes
        coords <- matrix(runif(N*2), ncol=2)
        # Adjust coordinates to avoid nodes too close or too far from each other
        for (i in 2:N) {
          while (min(sqrt(rowSums((coords[i,] - coords[1:(i-1),])^2))) < r/2 ||
                 max(sqrt(rowSums((coords[i,] - coords[1:(i-1),])^2))) > 2*r) {
            coords[i,] <- runif(2)
          }
        }
        # Calculate distance matrix
        dist_mat <- as.matrix(dist(coords))
        # Initialize adjacency matrix
        adj_mat <- matrix(0, nrow=N, ncol=N)
        # Initialize degree vector
        degree <- rep(0, N)
        # Add nodes one by one
        for (i in 2:N) {
          # Compute attachment probability vector for node i
          p <- rep(0, i-1)
          for (j in 1:(i-1)) {
            if (dist_mat[i,j] <= r) {
              p[j] <- (degree[j] + alpha) / (sum(degree) + alpha*(i-1)) * (dist_mat[i,j]/r)^(-beta)
            }
          }
          # Sample m neighbors to attach to from node i
          if (sum(p) > 0) {
            neighbors <- sample(1:(i-1), size=min(m, i-1), replace=FALSE, prob=p/sum(p))
            for (j in neighbors) {
              adj_mat[i,j] <- 1
              adj_mat[j,i] <- 1
              degree[i] <- degree[i] + 1
              degree[j] <- degree[j] + 1
            }
          }
          # Rewire edges with probability prewire
          for (j in 1:(i-1)) {
            if (adj_mat[i,j] == 1 && runif(1) < prewire) {
              # Select a random node within distance r to rewire to
              rewirable_nodes <- which(dist_mat[i,] <= r & degree > 0 & seq_along(degree) != i & seq_along(degree) != j)
              if (length(rewirable_nodes) > 0) {
                rewire_to <- sample(rewirable_nodes, size=1)
                adj_mat[i,j] <- 0
                adj_mat[j,i] <- 0
                adj_mat[i,rewire_to] <- 1
                adj_mat[rewire_to,i] <- 1
                degree[i] <- degree[i] - 1
                degree[j] <- degree[j] - 1
                degree[rewire_to] <- degree[rewire_to] + 1
              }
            }
          }
        }
        return(adj_mat)
      }
      
      
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      grow_spatial_scalefree_expander <- function(N, beta, alpha, mu, prewire, m, r) {
        
        # Set up initial network of m nodes with spatial structure
        nodes <- data.frame(id = 1:m, x = runif(m), y = runif(m))
        dist_mat <- as.matrix(dist(nodes[, 2:3]))
        adj_mat <- as.matrix(dist_mat <= r)
        
        # Initialize node degrees and communities
        degrees <- rep(m, m)
        communities <- rep(1:m, each = 1)
        num_communities <- m
        
        # Helper function to compute attachment probabilities
        compute_probabilities <- function(i, j) {
          d_ij <- dist(nodes[i, 2:3], nodes[j, 2:3])
          k_i <- degrees[i]
          k_j <- degrees[j]
          c_i <- communities[i]
          c_j <- communities[j]
          prob_degree <- (k_i + alpha) * (k_j + alpha) / (2 * m * alpha + sum(degrees) + alpha * N)
          prob_spatial <- exp(-beta * d_ij^2 / r^2)
          prob_community <- ifelse(c_i == c_j, 1 - mu, mu / (num_communities - 1))
          prob_degree * prob_spatial * prob_community
        }
        
        # Add remaining N-m nodes to network
        for (i in (m + 1):N) {
          
          # Select m nodes to connect to based on attachment probabilities
          selected <- sample(1:(i-1), size = m, replace = TRUE, prob = compute_probabilities(i, 1:(i-1)))
          
          # Connect new node to selected nodes and update degrees
          for (j in selected) {
            adj_mat[i, j] <- 1
            adj_mat[j, i] <- 1
            degrees[i] <- degrees[i] + 1
            degrees[j] <- degrees[j] + 1
          }
          
          # Assign new node to a community
          if (runif(1) < mu) {
            num_communities <- num_communities + 1
            communities[i] <- num_communities
          } else {
            neighbor_communities <- unique(communities[adj_mat[i,] == 1])
            community_sizes <- sapply(neighbor_communities, function(c) sum(communities == c))
            communities[i] <- neighbor_communities[which.min(community_sizes)]
          }
          
          # Randomly rewire edges with probability prewire
          for (j in 1:(i-1)) {
            if (adj_mat[i,j] == 1 && runif(1) < prewire) {
              neighbors <- which(adj_mat[j,] == 1 & dist_mat[j,i] < r)
              if (length(neighbors) > 0) {
                k <- sample(neighbors, size = 1)
                adj_mat[i,j] <- 0
                adj_mat[j,i] <- 0
                adj_mat[i,k] <- 1
                adj_mat[k,i] <- 1
              }
            }
          }
          
          # Update distance and community matrices
          nodes <- rbind(nodes, data.frame(id = i, x = runif(1), y = runif(1)))
          dist_mat <- cbind(dist_mat, dist(nodes[, 2:3], nodes[i, 2:3]))
          dist_mat <- rbind(dist_mat, c(dist(nodes
                                             
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
spatial_scale_free_graph <- function(N, beta, alpha, r, prewire) {
  # Generate initial nodes with spatial distribution
  coordinates <- matrix(runif(2*N), ncol = 2)
  while (min(dist(coordinates)) < 0.1) {
    coordinates <- matrix(runif(2*N), ncol = 2)
  }
  distances <- dist(coordinates)
  
  # Create adjacency matrix based on distance matrix
  adjacency_matrix <- as.matrix(distances <= r)
  
  # Assign nodes to communities
  community_memberships <- rep(NA, N)
  community_sizes <- numeric()
  for (i in 1:N) {
    if (i == 1) {
      community_memberships[i] <- 1
      community_sizes[1] <- 1
    } else {
      community_probabilities <- rep(0, max(community_memberships))
      for (j in 1:max(community_memberships)) {
        community_size <- sum(community_memberships == j)
        if (community_size > 0) {
          # Calculate attachment probability based on spatial and degree effects
          degree_prob <- (sum(adjacency_matrix[i, community_memberships == j]) + alpha) /
            (sum(degree(adjacency_matrix)[community_memberships == j]) + alpha*max(degree(adjacency_matrix)))
          spatial_prob <- exp(-beta * distances[i, community_memberships == j]^2)
          community_probabilities[j] <- degree_prob * spatial_prob
        }
      }
      # Assign node to a community with highest probability
      community_memberships[i] <- which.max(community_probabilities)
      if (community_memberships[i] > length(community_sizes)) {
        community_sizes[community_memberships[i]] <- 0
      }
      community_sizes[community_memberships[i]] <- community_sizes[community_memberships[i]] + 1
    }
  }
  
  # Rewire edges for small world effect
  for (i in 1:N) {
    for (j in (i+1):N) {
      if (adjacency_matrix[i,j] == 1) {
        if (runif(1) < prewire) {
          distances_to_j <- distances[i,] + distances[j,] - 2*distances[i,j]
          valid_nodes <- which((distances_to_j <= r) & (community_memberships != community_memberships[i]))
          if (length(valid_nodes) > 0) {
            new_j <- sample(valid_nodes, 1)
            adjacency_matrix[i,j] <- 0
            adjacency_matrix[j,i] <- 0
            adjacency_matrix[i,new_j] <- 1
            adjacency_matrix[new_j,i] <- 1
          }
        }
      }
    }
  }
  
  # Convert adjacency matrix to graph object
  g <- graph_from_adjacency_matrix(adjacency_matrix)
  return(g)
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
generate_graph <- function(N, beta, alpha, r, prewire) {
  # Generate initial nodes with spatial structure
  coords <- matrix(0, nrow = N, ncol = 2)
  coords[1, ] <- runif(2)
  for (i in 2:N) {
    repeat {
      candidate <- runif(2)
      if (min(sqrt(rowSums((coords[1:(i-1),] - candidate)^2))) > r) {
        coords[i,] <- candidate
        break
      }
    }
  }
  
  # Compute distance matrix and adjacency matrix
  dist_mat <- as.matrix(dist(coords))
  adj_mat <- as.matrix(dist_mat <= r)
  
  # Initialize node and community degrees
  node_degrees <- rep(0, N)
  comm_degrees <- rep(0, N)
  
  # Create initial community structure
  communities <- rep(1, N)
  comm_sizes <- rep(1, N)
  num_comms <- 1
  
  # Add nodes one by one
  for (i in 2:N) {
    # Compute attachment probabilities
    comm_probs <- rep(0, num_comms)
    for (j in 1:num_comms) {
      comm_nodes <- which(communities == j)
      comm_dists <- dist_mat[comm_nodes, i]
      comm_degrees[j] <- sum(node_degrees[comm_nodes])
      comm_sizes[j] <- length(comm_nodes)
      comm_probs[j] <- sum(comm_dists^(-beta) * (comm_degrees[j] + alpha) / 
                             (comm_sizes[j] * (N-1) + alpha * N))
    }
    comm_probs <- comm_probs / sum(comm_probs)
    
    # Assign node to a community
    comm_id <- sample(1:num_comms, size = 1, prob = comm_probs)
    communities[i] <- comm_id
    comm_nodes <- which(communities == comm_id)
    comm_degrees[comm_id] <- sum(node_degrees[comm_nodes])
    comm_sizes[comm_id] <- length(comm_nodes)
    
    # Compute attachment probabilities within community
    node_probs <- rep(0, length(comm_nodes))
    for (j in 1:length(comm_nodes)) {
      node_dists <- dist_mat[comm_nodes[j], comm_nodes]
      node_degrees[comm_nodes[j]] <- sum(adj_mat[comm_nodes[j],])
      node_probs[j] <- sum(node_dists^(-beta) * (node_degrees[comm_nodes] + alpha) / 
                             (comm_sizes[comm_id] - 1 + alpha))
    }
    node_probs <- node_probs / sum(node_probs)
    
    # Select node to connect to within community
    node_id <- sample(comm_nodes, size = 1, prob = node_probs)
    
    # Connect node to selected node
    adj_mat[i, node_id] <- 1
    adj_mat[node_id, i] <- 1
    
    # Rewire edges with small world probability
    if (prewire > 0 && runif(1) < prewire) {
      # Select random node within a certain distance
      dist_to_node <- dist_mat[i,]
      dist_to_node[i] <- Inf
      rewire_ids <- which(dist_to_node <= r & dist_to_node > 0)
      if (length(rewire_ids) > 0) {
        rewire_id <- sample(rewire_ids, size = 1)
        # Rewire edge
        adj_mat
        
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        generate_spatial_scale_free_expander_graph <- function(N, r, beta, m, alpha, prewire, plot = TRUE) {
          # Generate initial set of nodes
          nodes <- matrix(runif(N*2), ncol=2)
          # Create adjacency matrix based on distance matrix
          dist_mat <- as.matrix(dist(nodes))
          A <- ifelse(dist_mat <= r, 1, 0)
          diag(A) <- 0
          
          # Generate new nodes and attach to existing graph
          for (i in 1:(N-m)) {
            # Create community membership for new node
            comm_size <- sample(1:(m+1), 1, prob=rbeta(m+2, alpha))
            comm_sizes <- colSums(A)
            comm_sizes[length(comm_sizes)+1] <- 0
            comm_sizes[length(comm_sizes)-m+1] <- 0
            comm_probs <- (comm_sizes+1)^beta
            comm_probs[length(comm_probs)-m+1:length(comm_probs)] <- 0
            comm_probs <- comm_probs / sum(comm_probs)
            new_comm <- sample(1:length(comm_probs), 1, prob=comm_probs)
            comm <- rep(0, N)
            comm[which(A[new_comm,]==1)] <- 1:length(which(A[new_comm,]==1))
            comm[new_comm] <- length(which(A[new_comm,]==1)) + 1
            
            # Create attachment probability matrix
            attach_probs <- matrix(0, ncol=N, nrow=m+1)
            for (j in 1:(m+1)) {
              for (k in 1:N) {
                if (comm[k] == j) {
                  d <- sqrt((nodes[i,1]-nodes[k,1])^2 + (nodes[i,2]-nodes[k,2])^2)
                  if (d <= r) {
                    attach_probs[j,k] <- (d/r)^(-alpha)
                  }
                }
              }
            }
            attach_probs[length(attach_probs)-m+1:length(attach_probs),] <- 0
            attach_probs <- attach_probs / sum(attach_probs)
            
            # Attach new node to existing graph
            new_node <- rep(0, N)
            new_node[i] <- 1
            comm_prob <- attach_probs[new_comm,]
            comm_prob[i] <- 0
            comm_prob <- comm_prob / sum(comm_prob)
            attach_prob <- c(rep(0, i-1), comm_prob, rep(0, N-i))
            attach_prob[i] <- 0
            attach_prob <- attach_prob / sum(attach_prob)
            for (j in 1:m) {
              attach_node <- sample(1:N, 1, prob=attach_prob)
              attach_prob[attach_node] <- 0
              new_node[attach_node] <- 1
              A[i, attach_node] <- 1
              A[attach_node, i] <- 
                
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                library(Matrix)
              
              grow_spatial_scalefree_expander <- function(N, beta, alpha, r, m, prewire) {
                
                # Generate initial nodes
                coords <- matrix(runif(N*2), ncol=2)
                D <- as.matrix(dist(coords))
                diag(D) <- Inf
                
                # Create adjacency matrix based on distance matrix
                A <- (D <= r) + 0
                
                # Initialize node and edge lists
                nodes <- 1:N
                edges <- matrix(0, ncol=2)
                k <- rep(0, N)
                
                # Add new nodes to the network
                for (i in 2:N) {
                  
                  # Calculate attachment probability matrix
                  P <- matrix(0, ncol=N, nrow=N)
                  for (j in nodes) {
                    if (j == i) next
                    if (A[i,j] == 1) {
                      P[i,j] <- (k[j]+alpha) * (D[i,j]^(-beta))
                    } else {
                      P[i,j] <- (D[i,j]^(-beta))
                    }
                  }
                  
                  # Normalize probability matrix
                  P <- P / sum(P)
                  
                  # Choose m nodes to attach to
                  attach_nodes <- sample(nodes, size=m, replace=TRUE, prob=P[i,])
                  
                  # Add edges between new node and attach nodes
                  for (j in attach_nodes) {
                    edges <- rbind(edges, c(i, j))
                    k[i] <- k[i] + 1
                    k[j] <- k[j] + 1
                  }
                  
                  # Add new node to nodes list
                  nodes <- c(nodes, i)
                  
                  # Add new row and column to adjacency matrix
                  A <- rbind(A, rep(0, N))
                  A <- cbind(A, rep(0, N+1))
                  
                  # Update adjacency matrix with new edges
                  for (j in attach_nodes) {
                    A[i,j] <- 1
                    A[j,i] <- 1
                  }
                  
                  # Rewire edges with probability prewire
                  for (j in 1:nrow(edges)) {
                    if (runif(1) < prewire) {
                      d <- as.numeric(D[edges[j,1], edges[j,2]])
                      idx <- sample(which(D[edges[j,1],] > d & k > 0))
                      if (length(idx) > 0) {
                        edges[j,2] <- idx[1]
                        A[edges[j,1], edges[j,2]] <- 0
                        A[edges[j,2], edges[j,1]] <- 0
                        A[edges[j,1], idx[1]] <- 1
                        A[idx[1], edges[j,1]] <- 1
                      }
                    }
                  }
                  
                }
                
                # Create sparse adjacency matrix
                A_sparse <- as(Matrix(A, sparse=TRUE), "dgCMatrix")
                
                return(A_sparse)
                
              }
              
                
                
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              # Function to generate spatial scale-free expander graph
              spatial_scale_free_expander <- function(N, r, beta, m, alpha, prewire){
                
                # Generate initial nodes
                coords <- matrix(runif(N*2), ncol=2)
                
                # Define adjacency matrix
                adj <- matrix(0, nrow=N, ncol=N)
                
                # Generate distance matrix
                dist_matrix <- as.matrix(dist(coords))
                
                # Create adjacency matrix based on distance matrix
                adj[dist_matrix <= r] <- 1
                
                # Initialize node and edge lists
                nodes <- 1:N
                edges <- list()
                
                # Initialize community assignment vector
                comm_assign <- rep(0, N)
                
                # Assign initial communities based on spatial proximity
                comm <- 1
                for (i in 1:N){
                  if (comm_assign[i] == 0){
                    comm_assign[i] <- comm
                    neighbors <- which(adj[i,] == 1)
                    dist_to_neighbors <- dist_matrix[i, neighbors]
                    comm_assign[neighbors[dist_to_neighbors <= r]] <- comm
                    comm <- comm + 1
                  }
                }
                
                # Add new nodes to the network
                for (i in 1:(N-m)){
                  
                  # Calculate attachment probabilities for each community
                  degree <- colSums(adj)
                  degree[degree==0] <- 1 # Avoid divide-by-zero error
                  comm_sizes <- table(comm_assign)
                  comm_degrees <- tapply(degree, comm_assign, sum)
                  comm_probs <- ((comm_degrees/degree)^alpha) * ((comm_sizes^(beta-1))/(sum(comm_sizes^(beta-1))))
                  
                  # Choose community for new node to attach to
                  new_comm <- sample(names(comm_probs), size=1, prob=comm_probs)
                  
                  # Choose nodes to attach to in chosen community
                  comm_nodes <- nodes[comm_assign == new_comm]
                  if (length(comm_nodes) > m){
                    attach_nodes <- sample(comm_nodes, size=m, replace=FALSE)
                  } else {
                    attach_nodes <- comm_nodes
                  }
                  
                  # Attach new node to chosen nodes
                  for (j in attach_nodes){
                    if (runif(1) < prewire){
                      # Rewire edge to a random node within a certain distance
                      dist_to_all_nodes <- dist_matrix[i,]
                      valid_nodes <- which(dist_to_all_nodes > r & dist_to_all_nodes < r*2)
                      if (length(valid_nodes) > 0){
                        attach_node <- sample(valid_nodes, size=1)
                        edges[[length(edges)+1]] <- c(i, attach_node)
                      }
                    } else {
                      edges[[length(edges)+1]] <- c(i, j)
                    }
                  }
                  
                  # Assign new node to a community based on spatial proximity
                  neighbors <- which(adj[i,] == 1)
                  comm_counts <- table(comm_assign[neighbors])
                  if (length(comm_counts) > 0){
                    new_comm <- names(comm_counts)[which.max(comm_counts)]
                    comm_assign[i] <- new_comm
                  } else {
                    comm <- comm + 1
                    comm_assign[i] <- comm
                  }
                }
                
                # Convert edge list to adjacency matrix
                adj <- matrix(0, nrow=N, ncol=N)
                for (e in edges){
                  adj[e[1], e[2]] <- 1
                  adj
                  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  # Function to generate spatial scale-free expander graph
                  # Parameters: beta - power law parameter, alpha - strength of degree effect,
                  #             m - number of edges added for each new node, N - number of nodes,
                  #             r - cutoff distance for spatial distance effect,
                  #             prewire - probability of rewiring edges for small-world effect
                  # Returns: adjacency matrix of the generated graph
                  generate_spatial_scale_free_expander <- function(beta, alpha, m, N, r, prewire) {
                    # Initialize adjacency matrix and distance matrix
                    adj_matrix <- matrix(0, nrow = N, ncol = N)
                    dist_matrix <- matrix(0, nrow = N, ncol = N)
                    
                    # Initialize nodes randomly within a unit square with spatial structure
                    nodes <- data.frame(x = runif(N), y = runif(N))
                    
                    # Calculate distance matrix and set adjacency matrix based on distance matrix
                    for (i in 1:(N-1)) {
                      for (j in (i+1):N) {
                        dist <- sqrt((nodes$x[i] - nodes$x[j])^2 + (nodes$y[i] - nodes$y[j])^2)
                        dist_matrix[i,j] <- dist
                        dist_matrix[j,i] <- dist
                        if (dist <= r) {
                          adj_matrix[i,j] <- 1
                          adj_matrix[j,i] <- 1
                        }
                      }
                    }
                    
                    # Initialize community assignment vector randomly
                    comm <- sample(1:ceiling(sqrt(N)), N, replace = TRUE)
                    
                    # Initialize degree vector and total degree
                    deg <- rep(0, N)
                    total_deg <- 0
                    
                    # Add nodes to the graph one by one
                    for (i in 1:(N-m)) {
                      # Calculate attachment probabilities for each potential neighbor
                      prob <- rep(0, N)
                      for (j in 1:N) {
                        if (comm[j] == comm[i]) {
                          prob[j] <- (deg[j] + alpha) / (total_deg + alpha*(N-1))
                        } else {
                          prob[j] <- (1 - beta) * (deg[j] + alpha) / (total_deg + alpha*(N-1)) + beta/N
                        }
                      }
                      
                      # Select m neighbors based on attachment probabilities
                      neighbors <- sample(1:N, m, replace = FALSE, prob = prob)
                      
                      # Add edges to selected neighbors
                      for (j in neighbors) {
                        if (adj_matrix[i,j] == 0) {
                          adj_matrix[i,j] <- 1
                          adj_matrix[j,i] <- 1
                          deg[i] <- deg[i] + 1
                          deg[j] <- deg[j] + 1
                          total_deg <- total_deg + 2
                        }
                      }
                      
                      # Assign new node to a community randomly
                      comm[i+m] <- sample(1:ceiling(sqrt(N)), 1)
                    }
                    
                    # Add small-world effect by rewiring edges with probability prewire
                    for (i in 1:(N-1)) {
                      for (j in (i+1):N) {
                        if (adj_matrix[i,j] == 1) {
                          if (runif(1) < prewire) {
                            # Find potential rewiring nodes within a certain distance
                            potential_nodes <- which(dist_matrix[i,] > r & dist_matrix[i,] < r*10 & comm == comm[j])
                            if
                            
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            library(Matrix)
                            library(Matrix.utils)
                            
                            grow_spatial_scalefree_expander <- function(N, beta, alpha, r, m, prewire) {
                              # Initialize the graph with the first m nodes
                              G <- spMatrix(N, N, 0)
                              pos <- matrix(runif(N*2), N, 2)
                              dmat <- dist(pos)
                              dmat[dmat > r] <- Inf
                              dmat[1:(m+1), 1:(m+1)] <- Inf
                              adjmat <- (dmat <= r) * 1
                              
                              for (i in 1:m) {
                                adjmat[i, (i+1):m] <- 0
                              }
                              
                              G[1:m, 1:m] <- adjmat[1:m, 1:m]
                              
                              # Initialize the community structure
                              comm <- rep(1:m, each=floor(N/m))[1:N]
                              
                              # Grow the graph
                              for (i in (m+1):N) {
                                # Compute the attachment probabilities
                                d <- colSums(G) + 1
                                d[i] <- 0
                                p <- (d^beta) * ((1 + alpha)^m) / ((1 + alpha)^d + (1 + alpha)^m)
                                p[comm == comm[i]] <- 0
                                
                                # Add edges to new node
                                edges <- sample(1:(i-1), size=m, replace=FALSE, prob=p[1:(i-1)])
                                for (j in edges) {
                                  G[i,j] <- 1
                                  G[j,i] <- 1
                                }
                                
                                # Update community structure
                                comm[i] <- sample(1:m, size=1, prob=paste(table(comm), collapse="/"))
                                
                                # Rewire edges for small world effect
                                for (j in edges) {
                                  if (runif(1) < prewire) {
                                    distj <- sqrt((pos[i,1]-pos[j,1])^2 + (pos[i,2]-pos[j,2])^2)
                                    rewires <- which((colSums(G) - G[j,]) > 0 & dmat[i,] < distj)
                                    if (length(rewires) > 0) {
                                      k <- sample(rewires, size=1)
                                      G[j,k] <- 0
                                      G[k,j] <- 0
                                      G[i,k] <- 1
                                      G[k,i] <- 1
                                    }
                                  }
                                }
                                
                                # Update distance matrix and adjacency matrix
                                dmat[i,] <- sqrt((pos[i,1]-pos[,1])^2 + (pos[i,2]-pos[,2])^2)
                                adjmat[i,] <- (dmat[i,] <= r) * 1
                                
                                # Add new node to the graph
                                G[i,] <- adjmat[i,] * 1
                                
                              }
                              
                              return(G)
                            }
                            
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            generate_spatial_scale_free_graph <- function(N, beta, alpha, r, m, prewire) {
                              
                              # Create an initial set of nodes distributed spatially
                              coords <- matrix(runif(N * 2), ncol = 2)
                              while (min(dist(coords)) < 0.1) {
                                coords <- matrix(runif(N * 2), ncol = 2)
                              }
                              
                              # Initialize adjacency matrix and community membership vector
                              adj_mat <- matrix(0, nrow = N, ncol = N)
                              comm_mem <- rep(1, N)
                              
                              # Define the attachment probability matrix
                              P <- matrix(0, nrow = N, ncol = N)
                              for (i in 1:N) {
                                for (j in 1:N) {
                                  if (i != j) {
                                    dist_ij <- sqrt((coords[i,1]-coords[j,1])^2 + (coords[i,2]-coords[j,2])^2)
                                    if (dist_ij <= r) {
                                      P[i,j] <- (dist_ij^(-beta)) * (comm_mem[i]==comm_mem[j]) * (1 + alpha*degree(adj_mat,i))^alpha
                                    }
                                  }
                                }
                              }
                              
                              # Add nodes to the graph
                              for (i in (N+1):(2*N)) {
                                # Determine which community the new node will belong to
                                comm_probs <- colSums(P[,1:i-1]) + rowSums(P[1:i-1,]) + 1
                                comm_probs[1:i-1] <- 0
                                comm_probs <- comm_probs / sum(comm_probs)
                                comm_mem[i-N] <- sample.int(i-1, 1, prob = comm_probs)
                                
                                # Attach edges to the new node
                                for (j in 1:m) {
                                  attach_probs <- P[i-N,1:i-1] + P[1:i-1,i-N]
                                  attach_probs[1:i-1] <- 0
                                  attach_probs <- attach_probs / sum(attach_probs)
                                  attach_node <- sample.int(i-1, 1, prob = attach_probs)
                                  adj_mat[i-N,attach_node] <- 1
                                  adj_mat[attach_node,i-N] <- 1
                                }
                              }
                              
                              # Rewire edges for small-world effect
                              for (i in 1:N) {
                                for (j in (i+1):N) {
                                  if (adj_mat[i,j] == 1 && runif(1) < prewire) {
                                    # Find a random node within a certain distance to rewire the edge to
                                    dist_i <- sqrt((coords[i,1]-coords[,1])^2 + (coords[i,2]-coords[,2])^2)
                                    rewire_nodes <- which(dist_i <= 2*r)
                                    rewire_nodes <- rewire_nodes[rewire_nodes != i & adj_mat[rewire_nodes,j] == 0]
                                    if (length(rewire_nodes) > 0) {
                                      rewire_node <- sample(rewire_nodes, 1)
                                      adj_mat[i,j] <- 0
                                      adj_mat[j,i] <- 0
                                      adj_mat[i,rewire_node] <- 1
                                      adj_mat[rewire_node,i] <- 1
                                    }
                                  }
                                }
                              }
                              
                              # Return the adjacency matrix
                              return(adj_mat)
                            }
                            
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            #+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                            library(Matrix)
                            library(pracma)
                            
                            # Function to create a spatial scale-free expander graph
                            spatial_scale_free_expander <- function(N, beta, alpha, r, m, prewire) {
                              
                              # Generate spatial coordinates
                              coords <- matrix(runif(N*2, 0, 1), ncol = 2)
                              # Ensure no two nodes are too close to each other
                              while(min(dist(coords)) < r) {
                                coords <- matrix(runif(N*2, 0, 1), ncol = 2)
                              }
                              
                              # Create distance matrix
                              dist_mat <- dist(coords)
                              
                              # Initialize adjacency matrix
                              A <- matrix(0, nrow = N, ncol = N)
                              
                              # Add nodes to the network one by one
                              for (i in 1:N) {
                                # Compute attachment probabilities
                                p <- rep(0, N)
                                for (j in 1:(i-1)) {
                                  if (dist_mat[i, j] <= r) {
                                    p[j] <- (degree(A, j) + alpha) ^ beta
                                  }
                                }
                                # Assign nodes to communities
                                comm_size <- ceiling(N / 10) # 10 communities
                                comm_id <- ceiling(i / comm_size)
                                comm_nodes <- ((comm_id - 1) * comm_size + 1):(comm_id * comm_size)
                                # Generate attachment probability matrix
                                P <- matrix(0, nrow = comm_size, ncol = comm_size)
                                for (j in 1:comm_size) {
                                  for (k in 1:comm_size) {
                                    if (j != k) {
                                      if (comm_nodes[j] < i && comm_nodes[k] < i) {
                                        if (A[comm_nodes[j], comm_nodes[k]] == 0) {
                                          P[j, k] <- p[comm_nodes[j]] * p[comm_nodes[k]]
                                        }
                                      }
                                    }
                                  }
                                }
                                # Normalize attachment probability matrix
                                P <- P / sum(P)
                                # Add m edges to the new node
                                for (j in 1:m) {
                                  # Choose a node to connect to based on attachment probability
                                  neighbor <- sample(comm_nodes, 1, prob = P[i%%comm_size + 1, ])
                                  # Add edge
                                  A[i, neighbor] <- 1
                                  A[neighbor, i] <- 1
                                }
                              }
                              
                              # Add small-world effect
                              for (i in 1:N) {
                                for (j in (i+1):N) {
                                  if (A[i, j] == 1) {
                                    if (runif(1) < prewire) {
                                      # Choose a random node within a certain distance
                                      dist_ij <- dist_mat[i, j]
                                      candidates <- which(dist_mat[i, ] + dist_mat[j, ] <= 2*r)
                                      candidates <- candidates[candidates != i & candidates != j]
                                      if (length(candidates) > 0) {
                                        new_j <- sample(candidates, 1)
                                        # Rewire edge
                                        A[i, j] <- 0
                                        A[j, i] <- 0
                                        A[i, new_j] <- 1
                                        A[new_j, i] <- 1
                                      }
                                    }
                                  }
                                }
                              }
                              
                              # Return adjacency matrix
                              return(A)
                            }
                            
                            
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            #+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            spatial_sf_expander <- function(N, r, beta, alpha, m, prewire){
                              # Generate initial set of nodes on a unit square
                              pos <- matrix(runif(N*2, 0, 1), ncol=2)
                              
                              # Create an empty adjacency matrix
                              adj <- matrix(0, nrow=N, ncol=N)
                              
                              # Generate the distance matrix
                              dist_mat <- as.matrix(dist(pos))
                              
                              # Create the initial adjacency matrix based on distance
                              adj[dist_mat <= r] <- 1
                              
                              # Generate the stochastic block matrix
                              comm_sizes <- sample(x=1:(N/10), size=10, replace=TRUE)
                              comm <- rep(1:10, comm_sizes)
                              P <- matrix(rep(1, 100), ncol=10)
                              for(i in 1:10){
                                for(j in 1:10){
                                  if(i != j){
                                    P[i, j] <- runif(1, 0, 1)
                                    P[j, i] <- P[i, j]
                                  }
                                }
                              }
                              
                              # Grow the network
                              for(i in 1:(N-m)){
                                # Calculate the degree of each node
                                degree <- rowSums(adj)
                                
                                # Calculate the attachment probability matrix
                                P_att <- P
                                for(j in 1:N){
                                  if(degree[j] == 0){
                                    P_att[j,] <- rep(1/N, N)
                                  } else {
                                    d <- dist_mat[j,]
                                    d[degree == 0] <- max(d) + 1
                                    p <- (d^(-alpha)) * (degree^beta)
                                    p[degree == 0] <- 0
                                    p <- p/sum(p)
                                    P_att[j,] <- p
                                  }
                                }
                                
                                # Select the nodes to connect to
                                new_edges <- sample(1:N, m, replace=FALSE, prob=P_att[i,])
                                
                                # Connect the new nodes to the network
                                adj[i, new_edges] <- 1
                                adj[new_edges, i] <- 1
                                
                                # Rewire the edges with probability prewire
                                for(j in 1:m){
                                  if(runif(1) < prewire){
                                    k <- sample(1:N, 1)
                                    if(k != i & adj[i,k] == 0){
                                      adj[i, new_edges[j]] <- 0
                                      adj[new_edges[j], i] <- 0
                                      adj[i, k] <- 1
                                      adj[k, i] <- 1
                                      new_edges[j] <- k
                                    }
                                  }
                                }
                              }
                              
                              # Return the adjacency matrix and the positions of the nodes
                              list(adj=adj, pos=pos)
                            }
                            

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            spatial_scalefree_expander <- function(N, beta, alpha, r, m, prewire) {
                              # Create initial set of nodes randomly distributed on a unit square
                              coords <- matrix(runif(N*2), ncol = 2)
                              while (min(dist(coords)) < 0.05) {
                                coords <- matrix(runif(N*2), ncol = 2)
                              }
                              
                              # Initialize adjacency matrix
                              adj_matrix <- matrix(0, nrow = N, ncol = N)
                              
                              # Define attachment probability matrix
                              P <- matrix(0, nrow = N, ncol = N)
                              for (i in 1:N) {
                                for (j in 1:N) {
                                  if (i == j) {
                                    P[i, j] <- 0
                                  } else {
                                    dist_ij <- sqrt(sum((coords[i,] - coords[j,])^2))
                                    if (dist_ij <= r) {
                                      P[i, j] <- (dist_ij^alpha)/(1 + dist_ij^alpha)
                                    } else {
                                      P[i, j] <- 0
                                    }
                                  }
                                }
                              }
                              
                              # Define stochastic block matrix
                              k <- m - 1
                              B <- matrix(rep(0, N^2), nrow = N)
                              for (i in 1:N) {
                                for (j in 1:N) {
                                  if (i == j) {
                                    B[i,j] <- m
                                  } else {
                                    if (runif(1) < P[i,j]) {
                                      if (runif(1) < (k/N)^(1-beta)) {
                                        B[i,j] <- 1
                                        k <- k - 1
                                      }
                                    }
                                  }
                                }
                              }
                              
                              # Add nodes to the graph
                              for (i in (m+1):N) {
                                # Initialize adjacency row for new node
                                new_row <- rep(0, i-1)
                                
                                # Assign node to a community
                                node_comm <- sample(1:(i-1), 1, prob = colSums(B[1:(i-1),]))
                                
                                # Connect new node to existing nodes within the same community
                                comm_nodes <- which(B[node_comm,] > 0)
                                if (length(comm_nodes) == 0) {
                                  comm_nodes <- node_comm
                                }
                                comm_degrees <- colSums(adj_matrix[comm_nodes,])
                                comm_probs <- comm_degrees^beta
                                comm_probs[node_comm] <- 0
                                comm_probs <- comm_probs/sum(comm_probs)
                                for (j in 1:m) {
                                  new_neighbor <- sample(comm_nodes, 1, prob = comm_probs)
                                  new_row[new_neighbor] <- 1
                                  adj_matrix[new_neighbor, i] <- 1
                                }
                                
                                # Connect new node to nodes in other communities
                                other_nodes <- setdiff(1:(i-1), comm_nodes)
                                other_degrees <- colSums(adj_matrix[other_nodes,])
                                other_probs <- other_degrees^beta
                                other_probs <- other_probs/sum(other_probs)
                                for (j in 1:m) {
                                  new_neighbor <- sample(other_nodes, 1, prob = other_probs)
                                  new_row[new_neighbor] <- 1
                                  adj_matrix[new_neighbor, i] <- 1
                                }
                                
                                # Add rewiring probability
                                for (j in 1:(i-1)) {
                                  
                                  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                  spatial_scalefree_expander_graph <- function(N, r, m, beta, alpha, prewire) {
                                    
                                    # Function to calculate the distance between two nodes
                                    dist <- function(x1, y1, x2, y2) {
                                      sqrt((x2 - x1)^2 + (y2 - y1)^2)
                                    }
                                    
                                    # Generate initial set of nodes on a unit square with a minimum distance of 'r' between them
                                    nodes <- data.frame(x = runif(1), y = runif(1))
                                    while (nrow(nodes) < N) {
                                      x <- runif(1)
                                      y <- runif(1)
                                      if (all(dist(x, y, nodes$x, nodes$y) > r)) {
                                        nodes <- rbind(nodes, data.frame(x = x, y = y))
                                      }
                                    }
                                    
                                    # Create distance matrix
                                    dist_mat <- as.matrix(dist(nodes))
                                    
                                    # Create adjacency matrix based on distance matrix with a cutoff distance of 'r'
                                    adj_mat <- ifelse(dist_mat <= r, 1, 0)
                                    diag(adj_mat) <- 0
                                    
                                    # Initialize the degree vector and community assignment vector
                                    degree <- rep(0, N)
                                    community <- rep(0, N)
                                    
                                    # Assign each node to a community based on a stochastic block model with attachment probability matrix 'P_ij'
                                    # 'P_ij' should favor short spatial distances and have a scale-free degree distribution with power law parameter 'beta'
                                    for (i in 1:N) {
                                      if (i == 1) {
                                        community[i] <- 1
                                      } else {
                                        # Calculate attachment probability for each community
                                        p <- rep(0, max(community))
                                        for (j in 1:max(community)) {
                                          k <- sum(degree[community == j]) # degree sum of nodes in community j
                                          d <- mean(dist_mat[i, community == j]) # average distance between node i and nodes in community j
                                          p[j] <- (k + alpha) * (d^-beta)
                                        }
                                        p <- p / sum(p) # normalize probabilities
                                        
                                        # Assign node i to a community based on attachment probabilities
                                        community[i] <- sample(1:max(community), size = 1, prob = p)
                                      }
                                      
                                      # Add m edges for node i, based on preferential attachment to nodes in the same community
                                      for (j in 1:m) {
                                        # Calculate attachment probability for each node in the same community as node i
                                        p <- rep(0, sum(community == community[i]))
                                        for (k in which(community == community[i])) {
                                          p[k] <- (degree[k] + alpha) / sum(degree[community == community[i]]) # preferential attachment
                                        }
                                        p[i] <- 0 # exclude self-attachment
                                        p <- p / sum(p) # normalize probabilities
                                        
                                        # Select a node to attach to based on attachment probabilities
                                        attach_to <- sample(which(community == community[i]), size = 1, prob = p)
                                        
                                        # Add edge between nodes i and attach_to
                                        adj_mat[i, attach_to] <- 1
                                        adj_mat[attach_to, i] <- 1
                                        
                                        # Update degree vector
                                        degree[i] <- degree[i] + 1
                                        degree[attach_to] <- degree[attach_to] + 1
                                      }
                                    }
                                    
                                    # Rewire edges with probability 'prewire' for
                                    
                            
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                    # function to generate spatial scale-free expander graph
                                    # beta: power law parameter
                                    # N: number of nodes
                                    # alpha: degree strength parameter
                                    # prewire: rewiring probability for small world effect
                                    # r: distance cutoff
                                    # m: number of edges added for each new node
                                    # P: attachment probability matrix
                                    # community_sizes: sizes of communities
                                    # community_prob: probability of within-community connections
                                    # min_distance: minimum distance between nodes
                                    # max_distance: maximum distance between nodes
                                    # seed: random seed
                                    
                                    spatial_scale_free_graph <- function(beta, N, alpha, prewire, r, m, 
                                                                         P, community_sizes, community_prob, 
                                                                         min_distance = 0.1, max_distance = 0.3,
                                                                         seed = 123) {
                                      set.seed(seed)
                                      
                                      # generate initial set of nodes with spatial structure
                                      coords <- matrix(0, nrow = N, ncol = 2)
                                      coords[1, ] <- runif(2, 0, 1)
                                      for (i in 2:N) {
                                        while (TRUE) {
                                          new_coord <- runif(2, 0, 1)
                                          distances <- apply(coords[1:(i-1), ], 1, function(x) sqrt(sum((x-new_coord)^2)))
                                          if (all(distances > min_distance) && all(distances < max_distance)) {
                                            coords[i, ] <- new_coord
                                            break
                                          }
                                        }
                                      }
                                      
                                      # create distance matrix and adjacency matrix
                                      dist_mat <- as.matrix(dist(coords))
                                      adj_mat <- as.matrix(dist_mat <= r)
                                      
                                      # initialize block matrix
                                      num_blocks <- length(community_sizes)
                                      block_sizes <- rep(community_sizes, each = community_sizes)
                                      block_prob <- matrix(0, nrow = num_blocks, ncol = num_blocks)
                                      for (i in 1:num_blocks) {
                                        for (j in 1:num_blocks) {
                                          if (i == j) {
                                            block_prob[i, j] <- community_prob
                                          } else {
                                            block_prob[i, j] <- (1 - community_prob) / (num_blocks - 1)
                                          }
                                        }
                                      }
                                      
                                      # add nodes one by one
                                      for (i in 2:N) {
                                        # calculate attachment probabilities
                                        degrees <- colSums(adj_mat)
                                        distances <- dist_mat[i, ]
                                        probs <- P * (distances ^ (-beta)) * (degrees ^ alpha)
                                        probs[i] <- 0
                                        probs <- probs / sum(probs)
                                        
                                        # assign node to a community
                                        block_probs <- rep(0, num_blocks)
                                        for (j in 1:num_blocks) {
                                          block_probs[j] <- sum(probs[sum(community_sizes[1:(j-1)]) + 1:sum(community_sizes[1:j])])
                                        }
                                        block_probs <- block_probs / sum(block_probs)
                                        block_idx <- sample(1:num_blocks, size = 1, prob = block_probs)
                                        start_idx <- sum(community_sizes[1:(block_idx-1)]) + 1
                                        end_idx <- sum(community_sizes[1:block_idx])
                                        
                                        # add edges to existing nodes
                                        for (j in 1:m) {
                                          neighbor <- sample(1:(i-1), size = 1, prob = probs)
                                          adj_mat[i, neighbor] <- 1
                                          adj_mat[neighbor, i] <- 1
                                        }
                                        
                                        # add
                                        
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                        grow_spatial_scale_free_expander <- function(N, beta, alpha, r, m, prewire) {
                                          
                                          # Set up initial graph with 2 nodes randomly placed within unit square
                                          x <- runif(2, 0, 1)
                                          y <- runif(2, 0, 1)
                                          adjacency <- matrix(0, ncol = 2, nrow = 2)
                                          
                                          # Set up stochastic block matrix for community structure
                                          c <- sample(1:2, 2, replace = TRUE)
                                          P <- matrix(c(alpha, beta, beta, alpha), ncol = 2, nrow = 2)
                                          
                                          # Grow graph by adding N-2 nodes
                                          for (i in 3:N) {
                                            
                                            # Initialize node coordinates and community assignment
                                            xy <- cbind(runif(1, 0, 1), runif(1, 0, 1))
                                            comm <- sample(1:2, 1, prob = P[c, ])
                                            c[i] <- comm
                                            
                                            # Compute distance matrix
                                            dist <- as.matrix(dist(rbind(xy, cbind(x, y))))
                                            
                                            # Compute adjacency matrix based on distance matrix and r
                                            adj <- ifelse(dist[1, -1] <= r, 1, 0)
                                            
                                            # Apply stochastic block matrix to adjacency matrix
                                            prob <- P[comm, c[1:i-1]] * adj
                                            prob[lower.tri(prob)] <- 0
                                            prob <- prob / sum(prob)
                                            new_edges <- sample(1:(i-1), m, replace = FALSE, prob = prob)
                                            adjacency[i, new_edges] <- 1
                                            adjacency[new_edges, i] <- 1
                                            
                                            # Rewire edges with probability prewire
                                            for (j in new_edges) {
                                              if (runif(1) < prewire) {
                                                rewired <- sample(setdiff(1:N, c(j, which(adjacency[j, ] == 1))), 1)
                                                adjacency[i, j] <- 0
                                                adjacency[j, i] <- 0
                                                adjacency[i, rewired] <- 1
                                                adjacency[rewired, i] <- 1
                                              }
                                            }
                                            
                                            # Update coordinates and community assignment
                                            x[i] <- xy[1]
                                            y[i] <- xy[2]
                                            c <- cbind(c, comm)
                                            
                                          }
                                          
                                          # Convert adjacency matrix to graph object
                                          g <- graph.adjacency(adjacency, mode = "undirected")
                                          
                                          # Assign coordinates and community attributes to nodes
                                          V(g)$x <- x
                                          V(g)$y <- y
                                          V(g)$comm <- c[, N+1]
                                          
                                          return(g)
                                          
                                        }
                                        
                                                        
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                        generate_graph <- function(N, beta, alpha, r, m, prewire) {
                                          # Generate initial set of nodes with spatial structure
                                          coords <- matrix(runif(N*2, 0, 1), ncol = 2)
                                          while (TRUE) {
                                            dist_mat <- as.matrix(dist(coords))
                                            if (all(dist_mat[row(dist_mat) > col(dist_mat)] > 0.1 / sqrt(N))) {
                                              break
                                            }
                                            coords <- matrix(runif(N*2, 0, 1), ncol = 2)
                                          }
                                          
                                          # Initialize adjacency matrix and community assignment matrix
                                          A <- matrix(0, ncol = N, nrow = N)
                                          C <- rep(1, N)
                                          
                                          # Loop over nodes and add to the graph
                                          for (i in seq_len(N)[-1]) {
                                            # Compute attachment probabilities based on community and spatial distance
                                            p <- matrix(0, ncol = i-1, nrow = max(C))
                                            for (c in seq_len(max(C))) {
                                              within_c <- which(C[seq_len(i-1)] == c)
                                              between_c <- which(C[seq_len(i-1)] != c)
                                              within_dist <- dist_mat[i, within_c]
                                              between_dist <- dist_mat[i, between_c]
                                              p[c, within_c] <- (within_dist <= r) * (within_dist ^ beta) * (1 + sum(A[i, within_c]) / m) ^ alpha
                                              p[c, between_c] <- (between_dist <= r) * (between_dist ^ beta) * (1 + sum(A[i, between_c]) / m) ^ alpha
                                            }
                                            
                                            # Assign node to a community based on the highest attachment probability
                                            C[i] <- which.max(rowSums(p))
                                            
                                            # Perform preferential attachment within the assigned community
                                            within_c <- which(C[seq_len(i-1)] == C[i])
                                            A[i, sample(within_c, m, replace = TRUE)] <- 1
                                            
                                            # Rewire edges with small world effect
                                            for (j in which(A[i, seq_len(i-1)] == 1)) {
                                              if (runif(1) < prewire) {
                                                d <- dist(coords[c(i, j),])
                                                candidates <- which(dist_mat[i,] > r & dist_mat[j,] > r & dist_mat[i,] + dist_mat[j,] - 2*d <= 2*r)
                                                A[i,j] <- 0
                                                A[j,i] <- 0
                                                A[i,sample(candidates, 1)] <- 1
                                                A[sample(candidates, 1),j] <- 1
                                              }
                                            }
                                          }
                                          
                                          # Return adjacency matrix
                                          return(A)
                                        }
                                        
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            
generate_spatial_scalefree_graph <- function(N, beta, alpha, mu, pwithin, prewire, m, r, weighted = FALSE, node_attrs = FALSE){
  
  # Define initial graph with two nodes
  nodes <- data.frame(id = 1:2, x = runif(2), y = runif(2))
  distances <- dist(nodes[, 2:3])
  adjacency <- matrix(0, nrow = 2, ncol = 2)
  if(distances[2] <= r){
    adjacency[1, 2] <- 1
    adjacency[2, 1] <- 1
  }
  degrees <- c(1, 1)
  
  # Generate remaining nodes and edges
  for(i in 3:N){
    # Compute probabilities for attaching to existing nodes
    dists <- sqrt((nodes$x - nodes[i, "x"])^2 + (nodes$y - nodes[i, "y"])^2)
    dists[dists == 0] <- NA
    comm_probs <- pwithin^((nodes$comm == nodes$comm[i]) * alpha) * (1 - pwithin)^((nodes$comm != nodes$comm[i]) * alpha)
    dist_probs <- dnorm(dists, mean = 0, sd = r) # Gaussian kernel for distance
    degree_probs <- (degrees^(-beta)) * (sum(degrees)^beta) # Power-law degree distribution
    probs <- comm_probs * dist_probs * degree_probs
    probs[i] <- 0 # No self-connections
    
    # Select nodes to attach to
    attach_nodes <- sample(1:(i-1), size = m, replace = TRUE, prob = probs[1:(i-1)])
    if(length(attach_nodes) > 0){
      edges <- cbind(rep(i, m), attach_nodes)
      if(weighted){
        weights <- runif(m)
        edges <- cbind(edges, weights)
      }
      adjacency[edges] <- 1
      degrees[attach_nodes] <- degrees[attach_nodes] + 1
      degrees[i] <- m
    }
    
    # Assign node to community based on existing communities
    if(i %% 2 == 0){
      nodes[i, "comm"] <- 1
    } else{
      nodes[i, "comm"] <- 2
    }
    
    # Ensure nodes are not too close to each other
    while(TRUE){
      new_node <- data.frame(id = i, x = runif(1), y = runif(1))
      min_dists <- sqrt((nodes$x - new_node$x)^2 + (nodes$y - new_node$y)^2)
      if(all(min_dists > 0.05) && all(min_dists < 0.2)){
        break
      }
    }
    nodes <- rbind(nodes, new_node)
    
    # Add small world rewiring
    for(j in 1:m){
      if(runif(1) < prewire){
        k <- sample(1:(i-1), size = 1)
        if(k != j && adjacency[j, k] == 0 && dist(nodes[j, 2:3], nodes[k, 2:3]) <= r){
          adjacency[j, k] <- 1
          adjacency[k, j] <- 1
        }
      }
    }
  }
  
  # Create igraph object
  g <- graph_from_adjacency_matrix(adjacency, mode = "undirected")
  
  #+++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++
  generate_spatial_scale_free_expander <- function(N, beta, alpha, pwithin, prewire, r, m, node_attrs = FALSE, edge_weights = FALSE) {
    
    # Generate initial nodes with spatial structure
    nodes <- matrix(runif(N*2), ncol = 2)
    while (TRUE) {
      quadtree <- spatialwidget::quadtree(nodes)
      too_close <- apply(quadtree$within_radius(r/2), 1, function(x) length(x) > 1)
      if (sum(too_close) == 0) break
      nodes[too_close, ] <- matrix(runif(sum(too_close)*2), ncol = 2)
    }
    
    # Initialize adjacency matrix
    adj <- matrix(0, ncol = N, nrow = N)
    
    # Initialize node degrees and community assignments
    degrees <- rep(0, N)
    communities <- rep(1, N)
    num_communities <- 1
    
    # Grow the graph
    for (i in 2:N) {
      # Compute probability of connecting to each existing node
      prob <- rep(0, i-1)
      for (j in 1:(i-1)) {
        distance <- sqrt(sum((nodes[i, ] - nodes[j, ])^2))
        if (communities[i] == communities[j]) {
          prob[j] <- pwithin * (1 + alpha * degrees[j]) / (1 + alpha * sum(degrees[communities == communities[j]]))
        } else {
          prob[j] <- (1 - pwithin) * (1 + alpha * degrees[j]) / (1 + alpha * sum(degrees))
        }
        prob[j] <- prob[j] * (distance <= r) * (1 + degrees[j])^(-beta)
      }
      # Connect to m existing nodes with highest probability
      new_edges <- numeric(m)
      for (j in 1:m) {
        if (runif(1) < prewire) {
          # Rewire to a random node within a certain distance
          indices <- which(sqrt(rowSums((nodes - nodes[i, ])^2)) <= r & (1:N) != i)
          if (length(indices) > 0) {
            new_edges[j] <- sample(indices, 1)
          }
        } else {
          # Connect to existing node with highest probability
          if (sum(prob) == 0) break
          new_edges[j] <- sample((1:(i-1))[prob == max(prob)], 1)
        }
      }
      for (j in new_edges[new_edges > 0]) {
        adj[i, j] <- 1
        adj[j, i] <- 1
        degrees[i] <- degrees[i] + 1
        degrees[j] <- degrees[j] + 1
        # Check if nodes belong to same community
        if (communities[i] == communities[j]) {
          p <- pwithin
        } else {
          p <- 1 - pwithin
        }
        if (runif(1) < p) {
          communities[i] <- communities[j]
        } else {
          num_communities <- num_communities + 1
          communities[i] <- num_communities
        }
      }
    }
    
    # Convert adjacency matrix to igraph object
    g <- graph_from_adjacency_matrix(adj)
    
    # Add node attributes if requested
    if (node_attrs) {
      V(g)$x <- nodes[,
                      
                      #+++++++++++++++++++++++++++++++++++++++
                      grow_spatial_sf_expander <- function(N, beta, mu, prewire, m, r, p_within, p_between, node_attr = NULL, edge_weight = NULL) {
                        # Generate initial set of nodes with spatial structure
                        coords <- matrix(runif(N * 2), ncol = 2)
                        coords <- coords / max(coords) # Rescale to unit square
                        while (TRUE) {
                          dist_mat <- as.matrix(dist(coords))
                          if (min(dist_mat[dist_mat > 0]) >= 0.05) { # Ensure nodes are not too close
                            break
                          }
                          coords <- matrix(runif(N * 2), ncol = 2)
                          coords <- coords / max(coords)
                        }
                        
                        # Generate adjacency matrix based on spatial distance
                        adj_mat <- matrix(0, nrow = N, ncol = N)
                        for (i in 1:(N-1)) {
                          for (j in (i+1):N) {
                            if (dist_mat[i, j] <= r) {
                              adj_mat[i, j] <- 1
                              adj_mat[j, i] <- 1
                            }
                          }
                        }
                        
                        # Initialize vector of degrees
                        degree <- rep(0, N)
                        
                        # Loop to add nodes and edges
                        for (i in 2:m) {
                          # Choose community based on preferential attachment with beta parameter
                          community_probs <- (degree + mu) ^ beta
                          community_probs <- community_probs / sum(community_probs)
                          community <- sample(1:N, 1, prob = community_probs)
                          
                          # Choose location based on spatial proximity
                          if (i <= N) {
                            new_coord <- coords[i, ]
                          } else {
                            while (TRUE) {
                              new_coord <- runif(2)
                              dist_to_nodes <- dist(rbind(coords, new_coord))
                              if (min(dist_to_nodes[-(1:i)]) > 0.05) { # Ensure node is not too close to existing nodes
                                break
                              }
                            }
                            coords <- rbind(coords, new_coord)
                            adj_mat <- rbind(adj_mat, rep(0, i-1))
                            adj_mat <- cbind(adj_mat, rep(0, i))
                            degree <- c(degree, 0)
                            N <- i
                          }
                          
                          # Choose nodes to connect to based on community affiliation and distance
                          candidate_nodes <- which(adj_mat[community, ] == 0 & degree < N-1)
                          if (length(candidate_nodes) == 0) {
                            candidate_nodes <- which(degree < N-1)
                          }
                          dist_to_candidate_nodes <- dist_mat[community, candidate_nodes]
                          connection_probs <- p_within * (degree[candidate_nodes] + mu) ^ beta / (dist_to_candidate_nodes + 1) ^ p_between
                          connection_probs <- connection_probs / sum(connection_probs)
                          new_edge_to <- sample(candidate_nodes, 1, prob = connection_probs)
                          adj_mat[community, new_edge_to] <- 1
                          adj_mat[new_edge_to, community] <- 1
                          degree[community] <- degree[community] + 1
                          degree[new_edge_to] <- degree[new_edge_to] + 1
                          
                          # Rewire edges with probability prewire
                          if (prewire > 0) {
                            for (j in 1:N) {
                              if (adj_mat[community, j] == 0 && j != new_edge_to) {
                                if (
                                  #++++++++++++++++++++++++++++++++++++++
                                  library(Matrix)
                                  library(igraph)
                                  
                                  # function to generate a unit square with n nodes, distributed spatially with an underlying structure
                                  generate_nodes <- function(n, seed = 123){
                                    set.seed(seed)
                                    nodes <- matrix(runif(n*2), ncol = 2)
                                    while(length(unique(round(nodes, 4))) < n){
                                      nodes <- matrix(runif(n*2), ncol = 2)
                                    }
                                    return(nodes)
                                  }
                                  
                                  # function to create an adjacency matrix based on the distance matrix and cutoff distance r
                                  create_adjacency_matrix <- function(nodes, r){
                                    d <- as.matrix(dist(nodes))
                                    adj_matrix <- ifelse(d <= r, 1, 0)
                                    return(adj_matrix)
                                  }
                                  
                                  # function to grow a spatial scale-free expander graph with parameters beta, alpha, mu, prewire, m, p_within, p_between
                                  grow_graph <- function(N, beta, alpha, mu, prewire, m, p_within, p_between, r, seed = 123){
                                    
                                    # initialize nodes
                                    nodes <- generate_nodes(N, seed)
                                    num_nodes <- N
                                    
                                    # create initial adjacency matrix based on distance
                                    adj_matrix <- create_adjacency_matrix(nodes, r)
                                    
                                    # initialize community assignments
                                    community <- rep(1, N)
                                    
                                    # initialize degree vector
                                    degrees <- rep(0, N)
                                    
                                    # initialize Laplacian matrix
                                    L <- diag(N)
                                    
                                    # grow graph
                                    for(i in 2:N){
                                      # calculate probability of connecting to each existing node
                                      within_comm <- which(community == community[i])
                                      between_comm <- which(community != community[i])
                                      dists_within <- dist(nodes[i,], nodes[within_comm,])
                                      prob_within <- exp(-alpha*dists_within)/degrees[within_comm]^beta
                                      prob_within[i] <- 0
                                      prob_within <- prob_within/sum(prob_within)
                                      dists_between <- dist(nodes[i,], nodes[between_comm,])
                                      prob_between <- exp(-mu*dists_between)/degrees[between_comm]^beta
                                      prob_between <- prob_between/sum(prob_between)
                                      
                                      # choose m nodes to connect to
                                      neighbors <- c(sample(within_comm, size = m, replace = TRUE, prob = p_within*prob_within),
                                                     sample(between_comm, size = m, replace = TRUE, prob = p_between*prob_between))
                                      
                                      # update adjacency matrix and degree vector
                                      adj_matrix[i, neighbors] <- 1
                                      adj_matrix[neighbors, i] <- 1
                                      degrees[i] <- m*2
                                      degrees[neighbors] <- degrees[neighbors] + 1
                                      
                                      # update community assignments
                                      community[i] <- ifelse(runif(1) < p_within, community[i], sample(unique(community[neighbors]), size = 1))
                                      
                                      # update Laplacian matrix
                                      L[i, i] <- sum(adj_matrix[i,])
                                      L[neighbors, i] <- -1
                                      L[i, neighbors] <- -1
                                      
                                      # rewire edges with probability prewire
                                      for(j in neighbors){
                                        if(runif(1) < prewire){
                                          dists <- dist(nodes[i,], nodes[-c(i,j),])
                                          candidates <- which(dists <= r)
                                          candidates <- candidates[!neighbors %in% candidates]
                                          if(length(candidates) > 0){
                                            new_neighbor <- sample(candidates, size = 1)
                                            adj_matrix[i, j] <- 0
                                            adj
                                            #+++++++++++++++++++++++++++++++++++++++
                                            library(Matrix)
                                            library(igraph)
                                            
                                            # function to generate a unit square with n nodes, distributed spatially with an underlying structure
                                            generate_nodes <- function(n, seed = 123){
                                              set.seed(seed)
                                              nodes <- matrix(runif(n*2), ncol = 2)
                                              while(length(unique(round(nodes, 4))) < n){
                                                nodes <- matrix(runif(n*2), ncol = 2)
                                              }
                                              return(nodes)
                                            }
                                            
                                            # function to create an adjacency matrix based on the distance matrix and cutoff distance r
                                            create_adjacency_matrix <- function(nodes, r){
                                              d <- as.matrix(dist(nodes))
                                              adj_matrix <- ifelse(d <= r, 1, 0)
                                              return(adj_matrix)
                                            }
                                            
                                            # function to grow a spatial scale-free expander graph with parameters beta, alpha, mu, prewire, m, p_within, p_between
                                            grow_graph <- function(N, beta, alpha, mu, prewire, m, p_within, p_between, r, seed = 123){
                                              
                                              # initialize nodes
                                              nodes <- generate_nodes(N, seed)
                                              num_nodes <- N
                                              
                                              # create initial adjacency matrix based on distance
                                              adj_matrix <- create_adjacency_matrix(nodes, r)
                                              
                                              # initialize community assignments
                                              community <- rep(1, N)
                                              
                                              # initialize degree vector
                                              degrees <- rep(0, N)
                                              
                                              # initialize Laplacian matrix
                                              L <- diag(N)
                                              
                                              # grow graph
                                              for(i in 2:N){
                                                # calculate probability of connecting to each existing node
                                                within_comm <- which(community == community[i])
                                                between_comm <- which(community != community[i])
                                                dists_within <- dist(nodes[i,], nodes[within_comm,])
                                                prob_within <- exp(-alpha*dists_within)/degrees[within_comm]^beta
                                                prob_within[i] <- 0
                                                prob_within <- prob_within/sum(prob_within)
                                                dists_between <- dist(nodes[i,], nodes[between_comm,])
                                                prob_between <- exp(-mu*dists_between)/degrees[between_comm]^beta
                                                prob_between <- prob_between/sum(prob_between)
                                                
                                                # choose m nodes to connect to
                                                neighbors <- c(sample(within_comm, size = m, replace = TRUE, prob = p_within*prob_within),
                                                               sample(between_comm, size = m, replace = TRUE, prob = p_between*prob_between))
                                                
                                                # update adjacency matrix and degree vector
                                                adj_matrix[i, neighbors] <- 1
                                                adj_matrix[neighbors, i] <- 1
                                                degrees[i] <- m*2
                                                degrees[neighbors] <- degrees[neighbors] + 1
                                                
                                                # update community assignments
                                                community[i] <- ifelse(runif(1) < p_within, community[i], sample(unique(community[neighbors]), size = 1))
                                                
                                                # update Laplacian matrix
                                                L[i, i] <- sum(adj_matrix[i,])
                                                L[neighbors, i] <- -1
                                                L[i, neighbors] <- -1
                                                
                                                # rewire edges with probability prewire
                                                for(j in neighbors){
                                                  if(runif(1) < prewire){
                                                    dists <- dist(nodes[i,], nodes[-c(i,j),])
                                                    candidates <- which(dists <= r)
                                                    candidates <- candidates[!neighbors %in% candidates]
                                                    if(length(candidates) > 0){
                                                      new_neighbor <- sample(candidates, size = 1)
                                                      adj_matrix[i, j] <- 0
                                                      adj
                                                      #+++++++++++++++++++++++++++++++++++++++
                                                      library(Matrix)
                                                      
                                                      grow_spatial_scalefree_graph <- function(N, beta, alpha, mu, prewire, m, r, attribute_function = NULL, weight_function = NULL) {
                                                        # Define initial set of nodes with underlying spatial structure
                                                        coords <- matrix(runif(N * 2), ncol = 2)
                                                        while (TRUE) {
                                                          quad <- quadtree(coords)
                                                          if (max(quad$size) <= 1 / sqrt(N)) break # minimum distance between nodes
                                                          coords <- coords + rnorm(N * 2, sd = 1e-3) # jitter nodes until minimum distance is reached
                                                        }
                                                        
                                                        # Create adjacency matrix based on distance matrix
                                                        dmat <- as.matrix(dist(coords))
                                                        adj_mat <- +(dmat <= r)
                                                        
                                                        # Initialize graph with initial nodes and edges
                                                        graph <- graph_from_adjacency_matrix(adj_mat, mode = "directed")
                                                        V(graph)$coords <- coords # store node coordinates as attributes
                                                        
                                                        # Add nodes and edges using preferential attachment with spatial bias
                                                        for (i in N+1:m) {
                                                          # Choose target community for new node
                                                          community_probs <- degree(graph) ^ alpha
                                                          community_probs <- community_probs / sum(community_probs)
                                                          community <- sample(1:vcount(graph), size = 1, prob = community_probs)
                                                          
                                                          # Choose existing nodes to connect to within and outside target community
                                                          neighbors <- neighbors(graph, community, mode = "out")
                                                          non_neighbors <- setdiff(seq_len(vcount(graph)), c(community, neighbors))
                                                          
                                                          # Calculate probabilities for preferential attachment with spatial bias
                                                          dist_to_neighbors <- norm(coords[i,] - coords[neighbors,], type = "2")
                                                          dist_to_non_neighbors <- norm(coords[i,] - coords[non_neighbors,], type = "2")
                                                          pref_prob_neighbors <- (1 / (dist_to_neighbors + 1e-6))^beta
                                                          pref_prob_non_neighbors <- (1 / (dist_to_non_neighbors + 1e-6))^beta
                                                          pref_prob_neighbors <- pref_prob_neighbors / sum(pref_prob_neighbors)
                                                          pref_prob_non_neighbors <- pref_prob_non_neighbors / sum(pref_prob_non_neighbors)
                                                          
                                                          # Choose nodes to attach to based on probabilities
                                                          neighbors_to_attach <- sample(neighbors, size = m, replace = TRUE, prob = pref_prob_neighbors)
                                                          non_neighbors_to_attach <- sample(non_neighbors, size = m, replace = TRUE, prob = pref_prob_non_neighbors)
                                                          
                                                          # Add edges to target community and between target and non-target communities
                                                          edges_to_add <- cbind(rep(i, m), c(neighbors_to_attach, non_neighbors_to_attach))
                                                          graph <- add_edges(graph, edges_to_add)
                                                        }
                                                        
                                                        # Add rewiring for small world effect
                                                        for (i in 1:ecount(graph)) {
                                                          if (runif(1) < prewire) {
                                                            source <- as.integer(as.character(tail(ends(graph, i), n = 1)))
                                                            target <- as.integer(as.character(head(ends(graph, i), n = 1)))
                                                            dist_to_all_nodes <- norm(coords[source,] - coords, type = "2")
                                                            nodes_within_dist <- which(dist_to_all_nodes <= r & seq_along(dist_to_all_nodes) != source)
                                                            new_target <- sample(nodes_within_dist, size = 1)
                                                            graph <- delete_edges(graph, i)
                                                            graph <- add_edges(graph, cbind(rep(source, 1), new_target))
                                                          }
                                                        }
                                                        
                                                        # Add attributes and weights
                                                        if (!is.null(attribute_function)) {
                                                          V(graph)$
                                                            #+++++++++++++++++++++++++++++++++++++++
                                                            library(Matrix) # for sparse matrix operations
                                                          
                                                          grow_spatial_graph <- function(N, alpha, prewire, m, r, seed = 1, add_edge_weights = FALSE, add_node_attributes = FALSE) {
                                                            
                                                            set.seed(seed)
                                                            
                                                            # Generate initial set of nodes on a unit square
                                                            x <- runif(N)
                                                            y <- runif(N)
                                                            
                                                            # Define the adjacency matrix and the distance matrix
                                                            A <- Matrix(0, nrow = N, ncol = N, sparse = TRUE)
                                                            dist_mat <- matrix(0, nrow = N, ncol = N)
                                                            for (i in 1:(N-1)) {
                                                              for (j in (i+1):N) {
                                                                dx <- abs(x[i] - x[j])
                                                                dy <- abs(y[i] - y[j])
                                                                if (dx > 0.5) dx <- 1 - dx
                                                                if (dy > 0.5) dy <- 1 - dy
                                                                dist_mat[i,j] <- sqrt(dx^2 + dy^2)
                                                                if (dist_mat[i,j] <= r) {
                                                                  A[i,j] <- 1
                                                                  A[j,i] <- 1
                                                                }
                                                              }
                                                            }
                                                            
                                                            # Assign each node to a block/community
                                                            block_size <- floor(N / sqrt(N))
                                                            blocks <- rep(1:block_size, each = block_size)
                                                            blocks <- c(blocks, rep(blocks, each = block_size)[1:(N %% block_size)])
                                                            block_centers <- tapply(1:N, blocks, function(nodes) c(mean(x[nodes]), mean(y[nodes])))
                                                            
                                                            # Initialize the degree sequence with the initial set of nodes
                                                            degree_seq <- rep(m, N)
                                                            degrees <- sum(degree_seq)
                                                            
                                                            # Start growing the graph
                                                            for (i in (N+1):2*N) {
                                                              
                                                              # Pick a random block
                                                              block <- sample(1:block_size, 1)
                                                              
                                                              # Compute the attachment probabilities to each node within the same block and to each node in other blocks
                                                              block_nodes <- which(blocks == block)
                                                              in_block <- dist_mat[block_nodes, i-1] <= r
                                                              block_distances <- dist_mat[block_nodes, i-1]
                                                              block_degrees <- degree_seq[block_nodes]
                                                              block_probs <- (in_block * block_distances^alpha) / (sum(in_block * block_distances^alpha) + sum((!in_block) * block_degrees^alpha))
                                                              block_probs[is.na(block_probs)] <- 0
                                                              other_nodes <- which(blocks != block)
                                                              other_distances <- dist_mat[other_nodes, i-1]
                                                              other_degrees <- degree_seq[other_nodes]
                                                              other_probs <- (other_distances^alpha) / (sum(other_distances^alpha) + sum(other_degrees^alpha))
                                                              
                                                              #+++++++++++++++++++++++++++++++++++++++
                                                              # Function to generate a spatial scale-free expander graph
                                                              # with non-overlapping communities
                                                              # Arguments:
                                                              # - N: number of nodes
                                                              # - alpha: preferential attachment parameter
                                                              # - prewire: rewiring probability
                                                              # - m: number of edges to attach from a new node to existing nodes
                                                              # - r: cutoff distance for adjacency matrix
                                                              # - nblocks: number of communities
                                                              # Returns:
                                                              # - igraph object representing the generated graph
                                                              generate_spatial_scalefree_expander <- function(N, alpha, prewire, m, r, nblocks){
                                                                # Step 1: Generate spatial distribution of nodes
                                                                coords <- generate_spatial_coords(N, r)
                                                                
                                                                # Step 2: Compute adjacency matrix based on distance matrix
                                                                distmat <- dist(coords)
                                                                adjmat <- as.matrix(distmat <= r)
                                                                
                                                                # Step 3: Generate graph with non-overlapping communities
                                                                block_sizes <- sample(x = ceiling(N/nblocks):floor(N/nblocks), size = nblocks, replace = TRUE)
                                                                nodes <- rep(1:nblocks, times = block_sizes)
                                                                g <- igraph::make_empty_graph(N, directed = FALSE)
                                                                for (i in 1:N) {
                                                                  # preferential attachment
                                                                  if (i <= m + 1) {
                                                                    edges <- c(1:i-1)
                                                                  } else {
                                                                    degrees <- igraph::degree(g, mode = "out")
                                                                    p <- (degrees^alpha) / sum(degrees^alpha)
                                                                    edges <- sample(1:(i-1), size = m, replace = FALSE, prob = p[1:(i-1)])
                                                                  }
                                                                  
                                                                  # spatial attachment
                                                                  dists <- distmat[i, 1:(i-1)]
                                                                  p_spatial <- exp(-dists/mean(dists)) # exponential decay
                                                                  p_spatial[adjmat[i, 1:(i-1)] == 1] <- 0 # don't attach to nodes already connected
                                                                  p_spatial <- p_spatial / sum(p_spatial)
                                                                  edges_spatial <- sample(1:(i-1), size = m, replace = FALSE, prob = p_spatial)
                                                                  
                                                                  # assign nodes to blocks
                                                                  block <- nodes[i]
                                                                  for (j in 1:m) {
                                                                    if (runif(1) < prewire) {
                                                                      # rewiring to random node within certain distance
                                                                      dists_rewire <- distmat[edges[-j], i]
                                                                      dists_rewire[adjmat[edges[-j], i] == 1] <- Inf # don't rewire to nodes already connected
                                                                      rewire_idx <- which(dists_rewire <= r)
                                                                      if (length(rewire_idx) > 0) {
                                                                        edges[j] <- sample(edges[-j][rewire_idx], size = 1)
                                                                      }
                                                                    }
                                                                    # add edge to graph
                                                                    g <- igraph::add_edges(g, c(i, edges[j]))
                                                                    igraph::set_edge_attr(g, "block", e = E(g)[igraph::ecount(g)], value = block)
                                                                  }
                                                                }
                                                                
                                                                return(g)
                                                              }
                                                              
                                                              # Function to generate spatial coordinates for nodes on a unit square
                                                              # without allowing nodes to land too near each other
                                                              # Arguments:
                                                              # - N: number of nodes
                                                              # - r: minimum distance between nodes
                                                              # Returns:
                                                              # - matrix with (x,y) coordinates of nodes
                                                              generate_spatial_coords <- function(N, r) {
                                                                coords <- matrix(0, n
                                                                                 
                                                                                 #+++++++++++++++++++++++++++++++++++
                                                                                 library(Matrix)
                                                                                 
                                                                                 generate_spatial_graph <- function(N, alpha, prewire, m, r, seed = NULL) {
                                                                                   # Set seed
                                                                                   if (!is.null(seed)) {
                                                                                     set.seed(seed)
                                                                                   }
                                                                                   
                                                                                   # Define function to compute Euclidean distances between points
                                                                                   euclidean_distance <- function(x1, y1, x2, y2) {
                                                                                     sqrt((x1 - x2)^2 + (y1 - y2)^2)
                                                                                   }
                                                                                   
                                                                                   # Create initial set of nodes
                                                                                   nodes <- data.frame(x = runif(N), y = runif(N))
                                                                                   
                                                                                   # Use quadtree algorithm to avoid nodes being too close to each other
                                                                                   insert_into_quadtree <- function(qtree, x, y) {
                                                                                     if (is.null(qtree$nw)) {
                                                                                       if (is.null(qtree$x)) {
                                                                                         qtree$x <- x
                                                                                         qtree$y <- y
                                                                                       } else {
                                                                                         qtree$nw <- list(x = NULL, y = NULL, nw = NULL, ne = NULL, sw = NULL, se = NULL)
                                                                                         insert_into_quadtree(qtree, qtree$x, qtree$y)
                                                                                         insert_into_quadtree(qtree, x, y)
                                                                                       }
                                                                                     } else {
                                                                                       if (x <= (qtree$x + qtree$ne$x) / 2) {
                                                                                         if (y <= (qtree$y + qtree$sw$y) / 2) {
                                                                                           insert_into_quadtree(qtree$nw, x, y)
                                                                                         } else {
                                                                                           insert_into_quadtree(qtree$sw, x, y)
                                                                                         }
                                                                                       } else {
                                                                                         if (y <= (qtree$y + qtree$se$y) / 2) {
                                                                                           insert_into_quadtree(qtree$ne, x, y
                                                                                                                
                                                                                                                #+++++++++++++++++++++++++++++++++++++++
                                                                                                                # Define function to generate a spatial scale-free expander graph
                                                                                                                # Inputs:
                                                                                                                #   N: number of nodes
                                                                                                                #   alpha: parameter controlling the strength of the preferential attachment
                                                                                                                #   r: distance cutoff for adjacency matrix
                                                                                                                #   prewire: rewiring probability
                                                                                                                #   m: number of edges to attach from a new node
                                                                                                                # Outputs:
                                                                                                                #   adj_matrix: adjacency matrix of the generated graph
                                                                                                                #   pos: positions of the nodes
                                                                                                                #   communities: vector of community assignments for each node
                                                                                                                #   Laplacian: Laplacian matrix of the generated graph
                                                                                                                spatial_scalefree_expander <- function(N, alpha, r, prewire, m) {
                                                                                                                  # Set up the initial nodes
                                                                                                                  pos <- matrix(runif(N*2), ncol=2) # generate random positions
                                                                                                                  pos <- pos - apply(pos, 2, min) # shift to non-negative coordinates
                                                                                                                  pos <- pos / max(pos) # scale to unit square
                                                                                                                  
                                                                                                                  # Quadtree algorithm to generate non-overlapping initial nodes
                                                                                                                  quadtree <- function(xmin, ymin, xmax, ymax, depth, maxdepth, nodes) {
                                                                                                                    if (depth == maxdepth) {
                                                                                                                      return(nodes)
                                                                                                                    } else {
                                                                                                                      xmid <- (xmin + xmax)/2
                                                                                                                      ymid <- (ymin + ymax)/2
                                                                                                                      nodes_left <- nodes[, 1] < xmid
                                                                                                                      nodes_right <- nodes[, 1] >= xmid
                                                                                                                      nodes_bottom <- nodes[, 2] < ymid
                                                                                                                      nodes_top <- nodes[, 2] >= ymid
                                                                                                                      nodes_BL <- nodes[nodes_left & nodes_bottom, ]
                                                                                                                      nodes_BR <- nodes[nodes_right & nodes_bottom, ]
                                                                                                                      nodes_TL <- nodes[nodes_left & nodes_top, ]
                                                                                                                      nodes_TR <- nodes[nodes_right & nodes_top, ]
                                                                                                                      nodes_BL <- quadtree(xmin, ymin, xmid, ymid, depth+1, maxdepth, nodes_BL)
                                                                                                                      nodes_BR <- quadtree(xmid, ymin, xmax, ymid, depth+1, maxdepth, nodes_BR)
                                                                                                                      nodes_TL <- quadtree(xmin, ymid, xmid, ymax, depth+1, maxdepth, nodes_TL)
                                                                                                                      nodes_TR <- quadtree(xmid, ymid, xmax, ymax, depth+1, maxdepth, nodes_TR)
                                                                                                                      return(rbind(nodes_BL, nodes_BR, nodes_TL, nodes_TR))
                                                                                                                    }
                                                                                                                  }
                                                                                                                  nodes <- cbind(runif(N), runif(N))
                                                                                                                  nodes <- quadtree(0, 0, 1, 1, 0, 5, nodes)
                                                                                                                  pos <- nodes[1:N, ]
                                                                                                                  
                                                                                                                  # Create distance matrix and adjacency matrix
                                                                                                                  dist_matrix <- as.matrix(dist(pos))
                                                                                                                  adj_matrix <- as.matrix(dist_matrix <= r)
                                                                                                                  
                                                                                                                  # Set up initial communities
                                                                                                                  communities <- rep(1, N)
                                                                                                                  num_communities <- 1
                                                                                                                  
                                                                                                                  # Start adding nodes
                                                                                                                  for (i in 1:(N-1)) {
                                                                                                                    # Compute probabilities for each node based on spatial distance and degree
                                                                                                                    d <- dist(pos[i+1, ], pos[1:i, ])
                                                                                                                    deg <- rowSums(adj_matrix[1:i, 1:i])
                                                                                                                    prob <- d^(-alpha) + deg
                                                                                                                    prob[i+1] <- 0
                                                                                                                    prob <- prob / sum(prob)
                                                                                                                    
                                                                                                                    # Choose m nodes to attach to
                                                                                                                    attach_to <- sample(1:i, size=m, replace=TRUE, prob=prob)
                                                                                                                    adj_matrix[i+
                                                                                                                               
                                                                                                                               #+++++++++++++++++++++++++++++++++++++++
                                                                                                                               library(Matrix)
                                                                                                                               
                                                                                                                               # function to create spatial coordinates for nodes
                                                                                                                               # using a quadtree algorithm to avoid nodes landing too near or too far from each other
                                                                                                                               make_spatial_coords <- function(N) {
                                                                                                                                 coords <- matrix(0, nrow = N, ncol = 2)
                                                                                                                                 coords[1, ] <- runif(2, min = 0, max = 1)
                                                                                                                                 for (i in 2:N) {
                                                                                                                                   done <- FALSE
                                                                                                                                   while (!done) {
                                                                                                                                     quad <- sample(1:4, size = 1, replace = TRUE)
                                                                                                                                     if (quad == 1) {
                                                                                                                                       x <- runif(1, min = 0, max = coords[i-1,1])
                                                                                                                                       y <- runif(1, min = 0, max = coords[i-1,2])
                                                                                                                                     } else if (quad == 2) {
                                                                                                                                       x <- runif(1, min = coords[i-1,1], max = 1)
                                                                                                                                       y <- runif(1, min = 0, max = coords[i-1,2])
                                                                                                                                     } else if (quad == 3) {
                                                                                                                                       x <- runif(1, min = 0, max = coords[i-1,1])
                                                                                                                                       y <- runif(1, min = coords[i-1,2], max = 1)
                                                                                                                                     } else {
                                                                                                                                       x <- runif(1, min = coords[i-1,1], max = 1)
                                                                                                                                       y <- runif(1, min = coords[i-1,2], max = 1)
                                                                                                                                     }
                                                                                                                                     if (min(sqrt((coords[1:(i-1),1]-x)^2 + (coords[1:(i-1),2]-y)^2)) > 1/(4*sqrt(N))) {
                                                                                                                                       coords[i,] <- c(x, y)
                                                                                                                                       done <- TRUE
                                                                                                                                     }
                                                                                                                                   }
                                                                                                                                 }
                                                                                                                                 return(coords)
                                                                                                                               }
                                                                                                                               
                                                                                                                               # function to create an adjacency matrix based on the distance matrix
                                                                                                                               make_adjacency_matrix <- function(dist, r) {
                                                                                                                                 adj <- as.matrix(dist <= r)
                                                                                                                                 diag(adj) <- 0
                                                                                                                                 return(adj)
                                                                                                                               }
                                                                                                                               
                                                                                                                               # function to grow a spatial scale-free expander graph
                                                                                                                               # with non-overlapping communities, preferential attachment,
                                                                                                                               # and small world effect
                                                                                                                               spatial_scale_free_expander <- function(N, alpha, prewire, m, r, edge_weight = NULL, node_attribute = NULL) {
                                                                                                                                 # create spatial coordinates for nodes
                                                                                                                                 coords <- make_spatial_coords(N)
                                                                                                                                 
                                                                                                                                 # create distance matrix
                                                                                                                                 dist <- as.matrix(dist(coords))
                                                                                                                                 
                                                                                                                                 # create adjacency matrix based on distance matrix
                                                                                                                                 adj <- make_adjacency_matrix(dist, r)
                                                                                                                                 
                                                                                                                                 # initialize graph with m nodes and complete graph
                                                                                                                                 edges <- t(combn(1:m, 2))
                                                                                                                                 graph <- graph_from_edgelist(edges, directed = FALSE)
                                                                                                                                 degree <- degree(graph)
                                                                                                                                 
                                                                                                                                 # add new nodes with preferential attachment and small world effect
                                                                                                                                 for (i in (m+1):N) {
                                                                                                                                   # create probability distribution for preferential attachment
                                                                                                                                   pa_probs <- degree^(alpha)
                                                                                                                                   pa_probs <- pa_probs/sum(pa_probs)
                                                                                                                                   
                                                                                                                                   # choose m nodes to attach to preferentially
                                                                                                                                   attach_to <- sample(1:(i-1), size = m, replace = TRUE, prob = pa_probs)
                                                                                                                                   
                                                                                                                                   # create new edges to attach to chosen nodes
                                                                                                                                   new_edges <- matrix(c(rep(i, m), attach_to), ncol = 2)
                                                                                                                                   
                                                                                                                                   # rewire edges with small world probability
                                                                                                                                   for (j in 1:m)
                                                                                                                                     
                                                                                                                                     #+++++++++++++++++++++++++++++++++++
                                                                                                                                     # function to generate a spatial scale-free expander graph
                                                                                                                                     # N: number of nodes
                                                                                                                                     # alpha: parameter for degree distribution
                                                                                                                                     # prewire: rewiring probability for small world effect
                                                                                                                                     # m: number of edges to attach from a new node to existing nodes
                                                                                                                                     # r: distance cutoff for adjacency matrix
                                                                                                                                     # node_attr: node attributes
                                                                                                                                     # edge_weight: edge weights
                                                                                                                                     
                                                                                                                                   generate_spatial_scalefree_expander <- function(N, alpha, prewire, m, r, node_attr = NULL, edge_weight = NULL) {
                                                                                                                                     
                                                                                                                                     # initialize nodes with spatial structure
                                                                                                                                     coords <- matrix(0, nrow = N, ncol = 2)
                                                                                                                                     coords[1, ] <- runif(2)
                                                                                                                                     for (i in 2:N) {
                                                                                                                                       # use quadtree algorithm to find nearest neighbor and place new node
                                                                                                                                       qtree <- quadtree(coords[1:(i-1),], bbox = c(0, 0, 1, 1))
                                                                                                                                       repeat {
                                                                                                                                         new_coords <- runif(2)
                                                                                                                                         if (length(quadtree_find_within(qtree, new_coords, r/2)) == 0) {
                                                                                                                                           coords[i, ] <- new_coords
                                                                                                                                           break
                                                                                                                                         }
                                                                                                                                       }
                                                                                                                                     }
                                                                                                                                     
                                                                                                                                     # create adjacency matrix based on distance matrix
                                                                                                                                     dist_mat <- dist(coords)
                                                                                                                                     adj_mat <- ifelse(dist_mat <= r, 1, 0)
                                                                                                                                     diag(adj_mat) <- 0
                                                                                                                                     
                                                                                                                                     # initialize graph with m nodes and a fully connected subgraph
                                                                                                                                     edges <- matrix(0, nrow = m*(m-1), ncol = 2)
                                                                                                                                     for (i in 1:m) {
                                                                                                                                       for (j in 1:(i-1)) {
                                                                                                                                         edges[(i-1)*(m-1)+j,] <- c(i,j)
                                                                                                                                       }
                                                                                                                                     }
                                                                                                                                     adj_mat[1:m, 1:m] <- 0
                                                                                                                                     
                                                                                                                                     # grow the graph by adding new nodes and edges
                                                                                                                                     for (i in (m+1):N) {
                                                                                                                                       # calculate degree distribution for preferential attachment
                                                                                                                                       deg_dist <- colSums(adj_mat) + 1
                                                                                                                                       deg_dist[m:(i-1)] <- deg_dist[m:(i-1)] + alpha
                                                                                                                                       deg_prob <- deg_dist / sum(deg_dist)
                                                                                                                                       # preferential attachment for degree distribution
                                                                                                                                       neighbors <- sample(1:(i-1), m, prob = deg_prob[1:(i-1)], replace = TRUE)
                                                                                                                                       edges[((i-m-1)*m+1):(i-m)*m,] <- cbind(rep(i, m), neighbors)
                                                                                                                                       # preferential attachment for spatial distance
                                                                                                                                       dist_prob <- 1 / (dist_mat[i, neighbors] + 1e-6)
                                                                                                                                       dist_prob <- dist_prob / sum(dist_prob)
                                                                                                                                       attach_prob <- 0.5 * dist_prob + 0.5 * deg_prob[neighbors] / sum(deg_prob[neighbors])
                                                                                                                                       attach_prob <- attach_prob / sum(attach_prob)
                                                                                                                                       new_edges <- sample(neighbors, m, prob = attach_prob, replace = TRUE)
                                                                                                                                       edges[((i-m-1)*m+1):(i-m)*m,] <- cbind(rep(i, m), new_edges)
                                                                                                                                     }
                                                                                                                                     
                                                                                                                                     # add small world effect by rewiring edges
                                                                                                                                     for (i in 1:nrow(edges)) {
                                                                                                                                       if (runif(1) < prewire) {
                                                                                                                                         # find nodes within a certain distance for rewiring
                                                                                                                                         dist
                                                                                                                                         
                                                                                                                                         #+++++++++++++++++++++++++++++++++++++++
                                                                                                                                         generate_spatial_scale_free_graph <- function(N, alpha, prewire, m, r, weight_scale = 1, attr_scale = 1) {
                                                                                                                                           
                                                                                                                                           # Step 1: Generate initial set of nodes on a unit square with spatial structure
                                                                                                                                           nodes <- data.frame(x = runif(N), y = runif(N))
                                                                                                                                           
                                                                                                                                           # Create quadtree to keep track of spatial distances
                                                                                                                                           quadtree <- quadtree_create(nodes$x, nodes$y)
                                                                                                                                           
                                                                                                                                           # Step 2: Use the distance matrix to create an adjacency matrix for the graph
                                                                                                                                           dist_mat <- dist(nodes)
                                                                                                                                           adj_mat <- (dist_mat <= r) * 1
                                                                                                                                           
                                                                                                                                           # Step 3: Grow the network and organize nodes into communities with scale-free degree distribution
                                                                                                                                           # Initialize degree sequence with initial set of nodes
                                                                                                                                           degree_seq <- rep(m, N)
                                                                                                                                           # Initialize community assignment with all nodes in community 1
                                                                                                                                           comm_assign <- rep(1, N)
                                                                                                                                           # Initialize edge list with initial set of edges
                                                                                                                                           edge_list <- as.matrix(which(adj_mat == 1, arr.ind = TRUE))
                                                                                                                                           
                                                                                                                                           # Repeat until desired number of nodes is reached
                                                                                                                                           while (nrow(nodes) < N) {
                                                                                                                                             
                                                                                                                                             # Compute probability of attaching to existing nodes based on distance and degree
                                                                                                                                             dist_probs <- apply(nodes, 1, function(node) {
                                                                                                                                               comm_nodes <- nodes[comm_assign == comm_assign[which.min(dist(nodes, node))), ]
                                                                                                                                                                   degree_probs <- degree_seq[comm_assign == comm_assign[which.min(dist(nodes, node))), ]^alpha
                                                                                                                                                                                              sum((1 / (1 + dist(comm_nodes, node))) * degree_probs)
                                                                                                                                             })
                                                                attach_probs <- dist_probs / sum(dist_probs)
                                                                
                                                                # Choose m existing nodes to attach to based on probabilities
                                                                attach_nodes <- sample(1:N, m, replace = TRUE, prob = attach_probs)
                                                                attach_edges <- cbind(rep(nrow(nodes) + 1, m), attach_nodes)
                                                                
                                                                # Add new node to nodes and quadtree
                                                                new_node <- data.frame(x = runif(1), y = runif(1))
                                                                while (nrow(quadtree_query(quadtree, new_node$x, new_node$y, r / 2)) > 0) {
                                                                  new_node <- data.frame(x = runif(1), y = runif(1))
                                                                }
                                                                nodes <- rbind(nodes, new_node)
                                                                quadtree_insert(quadtree, new_node$x, new_node$y)
                                                                
                                                                # Assign new node to an existing community or create a new one
                                                                comm_probs <- table(comm_assign[attach_nodes])
                                                                if (length(comm_probs) == 1) {
                                                                  new_comm <- comm_probs
                                                                } else {
                                                                  comm_probs <- comm_probs^alpha
                                                                  new_comm <- sample(names(comm_probs), 1, prob = comm_probs / sum(comm_probs))
                                                                }
                                                                comm_assign <- c(comm_assign, new_comm)
                                                                degree_seq <- c(degree_seq, m)
                                                                
                                                                # Add new edges to edge list and update degree sequence
                                                                edge_list <- rbind(edge_list, attach_edges)
                                                                degree_seq[attach_nodes] <- degree_seq[
                                                                  
#+++++++++++++++++++++++++++++++++++++++
                                                               
#++++++++++++++++++++++++++++++++++++++
# Function to generate a spatial scale-free expander graph
# Arguments:
# - N: number of nodes
# - alpha: degree exponent
# - prewire: rewiring probability
# - m: number of edges to attach from a new node to existing nodes
# - r: cutoff distance for adjacency matrix
# - comm_size: average community size
# - max_attempts: maximum number of attempts to place a node without overlap
# Returns:
# - adj_matrix: adjacency matrix
# - node_attrs: node attributes (x, y coordinates)
# - edge_weights: edge weights
# - block_ids: block IDs
spatial_scale_free_expander <- function(N, alpha, prewire, m, r, comm_size, max_attempts = 100) {
  # Initialize data structures
  adj_matrix <- matrix(0, nrow = N, ncol = N)
  node_attrs <- matrix(0, nrow = N, ncol = 2)
  edge_weights <- matrix(0, nrow = N, ncol = N)
  block_ids <- integer(N)
  block_sizes <- integer(ceiling(N / comm_size))
  # Generate initial nodes with spatial structure
  node_attrs[1, ] <- runif(2)
  for (i in 2:N) {
    attempts <- 0
    repeat {
      node_attrs[i, ] <- runif(2)
      if (min(sqrt(rowSums((node_attrs[1:(i-1), ] - node_attrs[i, ])^2))) > r/2 ||
          attempts >= max_attempts) break
      attempts <- attempts + 1
    }
    if (attempts >= max_attempts) stop("Unable to place node without overlap")
  }
  # Create distance matrix and set adjacency matrix
  dist_matrix <- as.matrix(dist(node_attrs))
  adj_matrix[dist_matrix <= r] <- 1
  # Initialize communities
  for (i in 1:N) {
    block_sizes[block_ids[i] + 1] <- block_sizes[block_ids[i] + 1] + 1
  }
  # Add edges with preferential attachment
  for (i in 1:(N - m)) {
    # Compute probability of attaching to each node
    degrees <- colSums(adj_matrix[1:i, 1:i])
    block_degrees <- tapply(degrees, block_ids[1:i], sum)
    prob <- (degrees^alpha + sum(block_degrees)^alpha) / (i + sum(block_sizes)^alpha)
    prob[block_ids[1:i] == block_ids[i+1]] <- 0 # no self-loops or within-block edges
    prob <- prob / sum(prob)
    # Choose m nodes to attach to
    attach_nodes <- sample(1:i, size = m, replace = TRUE, prob = prob)
    # Add edges and update block IDs
    for (j in attach_nodes) {
      adj_matrix[i+1, j] <- 1
      adj_matrix[j, i+1] <- 1
      edge_weights[i+1, j] <- 1
      edge_weights[j, i+1] <- 1
      block_sizes[block_ids[j] + 1] <- block_sizes[block_ids[j] + 1] - 1
      if (runif(1) < prewire) {
        # Rewire edge to random node within distance r
        candidates <- which(dist_matrix[j, ] <= r & block_ids != block_ids[i+1] &
                              
                              #++++++++++++++++++++++++++++++++++++++
                              library(Matrix)
                            
                            # Function to generate a spatial scale-free expander graph
                            # n: number of nodes
                            # alpha: parameter that controls the strength of the degree effect
                            # prewire: probability of rewiring edges for small world effect
                            # m: number of edges to add for each new node
                            # r: distance cutoff for the spatial effect
                            # quadtreeDepth: depth of the quadtree for initial node placement
                            spatialScaleFreeExpander <- function(n, alpha, prewire, m, r, quadtreeDepth) {
                              # Function to compute the distance between two nodes in the unit square
                              distance <- function(x1, y1, x2, y2) {
                                sqrt((x1-x2)^2 + (y1-y2)^2)
                              }
                              
                              # Function to add a node to the graph with m edges to existing nodes
                              addNode <- function(graph, nodes, edges) {
                                # Compute the degree of each node
                                degrees <- rowSums(graph)
                                # Choose m existing nodes to connect to
                                targets <- sample(nodes, m, replace = TRUE, prob = degrees^alpha)
                                # Connect the new node to the chosen nodes
                                graph[nodes[length(nodes)], targets] <- 1
                                graph[targets, nodes[length(nodes)]] <- 1
                                # Add the new edges to the edge list
                                newEdges <- cbind(rep(nodes[length(nodes)], m), targets)
                                edges <- rbind(edges, newEdges)
                                return(list(graph = graph, nodes = nodes, edges = edges))
                              }
                              
                              # Function to compute the quadtree for node placement
                              quadtree <- function(x, y, depth) {
                                # Check if the quadtree is at the maximum depth
                                if (depth == 0) {
                                  return(list(x = x, y = y, nodes = NULL))
                                }
                                # Split the square into four quadrants
                                xmid <- (max(x) + min(x))/2
                                ymid <- (max(y) + min(y))/2
                                x1 <- x[x <= xmid]
                                x2 <- x[x > xmid]
                                y1 <- y[y <= ymid]
                                y2 <- y[y > ymid]
                                # Recursively compute the quadtree for each quadrant
                                q1 <- quadtree(x1, y1, depth-1)
                                q2 <- quadtree(x2, y1, depth-1)
                                q3 <- quadtree(x1, y2, depth-1)
                                q4 <- quadtree(x2, y2, depth-1)
                                # Combine the quadrants into a single quadtree
                                return(list(x = x, y = y, nodes = NULL, q1 = q1, q2 = q2, q3 = q3, q4 = q4))
                              }
                              
                              # Function to place a node in the quadtree and check for collisions with nearby nodes
                              placeNode <- function(qt, x, y, r) {
                                # Check if the current quadtree node has any nodes already
                                if (!is.null(qt$nodes)) {
                                  # Check the distance to
                                  
                                  #+++++++++++++++++++++++++++++++++++++++
                                  # Define function to generate spatial scale-free expander graph
                                  # N: number of nodes
                                  # alpha: parameter controlling degree effect (higher value means more scale-free)
                                  # prewire: rewiring probability for small-world effect
                                  # m: average degree of each node
                                  # r: cutoff distance for adjacency matrix
                                  generate_spatial_scale_free_expander_graph <- function(N, alpha, prewire, m, r) {
                                    # Step 1: Generate initial set of nodes with spatial structure
                                    coords <- matrix(nrow = N, ncol = 2)
                                    coords[1,] <- runif(2)
                                    for (i in 2:N) {
                                      repeat {
                                        candidate <- runif(2)
                                        distances <- sqrt(rowSums((coords[1:(i-1),] - candidate)^2))
                                        if (min(distances) > 1/N^(1/2)) {
                                          coords[i,] <- candidate
                                          break
                                        }
                                      }
                                    }
                                    
                                    # Step 2: Create distance matrix and adjacency matrix
                                    distances <- as.matrix(dist(coords))
                                    adj_matrix <- as.matrix(distances <= r)
                                    
                                    # Step 3: Grow graph by preferential attachment with spatial bias
                                    degrees <- rep(0, N)
                                    communities <- rep(0, N)
                                    max_community <- 1
                                    for (i in 1:m) {
                                      # Add new node with preferential attachment
                                      if (i == 1) {
                                        new_node <- sample(N, 1)
                                      } else {
                                        prob <- (degrees^alpha + 1) / sum(degrees^alpha + 1)
                                        new_node <- sample(N, 1, prob)
                                      }
                                      degrees[new_node] <- degrees[new_node] + 1
                                      
                                      # Assign node to community
                                      if (i <= max_community) {
                                        community <- i
                                      } else {
                                        prob <- (degrees[communities == 0]^alpha + 1) / sum(degrees[communities == 0]^alpha + 1)
                                        community <- sample(which(communities == 0), 1, prob)
                                      }
                                      communities[new_node] <- community
                                      
                                      # Connect node to existing nodes in community
                                      potential_nodes <- which(adj_matrix[new_node,] == 1 & communities == community)
                                      if (length(potential_nodes) > 0) {
                                        prob <- 1 / sqrt(rowSums((coords[new_node,] - coords[potential_nodes,])^2))
                                        prob <- prob / sum(prob)
                                        to_connect <- sample(potential_nodes, 1, prob)
                                        adj_matrix[new_node
                                                   
                                                   #+++++++++++++++++++++++++++++++++++++++
                                                   library(Matrix)
                                                   library(RANN)
                                                   
                                                   # function to generate a spatially distributed set of nodes
                                                   generate_nodes <- function(n, min_dist = 0.05, max_dist = 0.15) {
                                                     nodes <- matrix(runif(n*2), ncol=2)
                                                     quad <- quadtree(nodes)
                                                     while (TRUE) {
                                                       nn <- nn2(quad, nodes, k = 2)
                                                       dists <- nn$dists[,2]
                                                       if (all(dists >= min_dist) && all(dists <= max_dist)) break
                                                       nodes <- matrix(runif(n*2), ncol=2)
                                                       quad <- quadtree(nodes)
                                                     }
                                                     return(nodes)
                                                   }
                                                   
                                                   # function to create an adjacency matrix based on distance matrix and cutoff distance
                                                   create_adjacency_matrix <- function(dist_matrix, r) {
                                                     adj_matrix <- ifelse(dist_matrix <= r, 1, 0)
                                                     return(adj_matrix)
                                                   }
                                                   
                                                   # function to grow the spatial scale-free expander graph
                                                   grow_graph <- function(N, alpha, prewire, m) {
                                                     nodes <- generate_nodes(N)
                                                     dist_matrix <- as.matrix(dist(nodes))
                                                     adj_matrix <- create_adjacency_matrix(dist_matrix, r = mean(dist_matrix)/2)
                                                     degrees <- rowSums(adj_matrix)
                                                     for (i in 1:m) {
                                                       new_node <- runif(2)
                                                       new_dist <- sqrt(rowSums((nodes - new_node)^2))
                                                       new_degrees <- degrees
                                                       new_degrees <- c(new_degrees, 0)
                                                       for (j in 1:length(degrees)) {
                                                         p <- (1 + alpha) * degrees[j] * new_degrees[length(degrees) + 1] / ((1 + alpha) * degrees[j] + new_degrees[length(degrees) + 1])
                                                         if (runif(1) < p) {
                                                           adj_matrix[j, length(degrees) + 1] <- 1
                                                           adj_matrix[length(degrees) + 1, j] <- 1
                                                           degrees[j] <- degrees[j] + 1
                                                           new_degrees[length(degrees)] <- new_degrees[length(degrees)] + 1
                                                         }
                                                       }
                                                       if (runif(1) < prewire) {
                                                         p <- 0.01
                                                         while (runif(1) < p) {
                                                           node1 <- sample(1:N, 1)
                                                           node2 <- sample(1:N, 1)
                                                           if (node1 != node2 && adj_matrix[node1, node2] == 0 && dist_matrix[node1, node2] > r) {
                                                             adj_matrix[node1, node2] <- 1
                                                             adj_matrix[node2, node1] <- 1
                                                             break
                                                           }
                                                         }
                                                       }
                                                       nodes <- rbind(nodes, new_node)
                                                       degrees <- c(degrees, 0)
                                                       adj_matrix <- rbind(adj_matrix, rep(0, N))
                                                       adj_matrix <- cbind(adj_matrix, rep(0, N+1))
                                                       N <- N+1
                                                     }
                                                     return(adj_matrix)
                                                   }
                                                   
                                                   # example usage
                                                   adj_matrix <- grow_graph(N = 100, alpha = 1, prewire = 0.05, m = 50)
                                                   
                                                   #+++++++++++++++++++++++++++++++++++++++
                                                   generate_graph <- function(N, alpha, prewire, m) {
                                                     
                                                     # Step 1: Generate spatially distributed nodes using quadtree algorithm
                                                     coords <- generate_coords(N)
                                                     
                                                     # Step 2: Compute the distance matrix and create adjacency matrix based on cutoff distance r
                                                     r <- mean(dist(coords)) / 2 # cutoff distance based on mean distance
                                                     adj <- create_adjacency_matrix(coords, r)
                                                     
                                                     # Step 3: Organize nodes into communities and connect them with preferential attachment
                                                     communities <- find_communities(coords, m)
                                                     degrees <- calculate_degrees(adj)
                                                     for (i in 1:N) {
                                                       for (j in (i+1):N) {
                                                         if (communities[i] == communities[j]) { # within-community connection
                                                           p <- (degrees[i]^alpha + degrees[j]^alpha) / sum(degrees[communities == communities[i]]^alpha)
                                                           if (runif(1) < p) adj[i,j] <- adj[j,i] <- 1
                                                         } else { # between-community connection
                                                           d <- dist(coords[i], coords[j])
                                                           p <- (d^(-alpha) + degrees[j]) / sum((communities == communities[j]) * (d^(-alpha) + degrees))
                                                           if (runif(1) < p) adj[i,j] <- adj[j,i] <- 1
                                                         }
                                                       }
                                                     }
                                                     
                                                     # Step 4: Add rewiring probability for small world effect
                                                     for (i in 1:N) {
                                                       for (j in (i+1):N) {
                                                         if (adj[i,j] == 1 && runif(1) < prewire) {
                                                           d <- sample(seq(1,N), 1) # distance of random node to connect to
                                                           if (dist(coords[i], coords[d]) < r) { # rewiring within cutoff distance
                                                             adj[i,j] <- adj[j,i] <- 0
                                                             adj[i,d] <- adj[d,i] <- 1
                                                           }
                                                         }
                                                       }
                                                     }
                                                     
                                                     return(adj)
                                                   }
                                                   
                                                   generate_coords <- function(N) {
                                                     coords <- matrix(0, nrow = N, ncol = 2)
                                                     q <- quadtree(coords)
                                                     for (i in 1:N) {
                                                       repeat {
                                                         x <- runif(1)
                                                         y <- runif(1)
                                                         if (length(query(q, x, y)) == 0) break # make sure node doesn't land too close to others
                                                       }
                                                       coords[i,] <- c(x,y)
                                                       insert(q, coords[i,])
                                                     }
                                                     return(coords)
                                                   }
                                                   
                                                   create_adjacency_matrix <- function(coords, r) {
                                                     N <- nrow(coords)
                                                     adj <- matrix(0, nrow = N, ncol = N)
                                                     for (i in 1:N) {
                                                       for (j in (i+1):N) {
                                                         if (dist(coords[i], coords[j]) <= r) {
                                                           adj[i,j] <- adj[j,i] <- 1
                                                         }
                                                       }
                                                     }
                                                     return(adj)
                                                   }
                                                   
                                                   find_communities <- function(coords, m) {
                                                     N <- nrow(coords)
                                                     communities <- rep(1, N)
                                                     degrees <- calculate_degrees(create_adjacency_matrix(coords, Inf))
                                                     for (i in 1:m) { # grow m communities
                                                       center <- sample(seq(1,N), 1, prob = degrees/sum(degrees))
                                                       
                                                       #+++++++++++++++++++++++++++++++++++++





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

#  min_distance=r/2
#  n <- dim(points)[1]
#  quad <- list()
#  if (n > 1) {
#    mid_x <- median(points[,1])
#    mid_y <- median(points[,2])
#    quad[[1]] <- quadtree(points[points[,1] <= mid_x & points[,2] <= mid_y,], min_distance)
#    quad[[2]] <- quadtree(points[points[,1] <= mid_x & points[,2] > mid_y,], min_distance)
#    quad[[3]] <- quadtree(points[points[,1] > mid_x & points[,2] <= mid_y,], min_distance)
#    quad[[4]] <- quadtree(points[points[,1] > mid_x & points[,2] > mid_y,], min_distance)
# #   # Check for nodes that are too close and remove them
#    for (i in 1:4) {
#      if (is.null(quad[[i]])) {
#        next
#      }
#      if (dim(quad[[i]])[1] == 1) {
# #       continue
# #     }
# #     dist_matrix <- as.matrix(dist(quad[[i]]))
# #     close_nodes <- which(dist_matrix < min_distance & dist_matrix != 0, arr.ind = TRUE)
# #     if (length(close_nodes) > 0) {
# #       for (j in 1:nrow(close_nodes)) {
# #         quad[[i]] <- quadtree(quad[[i]][-close_nodes[j,2],], min_distance)
# #       }
# #     }
# #   }
# # } else {
# #   quad <- points
# # }
# 
# # Construct quadtree and set nodes as points that are not too close
# nodes <- do.call(rbind, quad)
