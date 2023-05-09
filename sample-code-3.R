  

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    generateSpatialGraph <- function(N, lambda, L, beta, p, probWithin, probBetween, add_weights = FALSE, add_attributes = FALSE) {
  
  # Generate points using 2D Poisson point process on a torus
  points <- matrix(stats::runif(N*2), ncol = 2)
  points <- L*points
  
  # Calculate distance matrix
  dist_mat <- as.matrix(dist(points, method = "euclidean"))
  
  # Create adjacency matrix
  adj_mat <- matrix(0, nrow = N, ncol = N)
  
  # Calculate preferential attachment weights
  pa_weights <- dist_mat^(-beta)
  
  # Connect nodes with preferential attachment
  for (i in 1:N) {
    neighbors <- sample(1:(i-1), size = floor(pa_weights[i, 1:(i-1)]/sum(pa_weights[i, 1:(i-1)])*N), replace = TRUE, prob = pa_weights[i, 1:(i-1)])
    adj_mat[i, neighbors] <- 1
    adj_mat[neighbors, i] <- 1
  }
  
  # Rewire edges for small world effect
  for (i in 1:N) {
    for (j in (i+1):N) {
      if (adj_mat[i,j] == 1 & runif(1) < p) {
        # Find a new neighbor for node i
        new_neighbor <- sample(setdiff(1:N, c(i,j)), size = 1)
        # Rewire edges
        adj_mat[i,j] <- 0
        adj_mat[j,i] <- 0
        adj_mat[i,new_neighbor] <- 1
        adj_mat[new_neighbor,i] <- 1
      }
    }
  }
  
  # Add community structure
  community_labels <- rep(1, N)
  num_communities <- max(2, floor(sqrt(N)))
  comm_size <- ceiling(N/num_communities)
  for (i in 1:num_communities) {
    community_labels[(i-1)*comm_size+1:min(i*comm_size,N)] <- i
  }
  within_probs <- matrix(probWithin, nrow = num_communities, ncol = num_communities)
  diag(within_probs) <- 1
  between_prob <- probBetween
  
  for (i in 1:N) {
    for (j in (i+1):N) {
      if (community_labels[i] == community_labels[j]) {
        if (runif(1) < within_probs[community_labels[i], community_labels[j]]) {
          adj_mat[i,j] <- 1
          adj_mat[j,i] <- 1
        }
      } else {
        if (runif(1) < between_prob) {
          adj_mat[i,j] <- 1
          adj_mat[j,i] <- 1
        }
      }
    }
  }
  
  # Create graph object
  graph <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected")
  
  # Add edge weights if specified
  if (add_weights) {
    weights <- pa_weights[adj_mat == 1]
    E(graph)$weight <- weights
  }
  
  # Add node attributes if specified
  if (add_attributes) {
    V(graph)$x <- points[,1]
    V(graph)$y <- points[,2]
    V(graph)$community <- community_labels
  }
  
  # Compute Laplacian matrix and eigenvalues

  
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Function to create a Spatial Expander Graph
  spatial_expander_graph <- function(N, lambda, L, beta, p, probWithin, probBetween, 
                                     create_edge_weights = FALSE, create_node_attributes = FALSE) {
    
    # Generate points using a 2D Poisson point process on an L-dimensional torus
    X <- matrix(runif(N*L), ncol=L) * 2*pi
    R <- sqrt(lambda/N)
    X <- R * qnorm(runif(N*L)) + X
    
    # Compute the distance matrix between all pairs of points
    dist_mat <- as.matrix(dist(X, method = "euclidean"))
    
    # Initialize adjacency matrix
    adj_mat <- matrix(0, N, N)
    
    # Connect nodes to one another with a probability attachment function
    for (i in 1:N) {
      neighbors <- which(dist_mat[i,] <= sort(dist_mat[i,])[min(beta, N)])
      degrees <- rowSums(adj_mat[neighbors,])
      degrees[neighbors == i] <- 0
      prob <- degrees^beta/sum(degrees^beta)
      prob[is.na(prob)] <- 0
      adj_mat[i,neighbors] <- rbinom(length(neighbors), 1, prob)
    }
    
    # Add rewiring probability for small world effect
    for (i in 1:N) {
      for (j in (i+1):N) {
        if (runif(1) < p) {
          adj_mat[i,j] <- 0
          adj_mat[j,i] <- 0
          k <- sample(1:N, 1)
          while (k == i || k == j || adj_mat[i,k] == 1 || adj_mat[j,k] == 1) {
            k <- sample(1:N, 1)
          }
          adj_mat[i,k] <- 1
          adj_mat[k,i] <- 1
        }
      }
    }
    
    # Add community structures to the graph
    comm_sizes <- sample(1:(N/10), N/10, replace=TRUE)
    comm_assign <- rep(1:(N/10), comm_sizes)
    comm_adj_mat <- matrix(0, N, N)
    for (i in 1:N) {
      for (j in (i+1):N) {
        if (comm_assign[i] == comm_assign[j]) {
          if (runif(1) < probWithin) {
            comm_adj_mat[i,j] <- 1
            comm_adj_mat[j,i] <- 1
          }
        } else {
          if (runif(1) < probBetween) {
            comm_adj_mat[i,j] <- 1
            comm_adj_mat[j,i] <- 1
          }
        }
      }
    }
    adj_mat <- adj_mat + comm_adj_mat
    
    # Create edge weights if desired
    if (create_edge_weights) {
      edge_weights <- matrix(runif(N^2), N, N)
      adj_mat <- adj_mat * edge_weights
    }
    
    # Create node attributes if desired
    if (create_node_attributes) {
      node_attributes <- data.frame(x=X[,1], y=X[,2], comm=comm_assign)
    }
    
    # Compute Laplacian matrix and eigenvalues
    D <- diag(rowSums(adj_mat))
    L <- D - adj_mat
    eigenvalues <- eigen(L)$values
    
   
            
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    spatial_expander_graph <- function(N, lambda, beta, p, probWithin, probBetween, ...){
      # generate points using 2D Poisson point process
      x <- runif(N)
      y <- runif(N)
      
      # calculate distance matrix
      dmat <- as.matrix(dist(cbind(x, y)))
      
      # create adjacency matrix
      adjmat <- matrix(0, N, N)
      for(i in 1:N){
        for(j in 1:i){
          if(runif(1) < exp(-lambda * dmat[i,j]^2)){
            adjmat[i,j] <- 1
            adjmat[j,i] <- 1
          }
        }
      }
      
      # preferential attachment with power law exponent beta
      degrees <- rowSums(adjmat)
      prob <- degrees^beta
      prob <- prob / sum(prob)
      
      for(i in 1:N){
        for(j in 1:(i-1)){
          if(runif(1) < p){
            # rewiring probability
            if(runif(1) < probWithin){
              if(degrees[i] > 0 && degrees[j] > 0){
                adjmat[i,j] <- 1
                adjmat[j,i] <- 1
                degrees[i] <- degrees[i] + 1
                degrees[j] <- degrees[j] + 1
                prob <- degrees^beta
                prob <- prob / sum(prob)
              }
            } else if(runif(1) < probBetween){
              if(degrees[i] > 0 && degrees[j] > 0){
                groups <- rep(1, N)
                groups[i] <- 2
                groups[j] <- 2
                group_prob <- matrix(0, 2, 2)
                for(k in 1:N){
                  for(l in 1:k){
                    if(groups[k] == groups[l]){
                      group_prob[groups[k], groups[l]] <- group_prob[groups[k], groups[l]] + adjmat[k,l]
                      group_prob[groups[l], groups[k]] <- group_prob[groups[k], groups[l]]
                    } else {
                      group_prob[groups[k], groups[l]] <- group_prob[groups[k], groups[l]] + adjmat[k,l]
                      group_prob[groups[l], groups[k]] <- group_prob[groups[k], groups[l]]
                    }
                  }
                }
                group_prob <- group_prob / sum(group_prob)
                if(runif(1) < group_prob[groups[i], groups[j]]){
                  adjmat[i,j] <- 1
                  adjmat[j,i] <- 1
                  degrees[i] <- degrees[i] + 1
                  degrees[j] <- degrees[j] + 1
                  prob <- degrees^beta
                  prob <- prob / sum(prob)
                }
              }
            }
          } else {
            if(runif(1) < prob[i] * prob[j]){
              adjmat[i,j] <- 1
              adjmat[j,i] <- 1
              degrees[i] <- degrees[i] + 1
              degrees[j] <- degrees[j] + 1
              prob <- degrees^beta
              prob <- prob / sum(prob)
            }
          }
        }
      }
      
      # create graph object
      graph <- list()
      graph$adjmat <- adjmat
      graph$x <- x
      graph$y <- y
      graph$degrees <- degrees
      graph$prob <- prob
      # add any additional arguments as needed
      
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      spatial_expander_graph <- function(N, d, lambda, L, beta,
attachment_type = "linear", p = 0, probWithin = 0, probBetween = 0, with_edge_weights = FALSE, with_node_attributes = FALSE) {
        # Generate nodes using a Poisson point process
        X <- matrix(rpois(N*d, lambda), ncol = d, byrow = TRUE) / lambda
        for (i in 1:d) {
          X[,i] <- X[,i] * L[i]
        }
        
        # Calculate distance matrix
        D <- as.matrix(dist(X))
        
        # Create adjacency matrix
        A <- matrix(0, N, N)
        for (i in 1:N) {
          # Calculate probability of attachment for each node
          if (attachment_type == "linear") {
            p_i <- 1 / (1 + beta * D[i,]^2)
          } else if (attachment_type == "exponential") {
            p_i <- exp(-beta * D[i,])
          } else if (attachment_type == "power-law") {
            p_i <- (D[i,] + 1) ^ (-beta)
          } else {
            stop("Invalid attachment type")
          }
          
          # Preferential attachment
          for (j in 1:N) {
            if (i == j) {
              next
            }
            
            if (runif(1) < p_i[j] + (1 - p_i[j]) * p) {
              A[i,j] <- 1
              A[j,i] <- 1
            }
          }
        }
        
        # Add community structure
        if (probWithin > 0 || probBetween > 0) {
          com <- rep(1, N)
          for (i in 2:(N/2)) {
            com[sample(N, i)] <- i
          }
          
          for (i in 1:(N-1)) {
            for (j in (i+1):N) {
              if (com[i] == com[j]) {
                if (runif(1) < probWithin) {
                  A[i,j] <- 1
                  A[j,i] <- 1
                }
              } else {
                if (runif(1) < probBetween) {
                  A[i,j] <- 1
                  A[j,i] <- 1
                }
              }
            }
          }
        }
        
        # Add edge weights
        if (with_edge_weights) {
          for (i in 1:N) {
            for (j in 1:N) {
              if (A[i,j] == 1) {
                A[i,j] <- exp(-D[i,j])
              }
            }
          }
        }
        
        # Add node attributes
        if (with_node_attributes) {
          attrs <- list()
          attrs$X <- X
          attrs$D <- D
          attrs$com <- com
          
          return(list(A = A, attrs = attrs))
        }
        
        # Compute Laplacian matrix
        deg <- rowSums(A)
        D <- diag(deg)
        L <- D - A
        
        # Check expansion properties
        lambda <- eigen(L)$values[2]
        mu <- eigen(L)$values[N]
        gap <- lambda / mu
        
        # Return graph
        return(list(A = A, gap = gap))
      }
      
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      # Function to generate Spatial Expander Graph
      # with given parameters
      
      generateSpatialExpanderGraph <- function(N, d, lambda, L, beta, attachment_type, p = 0, probWithin = 0, probBetween = 0, include_weights = FALSE, include_node_attributes = FALSE) {
        
        # Generate points using Poisson point process with intensity lambda
        points <- matrix(rpois(N*d, lambda), ncol = d, byrow = TRUE)
        
        # Adjust points to torus if L is specified
        if (!is.null(L)) {
          points <- points %% L
        }
        
        # Calculate pairwise distance matrix
        dist_matrix <- as.matrix(dist(points))
        
        # Create adjacency matrix using preferential attachment function
        adjacency_matrix <- matrix(0, nrow = N, ncol = N)
        for (i in 1:N) {
          for (j in 1:N) {
            if (i != j) {
              # Calculate probability of attachment for nodes i and j
              if (attachment_type == "linear") {
                prob_ij <- 1/(1+beta*dist_matrix[i,j])
              } else if (attachment_type == "exponential") {
                prob_ij <- exp(-beta*dist_matrix[i,j])
              } else if (attachment_type == "powerlaw") {
                prob_ij <- 1/(1+dist_matrix[i,j]^beta)
              } else {
                stop("Invalid attachment_type specified")
              }
              # Apply preferential attachment based on degree and distance
              if (runif(1) < p) {
                # Connect to random node
                k <- sample(1:N, 1)
                adjacency_matrix[i,k] <- adjacency_matrix[k,i] <- 1
              } else if (runif(1) < prob_ij) {
                adjacency_matrix[i,j] <- adjacency_matrix[j,i] <- 1
              }
            }
          }
        }
        
        # Add community structures if specified
        if (probWithin > 0 || probBetween > 0) {
          # Generate random community assignment
          community_assignment <- sample(1:(N/2), N, replace = TRUE)
          # Connect nodes within communities with probWithin and between communities with probBetween
          for (i in 1:N) {
            for (j in (i+1):N) {
              if (community_assignment[i] == community_assignment[j]) {
                if (runif(1) < probWithin) {
                  adjacency_matrix[i,j] <- adjacency_matrix[j,i] <- 1
                }
              } else {
                if (runif(1) < probBetween) {
                  adjacency_matrix[i,j] <- adjacency_matrix[j,i] <- 1
                }
              }
            }
          }
        }
        
        # Create edge weights if specified
        if (include_weights) {
          edge_weights <- matrix(runif(N*N), nrow = N, ncol = N)
          edge_weights[adjacency_matrix == 0] <- 0
        } else {
          edge_weights <- NULL
        }
        
        # Create node attributes if specified
        if (include_node_attributes) {
          node_attributes <- matrix(runif(N), nrow = N, ncol = 1)
        } else {
          node_attributes <- NULL
        }
        
        # Compute Laplacian matrix
        degrees <- rowSums(adjacency_matrix)
        diagonal <- matrix(0, nrow = N, ncol = N)
        diag(diagonal) <- degrees
        laplacian_matrix <- diagonal - adjacency_matrix
        

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
      
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      
      
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      spatial_expander_graph <- function(N, d, lambda, L, beta, p, probWithin = 0, probBetween = 0, weighted = FALSE, attribute = FALSE) {
        # N: number of nodes
        # d: dimension of torus
        # lambda: intensity parameter of Poisson point process
        # L: length of each side of torus
        # beta: power law parameter for preferential attachment
        # p: rewiring probability
        # probWithin: probability of connecting nodes within a community
        # probBetween: probability of connecting nodes between communities
        # weighted: whether to generate edge weights
        # attribute: whether to generate node attributes
        
        # Generate Poisson point process
        coords <- matrix(runif(N*d), ncol = d)
        coords <- L * coords # scale to torus
        
        # Calculate distance matrix
        dist_mat <- as.matrix(dist(coords, method = "euclidean", diag = TRUE, upper = TRUE))
        for (i in 1:N) {
          for (j in (i+1):N) {
            dist_mat[i,j] <- min(dist_mat[i,j], dist_mat[j,i])
          }
        }
        
        # Create adjacency matrix
        adj_mat <- matrix(0, nrow = N, ncol = N)
        for (i in 1:N) {
          for (j in (i+1):N) {
            if (runif(1) < probWithin) {
              if (floor(coords[i,1]/L) == floor(coords[j,1]/L)) { # same row
                if (floor(coords[i,2]/L) == floor(coords[j,2]/L)) { # same column
                  if (runif(1) < probBetween) {
                    adj_mat[i,j] <- 1
                    adj_mat[j,i] <- 1
                  }
                }
              }
            }
          }
        }
        
        # Preferential attachment
        for (i in 1:N) {
          neighbors <- which(adj_mat[i,] == 1)
          neighbor_degrees <- rowSums(adj_mat[neighbors,])
          for (j in 1:N) {
            if (j != i && adj_mat[i,j] == 0) {
              prob <- (dist_mat[i,j] ^ (-beta)) * (neighbor_degrees[j] + 1)
              if (runif(1) < prob) {
                adj_mat[i,j] <- 1
                adj_mat[j,i] <- 1
              } else if (runif(1) < p) {
                new_neighbor <- sample(neighbors, 1)
                adj_mat[i,new_neighbor] <- 0
                adj_mat[new_neighbor,i] <- 0
                adj_mat[i,j] <- 1
                adj_mat[j,i] <- 1
              }
            }
          }
        }
        
        # Generate edge weights if requested
        if (weighted) {
          edge_weights <- matrix(0, nrow = N, ncol = N)
          for (i in 1:N) {
            for (j in (i+1):N) {
              if (adj_mat[i,j] == 1) {
                edge_weights[i,j] <- runif(1)
                edge_weights[j,i] <- edge_weights[i

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Required packages
library(Matrix)
library(pracma)
library(spatstat)

# Function to generate Spatial Expander Graph
Spatial_Expander_Graph_3 <- function(N, lambda, L, beta, prewire, prob_within,
                                     prob_between, add_weights=FALSE, add_attr=FALSE,
                                     add_comm=F){
  
  # Generate Nodes using Poisson Point Process
  # points <- rpoispp(lambda, dimyx = c(L, L), win = owin(c(0, L), c(0, L)))
  # coords <- as.matrix(points$x)
  # Generate nodes using Poisson point process
  coords <- matrix(runif(2 * N), ncol = 2)
  coords <- coords * L
  
  # Calculate distance matrix
  dist_mat <- as.matrix(dist(coords, diag = TRUE, upper = TRUE))
  dist_mat <- pmin(dist_mat, L - dist_mat)
 
  # Create adjacency matrix using preferential attachment model
  adj_mat <- matrix(0, nrow = N, ncol = N)
  adj_mat[1, 2] <- 1
  adj_mat[2, 1] <- 1
  for(i in 3:N){
    deg <- colSums(adj_mat[1:(i-1), 1:(i-1)])
    for(j in 1:(i-1)){
      p1=(dist_mat[i,j]^(-beta))*(deg[j])
      p2=sum((dist_mat[i,j]^(-beta))*(deg[j]))
      prob_ij <- p1/p2
      if(runif(1) < prob_ij){
        adj_mat[i, j] <- 1
        adj_mat[j, i] <- 1
      }
    }
  }
  
  # Add rewiring probability for small world effect
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (adj_mat[i, j] == 1 && runif(1) < prewire) {
        adj_mat[i, j] <- 0
        adj_mat[j, i] <- 0
        
        # Choose random long-range neighbor to connect to
        k <- sample((1:N)[-c(i,j)], size = 1)
        adj_mat[i, k] <- 1
        adj_mat[k, i] <- 1
      }
    }
  }
  
  # Add communities if requested
 if (add_comm) {
    # Divide nodes into communities
    n_communities <- ceiling(sqrt(N))
    comm_size <- ceiling(N / n_communities)
    communities <- rep(1:n_communities, each = comm_size)[1:N]
    # Iterate through each node and connect to other nodes within and between communities
    for (i in 1:(N-1)) {
      for (j in (i+1):N) {
        
        # Check if nodes are in the same community
        if (communities[i] == communities[j]) {
          
          # Connect with probability prob_within
          if (runif(1) < prob_within) {
            adj_mat[i, j] <- 1
            adj_mat[j, i] <- 1
          }
          
        } else {
          # Connect with probability prob_between
          if (runif(1) < prob_between) {
            adj_mat[i, j] <- 1
            adj_mat[j,i] <-1
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

g=Spatial_Expander_Graph_3 (N=50,lambda=10, L=1,
                           beta=1, prewire=0.4, add_comm= F, 
                           prob_within = 0.5, prob_between = 0.2, add_weights=F,
                           add_attr=F)

plot(g$GraphObject,vertex.label=NA,vertex.size=2)

#   if(i == 1){
#     adj_mat[1, 2] <- 1
#     adj_mat[2, 1] <- 1
#   }else{
#     deg <- colSums(adj_mat[1:i-1, 1:i-1])
#     for(j in 1:i-1){
#       prob_ij <- ((dist_mat[i, j]^(-beta))*(deg[j]^p))/sum((dist_mat[i, 1:i-1]^(-beta))*(deg[1:i-1]^p))
#       if(runif(1) < prob_ij){
#         adj_mat[i, j] <- 1
#         adj_mat[j, i] <- 1
#       }
#     }
#   }
# }


# # Add community structures
  # community_size <- floor(N/4)
  # communities <- rep(1:4, each = community_size)
  # adj_mat_within <- matrix(0, nrow = N, ncol = N)
  # for(i in 1:4){
  #   community_nodes <- which(communities == i)
  #   for(j in community_nodes){
  #     for(k in community_nodes){
  #       if(j < k){
  #         if(runif(1) < probWithin){
  #           adj_mat_within[j, k] <- 1
  #           adj_mat_within[k, j] <- 1
  #         }
  #       }
  #     }
  #   }
  # }
  # 
  # adj_mat_between <- matrix(0, nrow = N, ncol = N)
  # for(i in 1:4){
  #   for(j in 1:4){
  #     if(i < j){
  #       community_i_nodes <- which(communities == i)
  #       community_j_nodes <- which(communities == j)
  #       for(k in community_i_nodes){
  #         for(l in community_j_nodes){
  #           if(runif(1) < probBetween){
  #             adj_mat_between[k, l] <- 1
  #             adj_mat_between[l, k] <- 1
  #           }
  #         }
  #       }
  #     }
  #   }
  # }
  # 
  # adj_mat <- adj_mat + adj_mat_within + adj_mat_between
  
  # Create weighted adjacency matrix and node attributes
  # if(edge_weight){
  #   adj_mat[adj_mat == 1] <- runif(sum(adj_mat == 1))
  # }
  # if(node_attribute){
    
    
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


Spatial_Expander_Graph <- function(N, beta, lambda, L, p, probWithin, probBetween) {
  # Generate points using a 2D Poisson point process on a torus
  points <- matrix(runif(N*L), ncol=L)
  points <- points %% 1 # torus
  
  # Calculate the distance matrix
  dist_mat <- as.matrix(dist(points))
  
  # Create an adjacency matrix using preferential attachment with power-law degree distribution
  adj_mat <- matrix(0, nrow=N, ncol=N)
  degree_seq <- rep(0, N)
  for (i in 1:N) {
    neighbors <- which(degree_seq > 0)
    dist_to_neighbors <- dist_mat[i, neighbors]
    attachment_probs <- (degree_seq[neighbors]^beta) * exp(-lambda*dist_to_neighbors)
    attachment_probs <- attachment_probs / sum(attachment_probs)
    neighbor_selected <- sample(neighbors, size=1, prob=attachment_probs)
    adj_mat[i, neighbor_selected] <- 1
    adj_mat[neighbor_selected, i] <- 1
    degree_seq[i] <- degree_seq[i] + 1
    degree_seq[neighbor_selected] <- degree_seq[neighbor_selected] + 1
  }
  
  # Rewire edges with probability p
  for (i in 1:N) {
    for (j in (i+1):N) {
      if (adj_mat[i,j] == 1) {
        if (runif(1) < p) {
          possible_neighbors <- which(degree_seq > 0 & dist_mat[i,] < dist_mat[i,j])
          if (length(possible_neighbors) > 0) {
            new_neighbor <- sample(possible_neighbors, size=1)
            adj_mat[i,j] <- 0
            adj_mat[j,i] <- 0
            adj_mat[i,new_neighbor] <- 1
            adj_mat[new_neighbor,i] <- 1
            degree_seq[j] <- degree_seq[j] - 1
            degree_seq[new_neighbor] <- degree_seq[new_neighbor] + 1
          }
        }
      }
    }
  }
  
  # Add community structures
  comm_sizes <- sample(2:(N/4), size=floor(N/10), replace=TRUE)
  comm_sizes <- c(comm_sizes, rep(max(comm_sizes), length.out=(N-length(comm_sizes))))
  comm_indices <- rep(1:length(comm_sizes), comm_sizes)
  comm_adj_mat <- matrix(0, nrow=N, ncol=N)
  for (i in 1:(length(comm_sizes)-1)) {
    for (j in ((i+1):length(comm_sizes))) {
      within_prob <- ifelse(comm_indices == i, probWithin, probBetween)
      between_prob <- ifelse(comm_indices == j, probWithin, probBetween)
      for (k in which(comm_indices == i)) {
        for (l in which(comm_indices == j)) {
          if (runif(1) < within_prob) {
            comm_adj_mat[k,l] <- 1
            comm_adj_mat[l,k] <- 1
          } else if (runif(1) < between_prob) {
            comm_adj_mat[k,l] <- 1
            comm_adj_mat[l,k] <- 1
          }
        }
      }
    }
  }
  adj_mat <- adj_mat * comm_adj_mat
  
  # Compute Laplacian matrix
  deg_mat <- diag(rowSums(adj_mat))
  laplacian_mat <- deg_mat
  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  spatial_expander_graph <- function(N, lambda, L, beta, p, probWithin, probBetween, add_weights=FALSE, add_attributes=FALSE) {
    
    # Generate points using Poisson point process
    points <- matrix(runif(N * L), ncol=L)
    lambda_vol <- lambda ^ L
    num_points <- rpois(lambda_vol, N)
    num_points <- ifelse(num_points > 0, num_points, 0)
    idx <- rep(1:N, num_points)
    points <- points[idx, ]
    
    # Calculate distance matrix
    dist_mat <- as.matrix(dist(points, method="euclidean"))
    torus_dist <- function(d) {
      return(sqrt(pmin(d^2, (1-d)^2)))
    }
    dist_mat <- torus_dist(dist_mat)
    
    # Create adjacency matrix using preferential attachment and small world effect
    adj_mat <- matrix(0, nrow=N, ncol=N)
    degree_seq <- rep(1, N)
    for (i in 2:N) {
      new_edges <- sample(1:(i-1), degree_seq[i] - 1, replace=TRUE, prob=dist_mat[i,]^(-beta))
      adj_mat[i, new_edges] <- 1
      adj_mat[new_edges, i] <- 1
      degree_seq[i] <- length(which(adj_mat[i, ] == 1))
      degree_seq[new_edges] <- degree_seq[new_edges] + 1
    }
    # Rewiring with probability p
    for (i in 1:(N-1)) {
      for (j in (i+1):N) {
        if (runif(1) < p) {
          adj_mat[i, j] <- 0
          adj_mat[j, i] <- 0
          candidates <- which(dist_mat[i,] + dist_mat[j,] > 0 & adj_mat[i,] == 0 & adj_mat[j,] == 0)
          if (length(candidates) > 0) {
            new_edge <- sample(candidates, 1)
            adj_mat[i, new_edge] <- 1
            adj_mat[new_edge, i] <- 1
          }
        }
      }
    }
    
    # Add community structure to the graph
    num_communities <- round(sqrt(N))
    community_sizes <- rep(ceiling(N/num_communities), num_communities)
    community_sizes[num_communities] <- community_sizes[num_communities] - (sum(community_sizes) - N)
    community_labels <- rep(1:num_communities, community_sizes)
    within_comm_prob <- probWithin
    between_comm_prob <- probBetween
    for (i in 1:(N-1)) {
      for (j in (i+1):N) {
        if (community_labels[i] == community_labels[j]) {
          if (runif(1) < within_comm_prob) {
            adj_mat[i, j] <- 1
            adj_mat[j, i] <- 1
          }
        } else {
          if (runif(1) < between_comm_prob) {
            adj_mat[i, j] <- 1
            adj_mat[j, i] <- 1
          }
        }
      }
    }
    
    # Add edge weights and node attributes if requested
    if (add_weights) {
      edge_weights <- dist_mat^(-beta)
      edge_weights[which(adj_mat == 0)] <- 0
    }
    if (add_attributes) {
      node_attributes <- matrix(runif(N))
    }
                                
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    spatial_expander_graph <- function(N, lambda, L, beta, p, probWithin, probBetween, add_weights = FALSE, add_attributes = FALSE) {
      # Generate a 2D Poisson point process with intensity parameter lambda on an L-dimensional torus
      X <- matrix(runif(2*N), ncol = 2)
      X <- X * matrix(rep(L, each = 2), ncol = 2, byrow = TRUE)
      points <- pmin(X, L-X+1e-10)
      
      # Calculate the distance matrix
      dist_matrix <- as.matrix(dist(points))
      
      # Create an adjacency matrix using preferential attachment
      adj_matrix <- matrix(0, nrow = N, ncol = N)
      adj_matrix[1:2, 1:2] <- 1
      for (i in 3:N) {
        # Calculate the probability of connecting to each existing node based on distance and degree
        dist_probs <- 1/(dist_matrix[i, 1:(i-1)]^beta)
        deg_probs <- colSums(adj_matrix[1:(i-1),])
        probs <- dist_probs * deg_probs
        probs <- probs / sum(probs)
        # Choose a random node to connect to with probability proportional to the calculated probabilities
        connected_to <- sample(1:(i-1), size = 1, prob = probs)
        adj_matrix[i, connected_to] <- 1
        adj_matrix[connected_to, i] <- 1
        # Add rewiring
        for (i in 1:N) {
          for (j in (i+1):N) {
            if (runif(1) < p) {
              if (runif(1) < 0.5) {
                adj_matrix[i,j] <- 1
                adj_matrix[j,i] <- 1
                adj_matrix[i, which(adj_matrix[j,] == 1)] <- 0
                adj_matrix[which(adj_matrix[i,] == 1), j] <- 0
              } else {
                adj_matrix[j,i] <- 1
                adj_matrix[i,j] <- 1
                adj_matrix[j, which(adj_matrix[i,] == 1)] <- 0
                adj_matrix[which(adj_matrix[j,] == 1), i] <- 0
              }
            }
          }
        }
        # Add community structure
        num_communities <- round(sqrt(N))
        community_assignment <- rep(1:num_communities, each = floor(N/num_communities))
        community_assignment <- c(community_assignment, rep(community_assignment[1:(N - length(community_assignment))], each = 1))
        comm_probs <- matrix(rep(probBetween, times = num_communities*num_communities), nrow = num_communities)
        diag(comm_probs) <- probWithin
        for (i in 1:N) {
          for (j in (i+1):N) {
            if (community_assignment[i] == community_assignment[j]) {
              if (runif(1) < comm_probs[community_assignment[i], community_assignment[j]]) {
                adj_matrix[i,j] <- 1
                adj_matrix[j,i] <- 1
              }
            } else {
              if (runif(1) < comm_probs[community_assignment[i], community_assignment[j]]) {
                adj_matrix[i,j] <- 1
                adj_matrix[j,i] <- 1
              }
            }
          }
        }
        
  }}      # Add edge weights if desired
        
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
spatial_expander_graph <- function(N, lambda, L, beta, p, probWithin, probBetween, use_edge_weights = TRUE) {
  # Generate N points using a 2D Poisson point process with intensity parameter lambda
  points <- matrix(runif(N * L), ncol = L)
  points <- points * sqrt(lambda / N)
  
  # Create distance matrix between all pairs of points
  dist_mat <- as.matrix(dist(points))
  
  # Create adjacency matrix based on preferential attachment
  adj_mat <- matrix(0, nrow = N, ncol = N)
  degrees <- rep(0, N)
  for (i in 2:N) {
    probs <- (dist_mat[1:(i-1), i]^(-beta)) * degrees[1:(i-1)]
    if (sum(probs) == 0) {
      probs <- rep(1/(i-1), i-1)
    } else {
      probs <- probs / sum(probs)
    }
    j <- sample(1:(i-1), size = 1, prob = probs)
    adj_mat[i, j] <- 1
    adj_mat[j, i] <- 1
    degrees[i] <- degrees[i] + 1
    degrees[j] <- degrees[j] + 1
  }
  
  # Add rewiring probability for small world effect
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (adj_mat[i, j] == 0 && runif(1) < p) {
        adj_mat[i, j] <- 1
        adj_mat[j, i] <- 1
      }
    }
  }

  # Add community structure to the graph
  num_communities <- round(sqrt(N))
  community_labels <- rep(1:num_communities, each = N/num_communities)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (community_labels[i] == community_labels[j]) {
        if (runif(1) < probWithin) {
          adj_mat[i, j] <- 1
          adj_mat[j, i] <- 1
        }
      } else {
        if (runif(1) < probBetween) {
          adj_mat[i, j] <- 1
          adj_mat[j, i] <- 1
        }
      }
    }
  }
  
  # Create edge weights based on distance
  if (use_edge_weights) {
    adj_mat[adj_mat == 1] <- dist_mat[adj_mat == 1]
  }
  
  # Compute Laplacian matrix to check expansion properties
  degree_mat <- diag(degrees)
  laplacian_mat <- degree_mat - adj_mat
  eigenvalues <- eigen(laplacian_mat, symmetric = TRUE)$values
  eigenvalues <- sort(eigenvalues)
  gap <- eigenvalues[2] - eigenvalues[1]
  
  # Return list of graph properties
  return(list(adj_mat = adj_mat, dist_mat = dist_mat, laplacian_mat = laplacian_mat,
              eigenvalues = eigenvalues, gap = gap))
}

set.seed(123)
graph <- spatial_expander_graph(N = 100, lambda = 10, L = 1, beta = 1, p = 0.1,
                                probWithin = 0.4,probBetween=0.01, use_edge_weights = F)

G <- graph.adjacency(as.matrix(graph$adj_mat), mode="undirected")
G=G%>%simplify(remove.loops = TRUE,remove.multiple = TRUE)
plot(G,vertex.label=NA,vertex.size=2)
                                
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to generate Spatial Expander Graph
generateSpatialExpanderGraph <- function(N, lambda, L, beta, p, probWithin, probBetween, edgeWeight=FALSE) {
  
  # Generate N points using 2D Poisson point process on L-dimensional torus
  points <- matrix(runif(N*L), N, L) %% 1  # wrap around torus
  
  # Compute pairwise distance matrix between points
  dist_mat <- as.matrix(dist(points, method = "euclidean", diag = TRUE, upper = TRUE))
  
  # Compute degree of each node
  degree <- rowSums(dist_mat < Inf)
  
  # Calculate attachment probabilities
  if (beta == 0) {
    prob_attach <- rep(1/N, N)
  } else {
    prob_attach <- (degree + 1)^(beta)
    prob_attach <- prob_attach / sum(prob_attach)
  }
  
  # Rewiring step
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (runif(1) < p) {
        dist_ij <- dist_mat[i, j]
        prob_ij <- prob_attach[i] * prob_attach[j]
        if (runif(1) < prob_ij) {
          dist_ik <- min(dist_mat[i,])  # closest neighbor to i
          dist_jk <- min(dist_mat[j,])  # closest neighbor to j
          if (dist_ij > max(dist_ik, dist_jk)) {
            adj_mat[i, j] <- dist_ij * runif(1)^(1/beta)
            adj_mat[j, i] <- dist_mat[i, j]
          }
        }
      }
    }
  }
  
 # Add community structures to graph
  # comm_size <- round(N/3)
  # community <- rep(1:3, each = comm_size)[1:N]
  # prob_within <- probBetween + probWithin
  # prob_between <- 1 - prob_within
  # adj_mat <- matrix(0, N, N)
  # for (i in 1:(N-1)) {
  #   for (j in (i+1):N) {
  #     if (is.na(community[i] == community[j])) {
  #       if (runif(1) < probWithin) {
  #         if (runif(1) < prob_attach[i]) {
  #           adj_mat[i, j] <- 1
  #           adj_mat[j, i] <- 1
  #         }
  #       }
  #     } else {
  #       if (runif(1) < probBetween) {
  #         if (runif(1) < prob_attach[i] * prob_attach[j]) {
  #           adj_mat[i, j] <- 1
  #           adj_mat[j, i] <- 1
  #         }
  #       }
  #     }
  #   }
  # }

  # Add edge weights if requested
  if (edgeWeight) {
    adj_mat[adj_mat == 1] <- dist_mat[adj_mat == 1]
  }
  
  # Compute Laplacian matrix
  diag_mat <- diag(rowSums(adj_mat))
  lap_mat <- diag_mat - adj_mat
  
  # Return adjacency matrix and Laplacian matrix
  return(list(adjacency = adj_mat, laplacian = lap_mat))
}

set.seed(123)
graph <- generateSpatialExpanderGraph(N = 100, lambda = 10, L = 1, beta = 1, p = 0.1,
                                probWithin = 0.4,probBetween=0.2, edgeWeight=FALSE)

G <- graph.adjacency(as.matrix(graph$adj_mat), mode="undirected")
G=G%>%simplify(remove.loops = TRUE,remove.multiple = TRUE)
plot(G,vertex.label=NA,vertex.size=2)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

generate_spatial_expander_graph <- function(N, lambda, L, beta, p, probWithin, probBetween, use_weights = FALSE) {
  
  # Generate points using a 2D Poisson point process
  points <- matrix(runif(N*L), ncol = L)
  points <- points * L^(1/L) # Adjust for torus shape
  
  # Calculate distance matrix
  dist_matrix <- as.matrix(dist(points))
  
  # Create adjacency matrix with preferential attachment
  adj_matrix <- matrix(0, nrow = N, ncol = N)
  for (i in 1:N) {
    for (j in (i+1):N) {
      prob_attachment <- (dist_matrix[i,j]^(-beta)) * (degree(adj_matrix)[i]^probWithin) * (degree(adj_matrix)[j]^probBetween)
      if (runif(1) < prob_attachment) {
        if (use_weights) {
          adj_matrix[i,j] <- dist_matrix[i,j]
          adj_matrix[j,i] <- dist_matrix[i,j]
        } else {
          adj_matrix[i,j] <- 1
          adj_matrix[j,i] <- 1
        }
      }
    }
  }
  
  # Rewire edges with probability p
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (adj_matrix[i,j] > 0) {
        if (runif(1) < p) {
          # Choose a random node to rewire to
          k <- sample(1:N, 1)
          while (k == i || k == j) {
            k <- sample(1:N, 1)
          }
          # Rewire edge
          adj_matrix[i,j] <- 0
          adj_matrix[j,i] <- 0
          prob_attachment <- (dist_matrix[i,k]^(-beta)) * (degree(adj_matrix)[i]^probWithin) * (degree(adj_matrix)[k]^probBetween)
          if (runif(1) < prob_attachment) {
            if (use_weights) {
              adj_matrix[i,k] <- dist_matrix[i,k]
              adj_matrix[k,i] <- dist_matrix[i,k]
            } else {
              adj_matrix[i,k] <- 1
              adj_matrix[k,i] <- 1
            }
          }
        }
      }
    }
  }
  
  # Add community structure
  num_communities <- floor(sqrt(N))
  community_ids <- rep(1:num_communities, each = floor(N/num_communities))
  community_ids <- c(community_ids, rep(community_ids[1], N - length(community_ids)))
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (community_ids[i] == community_ids[j]) {
        prob_attachment <- 0.5 + (degree(adj_matrix)[i]^probWithin) * (degree(adj_matrix)[j]^probWithin)
        if (runif(1) < prob_attachment) {
          if (use_weights) {
            adj_matrix[i,j] <- dist_matrix[i,j]
            adj_matrix[j,i] <- dist_matrix[i,j]
          } else {
            adj_matrix[i,j] <- 1
            adj_matrix[j,i] <- 1
          }
        }
      } else {
        prob_attachment <- (degree(adj_matrix)[i]^probBetween) * (degree(adj_matrix)[j]^probBetween)
        if (runif(1) < prob_attachment) {
          if (use_weights) {
            adj_matrix[i,j] <- dist_matrix[i,j]
            adj_matrix[j,i] <- dist_matrix[i,j]
          } else {
            adj_matrix[i,j] <- 1
            adj_matrix[j,i] <- 1
          }}}
            
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to create a spatial scale-free network with small-world and community effect
# N: total number of nodes
# gamma: degree distribution exponent
# k: number of nearest neighbors to connect to
# p: rewiring probability
# rewired_fraction: desired fraction of edges to be rewired

create_spatial_scale_free_network <- function(N, gamma, k, p, rewired_fraction) {
  
  # Define the normalization constant
  c <- k^(gamma-1) / sum(1:(N-1)^(-gamma))
  
  # Create a regular lattice
  adj_matrix <- matrix(0, nrow=N, ncol=N)
  for (i in 1:N) {
    for (j in 1:N) {
      if (abs(i-j) <= k/2 || abs(i-j) > N-k/2) {
        adj_matrix[i,j] <- 1
      }
    }
  }
  
  # Calculate degree for each node
  degrees <- c * (1:N)^(-gamma)
  
  # Create a list of all edges, sorted by the sum of their endpoints' degrees
  edges <- data.frame()
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (adj_matrix[i,j] == 1) {
        edges <- rbind(edges, data.frame(from=i, to=j, weight=degrees[i]+degrees[j]))
      }
    }
  }
  edges <- edges[order(-edges$weight),]
  
  # Rewire edges with probability p
  num_rewired <- 0
  for (i in 1:nrow(edges)) {
    # Check if edge should be rewired
    if (runif(1) < p) {
      # Choose one endpoint randomly
      endpoint <- sample(c(edges[i,]$from, edges[i,]$to), size=1)
      # Find nodes that are not its neighbors
      non_neighbors <- setdiff(1:N, c(endpoint, neighbors(
        graph_from_adjacency_matrix(adj_matrix,mode ="undirected"), endpoint)))
      # Choose a new endpoint uniformly at random from non-neighbors
      new_endpoint <- sample(non_neighbors, size=1)
      # Rewire the edge
      adj_matrix[edges[i,]$from, edges[i,]$to] <- 0
      adj_matrix[edges[i,]$to, edges[i,]$from] <- 0
      adj_matrix[endpoint, new_endpoint] <- 1
      adj_matrix[new_endpoint, endpoint] <- 1
      num_rewired <- num_rewired + 1
    }
    # Check if desired fraction of edges has been rewired
    if (num_rewired >= rewired_fraction * nrow(edges)) {
      break
    }
  }
  
  # Create graph object and return
  g <- graph.adjacency(adj_matrix, mode="undirected")
  graph= simplify(g,remove.multiple = T,remove.loops = T)
  return(graph)
}

set.seed(123)
N <- 300
gamma <- 4
k <- 2
p <- 0.1
rewired_fraction <- 0.2
g <- create_spatial_scale_free_network(N, gamma, k, p, rewired_fraction)
plot(g, vertex.size=2, vertex.label=NA)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(Matrix)
library(spatstat)

# Function to generate a spatial expander graph
spatial_expander_graph <- function(N, lambda, L, beta, r, p, probWithin, probBetween, use_edge_weights = FALSE) {
  
  # Generate points using a 2D Poisson point process on a torus
  points <- matrix(runif(N * L), ncol = L)#runifpoint(N, dim = L, win = torus(L))
  
  # Calculate the distance matrix
  dist_mat <- as.matrix(dist(points))
  
  # Calculate the adjacency matrix based on spatial distance and preferential attachment
  adj_mat <- matrix(0, nrow = N, ncol = N)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      prob <- ((dist_mat[i,j] <= r) + 1e-6)^(-beta)
      if (runif(1) < prob) {
        adj_mat[i,j] <- 1
        adj_mat[j,i] <- 1
      }
    }
  }
  
  # Add rewiring probability
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (adj_mat[i,j] == 1 && runif(1) < p) {
        # Find all nodes within distance r and not already connected
        available_nodes <- which(dist_mat[i,] <= r & adj_mat[i,] == 0)
        # If there are any available nodes, randomly select one and rewire
        if (length(available_nodes) > 0) {
          k <- sample(available_nodes, 1)
          adj_mat[i,j] <- 0
          adj_mat[j,i] <- 0
          adj_mat[i,k] <- 1
          adj_mat[k,i] <- 1
        }
      }
    }
  }
  
  # Add community structure
  community_sizes <- c(ceiling(N/2), floor(N/2))
  comm <- rep(1:length(community_sizes), community_sizes)
  prob_mat <- matrix(probBetween, nrow = length(community_sizes), ncol = length(community_sizes))
  diag(prob_mat) <- probWithin
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (comm[i] == comm[j] && runif(1) < probWithin) {
        adj_mat[i,j] <- 1
        adj_mat[j,i] <- 1
      } else if (comm[i] != comm[j] && runif(1) < prob_mat[comm[i], comm[j]]) {
        adj_mat[i,j] <- 1
        adj_mat[j,i] <- 1
      }
    }
  }
  
  # Create the graph object
  if (use_edge_weights) {
    graph <- graph.adjacency(adj_mat, weighted = TRUE,mode = "undirected")
    graph=graph%>%simplify(remove.loops = TRUE,remove.multiple = TRUE)
    E(graph)$weight <- dist_mat[adj_mat == 1]
  } else {
    graph <- graph.adjacency(adj_mat, weighted = NULL,mode = "undirected")
    graph=graph%>%simplify(remove.loops = TRUE,remove.multiple = TRUE)
  }
  
  # Calculate Laplacian matrix and eigenvalues
  laplacian <- laplacian_matrix(graph, normalized = TRUE)
  eigenvalues <- eigen(laplacian)$values
  
  # Return the graph object and eigenvalues
  list(graph = graph, eigenvalues = eigenvalues)
}

N=100
L=1
beta=1
r=0.1
p=0.1
lambda=2
probWithin=0
probBetween=0
use_edge_weights = FALSE

g=spatial_expander_graph(N, lambda, L, beta, r, p, probWithin, probBetween, use_edge_weights = FALSE)
plot(g$graph,vertex.label=NA,vertex.size=2)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
spatial_expander_graph <- function(N, lambda, L, r, beta, p, probWithin, probBetween, weighted=FALSE) {
  # Generate N points using a 2D Poisson point process on an L-dimensional torus
  coords <- matrix(runif(N*L), ncol=L) * rep(2*pi, N)
  points <- matrix(cos(coords), ncol=L) * sqrt(-log(runif(N)))
  if (L > 1) {
    points <- t(apply(points, 1, function(x) { apply(x, 2, cumsum) }))
  }
  
  # Calculate the distance matrix
  dist_mat <- as.matrix(dist(points, diag=TRUE, upper=TRUE))
  
  # Create the adjacency matrix based on the distance matrix and preferential attachment
  adj_mat <- matrix(0, nrow=N, ncol=N)
  degree <- rep(0, N)
  for (i in 1:N) {
    neighbors <- which(dist_mat[i,] <= r)
    neighbors <- neighbors[neighbors != i]
    if (length(neighbors) > 0) {
      degree[neighbors] <- degree[neighbors] + 1
      if (weighted) {
        w <- exp(-beta * dist_mat[i,neighbors])
        adj_mat[i,neighbors] <- w
        adj_mat[neighbors,i] <- w
      } else {
        adj_mat[i,neighbors] <- 1
        adj_mat[neighbors,i] <- 1
      }
    }
    if (runif(1) < p) {
      new_neighbor <- sample(setdiff(1:N, c(i, neighbors)), size=1)
      if (weighted) {
        w <- exp(-beta * dist_mat[i,new_neighbor])
        adj_mat[i,new_neighbor] <- w
        adj_mat[new_neighbor,i] <- w
      } else {
        adj_mat[i,new_neighbor] <- 1
        adj_mat[new_neighbor,i] <- 1
      }
      degree[c(i, new_neighbor)] <- degree[c(i, new_neighbor)] + 1
    }
  }
  
  # Add community structure to the graph
  num_communities <- ceiling(sqrt(N))
  community <- rep(1:num_communities, each=ceiling(N/num_communities))[1:N]
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (community[i] == community[j]) {
        if (runif(1) < probWithin) {
          if (weighted) {
            w <- exp(-beta * dist_mat[i,j])
            adj_mat[i,j] <- w
            adj_mat[j,i] <- w
          } else {
            adj_mat[i,j] <- 1
            adj_mat[j,i] <- 1
          }
          degree[c(i, j)] <- degree[c(i, j)] + 1
        }
      } else {
        if (runif(1) < probBetween) {
          if (weighted) {
            w <- exp(-beta * dist_mat[i,j])
            adj_mat[i,j] <- w
            adj_mat[j,i] <- w
          } else {
            adj_mat[i,j] <- 1
            adj_mat[j,i] <- 1
          }
          degree[c(i, j)] <- degree[c(i, j)] + 1
        }
      }
    }
    graph <- graph.adjacency(adj_mat,mode = "undirected")
    graph=graph%>%simplify(remove.loops = TRUE,remove.multiple = TRUE)
  }

  
  return(graph)
}  # Normalize the adjacency matrix to
  
N=100
L=2
beta=4
r=0.7
p=0.1
lambda=2
probWithin=0.1
probBetween=0.01
weighted = FALSE

g=spatial_expander_graph(N, lambda, L, r, beta, p, probWithin, probBetween, weighted=FALSE)
plot(g,vertex.label=NA,vertex.size=2)


plot(mg,vertex.label=NA,vertex.size=2)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
