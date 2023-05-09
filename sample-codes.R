#Function to generate a 2D Poisson point process
ppp <- function(N, lambda, L){
  x <- runif(N, 0, L)
  y <- runif(N, 0, L)
  r <- sqrt(runif(N))L/sqrt(lambdapi)
  x <- x + rcos(2pirunif(N))
  y <- y + rsin(2pirunif(N))
  return(data.frame(x=x, y=y))
}

#Function to calculate distance matrix
distMatrix <- function(points){
  dmat <- as.matrix(dist(points))
  return(dmat)
}

#Function to create an adjacency matrix
adjMatrix <- function(dist, r, beta, p){
  A <- ifelse(dist < r, exp(-beta*dist/r), 0)
  A <- ifelse(runif(length(A)) < p, 0, A)
  return(A)
}

#Function to add community structures
addCommunity <- function(A, probWithin, probBetween, numCommunities){
  N <- nrow(A)
  commSize <- ceiling(N/numCommunities)
  communities <- rep(1:numCommunities, each=commSize)
  communities[(numCommunities*commSize+1):N] <- numCommunities
  for(i in 1:numCommunities){
    for(j in 1:numCommunities){
      if(i == j){
        A[(communities == i), (communities == i)] <- ifelse(runif(sum(communities == i)^2) < probWithin, 1, 0)
      } else{
        A[(communities == i), (communities == j)] <- ifelse(runif(sum(communities == i)*sum(communities == j)) < probBetween, 1, 0)
      }
    }
  }
  return(A)
}

#Function to calculate Laplacian matrix
laplacian <- function(A){
  D <- diag(colSums(A))
  L <- D - A
  return(L)
}

generate_spatial_expander_graph <- function(N, lambda, L, beta, r, p, probWithin, probBetween, use_edge_weights = FALSE) {
  library(Matrix)
  library(spatialwarnings)
  
  # Generate N points on an L-dimensional torus using a Poisson point process with intensity parameter lambda
  points <- matrix(runif(N*L), ncol=L)
  for (i in 1:L) {
    points[,i] <- points[,i] * L
  }
  pp <- spatstat::ppp(points[,1], points[,2], c(0,L), c(0,L))
  pp <- spatstat::rpoispp(lambda, W=pp$window)
  points <- as.matrix(pp)
  
  # Calculate the distance matrix between all pairs of points
  dist_matrix <- dist(points, diag=TRUE, upper=TRUE)
  
  # Create an adjacency matrix using a function that connects nodes with a probability favoring short spatial distance and scale-free degree preferential attachment
  adj_matrix <- matrix(0, nrow=N, ncol=N)
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        # Calculate probability of connection based on distance and degree
        prob <- exp(-dist_matrix[i,j]/r) * (degree_vector[i]^beta * degree_vector[j]^beta) / sum(degree_vector^beta)
        # Add edge with probability prob
        if (runif(1) < prob) {
          if (use_edge_weights) {
            adj_matrix[i,j] <- dist_matrix[i,j]
          } else {
            adj_matrix[i,j] <- 1
          }
        }
      }
    }
  }
  
  # Add a rewiring probability p for small world effect to the graph
  for (i in 1:N) {
    for (j in (i+1):N) {
      if (adj_matrix[i,j] != 0 && runif(1) < p) {
        adj_matrix[i,j] <- 0
        adj_matrix[j,i] <- 0
        k <- sample(1:(N-1), 1)
        if (k >= i) {
          k <- k + 1
        }
        adj_matrix[i,k] <- 1
        adj_matrix[k,i] <- 1
      }
    }
  }
  
  # Add community structures to the graph by connecting nodes within communities with within node probability probWithin more than between communities with probability probBetween
  community_sizes <- floor(sqrt(N))
  community_assignments <- rep(1:community_sizes, each=community_sizes)
  for (i in 1:N) {
    for (j in (i+1):N) {
      if (community_assignments[i] == community_assignments[j]) {
        if (runif(1) < probWithin) {
          adj_matrix[i,j] <- 1
          adj_matrix[j,i] <- 1
        }
      } else {
        if (runif(1) < probBetween) {
          adj_matrix[i,j] <- 1
          adj_matrix[j,i] <- 1
        }
      }
    }
  }
  
  # Create a sparse matrix from the adjacency matrix
  adj_sparse <- Matrix(adj_matrix, sparse=TRUE)
  
  # Compute the Laplacian matrix to check the expansion properties of the graph
  laplacian <- laplacian(adj_sparse, normalized=TRUE)
  
  # Return the adjacency matrix and Laplac
  
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Main function to generate spatial expander graph
spatialGraph <- function(N, lambda, L, r, beta, p, probWithin, probBetween, numCommunities=1){
  points <- ppp(N, lambda, L)
  dist <- distMatrix(points)
  A <- adjMatrix(dist, r, beta, p)
  if(numCommunities > 1){
    A <- addCommunity(A, probWithin, probBetween, numCommunities)
  }
  A <- ifelse(A > 0, 1, 0)
  L <- laplacian(A)
  return(list(points=points, adj=A, laplacian=L))
}

#Example usage
graph <- spatialGraph(N=100, lambda=0.1, L=10, r=1, beta=1, p=0.1, probWithin=0.8, probBetween=0.2, numCommunities=2)
plot(graph$points, col=rep(1:2, each=50))
image(graph$adj)
image(graph$laplacian)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(Matrix)

# Function to generate a 2D Poisson point process with intensity parameter lambda on an L-dimensional torus
generate_points <- function(N, lambda, L) {
  coords <- matrix(runif(N*L), ncol = L)
  return(mod(coords + runif(L), 1) * sqrt(lambda))
}

# Function to calculate the distance matrix between all pairs of points
distance_matrix <- function(points, L) {
  N <- dim(points)[1]
  dist_mat <- matrix(0, nrow = N, ncol = N)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      dist_ij <- min(abs(points[i,] - points[j,]), abs(L - abs(points[i,] - points[j,])))
      dist_mat[i,j] <- sqrt(sum(dist_ij^2))
      dist_mat[j,i] <- dist_mat[i,j]
    }
  }
  return(dist_mat)
}

# Function to create an adjacency matrix for the graph with edge weights
create_adjacency_matrix <- function(dist_mat, r, beta, p) {
  N <- dim(dist_mat)[1]
  adj_mat <- matrix(0, nrow = N, ncol = N)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (dist_mat[i,j] < r) {
        weight <- exp(-beta*dist_mat[i,j])
        if (runif(1) < p) {
          # Small-world effect
          k <- sum(adj_mat[i,] != 0)
          candidates <- which(adj_mat[i,] == 0 & dist_mat[i,] < r)
          if (length(candidates) > 0) {
            probs <- weight*exp(-beta*dist_mat[i,candidates])*(1 + k*sqrt(log(N)/length(candidates)))
            probs <- probs/sum(probs)
            selected <- sample(candidates, 1, prob = probs)
            adj_mat[i,selected] <- weight
            adj_mat[selected,i] <- weight
          } else {
            adj_mat[i,j] <- weight
            adj_mat[j,i] <- weight
          }
        } else {
          adj_mat[i,j] <- weight
          adj_mat[j,i] <- weight
        }
      }
    }
  }
  return(adj_mat)
}

# Function to add community structures to the graph
add_communities <- function(adj_mat, probWithin, probBetween) {
  N <- dim(adj_mat)[1]
  labels <- sample(rep(1:ceiling(sqrt(N)), ceiling(N/ceiling(sqrt(N)))), N, replace = TRUE)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (labels[i] == labels[j]) {
        if (runif(1) < probWithin) {
          adj_mat[i,j] <- adj_mat[j,i] <- adj_mat[i,j]*2
        }
      } else {
        if (runif(1) < probBetween) {
          adj_mat[i,j] <- adj_mat[j,i] <- adj_mat[i,j]*2
        }
      }
    }
  }
  return(adj_mat)
}

# Function to compute the Laplacian matrix
laplacian_matrix <- function(adj_mat) {
  D <- diag(rowSums(adj_mat))
  return(D - adj_mat)
}

# Main function to generate the graph
generate_graph <- function(N, lambda,
                           
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Define a function to generate a spatial expander propagation graph
# on an L-dimensional torus
spatial_expander_propagation_graph <- function(n, L, beta, alpha, gamma, omega) {
  
  # Generate n points spatially distributed on the torus using a
  # 2D Poisson point process
  points <- matrix(runif(n * L), ncol = L)
  
  # Compute the distance matrix between all pairs of points on the torus
  distances <- as.matrix(dist(points, method = "euclidean", diag = TRUE, upper = TRUE))
  for (i in 1:n) {
    for (j in (i + 1):n) {
      d <- abs(points[i,] - points[j,])
      d <- min(d, 1 - d)
      distances[i,j] <- sqrt(sum(d^2))
      distances[j,i] <- distances[i,j]
    }
  }
  
  # Compute the adjacency matrix of the graph
  A <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in (i + 1):n) {
      p <- beta * exp(-alpha * distances[i,j]) + gamma * (degree(A, i) + degree(A, j)) / (2 * n)
      if (runif(1) < p) {
        A[i,j] <- 1
        A[j,i] <- 1
      }
    }
  }
  
  # Add community structure to the graph
  n_comms <- floor(sqrt(n))
  comm_size <- floor(n / n_comms)
  comm_indices <- rep(1:n_comms, each = comm_size)[1:n]
  for (i in 1:n) {
    for (j in (i + 1):n) {
      if (comm_indices[i] == comm_indices[j]) {
        if (runif(1) < omega) {
          A[i,j] <- 1
          A[j,i] <- 1
        }
      } else {
        if (runif(1) < (omega / n)) {
          A[i,j] <- 1
          A[j,i] <- 1
        }
      }
    }
  }
  
  # Add rewiring to the graph
  for (i in 1:n) {
    for (j in (i + 1):n) {
      if (A[i,j] == 1) {
        if (runif(1) < omega) {
          A[i,j] <- 0
          A[j,i] <- 0
          k <- sample(1:n, 1)
          while (k == i || A[i,k] == 1) {
            k <- sample(1:n, 1)
          }
          A[i,k] <- 1
          A[k,i] <- 1
        }
      }
    }
  }
  
  # Compute the degree matrix of the graph
  D <- diag(colSums(A))
  
  # Compute the Laplacian matrix of the graph
  L <- D - A
  
  # Check the expansion properties of the graph
  eigenvalues <- eigen(L)$values
  lambda_2 <- eigenvalues[2]
  if (lambda_2 > 1e-10) {
    cat("Graph is not an expander\n")
  } else {
    cat("Graph is an expander\n")
  }
  
  # Return the adjacency matrix, degree matrix,
  
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to create spatial expander propagation graph on an L-dimensional torus
# Arguments:
#   - n: number of nodes in the graph
#   - L: number of dimensions in the torus
#   - r: distance parameter for spatial distances
#   - p: probability of edge creation for adjacent nodes
#   - rewire_p: probability of rewiring edges
#   - comm_size: number of nodes per community
#   - attach_type: type of preferential attachment ('linear' or 'nonlinear')
#   - nonlinear_a: parameter for nonlinear attachment

create_torus_expander <- function(n, L, r, p, rewire_p, comm_size, attach_type, nonlinear_a) {
  
  # Create Poisson point process on L-dimensional torus
  X <- matrix(runif(n*L), ncol=L)
  for (i in 1:L) {
    X[,i] <- X[,i]*2*pi
  }
  X_torus <- apply(X, 2, function(x) {
    cos(x)
  })
  Y_torus <- apply(X, 2, function(x) {
    sin(x)
  })
  
  # Calculate distances between nodes
  D <- matrix(0, ncol=n, nrow=n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      dist <- 0
      for (k in 1:L) {
        d_k <- abs(X_torus[i,k] - X_torus[j,k])
        if (d_k > pi) {
          d_k <- 2*pi - d_k
        }
        dist <- dist + d_k^2
      }
      D[i,j] <- sqrt(dist)
      D[j,i] <- D[i,j]
    }
  }
  
  # Create adjacency matrix
  A <- matrix(0, ncol=n, nrow=n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (D[i,j] < r) {
        A[i,j] <- 1
        A[j,i] <- 1
      }
    }
  }
  
  # Preferential attachment for scale-free network
  if (attach_type == "linear") {
    degree <- rowSums(A)
    prob <- degree / sum(degree)
  } else if (attach_type == "nonlinear") {
    degree <- rowSums(A)
    prob <- (degree + nonlinear_a) / (sum(degree) + n*nonlinear_a)
  } else {
    stop("Invalid attachment type")
  }
  
  # Rewiring for small-world effect
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (A[i,j] == 1 && runif(1) < rewire_p) {
        # Find a node to rewire to
        node <- sample(1:n, 1, prob=prob)
        while (node == i || A[i,node] == 1) {
          node <- sample(1:n, 1, prob=prob)
        }
        # Rew
        
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        # Function to create a spatial expander propagation graph on an L-dimensional torus
        # Arguments:
        #   - n: the number of nodes in the graph
        #   - d: the degree of each node in the graph
        #   - p: the rewiring probability for small world effect
        #   - k: the number of communities in the graph
        #   - alpha: the preferential attachment parameter for scale-free
        #   - beta: the distance parameter for short spatial distances
        #   - weighted: a logical value indicating whether to use edge weights or not
        # Output:
        #   - A: the adjacency matrix of the graph
        
        create_spatial_expander <- function(n, d, p, k, alpha, beta, weighted = FALSE) {
          
          # Step 1: Create the 2D Poisson point process on an L-dimensional torus
          L <- 2 # we'll use a 2D Poisson point process
          lambda <- n^(1/L) # the Poisson parameter
          
          # Generate the points on the torus
          points <- matrix(runif(n*L), ncol=L)
          for (i in 1:L) {
            points[,i] <- points[,i] * 2*pi
          }
          
          # Map the points to the torus
          torus_points <- matrix(0, ncol=L, nrow=n)
          for (i in 1:L) {
            torus_points[,i] <- cos(points[,i])
            for (j in (i+1):L) {
              torus_points[,i] <- torus_points[,i] * cos(points[,j])
            }
            torus_points[,i] <- torus_points[,i] * sin(points[,i])
            for (j in (i+1):L) {
              torus_points[,i] <- torus_points[,i] * sin(points[,j])
            }
          }
          
          # Step 2: Create the adjacency matrix of the graph
          A <- matrix(0, nrow=n, ncol=n)
          
          # Connect nodes based on spatial distance and degree
          for (i in 1:n) {
            distances <- apply(torus_points - torus_points[i,], 1, function(x) {
              x[x > pi] <- x[x > pi] - 2*pi
              x[x < -pi] <- x[x < -pi] + 2*pi
              sum(x^2)
            })
            distances[i] <- Inf
            neighbors <- sample((1:n)[distances <= beta], d-1)
            A[i,neighbors] <- 1
            A[neighbors,i
              
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to create a spatial expander propagation graph on an L-dimensional torus
# Inputs:
#   n: number of nodes
#   L: dimension of torus
#   gamma: preferential attachment parameter
#   d0: distance parameter
#   p: rewiring probability
#   k: number of communities
#   alpha: community structures parameter
#   w: use edge weights (TRUE or FALSE)
# Output:
#   adjacency matrix of the graph

spatial_expander_propagation_graph <- function(n, L, gamma, d0, p, k, alpha, w) {
  
  # Generate 2D Poisson point process on torus
  xy <- matrix(runif(n*2, 0, L), ncol = 2)
  
  # Create distance matrix
  dist_mat <- as.matrix(dist(xy, method = "euclidean"))
  dist_mat[dist_mat > L/2] <- L - dist_mat[dist_mat > L/2]
  
  # Create adjacency matrix
  adj_mat <- matrix(0, nrow = n, ncol = n)
  
  # Preferential attachment
  for (i in 1:2) {
    for (j in (i+1):n) {
      d_ij <- dist_mat[i,j]
      p_ij <- gamma * (degree(adj_mat,i)^alpha) * (degree(adj_mat,j)^alpha) / (1 + d_ij/d0)^2
      if (runif(1) < p_ij) {
        if (w) {
          adj_mat[i,j] <- runif(1)
          adj_mat[j,i] <- adj_mat[i,j]
        } else {
          adj_mat[i,j] <- 1
          adj_mat[j,i] <- 1
        }
      }
    }
  }
  
  # Small world
  for (i in 1:n) {
    for (j in (i+1):n) {
      if (adj_mat[i,j] == 1) {
        if (runif(1) < p) {
          new_j <- sample((1:n)[-c(i,j)][dist_mat[j,-c(i,j)] < d0], 1)
          adj_mat[i,j] <- 0
          adj_mat[j,i] <- 0
          if (w) {
            adj_mat[i,new_j] <- runif(1)
            adj_mat[new_j,i] <- adj_mat[i,new_j]
            adj_mat[j,new_j] <- runif(1)
            adj_mat[new_j,j] <- adj_mat[j,new_j]
          } else {
            adj_mat[i,new_j] <- 1
            adj_mat[new_j,i] <- 1
            adj_mat[j,new_j] <- 1
            adj_mat[new_j,j] <- 1
          }
        }
      }
    }
  }
  
  # Community structures
  comm <- sample(1:k, n, replace = TRUE)
  for (i in 1:n) {
    for (j in (i+1):n) {
      if (comm[i] == comm[j]) {
        if (adj_mat[i,j] == 1) {
          if (runif(1) < p) {
            new_j <- sample(which(comm != comm[i] & dist_mat[j,] < d0), 1)
            adj_mat[i,j] <- 0
            adj_mat[j,i] <- 0
            if (w) {
              adj_mat[i,new_j] <-
                

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


library(Matrix)
library(spatstat)

# Function to create spatial expander propagation graph
# on an L-dimensional torus
#
# Arguments:
#   n: number of nodes
#   d: dimensionality of torus
#   r: radius of neighbourhood for preferential attachment
#   alpha: parameter for linear preferential attachment
#   gamma: parameter for nonlinear preferential attachment
#   sigma: parameter for spatial distance
#   p: rewiring probability for small world effect
#   k: number of communities
#   w: weight parameter for edges
#
# Returns:
#   A sparse adjacency matrix representing the graph
spatial_expander_propagation_graph <- function(n, d, r, alpha, gamma, sigma, p, k, w) {
  # Generate n nodes uniformly on the torus
  coords <- runifpoint(n, dim=d, win=owin(c(rep(0, d), rep(1, d))))
  coords <- cbind(coords, coords + 1) # wrap around torus
  
  # Compute pairwise distances
  dists <- as.matrix(dist(coords))
  dists[dists > 0.5] <- 1 - dists[dists > 0.5] # wrap around torus
  
  # Create adjacency matrix
  adj_mat <- Matrix(0, nrow=n, ncol=n, sparse=TRUE)
  
  # Create communities
  comm_sizes <- rep(floor(n/k), k) # equal sized communities
  comm_sizes[1:(n%%k)] <- comm_sizes[1:(n%%k)] + 1
  comm_starts <- cumsum(c(1, comm_sizes[-length(comm_sizes)]))
  comm_ends <- cumsum(comm_sizes)
  comm_nodes <- lapply(1:k, function(i) comm_starts[i]:comm_ends[i])
  
  # Loop over nodes and connect them
  for (i in 1:n) {
    # Compute distances and degrees of neighbouring nodes
    neighbours <- which(dists[i,] <= r & dists[i,] > 0)
    degrees <- degree(adj_mat, i)
    neighbour_dists <- dists[i, neighbours]
    neighbour_degrees <- degrees[neighbours]
    
    # Compute attachment probabilities
    if (length(neighbours) == 0) {
      probs <- rep(1/n, n)
    } else {
      if (gamma == 0) {
        probs <- neighbour_degrees
      } else {
        probs <- neighbour_degrees * (1 + alpha * neighbour_dists^(-gamma))
      }
    }
    
    # Rewire with probability p
    if (p > 0) {
      to_rewire <- sample(1:n, length(neighbours), replace=TRUE, prob=probs)
      for (j in neighbours) {
        if (runif(1) < p) {
          adj_mat[i, j] <- 0
          adj_mat[j, i] <- 0
          new_neighbour <- sample(setdiff(1:n, c(i, j)), 1, prob=probs)
          adj_mat[i, new_neighbour] <- w
          adj_mat[new_neighbour, i] <- w
        } else {
          adj_mat[i, j] <- w
          adj_mat[j, i] <- w
        }
      }
    } else {
      adj_mat[i, neighbours] <- w
      adj_mat[neighbours, i] <- w
    }
    
    # Preferential attachment
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

generate_spatial_expander <- function(n, L, use_weights = TRUE, beta = 1, r = 0.1, gamma = 1, mu = 0.1) {
  # n: number of nodes
  # L: number of dimensions
  # use_weights: boolean indicating whether or not to use edge weights
  # beta: preferential attachment parameter for scale-free behavior
  # r: rewiring probability for small-world behavior
  # gamma: distance parameter for spatial distances
  # mu: community structures parameter for within- and between-community connections
  
  # Generate random points on the torus using a 2D Poisson point process
  points <- matrix(runif(n*L, 0, 1), ncol = L)
  for (i in 1:L) {
    points[,i] <- points[,i] * 2 * pi
  }
  x <- cos(points[,1])
  for (i in 2:(L-1)) {
    x <- cbind(x, cos(points[,i]))
  }
  x <- cbind(x, sin(points[,L]))
  
  # Compute distance matrix
  d <- as.matrix(dist(points, diag = TRUE, upper = TRUE))
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) {
        d[i,j] <- 0
      } else {
        # Calculate minimum distance on the torus
        dists <- abs(points[i,] - points[j,])
        for (k in 1:L) {
          if (dists[k] > pi) {
            dists[k] <- 2*pi - dists[k]
          }
        }
        d[i,j] <- sum(dists^gamma)
      }
    }
  }
  
  # Initialize adjacency matrix
  adj_mat <- matrix(0, nrow = n, ncol = n)
  
  # Initialize degrees vector
  degrees <- rep(0, n)
  
  # Add edges between nodes
  for (i in 1:n) {
    # Preferential attachment
    if (i == 1) {
      choices <- c(2:n)
    } else {
      probs <- (degrees[1:(i-1)] + beta) / sum(degrees[1:(i-1)] + beta)
      choices <- sample(1:(i-1), 1, prob = probs)
    }
    
    for (j in choices) {
      # Small-world rewiring
      if (runif(1) < r) {
        new_choice <- sample(setdiff(1:n, c(i, j)), 1)
        adj_mat[i,new_choice] <- 1
        adj_mat[new_choice,i] <- 1
        adj_mat[j,i] <- 0
        adj_mat[i,j] <- 0
        degrees[i] <- degrees[i] + 1
        degrees[new_choice] <- degrees[new_choice] + 1
      } else {
        adj_mat[i,j] <- 1
        adj_mat[j,i] <- 1
        degrees[i] <- degrees[i] + 1
        degrees[j] <- degrees[j] + 1
      }
    }
  }
  
  # Add community structure
  num_communities <- ceiling(sqrt(n))
  community_sizes <- rep(floor(n/num_communities), num_communities)
  remainder <- n - sum(community_sizes)
  community_sizes[1:remainder] <- community_sizes[1
                                                  
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to generate a spatial expander propagation graph on an L-dimensional torus
# using a 2D poisson point process and preferential attachment
#
# Args:
#   n: number of nodes
#   L: number of dimensions
#   r: distance parameter for spatial distances
#   p: rewiring probability for small world effect
#   q: community structure parameter for connections within/between communities
#   use_weights: boolean indicating whether to use edge weights or not
#
# Returns:
#   adjacency matrix of the graph
#
spatial_expander_propagation_graph <- function(n, L=2, r=0.1, p=0.1, q=0.5, use_weights=TRUE) {
  
  # Generate random points on the torus using a 2D poisson point process
  coords <- matrix(runif(n*L), ncol=L) * 2 * pi
  coords_torus <- apply(coords, 2, function(x) cos(x)*(1+1e-3*cos(x))) # Torus coordinates
  
  # Calculate distance matrix between all pairs of points
  d <- as.matrix(dist(coords_torus, diag=TRUE, upper=TRUE))
  
  # Initialize the adjacency matrix with zeros
  adj_mat <- matrix(0, nrow=n, ncol=n)
  
  # Preferential attachment: connect each new node to an existing node with probability proportional to its degree
  for (i in 1:n) {
    # Calculate probability of connecting to each existing node
    prob <- d[i,]^(-r) * colSums(adj_mat) # degree-based attachment with distance parameter r
    # Normalize probabilities
    prob <- prob / sum(prob)
    # Select a node to connect to
    j <- sample.int(n, size=1, prob=prob)
    # Connect the nodes with an edge
    if (runif(1) > p) { # with probability 1-p, connect them directly
      if (use_weights) { # if using edge weights, assign weight as inverse of distance
        adj_mat[i,j] <- 1/d[i,j]
        adj_mat[j,i] <- 1/d[i,j]
      } else { # otherwise, assign binary value
        adj_mat[i,j] <- 1
        adj_mat[j,i] <- 1
      }
    } else { # with probability p, rewire the edge
      # Calculate probability of connecting to each node within/outside the same community
      comm <- ifelse(coords_torus[i,1] < pi, 1, 2) # divide into two communities along x-axis
      within_comm <- comm == ifelse(coords_torus[,1] < pi, 1, 2)
      prob_within <- within_comm * d[i,]^(-r) * colSums(adj_mat[within_comm,])
      prob_within <- prob_within / sum(prob_within)
      prob_between <- (!within_comm) * d[i,]^(-r) * colSums(adj_mat[!within_comm,])
      prob_between <- prob_between / sum(prob_between)
      # Select a node to connect to
      if (runif(1) < q) { # with probability q, connect within the same community
        j <- sample.int(n, size=1, prob=prob_within)
      } else { # with probability 1-q, connect between communities
        j <- sample.int(n, size=1, prob=prob_between)
      }
      #
      
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Define the parameters for adding edges
prob_within <- 0.5
prob_between <- 0.2

# Loop through each pair of nodes in the network
for (i in 1:nrow(adj_mat)) {
  for (j in 1:ncol(adj_mat)) {
    
    # Check if nodes i and j belong to the same community
    if (community[i] == community[j]) {
      
      # If the probability of a within-community edge is less than the threshold, add an edge
      if (runif(1) < prob_within && adj_mat[i, j] == 0) {
        adj_mat[i, j] <- 1
        adj_mat[j, i] <- 1
      }
    } else {
      
      # If the probability of a between-community edge is less than the threshold, add an edge
      if (runif(1) < prob_between && adj_mat[i, j] == 0) {
        adj_mat[i, j] <- 1
        adj_mat[j, i] <- 1
      }
    }
  }
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Function to create spatial expander propagation graph on L-dimensional torus
# Arguments:
#   n: number of nodes
#   L: dimension of the torus
#   edge_weights: logical value indicating whether to use edge weights or not
#   scale_free: boolean indicating whether to use scale-free preferential attachment or not
#   rewire_prob: rewiring probability for small world effect
#   comm_prob: probability of connection within communities
# Returns:
#   list containing the adjacency matrix, node attributes, and laplacian matrix of the graph
create_spatial_expander <- function(n, L, edge_weights = FALSE, scale_free = TRUE, rewire_prob = 0.1, comm_prob = 0.5) {
  # Generate points on torus using 2D poisson point process
  coords <- matrix(runif(n*L), ncol = L)
  points <- t(apply(coords, 1, function(x) L * x))
  
  # Initialize adjacency matrix
  adj_mat <- matrix(0, nrow = n, ncol = n)
  
  # Initialize node degrees
  deg <- rep(0, n)
  
  # Add edges between nodes with probability proportional to short spatial distance and high degree nodes
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      dist <- sqrt(sum((points[i,] - points[j,] - L*round((points[i,] - points[j,])/L))^2))
      prob <- exp(-dist) * (deg[i] + deg[j] + 1)
      if (runif(1) < prob) {
        adj_mat[i,j] <- 1
        adj_mat[j,i] <- 1
        deg[i] <- deg[i] + 1
        deg[j] <- deg[j] + 1
      }
    }
  }
  
  # Add scale-free preferential attachment if specified
  if (scale_free) {
    for (i in (n+1):n) {
      prob <- (deg + 1) / sum(deg + 1)
      neighbors <- sample(1:(i-1), size = 1, prob = prob[1:(i-1)])
      adj_mat[i,neighbors] <- 1
      adj_mat[neighbors,i] <- 1
      deg[i] <- deg[i] + 1
      deg[neighbors] <- deg[neighbors] + 1
    }
  }
  
  # Rewire edges with probability rewire_prob for small world effect
  for (i in 1:n) {
    for (j in (i+1):n) {
      if (adj_mat[i,j] == 1 & runif(1) < rewire_prob) {
        neighbors <- which(adj_mat[i,] == 1 & deg > 1)
        new_neighbor <- sample(neighbors, size = 1)
        adj_mat[i,j] <- 0
        adj_mat[j,i] <- 0
        adj_mat[i,new_neighbor] <- 1
        adj_mat[new_neighbor,i] <- 1
      }
    }
  }
  
  # Create community structure by adding more connections within communities
  comm_size <- round(n / 10)
  for (i in 1:(n/comm_size)) {
    comm_nodes <- ((i-1)*comm_size+1):(i*comm_size)
    for (j in comm_nodes) {
      for (k in comm_nodes) {
        
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


library(Matrix)

# Function to create a spatial expander propagation graph on an L-dimensional torus
# Arguments:
# n - number of nodes
# L - number of dimensions for torus
# use_weights - logical value indicating whether to use edge weights or not
# alpha - scale-free parameter for preferential attachment
# p_rewire - rewiring probability for small world effect
# q - proportion of edges within communities compared to between communities
# node_attributes - a matrix of size n x k containing node attributes
# Returns: a sparse adjacency matrix representing the graph
create_spatial_expander_propagation_graph <- function(n, L, use_weights = TRUE, alpha = 1, p_rewire = 0.1, q = 0.8, node_attributes = NULL) {
  
  # Generate nodes using a 2D Poisson point process on a torus
  points <- matrix(runif(n * L, min = 0, max = 1), ncol = L)
  torus_points <- abs(points - 0.5) * 2
  for (i in 1:L) {
    torus_points[,i] <- torus_points[,i] %% 1
  }
  
  # Create a fully connected graph
  graph <- spMatrix(n, n, 0)
  for (i in 1:n) {
    for (j in (i+1):n) {
      d <- min(abs(torus_points[i,] - torus_points[j,]), abs(torus_points[i,] - torus_points[j,] + 1))
      w <- 1 / (d + 1)
      graph[i,j] <- w
      graph[j,i] <- w
    }
  }
  
  # Add preferential attachment using a scale-free model
  degrees <- colSums(graph)
  for (i in (3:n)) {
    prob <- (degrees[1:(i-1)] + alpha) / (sum(degrees[1:(i-1)]) + alpha * (i-2))
    to <- sample(1:(i-1), size = 1, prob = prob)
    graph[i,to] <- 1
    graph[to,i] <- 1
  }
  
  # Rewire edges for small world effect
  for (i in 1:n) {
    for (j in (i+1):n) {
      if (graph[i,j] > 0 && runif(1) < p_rewire) {
        k <- sample(1:n, size = 1)
        graph[i,j] <- 0
        graph[j,i] <- 0
        graph[i,k] <- 1
        graph[k,i] <- 1
      }
    }
  }
  
  # Create community structure
  communities <- cutree(cluster::walktrap.community(as.undirected(graph)), k = 2)
  within_community_prob <- q / sum(communities == 1) + (1-q) / sum(communities == 2)
  between_community_prob <- (1-q) / sum(communities == 1) + q / sum(communities == 2)
  for (i in 1:n) {
    for (j in (i+1):n) {
      if (communities[i] == communities[j]) {
        if (runif(1) < within_community_prob) {
          graph[i,j] <- 1
          graph[j,i] <- 1
        }
      } else {
        # If the probability of a between-community edge is less than the threshold, add an edge
        if (runif(1) < between_community_prob && adj_mat[i, j] == 0) {
          graph[i, j] <- 1
          graph[j, i] <- 1
        }}  
        
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Required packages
library(spatstat)
library(Matrix)

# Function to generate spatial expander propagation graph
# Inputs:
#    - L: Number of dimensions
#    - n: Number of nodes
#    - radius: Radius of torus
#    - scale_free_param: Parameter for scale-free preferential attachment
#    - rewire_prob: Probability of rewiring edges
#    - community_prob: Probability of within-community edge connection
#    - weight_flag: TRUE for weighted edges, FALSE for unweighted edges
# Output:
#    - An adjacency matrix representing the spatial expander propagation graph
generate_spatial_expander <- function(L, n, radius, scale_free_param, rewire_prob, community_prob, weight_flag) {
  
  # Generate 2D Poisson point process on torus
  pp <- rPoisson(n, win=owin(c(rep(0, L), rep(radius, L))))
  pp <- unmark(pp) # Remove marks from spatstat format
  
  # Compute pairwise distance matrix
  dist_mat <- as.matrix(dist(pp, method = "euclidean", diag = FALSE, upper = TRUE))
  
  # Initialize adjacency matrix
  adj_mat <- Matrix(0, nrow=n, ncol=n, sparse=TRUE)
  
  # Construct edges with preferential attachment and within-community connection probability
  for (i in 1:n) {
    # Find nearest neighbors
    nn_idx <- order(dist_mat[i,])[2:(scale_free_param+1)]
    nn_degrees <- colSums(adj_mat[nn_idx,])
    # Preferentially attach to high degree neighbors
    new_edges <- sample(nn_idx, size=scale_free_param, replace=TRUE, prob=nn_degrees)
    # Connect to nodes in same community with higher probability
    for (j in new_edges) {
      if (runif(1) < community_prob || all(nn_degrees == 0)) {
        adj_mat[i,j] <- 1
        adj_mat[j,i] <- 1
      }
    }
  }
  
  # Rewire edges with small world probability
  for (i in 1:n) {
    for (j in (i+1):n) {
      if (adj_mat[i,j] == 1) {
        if (runif(1) < rewire_prob) {
          # Rewire edge
          k <- sample(1:n, size=1)
          while (k == i || k == j || adj_mat[i,k] == 1) {
            k <- sample(1:n, size=1)
          }
          adj_mat[i,j] <- 0
          adj_mat[j,i] <- 0
          adj_mat[i,k] <- 1
          adj_mat[k,i] <- 1
          if (weight_flag) {
            # Add edge weights proportional to inverse distance
            dist_ij <- dist_mat[i,j]
            dist_ik <- dist_mat[i,k]
            w_ij <- 1 / dist_ij
            w_ik <- 1 / dist_ik
            adj_mat[j,i] <- w_ij
            adj_mat[k,i] <- w_ik
          } else {
            adj_mat[j,k] <- 1
            adj_mat[k,j] <- 1
          }
        }
      }
    }
  }
}  
  # Add node attributes
 # node_attr <- data.frame(x=pp[,1],
                          

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
generate_spatial_expander <- function(n, L, r, k, beta, mu, communities, intra_prob, inter_prob, use_weights = TRUE, attributes = NULL) {
  
  # Generate points on torus using Poisson point process
  points <- matrix(runif(n*L), ncol = L)
  points <- r * ((2 * pi * points) %% 1)
  
  # Compute pairwise distances between points
  dist_mat <- as.matrix(dist(points))
  
  # Compute probability of connection based on distance and degree
  p_conn <- mu * dist_mat^(-L) + colSums(adjacency_matrix > 0) / k
  
  # Create adjacency matrix using preferential attachment and rewiring
  adjacency_matrix <- matrix(0, ncol = n, nrow = n)
  for (i in 1:n) {
    for (j in (i+1):n) {
      if (runif(1) < p_conn[i,j]) {
        if (runif(1) < intra_prob && abs(i-j) < n/communities) {
          adjacency_matrix[i,j] <- 1
          adjacency_matrix[j,i] <- 1
        } else if (runif(1) < inter_prob && abs(i-j) >= n/communities) {
          adjacency_matrix[i,j] <- 1
          adjacency_matrix[j,i] <- 1
        }
      }
    }
  }
  if (beta > 0) {
    for (i in 1:n) {
      for (j in (i+1):n) {
        if (adjacency_matrix[i,j] == 1 && runif(1) < beta) {
          new_j <- sample(setdiff(1:n, c(i,j)), 1)
          adjacency_matrix[i,j] <- 0
          adjacency_matrix[j,i] <- 0
          adjacency_matrix[i,new_j] <- 1
          adjacency_matrix[new_j,i] <- 1
        }
      }
    }
  }
  
  # Create edge weights based on distance
  if (use_weights) {
    weights <- dist_mat / mean(dist_mat[adjacency_matrix == 1])
    weights[adjacency_matrix == 0] <- 0
    adjacency_matrix <- weights
  }
  
  # Create Laplacian matrix
  degree_matrix <- diag(colSums(adjacency_matrix))
  laplacian_matrix <- degree_matrix - adjacency_matrix
  
  # Add node attributes
  if (!is.null(attributes)) {
    node_attributes <- attributes
  }
  
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
spatial_expander_propagation_graph <- function(num_points, torus_size, edge_weights = FALSE, rewiring_prob = 0.2, num_communities = 3, community_prob = 0.8, alpha = 1) {
  
  # Generate 2D Poisson point process
  points <- matrix(runif(num_points*2), ncol=2)
  points <- torus_size * points
  
  # Compute distances between points
  dists <- as.matrix(dist(points))
  
  # Create adjacency matrix based on distances and degree distribution
  adj_matrix <- matrix(0, nrow=num_points, ncol=num_points)
  for (i in 1:num_points) {
    # Sort distances to other nodes in increasing order
    sorted_dists <- sort(dists[i,])
    # Calculate probabilities of connection based on distances and degree distribution
    prob_connect <- alpha / (sorted_dists + 1)
    prob_connect[1] <- 0 # don't connect to self
    # Add edges to nodes with highest probabilities until degree is reached
    degree <- 0
    while (degree < sqrt(num_points)) {
      max_prob_index <- which.max(prob_connect)
      # Check if adding edge would create a cycle
      if (degree >= 2 && all(adj_matrix[max_prob_index,] == 0)) {
        # Randomly rewire edge with probability rewiring_prob
        if (runif(1) < rewiring_prob) {
          indices <- which(adj_matrix[max_prob_index,] == 1)
          new_index <- sample(setdiff(1:num_points, c(i,indices)), 1)
          adj_matrix[max_prob_index, indices] <- 0
          adj_matrix[indices, max_prob_index] <- 0
          adj_matrix[max_prob_index, new_index] <- 1
          adj_matrix[new_index, max_prob_index] <- 1
        } else {
          prob_connect[max_prob_index] <- 0
        }
      } else {
        adj_matrix[i, max_prob_index] <- 1
        adj_matrix[max_prob_index, i] <- 1
        prob_connect[max_prob_index] <- 0
        degree <- degree + 1
      }
    }
  }
  
  # Create community structure
  community_size <- num_points / num_communities
  community_labels <- rep(1:num_communities, each=community_size)
  community_adj_matrix <- matrix(0, nrow=num_points, ncol=num_points)
  for (i in 1:num_points) {
    for (j in 1:num_points) {
      if (community_labels[i] == community_labels[j]) {
        if (adj_matrix[i,j] == 1 && runif(1) < community_prob) {
          community_adj_matrix[i,j] <- 1
          community_adj_matrix[j,i] <- 1
        }
      }
    }
  }
  
  # Combine adjacency matrices with community structure
  adj_matrix <- adj_matrix * community_adj_matrix
  
  # Optionally add edge weights
  if (edge_weights) {
    weights <- runif(num_points^2)
    adj_matrix <- adj_matrix * weights
  }
  
  # Compute Laplacian matrix and check expansion properties
  d <- rowSums(adj_matrix)
  D <- diag(d)
  L <- D - adj_matrix
  eigenvalues <- eigen(L)$values
  lambda_2 <- sort(eigenvalues)[2]
  expansion_constant <- lambda_2 / d
  
  

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # N: the number of nodes in the graph
  # L: the dimensionality of the torus
  # lambda: the mean of the Poisson point process used to generate the nodes
  # alpha: the parameter that controls the number of edges per node
  # beta: the parameter that controls the strength of the spatial distance effect
  # mu: the parameter that controls the strength of the degree effect
  # sigma: the parameter that controls the decay of the degree effect with distance
  # rewiring_prob: the probability of rewiring an edge
  # use_weights: whether to use edge weights or not
  # 
  # Required packages
  library(Matrix)
  library(mvtnorm)
  library(Matrix.utils)
  
  # Function to generate a spatial expander propagation graph
  generate_spatial_expander <- function(N, L, lambda, alpha, beta, mu, sigma, rewiring_prob, use_weights) {
    
    # Step 1: Generate Poisson point process
    poisson_points <- rmvnorm(N, mean = rep(0, L), sigma = diag(rep(1/lambda, L)))
    poisson_points <- poisson_points %% 1 # L-dimensional torus
    
    # Step 2: Calculate distances and probabilities
    distances <- dist(poisson_points)
    adjacency_probabilities <- alpha * exp(-beta * distances^2) + mu / (1 + distances)^sigma
    
    # Step 3: Create adjacency matrix
    adjacency_matrix <- matrix(0, nrow = N, ncol = N)
    for (i in 1:N) {
      neighbors <- sample(1:N, size = round(alpha * N), prob = adjacency_probabilities[i, ], replace = TRUE)
      if (use_weights) {
        edge_weights <- rexp(length(neighbors), 1)
        adjacency_matrix[i, neighbors] <- edge_weights
      } else {
        adjacency_matrix[i, neighbors] <- 1
      }
    }
    
    # Step 4: Rewiring
    for (i in 1:N) {
      for (j in 1:N) {
        if (runif(1) < rewiring_prob) {
          k <- sample(1:N, size = 1)
          if (use_weights) {
            adjacency_matrix[i, j] <- 0
            adjacency_matrix[i, k] <- rexp(1, 1)
          } else {
            adjacency_matrix[i, j] <- 0
            adjacency_matrix[i, k] <- 1
          }
        }
      }
    }
    
    # Step 5: Laplacian matrix and expander check
    laplacian_matrix <- diag(rowSums(adjacency_matrix)) - adjacency_matrix
    eigenvalues <- eigen(laplacian_matrix)$values
    second_smallest_eigenvalue <- sort(eigenvalues)[2]
    expansion_gap <- second_smallest_eigenvalue / (2 * alpha)
    if (expansion_gap < 1) {
      print("Warning: Graph is not an expander")
    }
    
    # Step 6: Node attributes
    node_attributes <- runif(N)
    
    # Return the graph as a list
    return(list(adjacency_matrix = adjacency_matrix,
                node_attributes = node_attributes,
                expansion_gap = expansion_gap))
  }
  
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Install and load required packages
  if(!require("igraph")) install.packages("igraph")
  if(!require("spatstat")) install.packages("spatstat")
  library(igraph)
  library(spatstat)
  
  # Function to create spatial expander propagation graph on L-dimensional torus
  createSpatialExpanderGraph <- function(num_points, dim, radius, weight = FALSE, rewiring_prob = 0.1, community_prob = 0.3) {
    # Generate spatial point pattern
    pp <- rpoispp(num_points, win = hyperframe(replicate(dim, c(0, radius)), check=FALSE))
    
    # Create graph object
    g <- graph.empty(n = num_points, directed = FALSE)
    
    # Set node attributes
    node_attrs <- data.frame(id = 1:num_points)
    for (d in 1:dim) {
      node_attrs[paste0("coord", d)] <- pp$marks[[d]]
    }
    V(g)$attrs <- node_attrs
    
    # Compute pairwise distances
    dist_mat <- as.matrix(dist(node_attrs[, paste0("coord", 1:dim)]))
    
    # Create adjacency matrix
    adj_mat <- matrix(0, nrow = num_points, ncol = num_points)
    
    # Add edges with preferential attachment
    for (i in 1:num_points) {
      # Find short spatial distance neighbors
      neighbors <- which(dist_mat[i, ] <= mean(dist_mat[i, ]))
      
      # Find high degree nodes
      degrees <- degree(g)
      high_degree <- which(degrees > mean(degrees))
      
      # Combine neighbors and high degree nodes
      candidates <- unique(c(neighbors, high_degree))
      
      # Compute probability of connecting to each candidate node
      candidate_probs <- rep(1/length(candidates), length(candidates))
      candidate_probs[high_degree] <- candidate_probs[high_degree] * 2
      
      # Connect to nodes based on preferential attachment
      new_edges <- sample(candidates, size = round(sqrt(num_points)), replace = TRUE, prob = candidate_probs)
      for (j in new_edges) {
        if (i != j && adj_mat[i, j] == 0) {
          if (runif(1) > rewiring_prob && node_attrs[i, "coord1"] == node_attrs[j, "coord1"]) {
            # Connect nodes in same community
            adj_mat[i, j] <- 1
            adj_mat[j, i] <- 1
          } else {
            # Rewire with probability
            rewired_node <- sample(1:num_points, size = 1)
            adj_mat[i, rewired_node] <- 1
            adj_mat[rewired_node, i] <- 1
          }
        }
      }
    }
    
    # Set edge weights if desired
    if (weight) {
      edge_weights <- as.vector(dist_mat[adj_mat == 1])
      E(g)$weight <- edge_weights
    }
    
    # Set adjacency matrix
    g <- set.adjacency(g, adj_mat)
    
    # Ensure graph is connected
    g <- simplify(graph.edgelist(weakly.connected.components(g)$membership[[1]]))
    
    # Compute Laplacian matrix
    laplacian_mat <- laplacian_matrix(g, weight = weight)
    
    # Return graph object
    return(list(graph = g, laplacian = laplacian_mat))
  }
  
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to generate a spatial expander propagation graph on an L-dimensional torus
# Args:
#   n: number of nodes
#   L: number of dimensions in torus
#   radius: radius of 2D Poisson point process
#   k: number of nearest neighbors to connect each node to
#   p: probability of rewiring a connection
#   beta: parameter for tuning community structure
#   use_weights: boolean indicating whether to use edge weights
# Returns:
#   adjacency matrix of the graph
spatial_expander_propagation_graph <- function(n, L, radius, k, p, beta, use_weights) {
  # Generate nodes with 2D Poisson point process on L-dimensional torus
  points <- matrix(runif(n*L, min = 0, max = 1), ncol = L)
  points <- points * 2 * pi
  
  # Compute distances between nodes on torus
  dists <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      dists[i,j] <- sqrt(sum((points[i,] - points[j,])^2))
    }
  }
  
  # Compute probabilities for connecting nodes
  # proportional to short spatial distance and high degree nodes
  probs <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    degrees <- sort(dists[i,])[2:(k+1)]
    degrees <- degrees / degrees[k]
    for (j in 1:n) {
      if (i == j) {
        probs[i,j] <- 0
      } else {
        d_ij <- dists[i,j]
        probs[i,j] <- exp(-d_ij / radius) * degrees[j]
      }
    }
  }
  
  # Add community structure by rewiring connections
  for (i in 1:n) {
    for (j in (i+1):n) {
      if (runif(1) < p) {
        comm_i <- floor((i - 1) * beta) + 1
        comm_j <- floor((j - 1) * beta) + 1
        if (comm_i == comm_j) {
          # Rewire to a node within the same community
          nodes_in_comm <- which(floor((1:n - 1) * beta) + 1 == comm_i)
          probs_in_comm <- probs[i, nodes_in_comm]
          probs_in_comm[j - nodes_in_comm] <- 0
          probs_in_comm <- probs_in_comm / sum(probs_in_comm)
          new_j <- sample(nodes_in_comm, size = 1, prob = probs_in_comm)
          probs[i,j] <- 0
          probs[j,i] <- 0
          probs[i,new_j] <- exp(-dists[i,new_j] / radius) * degrees[new_j]
          probs[new_j,i] <- exp(-dists[i,new_j] / radius) * degrees[i]
        } else {
          # Rewire to a node in a different community
          nodes_in_comm_i <- which(floor((1:n - 1) * beta) + 1 == comm_i)
          nodes_in_comm_j <- which(floor((1:n - 1) * beta) + 1 == comm_j)
          probs_in_comm_i <- probs[i, nodes_in_comm_i]
          probs_in_comm_j <- probs[j, nodes_in_comm_j]
          
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to create spatial expander propagation graph
# on an L-dimensional torus
#
# Args:
#  - L: the dimension of the torus
#  - n: the number of nodes
#  - p: the probability of connecting two nodes with an edge
#  - q: the rewiring probability
#  - w: the weighting factor for edge weights
#  - community_structure: a vector of community assignments for each node
#  - node_attributes: a matrix of node attributes
#  - torus_radius: the radius of the torus (default: 1)
#
# Returns:
#  - A: the adjacency matrix of the graph
#  - D: the degree matrix of the graph
#  - L: the Laplacian matrix of the graph
#  - X: the node attributes matrix

create_spatial_expander_propagation_graph <- function(L, n, p, q, w, 
                                                      community_structure, node_attributes,
                                                      torus_radius=1) {
  # Generate n points uniformly distributed on the torus
  x <- matrix(runif(L * n, 0, torus_radius), n, L)
  
  # Compute the pairwise distances between points
  D <- as.matrix(dist(x, method="euclidean", diag=TRUE, upper=TRUE))
  
  # Compute the adjacency matrix
  A <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      # Compute the distance between nodes i and j on the torus
      dist_ij <- sqrt(sum(pmin(abs(x[i,] - x[j,]), 2*torus_radius - abs(x[i,] - x[j,]))^2))
      
      # Connect the nodes with probability proportional to short spatial distance
      if (runif(1) < p * exp(-dist_ij)) {
        A[i,j] <- 1
        A[j,i] <- 1
      }
    }
  }
  
  # Add edges to create a scale-free degree distribution
  for (i in 1:n) {
    k_i <- sum(A[i,])
    for (j in 1:n) {
      if (j != i && A[i,j] == 0) {
        k_j <- sum(A[j,])
        if (runif(1) < k_i * k_j / (n * (n-1)) * exp(-w * D[i,j])) {
          A[i,j] <- 1
          A[j,i] <- 1
        }
      }
    }
  }
  
  # Rewire edges to create small-world effect
  for (i in 1:n) {
    for (j in (i+1):n) {
      if (A[i,j] == 1 && runif(1) < q) {
        # Compute the distance between nodes i and j on the torus
        dist_ij <- sqrt(sum(pmin(abs(x[i,] - x[j,]), 2*torus_radius - abs(x[i,] - x[j,]))^2))
        
        # Find a random node k that is not i or j and is not already connected to i
        # or j, and connect i and j to k instead
        k <- sample(setdiff(1:n, c(i,j,which(A[i,] == 1 | A[j,] == 1))), 1)
        A[i,j
          
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Function to generate a spatial expander propagation graph
# Arguments:
#   - L: number of dimensions of the torus
#   - n: number of nodes
#   - r: radius of the torus
#   - p: probability of rewiring
#   - q: probability of intra-community edge
#   - w: use edge weights (TRUE/FALSE)
#   - beta: parameter controlling the strength of short-range connections
#   - alpha: parameter controlling the strength of high-degree connections
# Returns:
#   - adj: adjacency matrix
#   - X: node attributes
#   - lap: laplacian matrix
#   - degrees: degree of each node

generate_spatial_expander_propagation_graph <- function(L, n, r, p, q, w, beta, alpha) {
  
  # Generate n points on the torus using a 2D Poisson point process
  X <- matrix(runif(n*L, 0, r), ncol=L)
  while (nrow(X) < n) {
    X <- rbind(X, matrix(runif(n*L, 0, r), ncol=L))
    X <- X[!duplicated(X), ]
  }
  X <- X[1:n, ]
  
  # Compute pairwise distances between points
  dists <- as.matrix(dist(X))
  for (i in 1:n) {
    for (j in 1:i) {
      if (i != j) {
        d <- dists[i, j]
        for (k in 1:L) {
          if (d[k] > r/2) {
            d[k] <- r - d[k]
          }
        }
        dists[i, j] <- sqrt(sum(d^2))
        dists[j, i] <- dists[i, j]
      }
    }
  }
  
  # Compute probability of edge between nodes
  probs <- beta * exp(-dists/beta) + alpha * degree(X)^(-1)
  probs <- (1-q)*probs + q*diag(n)
  
  # Rewire edges
  for (i in 1:n) {
    for (j in (i+1):n) {
      if (runif(1) < p) {
        probs[i, j] <- q
        probs[j, i] <- q
      }
    }
  }
  
  # Generate adjacency matrix
  if (w) {
    adj <- matrix(0, nrow=n, ncol=n)
    for (i in 1:n) {
      for (j in (i+1):n) {
        if (runif(1) < probs[i, j]) {
          adj[i, j] <- adj[j, i] <- rnorm(1, mean=dists[i, j], sd=0.1)
        }
      }
    }
  } else {
    adj <- as.matrix(runif(n^2) < probs, nrow=n)
    adj <- lower.tri(adj) + t(lower.tri(adj))
  }
  
  # Compute laplacian matrix
  degrees <- rowSums(adj)
  lap <- diag(degrees) - adj
  
  return(list(adj=adj, X=X, lap=lap, degrees=degrees))
}

# Helper function to compute the degree of each node
degree <- function(X) {
  return(rowSums(as.matrix(dist(X)) < 0.1))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to generate a spatial expander propagation graph on an L-dimensional torus
# Arguments:
#   - n: number of nodes in the graph
#   - L: number of dimensions of the torus
#   - rho: density parameter for 2D Poisson point process
#   - edge_weights: logical indicating whether to use edge weights in the adjacency matrix
#   - rewire_prob: rewiring probability for small-world effect
#   - p_in: probability of connecting nodes within the same community
#   - p_out: probability of connecting nodes in different communities
#   - node_attr: matrix of node attributes, with one row per node and one column per attribute
# Returns:
#   - adjacency matrix of the graph
#   - matrix of node attributes
torus_expander_graph <- function(n, L, rho, edge_weights = FALSE, rewire_prob = 0.1, p_in = 0.7, p_out = 0.2, node_attr = NULL) {
  # Generate 2D Poisson point process on torus
  points <- matrix(runif(n * L), ncol = L) * 2 * pi
  points <- matrix(t(apply(points, 1, function(x) {cos(x)})), ncol = L)
  points <- sqrt(L / rho) * qnorm(points)
  
  # Compute pairwise distances between points
  distances <- as.matrix(dist(points))
  
  # Compute probabilities of edge formation based on distance and degree
  probs <- exp(-distances * sqrt(log(n) / n)) * (distances > 0)
  probs <- probs * (runif(n * n) < p_in)  # enforce community structure
  
  # Add rewiring for small-world effect
  rewire_mask <- runif(n * n) < rewire_prob
  probs[rewire_mask] <- exp(-abs(distances)[rewire_mask] * sqrt(log(n) / n)) * (distances[rewire_mask] > 0)
  probs[rewire_mask] <- probs[rewire_mask] * (runif(sum(rewire_mask)) < p_out)  # enforce community structure
  
  # Sample edges from probabilities
  edges <- matrix(runif(n * n) < probs, ncol = n)
  edges[lower.tri(edges)] <- t(edges)[lower.tri(edges)]
  
  # Construct adjacency matrix
  adj <- if (edge_weights) {
    distances * edges
  } else {
    edges
  }
  
  # Compute degree matrix
  deg <- diag(colSums(adj))
  
  # Compute Laplacian matrix
  lap <- deg - adj
  
  # Check expansion properties
  eigvals <- eigen(lap)$values
  min_eigval <- min(eigvals[-1])
  max_eigval <- max(eigvals)
  expansion_gap <- min_eigval / max_eigval
  if (expansion_gap < 0.01) {
    warning("Graph may not have strong expansion properties")
  }
  
  # Return adjacency matrix and node attributes
  return(list(adj, node_attr))
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to create a spatial expander propagation graph on an L-dimensional torus
# Arguments:
#   - n: number of nodes in the graph
#   - L: number of dimensions of the torus
#   - r: interaction radius (maximum spatial distance between connected nodes)
#   - p_rewire: rewiring probability for small world effect (0 = no rewiring, 1 = fully rewired graph)
#   - alpha: parameter for scale-free degree distribution (larger values lead to more nodes with high degree)
#   - n_comm: number of communities in the graph
#   - comm_size: vector of length n_comm specifying the size of each community
#   - node_attributes: data frame with attributes for each node (optional)
#   - weighted_edges: logical value indicating whether to use edge weights (default = FALSE)
# Returns:
#   - adjacency matrix of the graph
#   - node attributes with additional columns for community membership and degree

spatial_expander_graph <- function(n, L, r, p_rewire, alpha, n_comm, comm_size, node_attributes=NULL, weighted_edges=FALSE) {
  
  # Generate points uniformly on the torus
  coords <- matrix(runif(n*L), n, L)
  
  # Compute distances between points
  dists <- as.matrix(dist(coords, method="euclidean"))
  for (i in 1:n) {
    for (j in 1:i) {
      for (d in 1:L) {
        d1 <- abs(coords[i,d] - coords[j,d])
        d2 <- 1 - d1
        dists[i,j] <- dists[i,j] + min(d1, d2)^2
        dists[j,i] <- dists[i,j]
      }
    }
  }
  dists <- sqrt(dists)
  
  # Compute adjacency matrix based on spatial distance and degree
  adj_matrix <- matrix(0, n, n)
  for (i in 1:n) {
    # Select candidate neighbors within interaction radius
    candidates <- which(dists[i,] <= r)
    candidates <- candidates[candidates != i]
    if (length(candidates) == 0) {
      continue
    }
    # Compute probabilities proportional to spatial distance and degree
    probs <- (1 / (dists[i,candidates]^alpha)) * (degree(adj_matrix, candidates)^alpha)
    probs <- probs / sum(probs)
    # Sample neighbors according to probabilities
    neighbors <- sample(candidates, size=min(length(candidates), 5), replace=FALSE, prob=probs)
    # Add edges to adjacency matrix
    for (j in neighbors) {
      adj_matrix[i,j] <- 1
      adj_matrix[j,i] <- 1
    }
  }
  
  # Rewire edges for small world effect
  if (p_rewire > 0) {
    for (i in 1:n) {
      for (j in (i+1):n) {
        if (adj_matrix[i,j] == 1 && runif(1) < p_rewire) {
          # Remove edge
          adj_matrix[i,j] <- 0
          adj_matrix[j,i] <- 0
          # Choose new random neighbor
          candidates <- setdiff(1:n, c(i,j))
          probs <- rep(1/length(candidates), length(candidates))
          new_neighbor <- sample(candidates, size=1, replace=FALSE, prob=probs
                                 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+SAMPLE CODE
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#' Generate spatial expander propagation graph on an L-dimensional torus.
#'
#' @param n Number of nodes in the graph
#' @param L Number of dimensions for the torus
#' @param radius Spatial radius for 2D Poisson point process
#' @param edge_weights Logical indicating whether to use edge weights (default: TRUE)
#' @param p_rewire Probability of rewiring edges (default: 0.1)
#' @param p_intra_community Probability of connecting nodes within the same community (default: 0.5)
#' @param p_inter_community Probability of connecting nodes between different communities (default: 0.1)
#' @param q_degree Scaling factor for the degree-dependent connection probability (default: 0.5)
#' @param verbose Logical indicating whether to print progress messages (default: FALSE)
#' 
#' @return An adjacency matrix representing the graph
spatial_expander_propagation_graph <- function(n, L, radius, edge_weights = TRUE, p_rewire = 0.1, p_intra_community = 0.5, p_inter_community = 0.1, q_degree = 0.5, verbose = FALSE) {
  # Generate nodes using 2D Poisson point process on torus
  if (verbose) cat("Generating nodes...\n")
  coords <- matrix(runif(n * L), n, L)
  coords <- coords * 2 * pi # Map to [0, 2pi]^L
  coords <- apply(coords, 2, function(x) cos(x)) # Map to [-1, 1]^L
  coords <- apply(coords, 2, function(x) (x + 1) * radius) # Map to [0, radius]^L
  
  # Compute pairwise distances and probabilities of connection
  if (verbose) cat("Computing pairwise distances and connection probabilities...\n")
  dists <- as.matrix(dist(coords))
  probs <- exp(-q_degree * dists) # Degree-dependent connection probability
  diag(probs) <- 0 # No self-connections
  
  # Divide nodes into communities based on K-means clustering
  if (verbose) cat("Dividing nodes into communities...\n")
  k <- max(2, round(sqrt(n / 2))) # Choose K based on number of nodes
  comm_ids <- kmeans(coords, centers = k)$cluster
  
  # Connect nodes based on probabilities and community structure
  if (verbose) cat("Connecting nodes...\n")
  adj_matrix <- matrix(0, n, n)
  for (i in 1:n) {
    comm_i <- comm_ids[i]
    intra_comm <- sample.int(n, size = round(p_intra_community * (n - 1)), replace = FALSE)
    intra_comm <- intra_comm[comm_ids[intra_comm] == comm_i]
    inter_comm <- sample.int(n, size = round(p_inter_community * (n - 1)), replace = FALSE)
    inter_comm <- inter_comm[comm_ids[inter_comm] != comm_i]
    all_nodes <- c(intra_comm, inter_comm)
    all_probs <- probs[i, all_nodes]
    all_edges <- rbinom(length(all_probs), 1, all_probs)
    adj_matrix[i, all_nodes] <- all_edges
  }
  
  # Rewire edges with some probability
  if (p_rewire > 0) {
    if (verbose) cat("Rewiring edges...\n")
    for (i in 1:n) {
      for (j in which(adj_matrix[i,] == 1)) {
        if
        
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(Matrix)
library(Matrix.utils)
library(spatstat)

torus_graph <- function(L, N, radius, weight=TRUE, rewiring_prob=0.1, comm_prob=0.7, comm_size=5, seed=NULL) {
  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)
  
  # Generate Poisson point process
  pp <- rpoispp(N, win=owin(rep(0, L), rep(1, L)), cell=rep(radius, L))
  coords <- pp$coords
  
  # Create adjacency matrix
  adj_matrix <- matrix(0, ncol=N, nrow=N)
  
  # Calculate Euclidean distances
  dist_matrix <- dist(coords, method="euclidean")
  
  # Calculate weights
  if (weight) {
    sigma <- mean(dist_matrix)
    weights <- exp(-dist_matrix^2/(2*sigma^2))
  } else {
    weights <- matrix(1, ncol=N, nrow=N)
  }
  
  # Create graph with preferential attachment
  degrees <- rep(0, N)
  communities <- rep(0, N)
  for (i in 1:N) {
    # Connect to high degree nodes with probability proportional to degree
    if (i > 1) {
      weights[i,1:(i-1)] <- weights[i,1:(i-1)]*degrees[1:(i-1)]
      weights[i,1:(i-1)] <- weights[i,1:(i-1)]/sum(weights[i,1:(i-1)])
      for (j in 1:(comm_size-1)) {
        comm_indices <- which(communities == communities[i])
        comm_weights <- weights[i,comm_indices]
        comm_weights[comm_indices == i] <- 0
        comm_weights <- comm_weights*degrees[comm_indices]
        comm_weights <- comm_weights/sum(comm_weights)
        connect_to <- sample(comm_indices, size=1, prob=comm_weights)
        adj_matrix[i,connect_to] <- 1
        adj_matrix[connect_to,i] <- 1
        degrees[i] <- degrees[i] + 1
        degrees[connect_to] <- degrees[connect_to] + 1
      }
      for (j in comm_size:N) {
        if (runif(1) < rewiring_prob) {
          connect_to <- sample(N, size=1, prob=weights[i,])
          adj_matrix[i,connect_to] <- 1
          adj_matrix[connect_to,i] <- 1
          degrees[i] <- degrees[i] + 1
          degrees[connect_to] <- degrees[connect_to] + 1
        } else {
          connect_to <- sample(1:(i-1), size=1, prob=weights[i,1:(i-1)])
          adj_matrix[i,connect_to] <- 1
          adj_matrix[connect_to,i] <- 1
          degrees[i] <- degrees[i] + 1
          degrees[connect_to] <- degrees[connect_to] + 1
        }
      }
    } else {
      degrees[i] <- 1
    }
    
    # Assign community membership
    if (i <= comm_size) {
      communities[i] <- i
    } else {
      if (runif(1) < comm_prob) {
        comm_memberships <- table(communities[1:i-1])
        comm_weights <- comm_memberships^2
        comm_weights
        
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load required packages
library(MASS)
library(Matrix)
library(igraph)

# Function to generate spatial expander propagation graph
# Arguments:
#   L: Number of dimensions
#   n: Number of nodes
#   r: Radius of neighborhood
#   k: Mean degree of the graph
#   p: Probability of rewiring
#   comm: Vector of length n indicating the community membership of each node
#   use_weight: Logical indicating whether to use edge weights
# Output:
#   An adjacency matrix representing the graph
generate_spatial_expander <- function(L, n, r, k, p, comm = NULL, use_weight = TRUE) {
  # Generate Poisson point process on torus
  lambda <- n^(1/L)
  points <- matrix(rpois(n*L, lambda), ncol = L)
  points <- apply(points, 2, function(x) x/(max(x) + 1))
  
  # Compute pairwise distances
  dist_mat <- as.matrix(dist(points, method = "manhattan"))
  
  # Compute probability of connection based on distance and degree
  prob_mat <- exp(-(dist_mat/r)^2)
  if (is.null(comm)) {
    prob_mat <- prob_mat * (k/(n-1))
  } else {
    comm_mat <- outer(comm, comm, "==")
    prob_mat <- prob_mat * ((k*comm_mat + 1)/(sum(comm_mat * (n-1)/sum(comm_mat)) + n - 1))
  }
  
  # Rewire edges with probability p
  if (p > 0) {
    for (i in 1:n) {
      for (j in (i+1):n) {
        if (runif(1) < p) {
          if (is.null(comm) || comm[i] != comm[j]) {
            prob_mat[i,j] <- runif(1)
            prob_mat[j,i] <- prob_mat[i,j]
          }
        }
      }
    }
  }
  
  # Generate random graph from probability matrix
  g <- graph_from_adjacency_matrix(as.matrix(prob_mat), mode = "undirected")
  
  # Set node attributes
  V(g)$x <- points[,1]
  V(g)$y <- points[,2]
  
  # Add edge weights if desired
  if (use_weight) {
    E(g)$weight <- runif(ecount(g))
  }
  
  # Check expansion property
  expand_check(g)
  
  # Return adjacency matrix
  as_adjacency_matrix(g)
}

adj_mat <- generate_spatial_expander(L = 2, n = 100, r = 0.2, k = 6, p = 0.05, use_weight = TRUE)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
    neighbors <- sort(dist_mat[1,])[2:(beta+1)]
    for (j in 1:beta) {
      if (runif(1) < p_rewire) {
        # rewiring with probability p_rewire
        new_neighbor <- sample(1:n, size = 1)
        adj_mat[i, new_neighbor] <- 1
        adj_mat[new_neighbor, i] <- 1
      } else {
        # preferential attachment
        p1 <- neighbors[j] / sum(neighbors) #degee attachment
        connected_nodes <- which(adj_mat[i,] == 1)#connected nodes
        for (k in connected_nodes) {
          p2 <- p1 + beta * (dist_mat[i,k] / sum(dist_mat[i, connected_nodes]))
        } 
        p3 <- p2 / (beta + degree(graph_from_adjacency_matrix(adj_mat), i))
        if (runif(1) < p3) {
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

graph <- spatial_expander(L = 2, n = 100, r = 0.3, beta = 1, p_rewire = .1, attribute = rnorm(100))
G=graph_from_adjacency_matrix(graph$adj_mat,mode="undirected")
G=G%>%simplify(remove.loops = TRUE,remove.multiple = TRUE)
plot(G,vertex.size=2,vertex.label=NA)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Required libraries
library(Matrix)
library(Matrix.utils)
library(stats)

# Function to generate a spatial expander propagation graph on an L-dimensional torus
# Arguments:
# - L: the number of dimensions of the torus (integer)
# - n: the number of nodes in the graph (integer)
# - lambda: the density parameter of the 2D Poisson point process (numeric)
# - alpha: the exponent of the degree distribution (numeric)
# - rewire_prob: the probability of rewiring an edge (numeric)
# - weights: whether to use edge weights or not (logical)
# Returns:
# - A sparse adjacency matrix representing the graph
# - A list of node attributes (coordinates, degrees, communities)
torus_expander <- function(L, n, lambda, alpha, rewire_prob, weights = FALSE) {
  
  # Compute the total number of nodes in the torus
  N <- n^L
  
  # Generate the coordinates of the nodes using a 2D Poisson point process
  coords <- matrix(rpois(N*2, lambda), ncol = 2) / n
  
  # Compute the distance matrix between nodes on the torus
  D <- distanceMatrix(coords, p = 2, torus = rep(TRUE, L))
  
  # Compute the probability of connecting two nodes based on their distance and degree
  prob <- (1 + D)^(-alpha)
  
  # Set the diagonal of the probability matrix to zero
  diag(prob) <- 0
  
  # Apply the rewiring probability to the probability matrix
  prob[runif(N^2) < rewire_prob] <- runif(sum(runif(N^2) < rewire_prob))
  
  # Generate the adjacency matrix using the probability matrix
  adj <- rbinom(N^2, 1, prob)
  
  # Set the diagonal of the adjacency matrix to zero
  diag(adj) <- 0
  
  # Apply the symmetry of the torus to the adjacency matrix
  for (l in 1:L) {
    adj <- apply(adj, l, function(x) { pmax(x, t(x)) })
  }
  
  # Convert the adjacency matrix to a sparse matrix
  adj_sp <- as(adj, "sparseMatrix")
  
  # Compute the degrees of the nodes
  degrees <- colSums(adj_sp)
  
  # Compute the communities of the nodes using the Louvain algorithm
  communities <- fastgreedy.community(adj_sp)$membership
  
  # Normalize the coordinates of the nodes to the range [0, 1)
  coords_norm <- coords / n
  
  # Store the node attributes in a list
  node_attrs <- list(coords = coords_norm, degrees = degrees, communities = communities)
  
  # If edge weights are requested, compute them based on the distance matrix
  if (weights) {
    weights_mat <- exp(-D)
    weights_mat[!adj_sp] <- 0
    weights_sp <- as(weights_mat, "sparseMatrix")
    return(list(adj = weights_sp, node_attrs = node_attrs))
  } else {
    return(list(adj = adj_sp, node_attrs = node_attrs))
  }
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sample codes
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load necessary packages
library(Matrix)
library(spatstat)

# Define the function
construct_torus_expander <- function(n, L, r, lambda, use_weights = TRUE, p_rewire = 0.1, num_communities = 2, community_size = 20) {
  
  # Generate Poisson point process
  points <- rpoispp(lambda, win = owin(rep(0, L), rep(1, L)), nsim = n)
  
  # Compute pairwise distances between points
  dist_matrix <- as.matrix(nndist(points))
  
  # Compute weights based on pairwise distances (optional)
  if (use_weights) {
    weights <- exp(-dist_matrix^2 / r^2)
    adj_matrix <- as.matrix(weights)
  } else {
    adj_matrix <- as.matrix(dist_matrix <= r)
  }
  
  # Ensure graph is connected
  while (!is.connected(Matrix(adj_matrix))) {
    # Re-generate points until graph is connected
    points <- rpoispp(lambda, win = owin(rep(0, L), rep(1, L)), nsim = n)
    dist_matrix <- as.matrix(nndist(points))
    if (use_weights) {
      weights <- exp(-dist_matrix^2 / r^2)
      adj_matrix <- as.matrix(weights)
    } else {
      adj_matrix <- as.matrix(dist_matrix <= r)
    }
  }
  
  # Compute degree distribution
  degree_distribution <- rowSums(adj_matrix)
  
  # Compute probability of connection based on short spatial distance and high degree
  prob_conn <- outer(degree_distribution, degree_distribution, function(x, y) min(x, y) * exp(-dist_matrix^2 / r^2))
  
  # Add rewiring probability for small world effect
  prob_conn <- prob_conn + p_rewire / n
  
  # Add community structure
  comm_indices <- rep(1:num_communities, each = community_size)
  shuffle_indices <- sample.int(n)
  comm_indices <- comm_indices[shuffle_indices]
  comm_membership <- matrix(comm_indices, nrow = n, ncol = n)
  prob_conn <- prob_conn * (comm_membership != t(comm_membership))
  
  # Create adjacency matrix based on probabilities
  adj_matrix <- rbinom(n * n, 1, prob_conn)
  adj_matrix <- matrix(adj_matrix, ncol = n)
  
  # Compute Laplacian matrix
  degree_matrix <- diag(rowSums(adj_matrix))
  laplacian_matrix <- degree_matrix - adj_matrix
  
  # Compute eigenvalues and eigenvectors of Laplacian matrix
  eig <- eigen(Matrix(laplacian_matrix))
  eig_values <- Re(eig$values)
  eig_vectors <- eig$vectors
  
  # Compute the expansion ratio
  lambda_2 <- eig_values[2]
  lambda_n <- eig_values[nrow(laplacian_matrix)]
  expansion_ratio <- lambda_2 / lambda_n
  
  # Create node attributes
  node_attributes <- data.frame(x = points$x, y = points$y, community = comm_indices)
  
  # Return adjacency matrix, Laplacian matrix, eigenvectors, expansion ratio, and node attributes
  return(list(adj_matrix = adj_matrix,
              laplacian_matrix = laplacian_matrix,
              eig_vectors = eig_vectors,
              expansion_ratio = expansion_ratio,
              node_attributes = node_attributes))
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sample codes
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load necessary packages
library(Matrix)
library(spatstat)
library(igraph)

# Define the function
construct_torus_expander <- function(n, L, r, lambda, weighted = FALSE, p_rewire = 0.1, community_size = 5, community_prob = 0.5) {
  
  # Generate Poisson point process
  points <- rpoispp(lambda, win = owin(rep(0, L), rep(1, L)), nsim = n)
  
  # Compute pairwise distances between points
  dist_matrix <- as.matrix(nndist(points))
  
  # Compute degree of each node
  degree <- rowSums(dist_matrix <= r)
  
  # Compute probability of connecting two nodes based on distance and degree
  prob_matrix <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        d_ij <- dist_matrix[i, j]
        k_i <- degree[i]
        k_j <- degree[j]
        prob_matrix[i, j] <- (d_ij <= r) * k_i * k_j
      }
    }
  }
  
  # Generate scale-free graph based on the probability matrix
  g <- sample_pa(n, m = community_size, power = 1, directed = FALSE, zero.appeal = 1, 
          start.graph = NULL, algo = "psumtree", add.vertices = NULL,
          implementation = "psumtree", citation = FALSE, aging = FALSE, 
          aging.preserve = NULL, aging.decay = 1, duplicate.edges = FALSE,
          multiple = FALSE, allow.loops = TRUE, degree.sequence = NULL, 
          suppress.messages = FALSE)
  
  # Add edges with probability proportional to the computed probabilities
  for (i in 1:(n-1)) {
    for (j in (i + 1):n) {
      if (runif(1) < prob_matrix[i, j]) {
        if (!weighted) {
          add_edge(g, i, j)
        } else {
          add_edge(g, i, j, weight = dist_matrix[i, j])
        }
      }
    }
  }
  
  # Rewire edges with small probability
  for (i in 1:(n-1)) {
    for (j in (i + 1):n) {
      if (runif(1) < p_rewire) {
        if (runif(1) < community_prob) {
          # Rewire within community
          neighbors_i <- neighbors(g, i)
          neighbors_j <- neighbors(g, j)
          candidates_i <- setdiff(neighbors_j, neighbors_i)
          candidates_j <- setdiff(neighbors_i, neighbors_j)
          if (length(candidates_i) > 0 && length(candidates_j) > 0) {
            new_j <- sample(candidates_i, 1)
            new_i <- sample(candidates_j, 1)
            delete_edges(g, E(g, i, neighbors_i))
            delete_edges(g, E(g, j, neighbors_j))
            add_edge(g, i, new_i)
            add_edge(g, j, new_j)
          }
        } else {
          # Rewire between communities
          neighbors_i <- neighbors(g, i)
          neighbors_j <- neighbors(g, j)
          candidates_i <- setdiff(V(g), c(i, neighbors_i))
          candidates_j <- setdiff(V(g), c(j, neighbors_j))
          if (
            
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sample codes
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sample codes
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # Load necessary packages
            library(Matrix)
            library(spatstat)
            library(igraph)
            
            # Define the function
            construct_torus_expander <- function(n, L, r, lambda, edge_weight = FALSE, rewire_prob = 0, n_comm = 2, comm_strength = 0.2) {
              
              # Generate Poisson point process
              points <- rpoispp(lambda, win = owin(rep(0, L), rep(1, L)), nsim = n)
              
              # Compute pairwise distances between points
              dist_matrix <- as.matrix(nndist(points))
              
              # Compute degree distribution
              degree_dist <- degree_power_law(n, exponent = 2, gamma = 1)
              degree_dist <- degree_dist / sum(degree_dist)
              
              # Create adjacency matrix based on short spatial distance and degree distribution
              if (edge_weight) {
                adj_matrix <- matrix(0, n, n)
                for (i in 1:n) {
                  # Compute probabilities of connecting to other nodes based on distance and degree
                  probs <- dist_matrix[i, ] < r
                  probs <- probs * degree_dist
                  probs[i] <- 0
                  
                  # Connect to other nodes with probabilities proportional to distance and degree
                  neighbors <- sample.int(n, size = degree_dist[i] * n, replace = TRUE, prob = probs)
                  adj_matrix[i, neighbors] <- runif(length(neighbors)) < comm_strength
                }
              } else {
                adj_matrix <- (dist_matrix < r) * (runif(n^2) < comm_strength)
              }
              
              # Add rewiring probability for small-world effect
              if (rewire_prob > 0) {
                for (i in 1:n) {
                  for (j in (i+1):n) {
                    if (runif(1) < rewire_prob) {
                      adj_matrix[i, j] <- runif(1) < comm_strength
                      adj_matrix[j, i] <- adj_matrix[i, j]
                    }
                  }
                }
              }
              
              # Generate communities
              communities <- rep(1:n_comm, length.out = n)
              if (n_comm > 1) {
                for (i in 1:n) {
                  # Connect to nodes in the same community with higher probability
                  comm_neighbors <- sample(which(communities == communities[i]), size = degree_dist[i] * n_comm, replace = TRUE)
                  adj_matrix[i, comm_neighbors] <- runif(length(comm_neighbors)) < 0.8
                }
              }
              
              # Ensure graph is connected
              while (!is.connected(Matrix(adj_matrix))) {
                # Re-generate points until graph is connected
                points <- rpoispp(lambda, win = owin(rep(0, L), rep(1, L)), nsim = n)
                dist_matrix <- as.matrix(nndist(points))
                if (edge_weight) {
                  adj_matrix <- matrix(0, n, n)
                  for (i in 1:n) {
                    # Compute probabilities of connecting to other nodes based on distance and degree
                    probs <- dist_matrix[i, ] < r
                    probs <- probs * degree_dist
                    probs[i] <- 0
                    
                    # Connect to other nodes with probabilities proportional to distance and degree
                    neighbors <- sample.int(n, size = degree_dist[i] * n, replace = TRUE, prob = probs)
                    adj_matrix[i, neighbors] <- runif(length(neighbors)) < comm_strength
                    
                    

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sample codes
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sample codes
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sample codes
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sample codes
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sample codes
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Load necessary packages
library(Matrix)
library(spatstat)

# Define the function
construct_torus_expander <- function(n, L, r, lambda, p_rewire, n_communities, community_size, alpha, beta, node_attributes) {
  
  # Generate Poisson point process
  points <- rpoispp(lambda, win = owin(rep(0, L), rep(1, L)), nsim = n)
  
  # Compute pairwise distances between points
  dist_matrix <- as.matrix(nndist(points))
  
  # Compute edge weights based on distances
  if (!is.null(node_attributes$weights)) {
    edge_weights <- exp(-alpha * dist_matrix ^ beta)
    adj_matrix <- (edge_weights >= runif(n * (n - 1) / 2))
  } else {
    adj_matrix <- (dist_matrix <= r)
  }
  
  # Ensure graph is connected
  while (!is.connected(Matrix(adj_matrix))) {
    # Re-generate points until graph is connected
    points <- rpoispp(lambda, win = owin(rep(0, L), rep(1, L)), nsim = n)
    dist_matrix <- as.matrix(nndist(points))
    if (!is.null(node_attributes$weights)) {
      edge_weights <- exp(-alpha * dist_matrix ^ beta)
      adj_matrix <- (edge_weights >= runif(n * (n - 1) / 2))
    } else {
      adj_matrix <- (dist_matrix <= r)
    }
  }
  
  # Create scale-free degree distribution
  degree_sequence <- rpareto(n, alpha = alpha, m = 1, xi = 1)
  degree_sequence[degree_sequence > n - 1] <- n - 1
  degree_sequence <- sample(degree_sequence)
  
  # Create communities
  if (n_communities > 1) {
    community_ids <- rep(1:n_communities, each = community_size)
    remaining_nodes <- n - n_communities * community_size
    if (remaining_nodes > 0) {
      community_ids <- c(community_ids, rep(0, remaining_nodes))
    }
    community_ids <- sample(community_ids)
  } else {
    community_ids <- rep(1, n)
  }
  
  # Connect nodes based on degree and community
  for (i in 1:n) {
    degrees_remaining <- degree_sequence[i]
    while (degrees_remaining > 0) {
      # Find nodes in same community with fewest neighbors
      candidates <- which(community_ids == community_ids[i] & rowSums(adj_matrix) < degree_sequence)
      if (length(candidates) == 0) {
        # If no candidates in same community, look for candidates in other communities
        candidates <- which(community_ids != community_ids[i] & rowSums(adj_matrix) < degree_sequence)
      }
      if (length(candidates) > 0) {
        candidate <- sample(candidates, size = 1)
        adj_matrix[i, candidate] <- TRUE
        adj_matrix[candidate, i] <- TRUE
        degrees_remaining <- degrees_remaining - 1
      } else {
        break
      }
    }
  }
  
  # Rewire edges for small world effect
  if (p_rewire > 0) {
    for (i in 1:n) {
      for (j in (i+1):n) {
        if (adj_matrix[i,j]) {
          if (runif(1) < p_re
              
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sample codes
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(Matrix)
library(igraph)

spatial_expander_propagation_graph <- function(n, L, r, lambda, use_edge_weights=TRUE) {
  
  # Generate 2D poisson point process
  coords <- matrix(runif(n*L), n, L)
  radius_sq <- r^2
  
  # Create adjacency matrix
  A <- Matrix(0, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      # Compute shortest distance between points on torus
      dist_sq <- sum(pmin(abs(coords[i,]-coords[j,]), 1-abs(coords[i,]-coords[j,]))^2)
      if (dist_sq <= radius_sq) {
        if (use_edge_weights) {
          # Assign edge weight as inverse of distance
          A[i,j] <- 1/sqrt(dist_sq)
          A[j,i] <- A[i,j]
        } else {
          # Assign binary edge weight
          A[i,j] <- 1
          A[j,i] <- 1
        }
      }
    }
  }
  
  # Compute graph Laplacian
  D <- Matrix::Diagonal(rowSums(A))
  L <- D - A
  
  # Check for strong expansion properties
  spectral_gap <- eigen(Matrix::t(L) %*% L, symmetric=TRUE, only.values=TRUE)$values[n-1]
  if (spectral_gap <= 0) {
    warning("Graph may not exhibit strong expansion properties")
  }
  
  # Convert Laplacian to igraph object
  G <- graph.adjacency(as.matrix(L), mode="undirected")
  
  # Add node attributes
  V(G)$coords <- coords
  V(G)$community <- cutree(cluster_walktrap(G))
  
  return(G)
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(Matrix)
library(igraph)

spatial_expander_propagation_graph <- function(n, L, r, lambda) {
  
  # Generate 2D poisson point process
  coords <- matrix(runif(n*L), n, L)
   radius_sq <- r^2
  
  # Create adjacency matrix with edge weights
if (edge.weight==T){  
  W <- Matrix(0, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      # Compute shortest distance between points on torus
      dist_sq <- sum(pmin(abs(coords[i,]-coords[j,]), 1-abs(coords[i,]-coords[j,]))^2)
      if (dist_sq <= radius_sq) {
        # Assign edge weight as inverse of distance
        W[i,j] <- 1/sqrt(dist_sq)
        W[j,i] <- W[i,j]
      }
    }
  }
}
  #Construct adjacency matrix
  adj_mat <- matrix(0, n, n)
  for (i in 1:n) {
    neighbors <- which(dist_mat[i,] < r)
    neighbors <- neighbors[neighbors != i]
    if (length(neighbors) > deg) {
      neighbors <- sample(neighbors, deg)
    }
    adj_mat[i, neighbors] <- 1
  }
  # Convert adjacency matrix to igraph object
  G <- graph_from_adjacency_matrix(as.matrix(W),mode="undirected")
  
  # Add node attributes
  if (node.attribute==True){
  V(G)$coords <- coords
  V(G)$community <- cutree(cluster_walktrap(G))
  }
  
  # Check for strong expansion properties
  spectral_gap <- eigen(Matrix::t(W) %*% W, symmetric=TRUE, only.values=TRUE)$values[n-1]
  if (spectral_gap <= 0) {
    warning("Graph may not exhibit strong expansion properties")
  }
  
  return(G)
}


# Generate a spatial expander propagation graph with 500 nodes on a 4D torus
G <- spatial_expander_propagation_graph(500, 4, 0.1, 10)

# Plot the graph with node colors based on community structure
plot(G, vertex.color=V(G)$community)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

# Load necessary packages
library(Matrix)
library(spatstat)

# Define the function
construct_torus_expander <- function(n, L, r, lambda) {
  
  # Generate Poisson point process
  points <- rpoispp(lambda, win = owin(rep(0, L), rep(1, L)), nsim = n)
  
  # Compute pairwise distances between points
  dist_matrix <- as.matrix(nndist(points))
  
  # Create adjacency matrix based on short spatial distance
  adj_matrix <- as.matrix(dist_matrix <= r)
  
  # Ensure graph is connected
  while (!is.connected(Matrix(adj_matrix))) {
    # Re-generate points until graph is connected
    points <- rpoispp(lambda, win = owin(rep(0, L), rep(1, L)), nsim = n)
    dist_matrix <- as.matrix(nndist(points))
    adj_matrix <- as.matrix(dist_matrix <= r)
  }
  
  # Compute Laplacian matrix
  degree_matrix <- diag(rowSums(adj_matrix))
  laplacian_matrix <- degree_matrix - adj_matrix
  
  # Compute eigenvalues and eigenvectors of Laplacian matrix
  eig <- eigen(Matrix(laplacian_matrix))
  eig_values <- Re(eig$values)
  eig_vectors <- eig$vectors
  
  # Compute the expansion ratio
  lambda_2 <- eig_values[2]
  lambda_n <- eig_values[nrow(laplacian_matrix)]
  expansion_ratio <- lambda_2 / lambda_n
  
  # Return adjacency matrix, Laplacian matrix, eigenvectors, and expansion ratio
  return(list(adj_matrix = adj_matrix,
              laplacian_matrix = laplacian_matrix,
              eig_vectors = eig_vectors,
              expansion_ratio = expansion_ratio))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(spatstat)

# Function to generate spatial expander propagation graph on L-dimensional torus
# Inputs:
# n: number of nodes
# L: torus dimension
# r: radius for connecting nodes
# lambda: intensity parameter for Poisson point process
# Output: spatial expander propagation graph
generate_spatial_expander_propagation_graph <- function(n, L, r, lambda) {
  # Generate Poisson point process
  pp <- rpoispp(lambda, dim=L, n=n)
  
  # Create empty graph
  g <- make_empty_graph(n)
  i=1;j=2
  # Add edges between nodes within distance r
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      d <- sqrt(min(sum(abs(pp[i,]-pp[j,])^2), sum(abs(pp[i,]-(pp[j,]+L))^2), sum(abs(pp[i,]+L-pp[j,])^2)))
      if (d <= r) {
        g <- add_edges(g, i, j)
      }
    }
  }
  
  # Check if graph is connected
  if (!is.connected(g)) {
    # Find largest connected component
    cc <- largest_component(g)
    # Remove nodes outside largest connected component
    g <- delete_vertices(g, which(!is.element(seq_len(n), cc)))
  }
  
  # Return graph
  return(g)
}

g <- generate_spatial_expander_propagation_graph(n=100, L=3, r=0.1, lambda=10)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(RANN) # for nearest neighbors search

# Generate points on the torus using a 2D Poisson point process
torus_ppp <- function(lambda, L) {
  # generate points in [0,1)^L
  X.points <- matrix(runif(L*lambda), ncol=L)
  # map points to torus
  X.points.torus <- apply(X.points, 2, function(x) x %% 1)
  return(X.points.torus)
}

# Compute pairwise distances between points on torus
dist_torus <- function(X.points.torus) {
  #n <- nrow(X.points.torus)
  D <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      d <- sum(pmin(abs(X.points.torus[i,]-X.points.torus[j,]), 1-abs(X.points.torus[i,]-X.points.torus[j,]))^2)
      D[i,j] <- sqrt(d)
      D[j,i] <- D[i,j]
    }
  }
  return(D)
}

# Construct a spatial expander propagation graph from points on the torus
spatial_expander_propagation_torus <- function(lambda, L, r) {
  # Generate points on torus
  X <- torus_ppp(lambda, L)
  # Compute pairwise distances between points
  D <- dist_torus(X)
  # Create adjacency matrix
  A <- (D <= r)
  diag(A) <- FALSE
  # Check for strong expansion
  vol <- rowSums(A)
  if (min(vol) == 0) stop("Graph is not connected")
  if (any(vol < lambda/2)) stop("Graph does not have strong expansion properties")
  # Return graph
  G <- list(adj = A, n = nrow(A), X = X)
  return(G)
}

# Example usage
set.seed(123)
n=1000
G <- spatial_expander_propagation_torus(lambda = 50, L = 2, r = 0.1)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Spatial Expandergraph
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(RANN) # for nearest neighbors search

# Generate points on the torus using a 2D Poisson point process
torus_ppp <- function(lambda, d) {
  # generate points in [0,1)^d
  X <- matrix(runif(d*lambda), ncol=d)
  # map points to torus
  X <- apply(X, 2, function(x) x %% 1)
  return(X)
}

# Compute pairwise distances between points on torus
dist_torus <- function(X) {
  n <- nrow(X)
  D <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      d <- sum(pmin(abs(X[i,]-X[j,]), 1-abs(X[i,]-X[j,]))^2)
      D[i,j] <- sqrt(d)
      D[j,i] <- D[i,j]
    }
  }
  return(D)
}

# Construct a spatial expander graph from points on the torus
spatial_expander_torus <- function(lambda, d, r) {
  # Generate points on torus
  X <- torus_ppp(lambda, d)
  # Compute pairwise distances between points
  D <- dist_torus(X)
  # Create adjacency matrix
  A <- (D <= r)
  diag(A) <- FALSE
  # Check for strong expansion
  vol <- rowSums(A)
  if (min(vol) == 0) stop("Graph is not connected")
  if (any(vol < lambda/2)) stop("Graph does not have strong expansion properties")
  # Return graph
  G <- list(adj = A, n = nrow(A))
  return(G)
}

# Example usage
set.seed(123)
G <- spatial_expander_torus(lambda = 500, d = 4, r = 0.5)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Step 1: Generate Poisson point process in 2D
library(spatstat)
n <- 100 # number of points
r <- 0.1 # radius
lambda <- 10 # intensity parameter

# Generate Poisson point process in 2D
ppp <- rpoispp(lambda, win = square(1))

# Plot the Poisson point process
plot(ppp)

#Step 2: Map the 2D points onto an L-dimensional torus

# Define the L-dimensional torus
torus <- hyperframe(x = list(dim=rep(1,L), units=rep("mm",L)), type = "tore")

# Map the 2D points onto the L-dimensional torus
ppp.torus <- im.apply(ppp, function(x) { x <- x %% torus; x })

#Step 3: Construct the spatial expander propagation graph
# Compute the Delaunay triangulation of the mapped points
tri <- triangulate.ppp(ppp.torus)

# Compute the Gabriel graph from the Delaunay triangulation
gab <- gabriel.graph(tri)

# Construct the spatial expander propagation graph
spg <- propagation.graph(gab)

# Plot the spatial expander propagation graph
plot(spg, vertex.col = "blue", vertex.cex = 0.5, edge.col = "gray", main = "Spatial Expander Propagation Graph")
# Step 4: Tune the parameters to generate different spatial expander graphs
# You can tune the number of nodes by changing the value of n, the radius by
# changing the value of r, the intensity parameter by changing the value of lambda,
# and the torus dimension by changing the value of L. Note that the values of n, r, 
# and lambda may need to be adjusted depending on the value of L to generate 
# a well-connected spatial expander graph.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Spatial expander graph
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(deldir)

# Define parameters
d <- 2 # dimension
n <- 50 # number of points
radius <- 0.3 # radius of Poisson disk in torus
epsilon <- 0.5 # expansion parameter

# Generate Poisson point process on torus

# coords <- matrix(runif(n*d), n, d)
# Generate the 2D Poisson point process
set.seed(123)
lambda <- n^(1/d) #intensity parameter
points <- matrix(runif(n*d), ncol=d)
points <- points - floor(points) #wrap around torus
points <- points[runif(n) < lambda, ]

# Compute the euclidean (pairwise) distances between points
distances <- as.matrix(dist(points))
# coords <- coords %% 1 # wrap around torus
# coords <- coords * (1-2*radius) + radius # shift to avoid boundary effects
#del <- deldir(coords[,1], coords[,2], rw=c(0,1,0,1)) # Delaunay triangulation

# Compute adjacency matrix using Delaunay triangulation
# adj_mat <- matrix(0, n, n)
# for (i in 1:(n-1)) {
#   for (j in (i+1):n) {
#     if (length(which(del$summary$cdist[which(del$summary$nc == i & del$summary$ne == j)] == 1)) > 0) {
#       adj_mat[i,j] <- 1
#       adj_mat[j,i] <- 1
#     }
#   }
# }

# Generate the adjacency matrix and connects points to k-nearest neighbours
adjacency <- matrix(0,n,n)
for (i in 1:n) {
  neighbors <- which(distances[i,] < r)
  neighbors <- neighbors[neighbors != i]
  if (length(neighbors) > deg) {
    neighbors <- sample(neighbors, deg)
  }
  adjacency[i,neighbors] <- 1
}


# Compute Laplacian matrix
deg<- rowSums(adjacency)
laplacian <- diag(deg) - adjacency

# Compute eigenvectors and eigenvalues
eigen <- eigen(laplacian)
eigenvectors <- eigen$vectors[,2:(d+1)]
eigenvalues <- eigen$values[2:(d+1)]

# Apply expander mixing
expander_coords <- coords + epsilon * eigenvectors %*% diag(sqrt(eigenvalues))

# Wrap expander_coords back into the torus
expander_coords <- expander_coords %% 1

###---Compute adjacency matrix using Delaunay triangulation on expanded coordinates

# del <- deldir(expander_coords[,1], expander_coords[,2], rw=c(0,1,0,1))
# expander_adj_mat <- matrix(0, n, n)
# for (i in 1:(n-1)) {
#   for (j in (i+1):n) {
#     if (length(which(del$summary$cdist[which(del$summary$nc == i & del$summary$ne == j)] == 1)) > 0) {
#       expander_adj_mat[i,j] <- 1
#       expander_adj_mat[j,i] <- 1
#     }
#   }
# }

# Plot results
par(mfrow=c(1,2))
plot(coords[,1], coords[,2], pch=19, main="Original Poisson process on torus")
plot(expander_coords[,1], expander_coords[,2], pch=19, main="Spatial expander graph on torus")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Step 1: Generate a set of n random points on a d-dimensional torus using a 2D Poisson
#point process. To generate such a set of points, we can use the R package 'spatstat'. 
#First, we need to define the window of the torus in d-dimensions.
#For example, if we want to create a torus in 2D space, we can use the following code:

library(spatstat)
W <- torus(dimension=2, r.major=1, r.minor=0.5)
#Here, r.major and r.minor are the major and minor radii of the torus, respectively.
#We can adjust these parameters according to our preference.

#To generate a set of n random points in this window using a Poisson point process, 
#we can use the rpoispp function in 'spatstat':
  
n <- 100  # number of points
X <- rpoispp(W, n=n)

#Step 2: Compute the pairwise distances between all pairs of points. We can use the dist function in R to compute the Euclidean distance matrix.


D <- dist(X$coords)

#Step 3: Create the adjacency matrix of the graph. We can set a threshold value r and 
#connect two points if their Euclidean distance is less than or equal to r. For example, we can set r = 0.2.

r <- 0.2
A <- as.matrix(D <= r)

#Step 4: Make the graph connected. We can use the Breadth First Search (BFS) algorithm to find a connected component of the graph. We can start the BFS from a randomly chosen point.

set.seed(123)
start <- sample(1:n, 1)
visited <- rep(FALSE, n)
visited[start] <- TRUE
queue <- list(start)

while (length(queue) > 0) {
  v <- queue[[1]]
  queue <- queue[-1]
  neighbors <- which(A[v,])
  for (i in neighbors) {
    if (!visited[i]) {
      visited[i] <- TRUE
      queue <- c(queue, i)
    }
  }
}

if (sum(visited) < n) {
  # graph is not connected, find another component
  other_vertices <- which(!visited)
  start <- sample(other_vertices, 1)
  visited_other <- rep(FALSE, n)
  visited_other[start] <- TRUE
  queue <- list(start)
  
  while (length(queue) > 0) {
    v <- queue[[1]]
    queue <- queue[-1]
    neighbors <- which(A[v,])
    for (i in neighbors) {
      if (!visited_other[i]) {
        visited_other[i] <- TRUE
        queue <- c(queue, i)
      }
    }
  }
  
  # connect the two components
  A[start,] <- visited | visited_other
  A[,start] <- visited | visited_other
}

#Step 5: Check the expansion properties of the graph. A graph is said to have strong expansion properties
#if the second eigenvalue of its normalized Laplacian matrix is bounded away from zero.
#We can compute the second eigenvalue using the eigen function in R.
A=adjacency
D1 <- diag(rowSums(A))
L <- D1 - A
Lnorm <- solve(sqrt(D1)) %*% L %*% solve(sqrt(D1))
ev <- eigen$values

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(Matrix)

# Define the number of nodes and the dimension of the torus
n <- 100
d <- 2

# Generate the Poisson point process on the torus
lambda <- 10
X <- matrix(runif(n*d, 0, 1), ncol=d)
for (i in 1:d) {
  X[,i] <- X[,i] * (2*pi)
}
Y <- matrix(rpois(n, lambda), ncol=1)
nodes <- matrix(0, nrow=sum(Y), ncol=d)
index <- 1
for (i in 1:n) {
  for (j in 1:Y[i]) {
    nodes[index,] <- X[i,]
    index <- index + 1
  }
}

# Compute the distances between nodes on the torus
distances <- matrix(0, nrow=n, ncol=n)
for (i in 1:n-1) {
  for (j in (i+1):n) {
    d1 <- abs(nodes[i,] - nodes[j,])
    d2 <- 2*pi - d1
    distances[i,j] <- sqrt(sum(pmin(d1, d2)^2))
    distances[j,i] <- distances[i,j]
  }
}

# Compute the Cheeger constant
vol <- sum(Y)
vol_left <- 0
vol_right <- vol
perm <- order(Y, decreasing=TRUE)
h <- Inf
for (i in 1:(n-1)) {
  vol_left <- vol_left + Y[perm[i]]
  vol_right <- vol_right - Y[perm[i]]
  cut_size <- sum(distances[perm[1:i], perm[(i+1):n]])
  conductance <- cut_size / min(vol_left, vol_right)
  if (conductance < h) {
    h <- conductance
    perm_star <- perm[1:i]
  }
}

# Construct the adjacency matrix of the expander graph
A <- matrix(0, nrow=n, ncol=n)
for (i in 1:n) {
  for (j in perm_star) {
    if (i != j) {
      d1 <- abs(nodes[i,] - nodes[j,])
      d2 <- 2*pi - d1
      d <- sqrt(sum(pmin(d1, d2)^2))
      if (d <= h/2) {
        A[i,j] <- 1
        A[j,i] <- 1
      }
    }
  }
}

# Construct the combinatorial Laplacian matrix of the expander graph
D <- diag(rowSums(A))
L <- D - A

# Check the expansion properties of the expander graph
lambda_2 <- eigen(L, symmetric=TRUE, only.values=TRUE)$values[2]
if (lambda_2 > 0.1*h^2) {
  print("The expander graph has strong expansion properties.")
} else {
  print("The expander graph does not have strong expansion properties.")
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(Matrix)

# Define the number of nodes and the dimension of the torus
n <- 100
d <- 3

# Define the parameters of the Poisson point process
lambda <- 5

# Define the distance function on the torus
distance <- function(x, y) {
  abs(x - y)
}

# Initialize the list of edges and the number of attempts
edges <- list()
attempts <- 0

# Generate the nodes using a Poisson point process on the torus
nodes <- matrix(rpois(n*d, lambda), ncol=d)
for (i in 1:n) {
  for (j in 1:d) {
    nodes[i,j] <- nodes[i,j] %% lambda
  }
}

# Compute the initial Laplacian matrix
D <- matrix(0, nrow=n, ncol=n)
W <- matrix(0, nrow=n, ncol=n)
for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      d_ij <- apply(cbind(nodes[i,], nodes[j,]), 1, distance)
      W[i,j] <- exp(-sum(d_ij^2)/2)
      D[i,i] <- D[i,i] + W[i,j]
    }
  }
}
L <- Diagonal(n) - W

# Perform the Metropolis-Hastings algorithm to generate the edges
while (length(edges) < 3*n*log(n)) {
  # Select two random nodes
  i <- sample(1:n, 1)
  j <- sample(1:n, 1)
  if (i == j) {
    next
  }
  # Compute the difference in Laplacian energy
  dL <- 2*(W[i,j] - W[i,])*L[i,i] + 2*(W[j,i] - W[j,])*L[j,j] - 4*W[i,j]*L[i,j]
  # Accept or reject the move with probability min(1, exp(-dL))
  if (runif(1) < min(1, exp(-dL))) {
    if (!i %in% edges[[j]]) {
      edges[[i]] <- c(edges[[i]], j)
      edges[[j]] <- c(edges[[j]], i)
    }
  }
  # Increment the number of attempts
  attempts <- attempts + 1
}

# Convert the list of edges to a sparse matrix
row_indices <- rep(seq_len(n), sapply(edges, length))
col_indices <- unlist(edges)
adj_matrix <- sparseMatrix(i=row_indices, j=col_indices, dims=c(n, n))

# Compute the combinatorial Laplacian matrix
D_adj <- Matrix::Diagonal(rowSums(adj_matrix))
L_adj <- D_adj - adj_matrix

# Check the expansion properties
eigenvalues <- eigen(L_adj)$values
lambda2 <- eigenvalues[2]
print(paste("The second smallest eigenvalue is", lambda2))
expansion_ratio <- lambda2 / lambda
print(paste("The expansion ratio is", expansion_ratio))

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Define the number of nodes and the dimension of the torus
n <- 100
d <- 3

# Define the intensity of the Poisson point process
lambda <- n^(1/d)

# Define the Metropolis-Hastings proposal distribution
sigma <- 0.1

# Define the function to compute the acceptance probability
acceptance_prob <- function(x, y, eps) {
  # Compute the distance between the points
  d <- min(sqrt(sum((x - y + eps)^2)), sqrt(sum((x - y - eps)^2)), sqrt(sum((x - y)^2)))
  # Compute the acceptance probability
  if (d <= eps) {
    return(1)
  } else {
    return(exp(-d/eps))
  }
}

# Initialize the current state of the Markov chain
X <- matrix(runif(n*d), ncol=d)

# Initialize the adjacency matrix
W <- matrix(0, nrow=n, ncol=n)

# Run the Metropolis-Hastings algorithm
for (i in 1:1000) {
  # Propose a new state
  Y <- X + rnorm(n*d, 0, sigma)
  Y <- Y %% 1
  # Compute the acceptance probability
  a <- acceptance_prob(X[1,], Y[1,], 0.1)
  # Accept or reject the proposal
  if (runif(1) < a) {
    X <- Y
  }
  # Update the adjacency matrix
  for (j in 1:n) {
    for (k in 1:n) {
      if (j != k) {
        d <- min(sqrt(sum((X[j,] - X[k,] + 0.1)^2)), sqrt(sum((X[j,] - X[k,] - 0.1)^2)), sqrt(sum((X[j,] - X[k,])^2)))
        if (d <= 0.1) {
          W[j,k] <- 1
          W[k,j] <- 1
        }
      }
    }
  }
}

# Compute the degree matrix and the Laplacian matrix
D <- diag(rowSums(W))
L <- D - W

# Compute the eigenvalues and eigenvectors of the Laplacian matrix
e <- eigen(L)
lambda2 <- e$values[2]
v2 <- e$vectors[,2]

# Define the function to build the graph
build_graph <- function(X, v2, lambda2) {
  # Compute the projection of the nodes on the second eigenvector
  p <- X %*% v2
  # Compute the threshold for the projection
  eps <- sqrt(-log(lambda2/n))
  # Compute the adjacency matrix
  W <- matrix(0, nrow=n, ncol=n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        d <- abs(p[i] - p[j])
        if (d <= eps) {
          W[i,j] <- 1
        }
      }
    }
  }
  # Create an empty adjacency list
  G <- vector("list", n)
  for (i in 1:n) {
    G[[i]] <- integer(0)
  }
  # Add edges to the adjacency list
  for (i in 1:n) {
    for (j in 1:n) {
      if
      

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(stats)

# Define the number of nodes and the dimension of the torus
n <- 100
d <- 2
eps <- 0.1

# Define the intensity of the Poisson point process
lambda <- n^(1/d)

# Generate the Poisson point process
X <- matrix(rpois(lambda*n, 1), ncol=d)

# Rescale the points to the torus
X <- X - 1
X <- X / n
X <- X %% 1

# Compute the distance matrix
D <- as.matrix(dist(X, method = "euclidean"))

# Build the adjacency matrix
W <- (D <= eps)

# Ensure that the graph is connected
while (!is.connected(W)) {
  eps <- eps + 0.1
  W <- (D <= eps)
}

# Convert the adjacency matrix to an edgelist
E <- which(W, arr.ind = TRUE)
E <- E[E[, 1] < E[, 2], ]

# Define the function to build the graph
build_graph <- function(E) {
  # Create an empty adjacency list
  G <- vector("list", n)
  for (i in 1:n) {
    G[[i]] <- integer(0)
  }
  # Add edges to the adjacency list
  for (e in E) {
    i <- e[1]
    j <- e[2]
    G[[i]] <- c(G[[i]], j)
    G[[j]] <- c(G[[j]], i)
  }
  return(G)
}

# Build the graph
G <- build_graph(E)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Define the number of nodes and the dimension of the torus
n <- 100
d <- 3

# Define the intensity of the Poisson point process
lambda <- n^(1/d)

# Define the Metropolis-Hastings proposal distribution
sigma <- 0.1

# Define the function to compute the acceptance probability
acceptance_prob <- function(x, y, eps) {
  # Compute the distance between the points
  d <- min(sqrt(sum((x - y + eps)^2)), sqrt(sum((x - y - eps)^2)), sqrt(sum((x - y)^2)))
  # Compute the acceptance probability
  if (d <= eps) {
    return(1)
  } else {
    return(exp(-d/eps))
  }
}

# Initialize the current state of the Markov chain
X <- matrix(runif(n*d), ncol=d)

# Initialize the adjacency matrix
W <- matrix(0, nrow=n, ncol=n)

# Run the Metropolis-Hastings algorithm
for (i in 1:1000) {
  # Propose a new state
  Y <- X + rnorm(n*d, 0, sigma)
  Y <- Y %% 1
  # Compute the acceptance probability
  a <- acceptance_prob(X[1,], Y[1,], 0.1)
  # Accept or reject the proposal
  if (runif(1) < a) {
    X <- Y
  }
  # Update the adjacency matrix
  for (j in 1:n) {
    for (k in 1:n) {
      if (j != k) {
        d <- min(sqrt(sum((X[j,] - X[k,] + 0.1)^2)), sqrt(sum((X[j,] - X[k,] - 0.1)^2)), sqrt(sum((X[j,] - X[k,])^2)))
        if (d <= 0.1) {
          W[j,k] <- 1
          W[k,j] <- 1
        }
      }
    }
  }
}

# Convert the adjacency matrix to an edgelist
E <- which(W, arr.ind = TRUE)
E <- E[E[, 1] < E[, 2], ]

# Define the function to build the graph
build_graph <- function(E) {
  # Create an empty adjacency list
  G <- vector("list", n)
  for (i in 1:n) {
    G[[i]] <- integer(0)
  }
  # Add edges to the adjacency list
  for (e in E) {
    i <- e[1]
    j <- e[2]
    G[[i]] <- c(G[[i]], j)
    G[[j]] <- c(G[[j]], i)
  }
  return(G)
}

# Build the graph
G <- build_graph(E)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(igraph)

# Define parameters
N <- 1000  # number of vertices
m <- 4     # number of edges to attach to new vertex in Barabasi-Albert model
p <- 0.01  # probability of adding a random edge
K <- 5     # number of communities

# Step 1: Randomly distribute N points in a 2D space using a Poisson Point Process
coords <- matrix(rpois(N*2, 10), ncol=2)
g <- graph.empty(n=N, directed=FALSE)
V(g)$x <- coords[,1]
V(g)$y <- coords[,2]


# Step 2: Assign degrees to each vertex using the Barabasi-Albert model
g=igraph::barabasi.game(N, m, directed=FALSE)
g

# Step 3: Add random edges with a small probability p
g1=add_edges(g, get.edgelist(g)[sample(length(E(g)), round(p*length(E(g))))])


# Step 4: Create community structures using stochastic block model
# First, assign each vertex to a random community
community <- c(table(sample(1:K, N, replace=TRUE)))
# Then, create a block matrix with higher probability of edges within each community
B <- matrix(0.1, K, K) + 0.8*diag(K)
g2 <- sample_sbm(N, B, community,directed=FALSE)

# Step 5: Add edges between vertices in different communities to create an expander graph
# First, calculate the normalized Laplacian matrix of the graph
L <- laplacian_matrix(g2)
D <- diag(vcount(g2))
Dinv <- solve(D)
Lnorm <- Dinv %*% L
# Then, use spectral partitioning to find a cut that minimizes edge expansion
spectral_partition <- spectral_partition(Lnorm, K)
cut_edges <- get.cut_edges(g, spectral_partition$membership)
# Finally, add the cut edges to the graph
g <- add.edges(g, cut_edges)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sample code 13
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

generate_torus_expander <- function(n, m, r, p_short, p_long, beta) {
  # n: number of nodes to add
  # m: number of edges to add for each new node
  # r: radius of torus
  # p_short: probability of connecting to nearest neighbors
  # p_long: probability of connecting to nodes with distance greater than 1
  # beta: parameter controlling degree distribution (larger values give more power-law-like distributions)
  
  # Initialize graph
  x <- rnorm(n, 0, r)
  y <- rnorm(n, 0, r)
  adj <- matrix(0, n, n)
  degrees <- rep(0, n)
  
  # Helper function to compute shortest distance on torus
  torus_dist <- function(x1, y1, x2, y2, r) {
    dx <- abs(x1 - x2)
    dy <- abs(y1 - y2)
    dx <- min(dx, r - dx)
    dy <- min(dy, r - dy)
    sqrt(dx^2 + dy^2)
  }
  
  # Add nodes to graph
  for (i in 1:n) {
    # Compute distances to existing nodes
    dists <- sapply(1:(i-1), function(j) torus_dist(x[i], y[i], x[j], y[j], r))
    
    # Connect to nearest neighbors with probability p_short
    dists_short <- dists[dists <= 1]
    if (length(dists_short) > 0) {
      nearest <- which(dists <= 1)
      p <- p_short / sum(dists_short^beta)
      probs <- p / dists_short^beta
      probs[is.infinite(probs)] <- 0
      probs <- probs / sum(probs)
      to_connect <- sample(nearest, size = m, replace = TRUE, prob = probs)
      adj[i, to_connect] <- 1
      adj[to_connect, i] <- 1
      degrees[i] <- sum(adj[i,])
      degrees[to_connect] <- degrees[to_connect] + 1
    }
    
    # Connect to more distant nodes with probability p_long
    dists_long <- dists[dists > 1]
    if (length(dists_long) > 0) {
      distant <- which(dists > 1)
      p <- p_long / sum(dists_long^beta)
      probs <- p / dists_long^beta
      probs[is.infinite(probs)] <- 0
      probs <- probs / sum(probs)
      to_connect <- sample(distant, size = m, replace = TRUE, prob = probs)
      adj[i, to_connect] <- 1
      adj[to_connect, i] <- 1
      degrees[i] <- sum(adj[i,])
      degrees[to_connect] <- degrees[to_connect] + 1
    }
    
    x[i] <- runif(1, 0, r)
    y[i] <- runif(1, 0, r)
  }
  
  # Generate degree distribution with power-law tail
  deg_seq <- numeric(sum(degrees))
  k <- 1
  for (i in 1:n) {
    deg_seq[k:(k+degrees[i]-1)] <- rep(i, degrees[i])
    k <- k + degrees[i]
  }
  deg_seq <- sort(deg_seq, decreasing = TRUE)
  deg_seq <- deg
  
  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sample code 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
library(Matrix)
library(RANN)
  
  # Set parameters for the spatial point process
  n_points <- 1000 # number of nodes to add
  lambda <- 1 # density of the Poisson point process
  
  # Set parameters for the graph growth
  p_degree <- 0.1 # probability of connecting to existing vertices based on degree
  p_distance <- 0.9 # probability of connecting to existing vertices based on spatial distance
  r_distance <- 0.1 # radius of spatial neighborhood for connecting based on distance
  
  # Generate spatial positions of nodes using Poisson point process on torus
  x <- runif(n_points)
  y <- runif(n_points)
  xy <- cbind(x, y)
  theta <- 2*pi*(x+y)
  xy_torus <- cbind(cos(theta), sin(theta))
  
  # Initialize adjacency matrix and degrees vector
  adj_mat <- Matrix(0, n_points, n_points)
  degrees <- rep(0, n_points)
  
  # Add nodes to the graph one at a time
  for (i in 1:n_points) {
    # Connect to existing vertices based on degree and distance
    neighbors <- which(degrees > 0)
    if (length(neighbors) > 0) {
      distances <- RANN::nn2(xy[i,], xy[neighbors,], k = 1)$nn.dist
      p_distance_neigh <- exp(-distances/r_distance)
      p_degree_neigh <- degrees[neighbors]/sum(degrees[neighbors])
      p_neighbors <- p_degree*p_degree_neigh + p_distance*p_distance_neigh
      neighbors <- sample(neighbors, size = 5, prob = p_neighbors)
      adj_mat[i,neighbors] <- 1
      adj_mat[neighbors,i] <- 1
    }
    
    # Add node to graph and update degrees
    degrees[i] <- sum(adj_mat[i,])
  }
  
  # Convert adjacency matrix to a sparse matrix for efficiency
  adj_mat <- as(adj_mat, "dgCMatrix")
  
  # Calculate Laplacian and eigenvectors of Laplacian
  L <- Diagonal(n_points) - adj_mat
  eig <- eigen(L, symmetric = TRUE)
  
  # Extract eigenvectors corresponding to smallest eigenvalues
  n_eigenvectors <- 10 # number of eigenvectors to use
  eigenvectors <- eig$vectors[,1:n_eigenvectors]
  
  # Apply threshold to eigenvectors to get binary indicators of community membership
  community_threshold <- 0.1 # threshold for binary indicator
  communities <- apply(eigenvectors, 1, function(x) ifelse(x > community_threshold, 1, 0))
  
  # Calculate degree distribution and average clustering coefficient
  degree_dist <- table(degrees)/n_points
  avg_cc <- mean(transitivity(adj_mat))
  
  # Plot the graph with nodes colored by community membership
  plot(x, y, pch = 19, col = communities+1)
  
  
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Spatial Expander graphs
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # Load required packages
  library(Matrix)
  library(Matrix.utils)
  library(gtools)
  
  # Define parameters
  n <- 1000      # Number of nodes
  r <- 0.2       # Radius of torus
  k <- 10        # Number of neighbors to connect to
  p <- 0.2       # Probability of rewiring edges
  
  # Define function to generate torus coordinates
  generate_coordinates <- function(n, r) {
    # Generate polar coordinates
    theta <- runif(n, 0, 2 * pi)
    rho <- sqrt(runif(n, 0, r^2))
    
    # Convert to Cartesian coordinates
    x <- rho * cos(theta)
    y <- rho * sin(theta)
    
    return(cbind(x, y))
  }
  
  # Define function to calculate pairwise distances on torus
  torus_distance <- function(coords, r) {
    # Calculate Euclidean distances
    d <- dist(coords)
    
    # Wrap around torus edges
    d <- min(d, abs(d - 2 * r))
    
    return(as.matrix(d))
  }
  
  # Define function to connect nodes
  connect_nodes <- function(adj, coords, k, p) {
    # Calculate distance matrix
    d <- torus_distance(coords, r)
    
    # Loop over nodes
    for (i in 1:n) {
      # Find k nearest neighbors
      neighbors <- as.vector(order(d[i, ])[2:(k+1)])
      
      # Rewire edges with probability p
      for (j in neighbors) {
        if (runif(1) < p) {
          # Choose random node to connect to
          new_neighbor <- sample(setdiff(1:n, c(i, neighbors)), 1)
          
          # Connect nodes
          adj[i, j] <- 0
          adj[j, i] <- 0
          adj[i, new_neighbor] <- 1
          adj[new_neighbor, i] <- 1
        }
      }
    }
    
    return(adj)
  }
  
  # Define function to create expander graph on torus
  create_torus_graph <- function(n, r, k, p) {
    # Generate coordinates
    coords <- generate_coordinates(n, r)
    
    # Initialize adjacency matrix
    adj <- matrix(0, n, n)
    
    # Connect nodes
    adj <- connect_nodes(adj, coords, k, p)
    
    # Ensure graph is connected
    adj <- ensure.connected(adj)
    
    # Convert to sparse matrix format
    adj_sparse <- Matrix(adj, sparse=TRUE)
    
    return(adj_sparse)
  }
  
  
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Spatial Expander graphs
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  library(Matrix)
  library(stats)
  
  # Function to compute distance between two points on a torus
  torus_distance <- function(p1, p2, width, height) {
    dx <- min(abs(p1[1] - p2[1]), width - abs(p1[1] - p2[1]))
    dy <- min(abs(p1[2] - p2[2]), height - abs(p1[2] - p2[2]))
    return(sqrt(dx^2 + dy^2))
  }
  
  # Parameters
  n <- 1000   # number of nodes to add
  width <- 100   # width of torus
  height <- 100  # height of torus
  lambda <- 5    # parameter for Poisson point process
  alpha <- 0.5   # parameter for preferential attachment
  beta <- 0.5    # parameter for community structure
  
  # Generate Poisson point process for node locations
  n_positions <- rpois(n, lambda)
  x_positions <- runif(sum(n_positions), 0, width)
  y_positions <- runif(sum(n_positions), 0, height)
  positions <- cbind(x_positions, y_positions)
  
  # Initialize adjacency matrix
  adj_mat <- Matrix(0, n, n)
  
  # Add edges between nodes
  for (i in 1:n) {
    # Compute distances to all other nodes
    distances <- apply(positions, 1, function(p) torus_distance(p, positions[i,], width, height))
    # Compute probabilities of connecting to each node
    probs <- distances^(-alpha) * (degree(adj_mat) + 1)^beta
    # Normalize probabilities
    probs <- probs / sum(probs)
    # Sample edges based on probabilities
    edges <- sample(1:n, size=min(10, n), replace=FALSE, prob=probs)
    # Add edges to adjacency matrix
    adj_mat[i, edges] <- 1
  }
  
  # Add community structure
  community_size <- round(n / 10)
  community_members <- rep(1:n, each=community_size)
  community_members <- sample(community_members)
  for (i in 1:n) {
    for (j in i+1:n) {
      if (community_members[i] == community_members[j]) {
        # Add edge within community
        if (runif(1) < 0.5) {
          adj_mat[i,j] <- 1
          adj_mat[j,i] <- 1
        }
      } else {
        # Add edge between communities
        if (runif(1) < 0.05) {
          adj_mat[i,j] <- 1
          adj_mat[j,i] <- 1
        }
      }
    }
  }
  
  # Convert adjacency matrix to graph
  g <- graph.adjacency(adj_mat, mode="undirected")
  
  # Compute layout for visualization
  layout <- layout.fruchterman.reingold(g)
  
  # Plot graph
  plot(g, layout=layout)
  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Spatial Expander graphs
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  library(spatstat)
  library(igraph)
  
  # Parameters
  N <- 1000 # Number of points
  m <- 10 # Number of edges to attach to each new vertex in Barabasi-Albert model
  p <- 0.1 # Probability of adding a random edge
  ncommunities <- 7 # Number of communities
  community_size <- 200 # Number of vertices in each community
  
  # Spatial distribution of points
  x <- runifpoint(N)
  
  #add.vertices(N,g)
  
  ##---Make an empty graph and Assign degrees using Barabasi-Albert model
  ## and rewire with rewiring prop p for small world effect
 # g <- make_empty_graph(N,directed = FALSE)%>%
 # g <- make_empty_graph(N,directed = FALSE)
  g1=barabasi.game(N,power=2,m, directed = FALSE,algorithm="psumtree")%>%
  rewire(each_edge(p = .1, loops = FALSE)) 
  #degree <- degree(g1)
  
  # Add random edges(rewiring)
  #g2 <- add.edges(g1,m, p = p, verbose = FALSE)
   
  # Divide vertices into communities
  communities <- sample(rep(1:ncommunities, each = community_size), N, replace = TRUE)
  
  ##Connect vertices within each community
  within_community_prob <- 0.5
  for (i in 1:ncommunities) {
    community_vertices <- which(communities == 5)
    within_community_edges <- sample(community_vertices, size = length(community_vertices) * (length(community_vertices) - 1) * within_community_prob / 2, replace = T)
    g2 <- add.edges(g1, within_community_edges, verbose = FALSE)
  }
  
  # Create expander graph by adding edges between different communities
  cutsize <- 0.1 * N * (ncommunities - 1)
  adjacency_matrix <- as.matrix(as_adjacency_matrix(g2))
  spectral_partition <- spantree.partition(adjacency_matrix, cutsize)
  between_community_edges <- spectral_partition$cut
  g3 <- add.edges(g2, between_community_edges, verbose = FALSE)
  
  # Visualize the graph
  plot(g3, vertex.col = communities, vertex.cex = 0.5, main = "Expander Graph with Scale-Free Degree Distribution, Small World Properties, and Community Structures")
  
  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Spatial Expander graphs
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  

  # Algorithm:
  #   
  # 1.Initialize the graph with a small number of vertices randomly distributed on the torus.
  # 
  # 2.Add new vertices to the graph one by one,
  # following a 2D Poisson point process (this defines the spatial position).
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
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Spatial Expander graphs
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # # Set seed for reproducibility
  # set.seed(123)
  # 
  # # Parameters
  # n <- 1000  # number of nodes
  # r <- 0.05  # node density (i.e., expected number of nodes per unit area)
  # gamma <- 2.5  # exponent for power-law degree distribution
  # kmin <- 5  # minimum degree
  # kmax <- 50  # maximum degree
  # theta <- pi/2  # angle for edge probability kernel
  # beta <- 0.01  # scaling parameter for edge probability kernel
  # 
  # # Generate node locations using 2D Poisson point process
  # x <- runif(n)
  # y <- runif(n)
  # lambda <- r*n
  # m <- rpois(n, lambda)
  # while (sum(m) == 0) {
  #   # Ensure at least one node is generated
  #   m <- rpois(n, lambda)
  # }
  # pos <- data.frame(x=rep(x, m), y=rep(y, m))
  # 
  # # Calculate pairwise distances and create adjacency matrix
  # dmat <- as.matrix(dist(pos))
  # dmat[dmat > 0.5] <- 1 - dmat[dmat > 0.5]  # wraparound distance on torus
  # p <- exp(-beta*dmat^theta)  # edge probability kernel
  # p[lower.tri(p)] <- 0  # remove duplicate edges
  # g <- matrix(rbinom(n*(n-1)/2, 1, p), nrow=n, ncol=n)
  # g <- g + t(g)  # make undirected
  # diag(g) <- 0  # remove self-loops
  # 
  # # Ensure graph is connected
  # cc <- c()
  # while (length(cc) != 1) {
  #   cc <- clusters(g, mode="weak")$membership
  #   if (length(cc) != 1) {
  #     # Remove all but the largest connected component
  #     idx <- which(cc != which.max(table(cc)))
  #     g[idx, ] <- 0
  #     g[, idx] <- 0
  #   }
  # }
  # 
  # # Generate degree sequence with power-law distribution
  # dseq <- rpl(n, gamma, kmin, kmax)
  # dseq <- dseq[order(dseq, decreasing=TRUE)]
  # while (sum(dseq) %% 2 != 0) {
  #   # Ensure even number of edges for each node
  #   dseq[which.max(dseq)] <- dseq[which.max(dseq)] - 1
  # }
  # 
  # # Generate graph by connecting nodes based on degree sequence
  # v <- 1
  # while (sum(dseq) > 0) {
  #   if (dseq[v] > 0) {
  #     # Connect node to random nodes with high degree and short spatial distance
  #     d <- sum(g[v, ])
  #     while (d < dseq[v]) {
  #       nbrs <- which(g[v, ] == 1)
  #       nbrs <- nbrs[which.max(dseq[nbrs])]  # select neighbor with highest degree
  #       if (length(nbrs) > 0) {
  #         # Choose neighbor with highest edge probability
  #         p <- exp(-beta*dmat[v, nbrs]^theta)
  #         p[nbrs] <- 0  # remove existing edges
  #         if (sum(p) > 0) {
  #           nbr <- sample(1:n, size=1, prob=p/sum(p))
  #           g[v, nbr] <- 1
  #           g[nbr, v]
  
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Spatial Expander graphs
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Set parameters
  n <- 1000 # number of nodes
  r <- 0.1 # radius of the torus
  k <- 10 # average degree
  p_rewire <- 0.1 # probability of rewiring
  n_communities <- 5 # number of communities
  p_within_community <- 0.3 # probability of connection within a community
  
  # Generate nodes based on a 2D Poisson point process on the torus
  x <- runif(n)
  y <- runif(n)
  x <- r*cos(2*pi*x)
  y <- r*sin(2*pi*y)
  
  # Calculate pairwise distances between nodes
  dist_mat <- as.matrix(dist(cbind(x, y)))
  dist_mat[dist_mat > r/2] <- r - dist_mat[dist_mat > r/2] # wrap around torus
  
  # Connect each node to its k nearest neighbors
  adj_mat <- matrix(0, n, n)
  for (i in 1:n) {
    neighbors <- sort(dist_mat[i, ])[2:(k+1)]
    for (j in which(dist_mat[i, ] <= max(neighbors))) {
      adj_mat[i, j] <- 1
    }
  }
  
  # Rewire edges to achieve small-world property
  for (i in 1:n) {
    for (j in which(adj_mat[i, ] == 1)) {
      if (runif(1) < p_rewire) {
        # Choose a random node to rewire to
        potential_neighbors <- which(adj_mat[i, ] == 0 & dist_mat[i, ] <= r/2)
        rewire_to <- sample(potential_neighbors, 1)
        # Rewire the edge
        adj_mat[i, j] <- 0
        adj_mat[j, i] <- 0
        adj_mat[i, rewire_to] <- 1
        adj_mat[rewire_to, i] <- 1
      }
    }
  }
  
  # Group nodes into communities
  communities <- cutree(cluster_walktrap(graph.adjacency(adj_mat)), n_communities)
  
  # Connect nodes within each community with higher probability
  for (i in unique(communities)) {
    nodes <- which(communities == i)
    n_nodes <- length(nodes)
    for (j in 1:n_nodes) {
      for (k in (j+1):n_nodes) {
        if (runif(1) < p_within_community) {
          adj_mat[nodes[j], nodes[k]] <- 1
          adj_mat[nodes[k], nodes[j]] <- 1
        }
      }
    }
  }
  
  # Check the expansion property of the graph
  eigenvalues <- eigen(adj_mat)$values
  second_eigenvalue <- sort(eigenvalues, decreasing = TRUE)[2]
  lambda <- min(second_eigenvalue, 2*sqrt(k*log(n)/n))
  
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Spatial Expander graphs
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  library(Matrix)
  library(RANN)
  
  #Step 1: Set up the parameters for the graph
  
  n <- 100 # number of nodes
  k <- 4   # average degree
  m <- 2   # number of communities
  p <- 0.3 # rewiring probability
  l <- 0.1 # probability of connecting within communities
  t <- 0.1 # probability of connecting between communities
  
  # Here, we specify the number of nodes n, the average degree k, the number of communities m, the rewiring probability p, and the probabilities of connecting within and between communities l and t.
  # 
  # Step 2: Generate spatial positions for the nodes
  
  x <- runif(n, 0, 1)
  y <- runif(n, 0, 1)
  pos <- cbind(x, y)
  
  # We generate n random values between 0 and 1 for the x and y coordinates of the nodes, and store them in the pos matrix.
  # Step 3: Construct the torus
  dist_x <- abs(outer(x, x, "-"))
  dist_y <- abs(outer(y, y, "-"))
  dist_x <- pmin(dist_x, 1 - dist_x)
  dist_y <- pmin(dist_y, 1 - dist_y)
  dist <- sqrt(dist_x^2 + dist_y^2)
  
  # We calculate the distance between every pair of nodes on the torus, taking into account the periodic boundary conditions. We store the distances in the dist matrix.
  # Step 4: Connect the nodes
  adj <- matrix(0, n, n)
  for (i in 1:n) {
    # Find k nearest neighbors
    knn <- knnsearch(pos[-i,], pos[i,], k)
    # Add edges to nearest neighbors
    for (j in knn) {
      if (runif(1) < p) {
        adj[i, j] <- 1
        adj[j, i] <- 1
      }
    }
  }
  
  
  # We loop over every node and connect it to its k nearest neighbors. We use a k-d tree to efficiently find the nearest neighbors. We add an edge between the nodes with a probability of p, allowing for rewiring. We store the resulting adjacency matrix in the adj matrix.
  # 
  # Step 5: Add community structure
  
  communities <- sample(rep(1:m, each = n/m))
  for (i in 1:n) {
    for (j in (i+1):n) {
      if (communities[i] == communities[j]) {
        if (runif(1) < l * exp(-dist[i, j]^2 / (2 * t^2))) {
          adj[i, j] <- 1
          adj[j, i] <- 1
        }
      } else {
        if (runif(1) < t * exp(-dist[i, j]^2 / (2 * t^2))) {
          adj[i, j] <- 1
          adj[j, i] <- 1
        }
      }
    }
  }
  
  
  # We first assign each node to a community using a random assignment. We then loop over every pair of nodes and connect them with a probability based on their community membership and their spatial distance. We use the l and t parameters to specify the probabilities of connecting 
  # within and between communities, respectively. We use a Gaussian kernel to favor
  
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Spatial Expander graphs
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  torus_model <- function(n, d, gamma, mu, k, c) {
    # n = number of nodes to generate
    # d = number of dimensions for the torus
    # gamma = power-law exponent for degree distribution
    # mu = spatial decay parameter
    # k = number of nearest neighbors to connect
    # c = probability of rewiring
    
    # Generate initial lattice
    coords <- matrix(rep(0, n*d), ncol = d)
    for (i in 1:n) {
      coords[i, ] <- sample(1:n^(1/d), d, replace = TRUE)
    }
    
    # Generate edges
    adj_matrix <- matrix(0, nrow = n, ncol = n)
    for (i in 1:n) {
      # Generate random points within distance mu of node i
      dists <- apply(coords - coords[i, ], 1, function(x) min(abs(x), n - abs(x)))
      nearby <- which(rowSums(dists <= mu) > 0)
      # Connect node i to its k nearest neighbors
      if (length(nearby) > k) {
        dists <- dists[nearby]
        nearest <- nearby[order(dists)[1:k]]
        adj_matrix[i, nearest] <- 1
      } else {
        adj_matrix[i, nearby] <- 1
      }
    }
    
    # Rewire edges with probability c
    for (i in 1:n) {
      for (j in which(adj_matrix[i, ] == 1)) {
        if (runif(1) < c) {
          # Choose random node to connect to
          possible_nodes <- setdiff(1:n, c(i, j, which(adj_matrix[, j] == 1)))
          if (length(possible_nodes) > 0) {
            new_j <- sample(possible_nodes, 1)
            adj_matrix[i, j] <- 0
            adj_matrix[j, i] <- 0
            adj_matrix[i, new_j] <- 1
            adj_matrix[new_j, i] <- 1
          }
        }
      }
    }
    
    # Assign degrees according to power-law distribution
    degrees <- rpl(n, gamma)
    
    # Add new nodes with Poisson point process
    new_coords <- matrix(rpois(n*d, 1/mu), ncol = d)
    new_degrees <- rpl(n, gamma)
    coords <- rbind(coords, new_coords)
    adj_matrix <- cbind(adj_matrix, matrix(0, nrow = n, ncol = n))
    adj_matrix <- rbind(adj_matrix, matrix(0, nrow = n, ncol = n + n))
    for (i in (n + 1):(n + nrow(new_coords))) {
      # Connect to existing nodes proportional to their degree
      probs <- degrees / sum(degrees)
      connect_to <- sample(1:n, new_degrees[i - n], replace = TRUE, prob = probs)
      adj_matrix[i, connect_to] <- 1
      adj_matrix[connect_to, i] <- 1
    }
    
    # Connect new nodes to torus with nearest neighbor edges
    # for (i in (n + 1):(n + nrow(new_coords))) {
    #   dists <- apply(coords - coords[i, ], 1, function(x) min(abs(x), n - abs(x)))
    #   nearby <- which(rowSums(dists <= mu) > 
    #                     
    
    
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Spatial Expander graphs
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    # Function to create a d-dimensional torus with a scale-free degree distribution, small-world structure, and community structure
    generate_torus <- function(d, n, m, p, q) {
      # d: dimension of the torus
      # n: number of nodes in the torus
      # m: number of edges to attach from a new node to existing nodes
      # p: probability of rewiring each edge
      # q: number of communities
      
      # Create a torus lattice with n nodes
      coords <- matrix(0, n, d)
      for (i in 1:d) {
        coords[,i] <- rep(seq(0, 1, length.out = sqrt(n)), each = sqrt(n))
      }
      
      # Initialize the adjacency matrix
      adj <- matrix(0, n, n)
      
      # Create a vector to store the degrees of nodes
      degrees <- rep(0, n)
      
      # Add nodes to the torus one by one
      for (i in 1:n) {
        # Generate a spatial position for the new node using a 2D Poisson point process
        new_coord <- runif(d, 0, 1)
        coords <- rbind(coords, new_coord)
        
        # Calculate the distances between the new node and existing nodes
        distances <- apply(coords, 1, function(x) sqrt(sum((new_coord - x)^2)))
        distances[n+1] <- Inf
        
        # Choose the m nodes with the shortest distances and connect the new node to them with a probability proportional to their degree
        neighbors <- order(degrees)[1:m]
        probs <- degrees[neighbors] / sum(degrees[neighbors])
        for (j in neighbors) {
          if (runif(1) < probs[j - n]) {
            adj[n+1, j] <- 1
            adj[j, n+1] <- 1
            degrees[n+1] <- degrees[n+1] + 1
            degrees[j] <- degrees[j] + 1
          }
        }
        
        # Rewire each edge with probability p
        for (j in (n+1):(n+m)) {
          if (runif(1) < p) {
            # Choose a random node to rewire the edge to
            old_neighbor <- sample(1:n, 1)
            # Check that the new neighbor is not already connected to the node
            if (adj[j, old_neighbor] == 0) {
              adj[j, neighbors[j-(n+1)+1]] <- 0
              adj[neighbors[j-(n+1)+1], j] <- 0
              adj[j, old_neighbor] <- 1
              adj[old_neighbor, j] <- 1
            }
          }
        }
      }
      
      # Create communities by randomly assigning nodes to q groups
      groups <- sample(1:q, n, replace = TRUE)
      for (i in 1:n) {
        for (j in 1:n) {
          if (groups[i] == groups[j] && adj[i,j] == 1) {
            adj[i,j] <- 0
            adj[j,i] <- 0
          }
        }
      }
      
      # Return the adjacency matrix and the coordinates of nodes
      return(list(adj = adj, coords = coords))
    }
    
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Spatial Expander graph
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    expandertorus <- function(n, d, p, q) {
      # Generate a torus graph
      adj_mat <- matrix(0, n, n)
      for (i in 1:n) {
        neighbors <- c()
        for (j in 1:d) {
          neighbors <- c(neighbors, (i + c(rep(1, j-1), 0, rep(-1, d-j))) %% n + 1)
          neighbors <- c(neighbors, (i - c(rep(1, j-1), 0, rep(-1, d-j))) %% n + 1)
        }
        adj_mat[i, neighbors] <- 1
      }
      degrees <- rowSums(adj_mat)
      
      # Generate a random point process
      points <- matrix(runif(n*2), n, 2)
      
      # Add new nodes to the graph
      for (i in 1:(n-1)) {
        # Compute distances to existing nodes
        distances <- apply(points[1:i, ], 1, function(x) sqrt((points[i+1, 1]-x[1])^2 + (points[i+1, 2]-x[2])^2))
        # Compute probabilities proportional to the degrees of existing nodes
        probs <- p * degrees[1:i] / sum(degrees[1:i]) + q * exp(-distances)
        # Normalize probabilities
        probs <- probs / sum(probs)
        # Choose a node to connect to
        new_neighbor <- sample(1:i, 1, prob=probs)
        adj_mat[i+1, new_neighbor] <- 1
        adj_mat[new_neighbor, i+1] <- 1
        degrees[c(i+1, new_neighbor)] <- degrees[c(i+1, new_neighbor)] + 1
      }
      
      # Scale-free degree distribution
      deg_seq <- sample(degrees, n, replace=TRUE)
      
      # Small world properties
      for (i in 1:n) {
        for (j in (i+1):n) {
          if (adj_mat[i, j] == 0) {
            if (runif(1) < p * deg_seq[i] * deg_seq[j] / n) {
              adj_mat[i, j] <- 1
              adj_mat[j, i] <- 1
            }
          }
        }
      }
      
      # Community structure
      for (i in 1:n) {
        for (j in (i+1):n) {
          if (adj_mat[i, j] == 0) {
            if (runif(1) < q) {
              if ((points[i, 1] < 0.5 && points[j, 1] < 0.5) || (points[i, 1] >= 0.5 && points[j, 1] >= 0.5)) {
                adj_mat[i, j] <- 1
                adj_mat[j, i] <- 1
              }
            }
          }
        }
      }
      
      return(adj_mat)
    }
    
    expandertorus(10, 2, 0.2, 0.1)
    
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Spatial Expander graph
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    expandTorus <- function(n, d, r, p, q) {
      # n: number of nodes
      # d: number of dimensions in torus
      # r: radius of neighborhood for small world property
      # p: probability of rewiring for small world property
      # q: fraction of edges within communities
      
      # Initialize the adjacency matrix and list of communities
      adj_matrix <- matrix(0, n, n)
      communities <- list()
      
      # Generate the first node and add to community 1
      nodes <- data.frame(matrix(runif(d), nrow = 1))
      communities[[1]] <- 1
      
      # Add remaining nodes
      for (i in 2:n) {
        # Generate new node according to 2D poison point process
        new_node <- matrix(rpois(d, lambda = i), nrow = 1)
        
        # Calculate distances to all existing nodes in torus
        distances <- apply(abs(nodes - new_node), 1, function(x) min(x, i - x))
        
        # Determine neighbors for expander property
        neighbors <- order(distances)[2:(1 + log2(i))]
        
        # Connect new node to neighbors proportional to their degree for scale-free property
        neighbor_degrees <- colSums(adj_matrix[neighbors, ])
        neighbor_probs <- neighbor_degrees / sum(neighbor_degrees)
        new_edges <- sample(neighbors, 1, prob = neighbor_probs)
        
        # Connect new node to a fraction of edges within its community
        community <- communities[[which.max(distances)]]
        community_nodes <- which(communities == community)
        community_edges <- adj_matrix[community_nodes, community_nodes]
        community_degrees <- colSums(community_edges)
        community_probs <- community_degrees / sum(community_degrees)
        community_edges_to_add <- sample(which(community_probs > 0), round(q * sum(community_probs)), replace = FALSE)
        new_edges <- c(new_edges, community_nodes[community_edges_to_add])
        
        # Rewire edges for small world property
        for (j in new_edges) {
          if (runif(1) < p) {
            non_neighbors <- setdiff(1:n, c(j, new_edges))
            new_edge <- sample(non_neighbors, 1)
            adj_matrix[j, new_edge] <- 1
            adj_matrix[new_edge, j] <- 1
            adj_matrix[j, which(new_edges == j)] <- 0
            adj_matrix[which(new_edges == j), j] <- 0
            new_edges[which(new_edges == j)] <- new_edge
          } else {
            adj_matrix[j, i] <- 1
            adj_matrix[i, j] <- 1
          }
        }
        
        # Add new node to nodes data frame and community list
        nodes <- rbind(nodes, new_node)
        communities[[i]] <- community
      }
      
      # Return the adjacency matrix and list of communities
      return(list(adj_matrix = adj_matrix, communities = communities))
    }
    
    
    result <- expandTorus(n = 100, d = 2, r = 5, p = 0.1, q = 0.1)
    
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Spatial Expander graph
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # ###[1]---Step 1: Set up the environment
    # #Load required libraries
    # library(Matrix)
    # library(scales)
    # library(ggplot2)
    # library(gridExtra)
    # 
    # # Set the seed for reproducibility
    # set.seed(123)
    # 
    # ###[2]--Step 2: Define a function to generate spatial points on a d-dimensional torus using a 2D 
    # #Poisson point process. This function will generate a matrix of n points in d dimensions,
    # #such that each point lies within the torus.
    # 
    # generate_points <- function(n, d, torus_radius) {
    #   # Generate n points in d dimensions within the torus of radius torus_radius
    #   u <- matrix(runif(n * d, -torus_radius, torus_radius), ncol = d)
    #   norms <- apply(u, 1, function(x) sqrt(sum(x^2)))
    #   return(u * torus_radius / matrix(norms, ncol = d, nrow = n))
    # }
    # 
    # ###[3]--Step 3:Define a function to calculate the distance matrix between two sets of points
    # #on a d-dimensional torus. This function will use the distance formula for the torus.
    # 
    # distance_matrix <- function(X, Y, torus_radius) {
    #   # Calculate the distance matrix between X and Y on the torus
    #   X <- X %% torus_radius
    #   Y <- Y %% torus_radius
    #   d <- ncol(X)
    #   n <- nrow(X)
    #   m <- nrow(Y)
    #   XY <- matrix(0, n, m)
    #   for (i in 1:d) {
    #     dxy <- abs(X[,i][,drop = TRUE] - t(Y[,i]))
    #     dxy <- pmin(dxy, torus_radius - dxy)
    #     XY <- XY + dxy^2
    #   }
    #   return(sqrt(XY))
    # }
    
    
    ###[3]--Step 4: Define a function to generate an expander graph with the mentioned properties. This function will take the following arguments: n (number of nodes), d (dimensionality of the torus), torus_radius (radius of the torus),
    #m (number of edges per node), p (rewiring probability), and q (community parameter)
    # expander_graph <- function(n, d, torus_radius, m, p, q) {
    #   # Generate an expander graph with n nodes, d-dimensional torus, m edges per node,
    #   # rewiring probability p, and community parameter q.
    #   
    #   # Generate n points on the d-dimensional torus using a 2D Poisson point process
    #   X <- generate_points(n, d, torus_radius)
    #   
    #   # Compute the distance matrix between the points
    #   D <- distance_matrix(X, X, torus_radius)
    #   
    #   # Initialize the adjacency matrix and degree vector
    #   A <- Matrix(0, nrow = n, ncol = n)
    #   degree <- numeric(n)
    #   
    #   # Add m edges to each node
    #   for (i in 1:n) {
    #     # Get the indices of the m nearest neighbors
    #     neighbors <- order(D[i,])[2:(m+1)]
    #     
    #     # Add edges to the nearest neighbors
    #     for (j in neighbors) {
    #       # Check if edge already exists
    #       if (A[i,j] == 0) {
    #         # Add edge with probability proportional to degree
    #         A[i,j] <- rbinom(1, 1, min(1, q * degree[j] / sum(degree)))
    #         A[j,i] <- A[i,j]
    #         if (A[i,j] == 1) {
    #           
    
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Spatial Expander graph
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # library(powerlaw)
    # 
    # expandGraph <- function(n, d, p) {
    #   # Step 1: initialize the graph with one node
    #   graph <- matrix(0, nrow = 1, ncol = d)
    #   
    #   for (i in 2:n) {
    #     # Step 2: generate new node coordinates using 2D Poisson point process
    #     new_node <- matrix(rpois(d, 10), nrow = 1, ncol = d)
    #     
    #     # Step 3: calculate distances to existing nodes
    #     distances <- apply(graph, 1, function(x) sum((x - new_node)^2))
    #     
    #     # Step 4: connect new node to k existing nodes with smallest distances
    #     k <- rdiscrete1(length(distances), power.law.dispersion = 1)
    #     neighbors <- order(distances)[1:k]
    #     graph <- rbind(graph, new_node)
    #     for (j in neighbors) {
    #       graph[j, ] <- (graph[j, ] + new_node)/2
    #     }
    #     
    #     # Step 6: rewire a connection with probability p
    #     if (runif(1) < p) {
    #       node <- sample(1:n, 1)
    #       neighbor <- sample(1:k, 1)
    #       graph[node, neighbor] <- sample(1:n, 1)
    #     }
    #   }
    #   
    #   # Step 7: introduce community structure using Girvan-Newman algorithm
    #   g <- graph.adjacency(graph, mode = "undirected")
    #   c <- clusters(g)
    #   num_communities <- length(unique(c$membership))
    #   for (i in 1:num_communities) {
    #     subgraph <- induced.subgraph(g, which(c$membership == i))
    #     num_nodes <- vcount
    #   }
    
    
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Spatial Expander graph
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    #Step 1: Create a function to generate 2D Poisson point process with a given
    #intensity parameter
    
    ppp <- function(lambda, x_max, y_max) {
      n <- rpois(1, lambda * x_max * y_max)
      x <- runif(n, 0, x_max)
      y <- runif(n, 0, y_max)
      return(data.frame(x, y))
    }
    
    #Step 2: Create a function to generate a d-dimensional torus with given dimensions 
    #and periodic boundary conditions
    
    torus <- function(dims) {
      n <- prod(dims)
      adj_mat <- matrix(0, n, n)
      for (i in 1:n) {
        coord_i <- as.numeric(apply(ind2sub(dims, i), 2, function(x) x - 1))
        for (j in 1:n) {
          if (i == j) {
            next
          }
          coord_j <- as.numeric(apply(ind2sub(dims, j), 2, function(x) x - 1))
          dist <- pmin(pmin(abs(coord_i - coord_j), dims - abs(coord_i - coord_j)))
          if (sum(dist > 1) == 0) {
            adj_mat[i, j] <- 1
          }
        }
      }
      return(adj_mat)
    }
    
    #Step 3: Create a function to generate an expander graph from a given adjacency matrix using the spectral method
    
    expander <- function(adj_mat) {
      d <- colSums(adj_mat)
      laplacian <- diag(d) - adj_mat
      eigen_vals <- eigen(laplacian)$values
      eig_gap <- min(diff(sort(eigen_vals)))
      scaling_factor <- sqrt(eig_gap / 2)
      norm_laplacian <- (2 / eig_gap) * laplacian
      spectral_gap <- 1 - max(abs(eigen(norm_laplacian, only.values = TRUE)$values))
      norm_adj_mat <- norm_laplacian / spectral_gap
      return(norm_adj_mat)
    }
    
    #Step 4: Create a function to generate a scale-free network with small-world and community structures using the Barabasi-Albert model
    
    barabasi_albert <- function(n, m, k) {
      #  Initialize graph with m nodes
      G <- matrix(0, m, m)
      for (i in 1:m) {
        for (j in (i + 1):m) {
          G[i, j] <- 1
          G[j, i] <- 1
        }
      }
      degrees <- rep(m, m)
      
      #  Add n-m nodes to the graph using preferential attachment
      
      for (i in (m + 1):n) {
        G <- rbind(G, rep(0, i))
        G <- cbind(G, rep(0, i))
        prob <- degrees / sum(degrees)
        neighbors <- sample(1:(i-1), m, replace = FALSE, prob = prob)
        for (j in neighbors) {
          G[i, j] <- 1
          G[j, i] <- 1
          degrees[i] <- degrees[i] + 1
          degrees[j] <- degrees[j] + 1
        }
      }
      
      #  Add small-world and community structures
      for (i in 1:n) {
        for (j in 1:n) {
          if (i == j) {
            next
          }
          if (G[i, j] == 0 && rbinom(1, 1, k / degrees[i]) == 1) {
            G
          }}}}
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Spatial Expander graph
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # This function takes four parameters:
    # n: the number of nodes to generate
    # d: the number of dimensions of the torus
    # lambda: the Poisson process parameter, which controls the average number of connections 
    # each new node makes
    # p: the preferential attachment parameter, which controls the probability of connecting to
    # existing nodes based on their degree
    # The function also has two internal parameters:
    # q: the rewiring parameter, which controls how strongly the preferential attachment is skewed
    # towards high-degree nodes. A value of 0 corresponds to uniform attachment, and higher values 
    # correspond to more strongly skewed attachment.
    # small_world: a Boolean flag that controls whether to add small-world and community
    # structure to the graph. If TRUE, the function adds edges to nodes that are within a 
    # certain distance of each other, with probability decaying exponentially with distance. 
    # This promotes clustering and short path lengths, and creates small-world and community 
    # structure in the graph.
    # To generate a torus expander graph with all the specified properties, you can call the 
    # torus_expander function with appropriate values of the parameters. For example:
    
    
    
    
    
    
    
    
    torus_expander <- function(n, d, lambda, p, q){
      # n: number of nodes
      # d: number of dimensions
      # lambda: Poisson process parameter
      # p: preferential attachment parameter
      # q: rewiring parameter
      
      # Generate a random set of points on the d-dimensional torus
      points <- matrix(runif(n * d, 0, 1), ncol = d)
      
      # Compute pairwise distances between points
      dists <- as.matrix(dist(points))
      
      # Generate an empty adjacency matrix
      adj <- matrix(0, n, n)
      
      # Connect nodes according to Poisson process and preferential attachment
      for (i in 1:n) {
        # Compute probabilities of connecting to each node
        probs <- p * (colSums(adj) / sum(colSums(adj))) + (1 - p) * (1 / (1 + dists[i,]^q))
        # Sample from a Poisson distribution to determine number of connections
        k <- rpois(1, lambda)
        # Sample nodes to connect to based on probabilities
        to_connect <- sample(1:n, k, replace = TRUE, prob = probs)
        # Connect nodes
        adj[i, to_connect] <- 1
        adj[to_connect, i] <- 1
      }
      
      # Add small-world and community structure
      for (i in 1:n) {
        for (j in (i+1):n) {
          if (adj[i,j] == 0) {
            # Compute distance between nodes
            dist_ij <- dists[i,j]
            # Compute probability of connecting
            if (dist_ij <= 2) {
              prob_ij <- exp(-dist_ij / 0.1)
            } else if (dist_ij <= 4) {
              prob_ij <- exp(-dist_ij / 0.2)
            } else {
              prob_ij <- exp(-dist_ij / 0.3)
            }
            # Connect nodes with probability prob_ij
            if (runif(1) < prob_ij) {
              adj[i,j] <- 1
              adj[j,i] <- 1
            }
          }
        }
      }
      
      # Return adjacency matrix
      return(adj)
    }
    
    adj <- torus_expander(n = 100, d = 2, lambda = 2,p=2,q=0.1)
    
    
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Spatial Expander graph
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    #To generate an expander graph on a torus with the properties you mentioned, we can follow the following steps:
    #--Step 1: Initialize the graph with a seed node.
    
    graph <- list()
    graph$edges <- matrix(0, ncol = 2, nrow = 1)
    graph$edges[1, ] <- c(1, 1)
    graph$degrees <- rep(1, 1)
    
    
    ##--Step 2: Define a function to add a new node to the graph.
    
    add_node <- function(graph) {
      # Add a new node
      new_node <- nrow(graph$edges) + 1
      
      # Generate the Poisson process for the new node's position
      lambda <- 5 # Mean number of points per unit area
      num_points <- rpois(1, lambda)
      points <- matrix(runif(num_points*2, 0, 1), ncol=2, byrow=TRUE)
      
      # Compute the distances between the new node and existing nodes
      distances <- apply(graph$edges, 1, function(edge) {
        p1 <- graph$positions[edge[1], ]
        p2 <- graph$positions[edge[2], ]
        min(c(norm(p1 - p2), norm(p1 - (p2 + c(1, 0))), norm(p1 - (p2 + c(0, 1))), norm(p1 - (p2 + c(1, 1)))))
      })
      
      # Compute the probabilities of connecting to each existing node
      degrees <- graph$degrees
      degrees[degrees == 0] <- 1 # Avoid dividing by zero
      probs <- degrees / sum(degrees)
      probs <- probs^2 # Preferential attachment
      
      # Add edges to existing nodes
      for (i in 1:num_points) {
        # Choose an existing node to connect to
        node <- sample(length(probs), size=1, prob=probs)
        
        # Add an edge to the chosen node
        graph$edges <- rbind(graph$edges, c(new_node, node))
        graph$degrees[node] <- graph$degrees[node] + 1
      }
      
      # Add the new node to the graph
      graph$edges <- rbind(graph$edges, c(new_node, new_node))
      graph$degrees <- c(graph$degrees, num_points + 1)
      graph$positions <- rbind(graph$positions, points)
      
      return(graph)
    }
    
    ##--Step 3: Generate the graph by iteratively adding nodes.
    
    num_nodes <- 100 # Number of nodes in the graph
    graph <- add_node(graph)
    for (i in 2:num_nodes) {
      graph <- add_node(graph)
    }
    
    ##--Step 4: Visualize the graph on a torus using the plot function
    
    # Compute the edges on the torus
    edges <- graph$edges
    num_edges <- nrow(edges)
    for (i in 1:num_edges) {
      p1 <- graph$positions[edges[i, 1], ]
      p2 <- graph$positions[edges[i, 2], ]
      edges[i, ] <- c(project(p1), project(p2))
    }
    
    # Plot the graph
    plot(0, 0, type="n", xlim=c(0, 1), ylim=c(0, 1))
    segments(edges[, 1], edges[, 2], edges[, 3], edges[, 4], col="gray", lwd=2)
    
    
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Spatial Expander graph
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # Function to generate an expander graph on a torus
    # Inputs:
    #   n: number of vertices in the graph
    #   k: number of edges each vertex should have
    #   r: radius of the torus
    #   beta: parameter for small world properties
    #   mu: parameter for community structure
    # Output: adjacency matrix of the generated graph
    generate_torus_expander <- function(n, k, r, beta, mu) {
      
      # Create a matrix to store the positions of the nodes
      positions <- matrix(0, n, 2)
      
      # Generate node positions using a 2D Poisson process
      lambda <- n/(pi*r^2) # Poisson rate
      for (i in 1:n) {
        theta <- runif(1, 0, 2*pi)
        rho <- sqrt(r^2*runif(1))
        positions[i,] <- c(rho*cos(theta), rho*sin(theta))
      }
      
      # Calculate distances between each pair of nodes
      distances <- as.matrix(dist(positions))
      
      # Create an adjacency matrix
      adj_matrix <- matrix(0, n, n)
      
      # Add edges to each node using a scale-free degree distribution
      degrees <- numeric(n)
      degrees[1:k] <- k # Initial connected graph
      for (i in (k+1):n) {
        # Calculate probabilities for attaching to each existing node
        probs <- rep(0, n)
        for (j in 1:(i-1)) {
          d_ij <- distances[i,j]
          probs[j] <- (degrees[j]+beta)/(2*k) * exp(-mu*d_ij)
        }
        probs <- probs/sum(probs)
        # Choose k nodes to connect to
        connections <- sample(1:(i-1), k, replace = FALSE, prob = probs)
        for (j in connections) {
          adj_matrix[i,j] <- 1
          adj_matrix[j,i] <- 1
          degrees[i] <- degrees[i] + 1
          degrees[j] <- degrees[j] + 1
        }
      }
      
      return(adj_matrix)
    }
    
    
    # Generate and plot the graph
    set.seed(123) # for reproducibility
    adj_matrix <- generate_torus_expander(100, 6, 5, 0.2, 0.1)
    plot(graph.adjacency(adj_matrix))
    
    
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Spatial Expander graph
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    poisson_torus <- function(lambda, torus_width, torus_height) {
      # generate number of points using Poisson distribution
      num_points <- rpois(1, lambda)
      # generate random coordinates on torus
      x <- runif(num_points, 0, torus_width)
      y <- runif(num_points, 0, torus_height)
      # adjust coordinates to account for torus topology
      x <- x %% torus_width
      y <- y %% torus_height
      # return coordinates as a data frame
      data.frame(x = x, y = y)
    }
    
    
    torus_distance <- function(x1, y1, x2, y2, torus_width, torus_height) {
      dx <- abs(x1 - x2)
      dy <- abs(y1 - y2)
      dx <- min(dx, torus_width - dx)
      dy <- min(dy, torus_height - dy)
      return(sqrt(dx^2 + dy^2))
    }
    
    
    torus_expander <- function(num_vertices, lambda, torus_width, torus_height, alpha, beta, num_communities) {
      # create initial vertex set with one random point on torus
      vertices <- data.frame(x = runif(1, 0, torus_width), y = runif(1, 0, torus_height))
      # create initial edge set
      edges <- data.frame(from = integer(), to = integer())
      # grow graph by adding new vertices
      for (i in 2:num_vertices) {
        # generate Poisson point process
        points <- poisson_torus(lambda, torus_width, torus_height)
        # calculate distances between new points and existing vertices
        distances <- matrix(0, nrow = nrow(points), ncol = nrow(vertices))
        for (j in 1:nrow(points)) {
          for (k in 1:nrow(vertices)) {
            distances[j,k] <- torus_distance(points[j,1], points[j,2], vertices[k,1], vertices[k,2], torus_width, torus_height)
          }
        }
        # calculate probabilities of attaching to existing vertices
        probabilities <- (distances)^(-alpha)
        rowSums(probabilities) <- 0
        probabilities <- probabilities / rowSums(probabilities)
        # add edges to graph with preferential attachment
        for (j in 1:nrow(points)) {
          selected <- sample.int(nrow(vertices), 1, prob = probabilities[j,])
          edges <- rbind(edges, data.frame(from = i, to = selected))
        }
        # add new vertex to vertex set
        vertices <- rbind(vertices, points)
      }
      # add small world edges
      for (i in 1:nrow(vertices)) {
        for (j in (i+1):nrow(vertices)) {
          if (runif(1) < beta) {
            edges <- rbind(edges, data.frame(from = i, to = j))
          }
        }
      }
      # add community structure
      community_size <- floor(num_vertices / num_communities)
      community_members <- rep(1:num_communities,x)
      
      
      
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Spatial Expander graph
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      generate_graph <- function(n, k, p, q, num_communities) {
        # Initialize adjacency matrix
        A <- matrix(0, nrow = n, ncol = n)
        
        # Generate initial graph using a 2-D grid
        side <- round(sqrt(n))
        for (i in 1:side) {
          for (j in 1:side) {
            v1 <- (i - 1) * side + j
            if (i > 1) {
              v2 <- v1 - side
              A[v1, v2] <- A[v2, v1] <- 1
            }
            if (j > 1) {
              v2 <- v1 - 1
              A[v1, v2] <- A[v2, v1] <- 1
            }
          }
        }
        
        # Compute degree of each vertex
        deg <- rowSums(A)
        
        # Grow graph using a 2-D Poisson process and preferential attachment
        for (i in (side^2 + 1):n) {
          # Generate spatial distribution of points on torus
          angle <- runif(k, 0, 2 * pi)
          radius <- rexp(k, p)
          x <- round(cos(angle) * radius, 2)
          y <- round(sin(angle) * radius, 2)
          coords <- cbind(x, y)
          
          # Compute distance from each existing vertex
          dists <- matrix(0, nrow = i-1, ncol = k)
          for (j in 1:(i-1)) {
            dx <- coords[,1] - x[j]
            dy <- coords[,2] - y[j]
            dx[abs(dx) > 0.5] <- dx[abs(dx) > 0.5] - sign(dx[abs(dx) > 0.5])
            dy[abs(dy) > 0.5] <- dy[abs(dy) > 0.5] - sign(dy[abs(dy) > 0.5])
            dists[j,] <- sqrt(dx^2 + dy^2)
          }
          
          # Compute probabilities for preferential attachment
          if (sum(deg) == 0) {
            probs <- rep(1/num_communities, num_communities)
          } else {
            comm_deg <- tapply(deg, rep(1:num_communities, each = n/num_communities), sum)
            comm_dists <- tapply(dists, rep(1:num_communities, each = n/num_communities), mean)
            comm_probs <- comm_deg^q / sum(comm_deg^q)
            comm_probs <- comm_probs * (1 - num_communities * 0.01) + 0.01
            comm_probs <- comm_probs * exp(-comm_dists/p)
            probs <- comm_probs / sum(comm_probs)
          }
          
          # Attach vertex to existing vertices
          for (j in 1:k) {
            attach_prob <- sample(probs, 1)
            attach_comm <- which(probs == attach_prob)
            attach_vertex <- sample(which(rep(1:num_communities, each = n/num_communities) == attach_comm))
            A[i, attach_vertex] <- A[attach_vertex, i] <- 1
            deg[c(i, attach_vertex)] <- deg[c(i, attach_vertex)] + 1
          }
        }
        
        # Convert adjacency matrix to igraph
        
      }  
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Spatial Expander graph
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      generate_graph <- function(num_vertices, lambda, beta, mu) {
        # Set up empty adjacency matrix
        adjacency_matrix <- matrix(0, ncol=num_vertices, nrow=num_vertices)
        
        # Set up empty list to hold the communities
        communities <- list()
        
        # Set up initial vertex randomly in 2-D space
        vertex_locations <- matrix(runif(2), ncol=2)
        
        # Add vertices one-by-one and connect them to existing vertices
        for (i in 2:num_vertices) {
          # Sample location for new vertex from 2-D Poisson process
          new_vertex_location <- matrix(rpois(2, lambda), ncol=2)
          new_vertex_location <- new_vertex_location + vertex_locations[i-1,]
          new_vertex_location <- new_vertex_location %% 1 # Wrap around to torus
          
          # Compute distances to existing vertices
          distances <- apply(vertex_locations, 1, function(v) sqrt(sum((v - new_vertex_location)^2)))
          
          # Compute probabilities of connecting to existing vertices based on distances
          probabilities <- beta * exp(-mu * distances)
          probabilities[i-1] <- 0 # Don't connect to self
          
          # Sample which vertices to connect to based on probabilities
          connection_indices <- sample(1:(i-1), size=rpois(1, lambda), replace=TRUE, prob=probabilities)
          
          # Add connections to adjacency matrix
          adjacency_matrix[i-1, connection_indices] <- 1
          adjacency_matrix[connection_indices, i-1] <- 1
          
          # Add vertex to list of communities
          community_index <- sample(length(communities), size=1, prob=sapply(communities, function(c) length(c)/(num_vertices-1)))
          communities[[community_index]] <- c(communities[[community_index]], i-1)
          
          # Update vertex locations
          vertex_locations <- rbind(vertex_locations, new_vertex_location)
        }
        
        # Compute degree distribution and small world clustering coefficient
        degree_distribution <- colSums(adjacency_matrix)
        clustering_coefficient <- mean(sapply(1:num_vertices, function(i) {
          neighbors <- which(adjacency_matrix[i,])
          num_neighbors <- length(neighbors)
          if (num_neighbors < 2) {
            return(0)
          }
          possible_edges <- combn(neighbors, 2)
          num_triangles <- sum(apply(possible_edges, 2, function(e) adjacency_matrix[e[1], e[2]]))
          num_triangles / choose(num_neighbors, 2)
        }))
        
        # Return results
        list(adjacency_matrix=adjacency_matrix, communities=communities, degree_distribution=degree_distribution, clustering_coefficient=clustering_coefficient)
      }
      
      
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Spatial Expander graph
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      # # Function to create a generative model
      # generative_model <- function(n, m, k, p, q, L) {
      #   # n: number of nodes
      #   # m: number of initial nodes to start the graph
      #   # k: average degree of each node
      #   # p: probability of adding a new edge between two existing nodes
      #   # q: probability of rewiring an existing edge
      #   # L: side length of the torus
      #   
      #   # Initialize the graph
      #   adj_mat <- matrix(0, n, n)
      #   pos <- matrix(runif(n*2, 0, L), n, 2)
      #   
      #   # Add the initial nodes
      #   for (i in 1:m) {
      #     for (j in 1:m) {
      #       if (i != j) {
      #         adj_mat[i, j] <- 1
      #       }
      #     }
      #   }
      #   
      #   # Grow the graph
      #   for (i in (m+1):n) {
      #     # Compute the distances to existing nodes
      #     dist <- sqrt((pos[1:(i-1), 1] - pos[i, 1])^2 + (pos[1:(i-1), 2] - pos[i, 2])^2)
      #     
      #     # Compute the degree distribution of existing nodes
      #     deg_dist <- rowSums(adj_mat[1:(i-1), 1:(i-1)])
      #     
      #     # Add k edges to the new node
      #     for (j in 1:k) {
      #       # Choose an existing node with probability proportional to its degree and distance
      #       prob <- deg_dist * exp(-dist^2/(2*k^2))
      #       prob[i] <- 0
      #       prob <- prob/sum(prob)
      #       j_new <- sample((1:(i-1)), size=1, prob=prob)
      #       
      #       # Add an edge between the new node and the chosen node with probability p
      #       if (runif(1) < p) {
      #         adj_mat[i, j_new] <- 1
      #         adj_mat[j_new, i] <- 1
      #       } else {
      #         # Rewire an existing edge with probability q
      #         adj_new <- adj_mat[i, ]
      #         adj_new[j_new] <- 0
      #         adj_new[i] <- 0
      #         prob_new <- deg_dist * (1-exp(-dist^2/(2*k^2)))
      #         prob_new[i] <- 0
      #         prob_new[j_new] <- 0
      #         prob_new <- prob_new/sum(prob_new)
      #         j_new <- sample((1:(i-1)), size=1, prob=prob_new)
      #         adj_new[j_new] <- 1
      #         adj_new[i] <- 1
      #         adj_mat[i, ] <- adj_new
      #         adj_mat[, i] <- adj_new
      #       }
      #     }
      #   }
      #   
      #   # Compute the clustering coefficient of each node
      #   clustering_coef <- apply(adj_mat, 1, function(x) {
      #     k <- sum(x)
      #     if (k <= 1) {
      #       return(0)
      #     } else {
      #       tri <- sum(x %*% t(x) * adj_mat)
      #       return(2*tri/(k*(k-1)))
      #     }
      #   })
      #   
      #   # Partition the nodes into communities
      #   modularity <- function(g, communities) {
      #     # Compute the modularity of a graph given a partition of its nodes into communities
      #     m <- sum(g)/2
      #     e <- sum(g %*% as.matrix(communities))
      #     a <- rowSums
      #   }
      # 
      # }
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Spatial Expander graph
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      # # Function to generate 2-D Poisson process on torus
      # poisson_torus <- function(n, lambda, torus_width, torus_height) {
      #   # Determine expected number of points
      #   expected_points <- round(lambda * torus_width * torus_height)
      #   
      #   # Generate Poisson process in rectangle
      #   x <- runif(expected_points, min = 0, max = torus_width)
      #   y <- runif(expected_points, min = 0, max = torus_height)
      #   
      #   # Apply periodic boundary conditions to wrap around torus
      #   x[x < n] <- x[x < n] + torus_width
      #   x[x > torus_width - n] <- x[x > torus_width - n] - torus_width
      #   y[y < n] <- y[y < n] + torus_height
      #   y[y > torus_height - n] <- y[y > torus_height - n] - torus_height
      #   
      #   # Return coordinates as matrix
      #   return(matrix(c(x, y), ncol = 2))
      # }
      # 
      # 
      # # Function to generate Barab?si-Albert scale-free graph
      # barabasi_albert <- function(n, m) {
      #   # Generate initial fully connected graph with m nodes
      #   edges <- matrix(1:m, ncol = 2, byrow = TRUE)
      #   for (i in m+1:n) {
      #     # Attach m edges to existing nodes using preferential attachment
      #     p <- colSums(table(as.vector(edges)))
      #     p <- p / sum(p)
      #     targets <- sample(1:(i-1), m, replace = TRUE, prob = p)
      #     edges <- rbind(edges, cbind(rep(i, m), targets))
      #   }
      #   # Convert edge list to adjacency matrix
      #   adj <- matrix(0, nrow = n, ncol = n)
      #   for (i in 1:nrow(edges)) {
      #     adj[edges[i,1], edges[i,2]] <- 1
      #     adj[edges[i,2], edges[i,1]] <- 1
      #   }
      #   # Return adjacency matrix
      #   return(adj)
      # }
      # 
      # 
      # # Function to generate Watts-Strogatz small-world graph
      # watts_strogatz <- function(n, k, p) {
      #   # Generate initial ring lattice
      #   adj <- matrix(0, nrow = n, ncol = n)
      #   for (i in 1:n) {
      #     for (j in 1:k) {
      #       j <- (j + i - 1) %% n + 1
      #       adj[i,j] <- 1
      #       
      
      
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Spatial Expander graph
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      
      
      # torus_graph <- function(n, k, p, q) {
      #   # n: number of nodes
      #   # k: degree of nodes
      #   # p: probability of rewiring edges
      #   # q: number of communities
      #   
      #   # Create an initial regular lattice with k neighbors
      #   coordinates <- matrix(0, nrow=n, ncol=2)
      #   for (i in 1:n) {
      #     coordinates[i,] <- c(cos(2*pi*(i-1)/n), sin(2*pi*(i-1)/n))
      #   }
      #   adj_matrix <- matrix(0, nrow=n, ncol=n)
      #   for (i in 1:n) {
      #     for (j in 1:k/2) {
      #       adj_matrix[i, ((i+j-1) %% n) + 1] <- 1
      #       adj_matrix[((i+j-1) %% n) + 1, i] <- 1
      #       adj_matrix[i, ((i-j-1+n) %% n) + 1] <- 1
      #       adj_matrix[((i-j-1+n) %% n) + 1, i] <- 1
      #     }
      #   }
      #   
      #   # Rewire edges with probability p
      #   for (i in 1:n) {
      #     for (j in 1:n) {
      #       if (i < j && adj_matrix[i,j] == 1 && runif(1) < p) {
      #         # Find a new random neighbor for i
      #         new_neighbor <- sample(1:n, 1)
      #         while (new_neighbor == i || adj_matrix[i, new_neighbor] == 1) {
      #           new_neighbor <- sample(1:n, 1)
      #         }
      #         # Rewire the edge
      #         adj_matrix[i,j] <- 0
      #         adj_matrix[j,i] <- 0
      #         adj_matrix[i,new_neighbor] <- 1
      #         adj_matrix[new_neighbor,i] <- 1
      #       }
      #     }
      #   }
      #   
      #   # Assign nodes to communities
      #   community_assignment <- rep(1:q, each=floor(n/q))
      #   community_assignment <- c(community_assignment, sample(1:q, n %% q, replace=TRUE))
      #   
      #   # Create a spatial distribution of points for each community
      #   community_centers <- matrix(0, nrow=q, ncol=2)
      #   for (i in 1:q) {
      #     community_centers[i,] <- c(runif(1), runif(1))
      #   }
      #   community_points <- list()
      #   for (i in 1:q) {
      #     community_points[[i]] <- matrix(rpois(n=k, lambda=1)*0.01, ncol=2) + community_centers[i,]
      #   }
      #   
      #   # Calculate distances between nodes and assign edges based on proximity and community
      #   edge_list <- list()
      #   for (i in 1:n) {
      #     distances <- matrix(0, nrow=n, ncol=2)
      #     for (j in 1:n) {
      #       distances[j,] <- coordinates[j,] - coordinates[i,]
      #       distances[j,] <- ifelse(distances[j,] > 0.5, distances[j,] - 1, distances[j,])
      #       distances[j,] <- ifelse(distances[j,] < -0.5, distances[j,] + 1, distances[j,])
      #     }
      #     distances <- apply(distances^2, 1, sum)
      #     for (j in 1:n) {
      #       if (i < j
      #           
      
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Spatial Expander graph
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      # # Set parameters
      # n <- 1000 # number of nodes
      # r <- 0.05 # radius of point process
      # alpha <- 2.5 # exponent for degree distribution
      # k <- 10 # target average degree
      # p <- 0.1 # rewiring probability
      # c <- 0.2 # community size proportion
      # 
      # # Generate spatial points
      # theta <- runif(n, 0, 2*pi)
      # rho <- sqrt(runif(n, 0, r^2))
      # x <- rho * cos(theta)
      # y <- rho * sin(theta)
      # 
      # # Calculate distance matrix on torus
      # dist_mat <- matrix(0, nrow=n, ncol=n)
      # for (i in 1:n) {
      #   for (j in 1:n) {
      #     dx <- abs(x[i] - x[j])
      #     dy <- abs(y[i] - y[j])
      #     dx <- min(dx, 1-dx)
      #     dy <- min(dy, 1-dy)
      #     dist_mat[i,j] <- sqrt(dx^2 + dy^2)
      #   }
      # }
      # 
      # # Generate edges
      # adj_mat <- matrix(0, nrow=n, ncol=n)
      # deg_seq <- rpois(n, r*k) # degree sequence
      # for (i in 1:n) {
      #   # Select neighbors within distance r
      #   neighbor_idxs <- which(dist_mat[i,] < r & dist_mat[i,] > 0)
      #   # Sample from potential neighbors based on distance
      #   prob <- exp(-dist_mat[i, neighbor_idxs])
      #   prob <- prob/sum(prob)
      #   new_neighbors <- sample(neighbor_idxs, deg_seq[i], replace=TRUE, prob=prob)
      #   # Connect to new neighbors
      #   adj_mat[i, new_neighbors] <- 1
      #   adj_mat[new_neighbors, i] <- 1
      # }
      # 
      # # Rewire edges
      # for (i in 1:n) {
      #   for (j in (i+1):n) {
      #     if (adj_mat[i,j] == 1 & runif(1) < p) {
      #       # Rewire edge
      #       candidate_idxs <- which(adj_mat[i,] == 0 & dist_mat[i,] < 1/2)
      #       prob <- dpois(sum(adj_mat[i,]), candidate_idxs, log=TRUE)
      #       prob[candidate_idxs==j] <- 0
      #       prob <- prob - max(prob)
      #       prob <- exp(prob)
      #       prob[candidate_idxs==j] <- 0
      #       prob <- prob/sum(prob)
      #       new_j <- sample(candidate_idxs, 1, prob=prob)
      #       adj_mat[i,j] <- 0
      #       adj_mat[j,i] <- 0
      #       adj_mat[i,new_j] <- 1
      #       adj_mat[new_j,i] <- 1
      #     }
      #   }
      # }
      # 
      # # Create communities
      # num_communities <- floor(n*c)
      # community_sizes <- rep(ceiling(n/num_communities), num_communities)
      # community_sizes[num_communities] <- n - sum(community_sizes[1:(num_communities-1)])
      # community_idxs <- rep(1:num_communities, community_sizes)
      # community_adj_mat <- matrix(0, nrow=n, ncol=n)
      # for (i in 1:n) {
      #   for (j in (i+1):n) {
      #     if (community_idxs[i] == community_idxs[j] & adj_mat[i,j] == 1) {
      #       community_adj_mat[i,j] <- 1
      #       community_adj_mat[j,i] <- 1
      #     }
      #   }
      # }
      # 
      # # Create degree sequence
      # degree_seq <-
      #   
      
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Spatial Expander graph
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # 
      # generate_torus_spatial_expander <- function(n, k, p, r, q) {
      #   # n: number of nodes
      #   # k: number of initial connections for each new node
      #   # p: probability of connecting to a spatially nearby node
      #   # r: radius of the torus
      #   # q: number of communities
      #   
      #   # create the initial network with k nodes
      #   nodes <- 1:k
      #   edges <- matrix(c(rep(nodes, each=k-1), 
      #                     unlist(sapply(nodes, function(x) setdiff(nodes, x)))), 
      #                   ncol=2)
      #   g <- list(nodes=nodes, edges=edges)
      #   
      #   # add nodes to the network
      #   for (i in (k+1):n) {
      #     # generate a 2D Poisson point process with intensity parameter proportional to the degree of the nodes
      #     intensities <- sapply(g$nodes, function(x) sum(g$edges %in% x))
      #     xy <- cbind(runif(i-1, 0, r), runif(i-1, 0, r))
      #     poi <- rpois(i-1, intensities*mean(p))
      #     poi <- split(xy, rep(1:(i-1), poi))
      #     
      #     # connect the new node to k spatially nearby nodes with probability p
      #     distances <- as.matrix(dist(rbind(xy, c(0,0)), method="euclidean")[1:(i-1),])
      #     for (j in 1:k) {
      #       neighbors <- which(distances[j,] < r)
      #       neighbors <- neighbors[neighbors != i]
      #       neighbors <- neighbors[sample(length(neighbors), 1)]
      #       if (runif(1) < p) {
      #         g$edges <- rbind(g$edges, c(i, neighbors))
      #       }
      #     }
      #     
      #     # connect the new node to m random nodes
      #     m <- max(1, round(rpois(1, k-1)))
      #     if (m > 0) {
      #       neighbors <- sample(setdiff(g$nodes, i), m)
      #       g$edges <- rbind(g$edges, c(rep(i, m), neighbors))
      #     }
      #     
      #     # assign the new node to a community
      #     if (q > 1) {
      #       g$nodes <- c(g$nodes, i)
      #       communities <- sample(1:q, i, replace=TRUE)
      #       g$communities <- communities
      #     }
      #   }
      #   
      #   # reindex the nodes
      #   if (q == 1) {
      #     g$nodes <- 1:n
      #   }
      #   
      #   # calculate the clustering coefficient for each node
      #   g$cc <- sapply(g$nodes, function(x) {
      #     neighbors <- g$edges[g$edges[,1] == x,2]
      #     if (length(neighbors) < 2) {
      #       return(0)
      #     }
      #     else {
      #       subgraph <- g$edges[g$edges[,1] %in% neighbors,]
      #       subgraph <- rbind(subgraph, subgraph[,2:1])
      #       subgraph <- unique(subgraph)
      #       n_triangles <- sum(apply(subgraph, 1, function(x) length(intersect(neighbors, g$edges[g$edges[,1] == x[2],2]))))
      #       return(n_triangles / (length(neighbors) * (length(neighbors)-1) / 2))
      #     }
      #   })
      #   
      #   # calculate the degree distribution
      #   g$degree <- sapply(g$nodes, function(x)
      #     
      # 





#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sample code 1
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # Load required packages
# library(igraph)
# library(caret)
# 
# # Load data
# data(karate)
# 
# # Generate synthetic graphs using Barabasi-Albert model and Erdos-Renyi model
# ba_graph <- sample_pa(30, m=3, directed=FALSE)
# er_graph <- erdos.renyi.game(30, p.or.m=0.5, type="gnp", directed=FALSE)
# 
# # Extract graph features
# ba_features <- as.data.frame(t(sapply(degree(ba_graph), function(x) hist(x, breaks=0:30, plot=FALSE)$counts)))
# er_features <- as.data.frame(t(sapply(degree(er_graph), function(x) hist(x, breaks=0:30, plot=FALSE)$counts)))
# 
# # Combine the features of synthetic graphs into a single data frame
# syn_features <- rbind(ba_features, er_features)
# syn_labels <- c(rep("Barabasi-Albert", nrow(ba_features)), rep("Erdos-Renyi", nrow(er_features)))
# 
# # Extract features from the empirical graph
# emp_features <- as.data.frame(t(sapply(degree(karate), function(x) hist(x, breaks=0:30, plot=FALSE)$counts)))
# emp_labels <- rep("Empirical", nrow(emp_features))
# 
# # Combine the features of synthetic and empirical graphs into a single data frame
# all_features <- rbind(syn_features, emp_features)
# all_labels <- c(syn_labels, emp_labels)
# 
# # Split the data into training and testing sets
# set.seed(123)
# train_index <- createDataPartition(all_labels, p=0.7, list=FALSE)
# train_features <- all_features[train_index,]
# train_labels <- all_labels[train_index]
# test_features <- all_features[-train_index,]
# test_labels <- all_labels[-train_index]
# 
# # Train a Random Forest classifier
# model <- train(x=train_features, y=train_labels, method="rf")
# 
# # Predict on the testing set
# predictions <- predict(model, newdata=test_features)
# 
# # Evaluate the model accuracy
# confusionMatrix(predictions, test_labels)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sample code 2
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load required libraries
# library(igraph)
# library(caret)
# 
# # Load the data
# data(polbooks)
# 
# # Create a matrix of graph properties for each graph
# graph_props <- sapply(polbooks, function(x) {
#   g <- graph.adjacency(as.matrix(x), weighted=TRUE)
#   return(c(degree(g),betweenness(g),closeness(g)))
# })
# 
# # Create the training data
# train_data <- graph_props[,1:10]
# train_labels <- c("Model A", "Model A", "Model A", "Model A", "Model A", "Model B", "Model B", "Model B", "Model B", "Model B")
# 
# # Train the classification model using a Random Forest
# model <- train(x = train_data, y = train_labels, method = "rf", trControl = trainControl(method = "cv", number = 10))
# 
# # Predict on new graph properties
# new_graph_props <- graph_props[,11:15]
# predictions <- predict(model, new_graph_props)
# 
# # Evaluate the model performance
# confusionMatrix(predictions, c("Model A", "Model A", "Model B", "Model B", "Model B"))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sample code 3
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# library(igraph)
# library(randomForest)
# 
# # Generate Synthetic Graphs
# set.seed(1234)
# g1 <- barabasi.game(100, m=2)
# g2 <- erdos.renyi.game(100, p.or.m=0.1)
# g3 <- watts.strogatz.game(100, dim=1, size=3, p=0.1)
# 
# # Add labels to the graphs
# V(g1)$label <- "Barabasi"
# V(g2)$label <- "Erdos Renyi"
# V(g3)$label <- "Watts Strogatz"
# 
# # Combine the graphs into a list
# graphs <- list(g1, g2, g3)
# 
# # Calculate graph features for each graph
# features <- data.frame()
# for (i in 1:length(graphs)) {
#   g <- graphs[[i]]
#   temp_features <- c(degree(g), clustering(g), transitivity(g))
#   temp_df <- data.frame(temp_features)
#   temp_df$label <- V(g)$label[1]
#   features <- rbind(features, temp_df)
# }
# 
# # Split data into training and test sets
# trainIndex <- sample(1:nrow(features), size = 0.7*nrow(features))
# train_data <- features[trainIndex,]
# test_data <- features[-trainIndex,]
# 
# # Train the classification model
# model <- randomForest(label ~ degree + clustering + transitivity, data = train_data, ntree=1000)
# 
# # Predict the class labels on test data
# predictions <- predict(model, newdata = test_data)
# 
# # Evaluate the model accuracy
# confusionMatrix(predictions, test_data$label)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sample code 4
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # Load required libraries
# library(igraph)
# library(caret)
# 
# # Generate example graphs using different generative models
# set.seed(123)
# g1 <- sample_pa(1000, m=5)
# g2 <- sample_ws(1000, m=5, k=5)
# g3 <- sample_gnp(1000, m=5, p=0.5)
# 
# # Extract graph features for each graph
# g1_features <- data.frame(degree_centrality(g1), clustering_coef(g1), t(apply(get.adjacency(g1), 1, sort)))
# g2_features <- data.frame(degree_centrality(g2), clustering_coef(g2), t(apply(get.adjacency(g2), 1, sort)))
# g3_features <- data.frame(degree_centrality(g3), clustering_coef(g3), t(apply(get.adjacency(g3), 1, sort)))
# 
# # Combine graph features into one data frame
# df <- rbind(cbind(g1_features, model = "PA"), cbind(g2_features, model = "WS"), cbind(g3_features, model = "GNP"))
# 
# # Split data into training and testing sets
# set.seed(456)
# split_index <- createDataPartition(df$model, p = 0.7, list = FALSE)
# training_data <- df[split_index,]
# testing_data <- df[-split_index,]
# 
# # Train classification model
# model <- train(model ~ ., data = training_data, method = "rf")
# 
# # Predict model on testing data
# predictions <- predict(model, testing_data)
# 
# # Evaluate accuracy of model
# confusionMatrix(predictions, testing_data$model)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sample code 5
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# library(igraph)
# library(caret)
# 
# # Generate data from three different generative models: Barabasi-Albert, Watts-Strogatz and Random graph
# g1 <- barabasi.game(1000, m = 2, directed = FALSE)
# g2 <- watts.strogatz.game(1000, dim = 2, size = 6, rewire = 0.05, circular = FALSE)
# g3 <- sample_gnp(1000, 0.1, directed = FALSE, loops = FALSE)
# 
# # Extract graph features to use as predictors
# features1 <- c(degree.distribution(g1)$x, g1$transitivity)
# features2 <- c(degree.distribution(g2)$x, g2$transitivity)
# features3 <- c(degree.distribution(g3)$x, g3$transitivity)
# 
# # Combine features into a data frame and add target variable
# features <- rbind(features1, features2, features3)
# labels <- c(rep("Barabasi-Albert", 1000), rep("Watts-Strogatz", 1000), rep("Random Graph", 1000))
# data <- data.frame(features, labels)
# 
# # Split data into training and testing sets
# set.seed(123)
# index <- createDataPartition(data$labels, p = 0.8, list = FALSE)
# train_data <- data[index, ]
# test_data <- data[-index, ]
# 
# # Train the model using Random Forest
# model <- train(labels ~ ., data = train_data, method = "rf")
# 
# # Test the model on the testing data
# predictions <- predict(model, newdata = test_data[, -ncol(test_data)])
# 
# # Evaluate model accuracy
# confusionMatrix(predictions, test_data$labels)
# # #


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sample code 6
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# library(igraph)
# library(randomForest)
# 
# # Define the generative models
# 
# # Model 1: Barabasi-Albert Model
# ba_model <- function(n, m){
#   return(barabasi.game(n, m, directed = FALSE, power = 1))
# }
# 
# # Model 2: Small World Model
# sw_model <- function(n, k, p){
#   return(watts.strogatz.game(n, k, p, loops = FALSE))
# }
# 
# # Model 3: Regular Graph Model
# rg_model <- function(n, k){
#   return(make_ring(n, directed = FALSE, mutual = FALSE, k = k))
# }
# 
# # Generate 100 graphs for each model
# ba_graphs <- lapply(1:100, function(x){ba_model(50, 2)})
# sw_graphs <- lapply(1:100, function(x){sw_model(50, 10, 0.05)})
# rg_graphs <- lapply(1:100, function(x){rg_model(50, 4)})
# 
# # Extract graph features
# ba_features <- lapply(ba_graphs, function(x){
#   c(degree.distribution(x), transitivity(x), average.path.length(x))
# })
# sw_features <- lapply(sw_graphs, function(x){
#   c(degree.distribution(x), transitivity(x), average.path.length(x))
# })
# rg_features <- lapply(rg_graphs, function(x){
#   c(degree.distribution(x), transitivity(x), average.path.length(x))
# })
# 
# # Combine the features into a single data frame
# ba_df <- data.frame(model = rep("ba", 100), do.call(rbind, ba_features))
# sw_df <- data.frame(model = rep("sw", 100), do.call(rbind, sw_features))
# rg_df <- data.frame(model = rep("rg", 100), do.call(rbind, rg_features))
# df <- rbind(ba_df, sw_df, rg_df)
# 
# # Train the random forest model
# model <- randomForest(model ~ ., data = df, ntree = 500, importance = TRUE)
# 
# # Predict the model for the empirical graph
# empirical_graph <- read.graph("empirical_graph.gml", format = "gml")
# empirical_features <- c(degree.distribution(empirical_graph), 
#                         transitivity(empirical_graph), 
#                         average.path.length(empirical_graph))
# prediction <- predict(model, newdata = data.frame(empirical_features))
# 
# # Print the prediction
# print(prediction)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sample code 7
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # Load the igraph package
# library(igraph)
# 
# # Generate a random normal distribution of points
# points <- rnorm(100, mean = 0, sd = 1)
# 
# # Generate a scale-free degree distribution
# degrees <- degree.sequence(100, power.law.fit(100))
# graph <- degree.sequence.game(degrees)
# 
# # Implement the small world effect
# graph <- rewire(graph, n = 10)
# 
# # Plot the graph
# plot(graph, vertex.size = 5, vertex.color = "blue")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sample code 8
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# library(igraph)
# 
# # Function to generate a spatial distribution of points
# generate_points <- function(n, xlim, ylim){
#   points <- matrix(runif(2 * n, min = 0, max = 1), ncol = 2)
#   points[, 1] <- points[, 1] * xlim
#   points[, 2] <- points[, 2] * ylim
#   return(points)
# }
# 
# # Set parameters
# n = 1000 # number of nodes
# m = 5 # number of neighbors for each node
# p = 0.01 # probability of rewiring
# 
# # Generate a random spatial distribution of points
# points = matrix(runif(2*n), n, 2)
# 
# # Generate a small-world network using the Watts-Strogatz model
# g = watts.strogatz.game(1, n, m, p)
# 
# # Plot the network with points
# plot(g, layout=points, vertex.size=3, vertex.label=NA)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sample code 9
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Generative Model on a Torus
# 
# #Set the number of nodes in the network
# n <- 1000
# 
# #Generate a spatial distribution of points
# x <- runif(n, 0, 1) #x-coordinate
# y <- runif(n, 0, 1) #y-coordinate
# 
# #Create a 2D Poisson process for the nodes
# lambda <- 0.1 #intensity parameter
# 
# #Calculate the expected number of points in each unit square
# k <- lambda * (1-0)^2 #area of unit square
# 
# #Generate the number of points in each unit square using a Poisson distribution
# n_points <- rpois(k, lambda)
# 
# #Set the scale-free degree distribution of the nodes
# d <- rgamma(n, shape = 3, rate = 1/2)
# 
# #Calculate the probability of a node being connected to another node
# p <- d/sum(d)
# 
# #Randomly connect nodes to create a small world effect
# for (i in 1:n) {
#   for (j in (i+1):n) {
#     if (runif(1) < p[i] * p[j]) {
#       #create an edge between nodes i and j
#     }
#   }
# }
# 
# #Plot the network on a torus
# library(ggplot2)
# 
# #Create a data frame with the x and y coordinates of the nodes
# df <- data.frame(x = x, y = y)
# 
# #Create a torus shaped plot
# ggplot(df) + 
#   geom_point(aes(x, y)) + 
#   coord_fixed(ratio = 1) + 
#   xlim(0, 1) + 
#   ylim(0, 1) + 
#   theme_void()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sample code 10
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Required Libraries
# library(spatstat)
# library(pracma)
# 
# # Parameter Settings
# n_nodes <- 500 # number of nodes
# mean_points <- 10 # mean number of points per node
# radius <- 1 # radius of the torus
# 
# # Generate a Poisson Process on a Torus
# poisson_process <- rpoispp(n_nodes, mean_points, torus(radius))
# 
# # Generate a 2D Graph
# adj_matrix <- matrix(0, nrow = n_nodes, ncol = n_nodes)
# for (i in 1:n_nodes) {
#   for (j in (i + 1):n_nodes) {
#     dist_ij <- min_dist(poisson_process$x[i], poisson_process$y[i], poisson_process$x[j], poisson_process$y[j])
#     if (dist_ij <= radius) {
#       adj_matrix[i, j] <- 1
#       adj_matrix[j, i] <- 1
#     }
#   }
# }
# 
# # Random Connection of Nodes
# random_adj_matrix <- sample_random_graph(adj_matrix, method = "ER", p = 0.1)
# 
# # Small World Effect
# small_world_adj_matrix <- rewire(random_adj_matrix, method = "Watts-Strogatz", p = 0.1)
# 
# # Scale-Free Degree Distribution
# degree_sequence <- degree_sequence(small_world_adj_matrix, cumulative = FALSE)
# scale_free_adj_matrix <- configure_model(n = n_nodes, exponent = 2, method = "BA", m = 3,
#                                          cumdegree = degree_sequence)
# 
# ################%%%%%%%##########
# # Required Libraries
# library(ggplot2)
# 
# # Function to create a torus
# create_torus <- function(radius_a, radius_b) {
#   torus <- function(theta, phi) {
#     x <- (radius_a + radius_b * cos(theta)) * cos(phi)
#     y <- (radius_a + radius_b * cos(theta)) * sin(phi)
#     z <- radius_b * sin(theta)
#     return(c(x, y, z))
#   }
#   return(torus)
# }
# 
# # Function to generate points on a torus
# generate_points_on_torus <- function(n, radius_a, radius_b) {
#   torus <- create_torus(radius_a, radius_b)
#   theta <- runif(n, 0, 2 * pi)
#   phi <- runif(n, 0, 2 * pi)
#   points <- cbind(theta, phi)
#   points <- t(apply(points, 1, torus))
#   return(points)
# }
# 
# # Function to generate connections between nodes
# generate_connections <- function(n, m, points) {
#   connections <- matrix(0, nrow = n, ncol = n)
#   for (i in 1:n) {
#     distances <- apply(points, 1, function(x) sqrt(sum((x - points[i, ])^2)))
#     closest_indices <- order(distances, index.return = TRUE)$ix[1:(m + 1)]
#     for (j in closest_indices[2:length(closest_indices)]) {
#       connections[i, j] <- 1
#       connections[j, i] <- 1
#     }
#   }
#   return(connections)
# }
# 
# # Function to generate a scale-free degree distribution
# generate_scale_free_degree_distribution <- function(n, m) {
#   degrees <- numeric(n)
#   for (i in 1:n) {
#     degrees[i] <- rpois(1, m)
#   }
#   return(degrees)
# }
# 
# # Function to generate small world effects
# generate_small_world_effects <- function(n, m, connections, p) {
#   for (i in 1:n) {
#     for (j in (i + 1):n) {
#       if (connections[i, j] == 1) {
#         if (runif(1) < p) {
#           new_connection <- sample(setdiff(1:n, c(i, j)), 1)
#           connections[i, new_connection] <- 1
#           connections[new_connection, i] <- 1
#           connections[j, new_connection] <- 0
#           connections[new_connection, j] <- 0
#         }
#       }
#     }
#   }
#   return(connections)
# }
# 
# # Generate the model
# n <- 1000 # number of nodes
# m <- 10 # number of nearest neighbors
# radius_a <- 30 # radius of the torus
# radius_b <- 10 # radius
# 
# #########%%%%####################
# library(spatstat)
# 
# # Create a torus of a specified size
# torus_size <- 50
# torus <- torus(torus_size)
# 
# # Generate the initial points on the torus using a 2D Poisson process
# lambda <- 0.1
# initial_points <- rpoispp(lambda, torus)
# 
# # Create a function for adding new nodes to the network
# add_node <- function(points, torus) {
#   # Generate new points on the torus following a 2D Poisson process
#   new_points <- rpoispp(lambda, torus)
#   # Combine the existing and new points
#   combined_points <- c(points, new_points)
#   # Randomly connect the nodes
#   connections <- sample(combined_points, size = 2, replace = TRUE)
#   return(list(points = combined_points, connections = connections))
# }
# 
# # Define the number of nodes to add
# num_nodes <- 100
# 
# # Create an empty list to store the nodes and connections
# network <- list()
# 
# # Add nodes to the network
# for (i in 1:num_nodes) {
#   result <- add_node(network$points, torus)
#   network$points <- result$points
#   network$connections <- c(network$connections, result$connections)
# }
# 
# # Plot the network
# plot(network$points, main = "Generative model on a torus")
# 
# ################%%%%%%#############
# library(spatstat)
# 
# # Define the torus
# torus <- owin(c(-pi, -pi), c(pi, pi), shape="square")
# 
# # Generate points on the torus
# points <- rpoispp(100, torus)
# 
# # Generate a random connection matrix based on a 2D Poisson process
# distance_matrix <- as.matrix(dist(points))
# connection_matrix <- as.matrix(exp(-distance_matrix))
# 
# # Generate scale-free degree distribution
# scale_free_matrix <- degree.sequence.game(connection_matrix, method="psumtree")
# 
# # Generate small-world effects
# small_world_matrix <- rewire(scale_free_matrix, method="sw")
# 
# # Plot the generated network
# plot(points, xlim=c(-pi, pi), ylim=c(-pi, pi), pch=20)
# lines(graph.adjacency(small_world_matrix, mode="directed", weighted=NULL, diag=FALSE), col="red")
# 
# ###########%%%%%%%%#################
# library(spatstat)
# 
# # Define the torus
# torus_fun <- function(x, y) {
#   return(c(cos(2*pi*x), sin(2*pi*y)))
# }
# 
# torus_data <- t(sapply(runif(100, 0, 1), torus_fun))
# 
# # Create a 2D Poisson process
# ppp <- ppp(torus_data, c(0,0), 1, torus_data)
# 
# # Define a function for creating connections between nodes
# create_connections <- function(ppp, prob=0.5) {
#   n <- npoints(ppp)
#   adjacency_matrix <- matrix(0, nrow=n, ncol=n)
#   for (i in 1:(n-1)) {
#     for (j in (i+1):n) {
#       if (runif(1) < prob) {
#         adjacency_matrix[i,j] <- 1
#         adjacency_matrix[j,i] <- 1
#       }
#     }
#   }
#   return(adjacency_matrix)
# }
# 
# # Generate the adjacency matrix
# adjacency_matrix <- create_connections(ppp, prob=0.3)
# 
# # Define a function for creating a scale-free network
# create_scale_free_network <- function(adjacency_matrix, n_steps=100) {
#   degree_sequence <- rowSums(adjacency_matrix)
#   while (n_steps > 0) {
#     node_index <- sample(1:length(degree_sequence), 1, prob=degree_sequence)
#     degree_sequence[node_index] <- degree_sequence[node_index] + 1
#     n_steps <- n_steps - 1
#   }
#   return(degree_sequence)
# }
# 
# # Generate the scale-free network
# degree_sequence <- create_scale_free_network(adjacency_matrix, n_steps=100)
# 
# # Define a function for creating small-world network
# create_small_world_network <- function(adjacency_matrix, rewiring_prob=0.1) {
#   n <- nrow(adjacency_matrix)
#   for (i in 1:(n-1)) {
#     for (j in (i+1):n) {
#       if (runif(1) < rewiring_prob) {
#         k <- sample(setdiff(1:n, c(i,j)), 1)
#         adjacency_matrix[i,j] <- 0
#         adjacency_matrix[j,i] <- 0
#         adjacency_matrix[i,k] <- 1
#         adjacency_matrix[k,i] <- 1
#       }
#     }
#   }
#   return(adjacency_matrix)
# }
# 
# # Generate the small-world network
# adjacency_matrix <- create_small_world_network(adjacency_matrix, rewiring_prob=0.1)
# 
# ########%%%%%%####################
# # Define the number of nodes in the network
# nodes <- 1000
# 
# # Define the torus size
# torus_size <- 100
# 
# # Generate the initial set of nodes on the torus using a 2D Poisson process
# nodes_2d <- rpois2d(nodes, lambda = 1, xlim = c(-torus_size, torus_size), ylim = c(-torus_size, torus_size))
# 
# # Create a matrix to store the connections between nodes
# connections <- matrix(0, nrow = nodes, ncol = nodes)
# 
# # Randomly connect nodes with a probability that can be adjusted
# p_connect <- 0.01
# for(i in 1:nodes) {
#   for(j in (i+1):nodes) {
#     if(runif(1) < p_connect) {
#       connections[i, j] <- 1
#       connections[j, i] <- 1
#     }
#   }
# }
# 
# # Define a function to calculate the degree of each node
# degree_of_node <- function(matrix) {
#   rowSums(matrix)
# }
# 
# # Calculate the degree of each node
# degrees <- degree_of_node(connections)
# 
# # Create a function to calculate the clustering coefficient of each node
# clustering_coeff <- function(matrix, node) {
#   neighbors <- which(matrix[node, ] == 1)
#   num_neighbors <- length(neighbors)
#   if(num_neighbors < 2) {
#     return(0)
#   }
#   connections_between_neighbors <- matrix[neighbors, neighbors]
#   triads <- sum(connections_between_neighbors) / 2
#   return(triads / (num_neighbors * (num_neighbors - 1) / 2))
# }
# 
# # Calculate the clustering coefficient for each node
# clustering_coefficients <- sapply(1:nodes, function(x) clustering_coeff(connections, x))
# 
# # Create a function to calculate the average shortest path length between nodes
# average_shortest_path_length <- function(matrix) {
#   distances <- graph.shortest.paths(graph.adjacency(matrix, mode = "undirected", weighted = NULL), mode = "all")
#   return(mean(distances, na.rm = TRUE))
# }
# 
# # Calculate the average shortest path length
# avg_shortest_path_length <- average_shortest_path_length(connections)
# 
# # Plot the distribution of node degrees to check for scale-free behavior
# hist(degrees, main = "Distribution of Node Degrees", xlab = "Degree", col = "blue")
# 
# # Plot the clustering coefficients to check for small world behavior
# plot(clustering_coefficients, main = "Clustering Coefficient of Nodes", xlab = "Node", ylab = "Clustering Coefficient", col = "red")
# 
# ##########%%%%%############%
# # Required Libraries
# library(spatstat)
# library(spatstat.utils)
# 
# # Set Parameters
# nodes <- 1000 # number of nodes
# d <- 20 # dimension of the torus
# radius <- 5 # radius of the connection function
# p_conn <- 0.1 # probability of connection between two nodes
# 
# # Create a Torus
# torus <- torus(d)
# 
# # Create Points on the Torus
# points <- rpoispp(nodes, torus)
# 
# # Compute Distances between Nodes
# distances <- spatstat.utils::rdist(points)
# 
# # Create Adjacency Matrix
# adj_matrix <- as.matrix(distances <= radius)
# 
# # Apply Random Connection
# random_connection <- runif(nodes^2)
# adj_matrix[random_connection > p_conn] <- 0
# 
# # Symmetrize Adjacency Matrix
# adj_matrix <- t(adj_matrix) + adj_matrix
# 
# # Create Network Object
# network <- as.network(adj_matrix)
# 
# # Plot the Network
# plot(network, vertex.cex=0.3)
# 
# ##############%%%%%#############
# # Generate a 2D Poisson Process on a Torus
# generate_poisson_process <- function(lambda, L) {
#   # lambda: intensity of the Poisson process
#   # L: size of the torus
#   n <- rpois(1, lambda * L^2) # number of points generated
#   x <- runif(n, 0, L) # x-coordinates of the points
#   y <- runif(n, 0, L) # y-coordinates of the points
#   # periodic boundary conditions on the torus
#   x <- x - floor(x / L) * L
#   y <- y - floor(y / L) * L
#   return(cbind(x, y))
# }
# 
# # Generate random connections between nodes
# generate_random_connections <- function(n, p) {
#   # n: number of nodes
#   # p: probability of connections
#   adjacency_matrix <- matrix(0, n, n)
#   for (i in 1:(n - 1)) {
#     for (j in (i + 1):n) {
#       if (runif(1) < p) {
#         adjacency_matrix[i, j] <- 1
#         adjacency_matrix[j, i] <- 1
#       }
#     }
#   }
#   return(adjacency_matrix)
# }
# 
# # Generate scale-free degree distribution
# generate_scale_free_degree_distribution <- function(n, m0, m) {
#   # n: number of nodes
#   # m0: initial number of nodes with degree 0
#   # m: number of edges added at each iteration
#   degree_sequence <- numeric(n)
#   degree_sequence[1:m0] <- 0
#   for (i in (m0 + 1):n) {
#     prob <- degree_sequence[1:(i - 1)]^(-1)
#     prob <- prob / sum(prob)
#     new_edges <- sample(1:(i - 1), m, replace = TRUE, prob = prob)
#     degree_sequence[new_edges] <- degree_sequence[new_edges] + 1
#     degree_sequence[i] <- m
#   }
#   adjacency_matrix <- matrix(0, n, n)
#   for (i in 1:n) {
#     if (degree_sequence[i] > 0) {
#       neighbors <- sample(1:n, degree_sequence[i], replace = TRUE)
#       adjacency_matrix[i, neighbors] <- 1
#     }
#   }
#   return(adjacency_matrix)
# }
# 
# # Implement small-world effects
# generate_small_world_network <- function(n, k, beta) {
#   # n: number of nodes
#   # k: mean degree
#   # beta: rewiring probability
#   adjacency_matrix <- matrix(0, n, n)
#   for (i in 1:(n - 1)) {
#     for (j in (i + 1):n) {
#       if (runif(1) < k / (n - 1)) {
#         adjacency_matrix[i, j] <- 1
#         adjacency


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sample code 11
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Generate a generative model on a torus
# 
# # Required Libraries
# library(deldir)
# 
# # Set the seed for reproducibility
# set.seed(123)
# 
# # Define the number of nodes in the network
# n <- 100
# 
# # Define the spatial distribution of points on a torus
# x <- (1 + cos(seq(0, 2 * pi, length.out = n))) * cos(seq(0, 2 * pi, length.out = n))
# y <- (1 + cos(seq(0, 2 * pi, length.out = n))) * sin(seq(0, 2 * pi, length.out = n))
# 
# # Generate a 2D Poisson process
# points <- cbind(x, y)
# 
# # Create a Delaunay triangulation
# triangulation <- deldir(points[,1], points[,2])
# 
# # Extract the edges of the triangulation
# edges <- triangulation$tile.list
# 
# # Define the average degree of the network
# avg_degree <- 8
# 
# # Define the probability of rewiring an edge
# rewiring_prob <- 0.1
# 
# # Generate a random network with small-world properties
# adj_matrix <- matrix(0, n, n)
# for (i in 1:n) {
#   for (j in i:n) {
#     if (i != j) {
#       if (length(grep(paste(c(i, j), collapse = "-"), edges)) > 0) {
#         adj_matrix[i, j] <- 1
#         adj_matrix[j, i] <- 1
#       } else if (runif(1) < rewiring_prob) {
#         adj_matrix[i, j] <- 1
#         adj_matrix[j, i] <- 1
#       }
#     }
#   }
# }
# 
# # Generate a scale-free network
# degree_sequence <- degree(adj_matrix)
# while (mean(degree_sequence) < avg_degree) {
#   new_node <- n + 1
#   n <- n + 1
#   prob <- degree_sequence / sum(degree_sequence)
#   target_node <- sample(1:n, 1, prob = prob)
#   adj_matrix[new_node, target_node] <- 1
#   adj_matrix[target_node, new_node] <- 1
#   degree_sequence <- degree(adj_matrix)
# }
# 
# # Plot the generated network on a torus
# plot(x, y, pch = 19, col = "blue", xlim = c(-2, 2), ylim = c(-2, 2))
# for (i in 1:n) {
#   for (j in i:n) {
#     if (adj_matrix[i, j] == 1) {
#       lines(x[c(i, j)], y[c(i, j)], col = "red", lwd = 2)
#       ##########%%%%%######################
#       # Set parameters
#       n_points <- 1000 # number of nodes
#       lambda <- 10 # mean number of points per unit area
#       torus_size <- 20 # size of the torus
#       p_rewire <- 0.1 # probability of rewiring
#       
#       # Generate a spatial distribution of points on a torus
#       x <- (torus_size * (runif(n_points) - 0.5)) + (torus_size/2)
#       y <- (torus_size * (runif(n_points) - 0.5)) + (torus_size/2)
#       
#       # Use a 2D Poisson process to determine the number of connections between nodes
#       distances <- as.matrix(dist(cbind(x, y)))
#       probabilities <- exp(-lambda * distances)
#       
#       # Randomly connect nodes based on probabilities
#       connections <- rbinom(n_points^2, 1, probabilities)
#       connections[lower.tri(connections)] <- 0
#       
#       # Convert connections to an adjacency matrix
#       adj_matrix <- as.matrix(connections)
#       diag(adj_matrix) <- 0
#       
#       # Implement small-world effects by rewiring edges with probability p_rewire
#       for (i in 1:n_points) {
#         for (j in (i+1):n_points) {
#           if (adj_matrix[i, j] == 1 & runif(1) < p_rewire) {
#             new_connection <- sample(1:n_points, 1, replace = TRUE)
#             adj_matrix[i, j] <- 0
#             adj_matrix[i, new_connection] <- 1
#           }
#         }
#       }
#       
#       # Plot the network
#       plot(x, y, xlab = "X", ylab = "Y", main = "Generated Network on a Torus")
#       for (i in 1:n_points) {
#         for (j in (i+1):n_points) {
#           if (adj_matrix[i, j] == 1) {
#             lines(c(x[i], x[j]), c(y[i], y[j]), col = "blue")
#           }
#         }
#       }
#       
#       #######%%%%%##############
#       # Define a function to generate a torus
#       generate_torus <- function(N, r, R) {
#         # N: number of nodes
#         # r: radius of small circle
#         # R: radius of big circle
#         
#         # Generate theta values
#         theta <- runif(N, 0, 2 * pi)
#         
#         # Generate x and y values
#         x <- (R + r * cos(theta)) * cos(theta)
#         y <- (R + r * cos(theta)) * sin(theta)
#         
#         return(data.frame(x, y))
#       }
#       
#       # Generate nodes on a torus
#       nodes <- generate_torus(1000, 0.1, 1)
#       
#       # Define a function to calculate the Euclidean distance between two nodes
#       distance <- function(x1, y1, x2, y2) {
#         return(sqrt((x1 - x2)^2 + (y1 - y2)^2))
#       }
#       
#       # Define a function to create a random connection between nodes
#       connect_nodes <- function(nodes, p) {
#         # nodes: data frame containing node positions
#         # p: probability of connection between two nodes
#         
#         N <- nrow(nodes)
#         edges <- c()
#         
#         for (i in 1:(N - 1)) {
#           for (j in (i + 1):N) {
#             if (runif(1) < p) {
#               edge <- c(i, j)
#               edges <- rbind(edges, edge)
#             }
#           }
#         }
#         
#         return(edges)
#       }
#       
#       # Create random connections between nodes
#       edges <- connect_nodes(nodes, 0.05)
#       
#       # Define a function to calculate the degree of each node
#       degree <- function(edges, N) {
#         # edges: matrix containing edge connections
#         # N: number of nodes
#         
#         degrees <- numeric(N)
#         
#         for (i in 1:N) {
#           degrees[i] <- sum(edges[, 1] == i | edges[, 2] == i)
#         }
#         
#         return(degrees)
#       }
#       
#       # Calculate the degree of each node
#       degrees <- degree(edges, nrow(nodes))
#       
#       # Define a function to rewire the network to generate small-world effects
#       rewire_network <- function(edges, degrees, p) {
#         # edges: matrix containing edge connections
#         # degrees: vector containing the degree of each node
#         # p: probability of rewiring an edge
#         
#         N <- length(degrees)
#         new_edges <- edges
#         
#         for (i in 1:nrow(edges)) {
#           if (runif(1) < p) {
#             n1 <- edges[i, 1]
#             n2 <- edges[i, 2]
#             while (TRUE) {
#               n3 <- sample(1:N
#                            #####################%%%%###########
#                            # Generative model on a torus
#                            
#                            # Set the seed for reproducibility
#                            set.seed(123)
#                            
#                            # Parameters
#                            n_nodes <- 1000 # number of nodes
#                            lambda <- 2 # rate parameter for Poisson process
#                            radius <- 100 # radius of the torus
#                            p_rewiring <- 0.1 # probability of rewiring in small world effects
#                            
#                            # Create a matrix to store the node positions
#                            positions <- matrix(NA, ncol = 2, nrow = n_nodes)
#                            
#                            # Create the first node at the origin
#                            positions[1, ] <- c(0, 0)
#                            
#                            # Create the rest of the nodes
#                            for (i in 2:n_nodes) {
#                              # Generate a random angle
#                              theta <- runif(1, 0, 2 * pi)
#                              
#                              # Generate a random radius
#                              r <- sqrt(-2 * log(runif(1))) / sqrt(lambda)
#                              
#                              # Calculate the x and y positions
#                              x <- r * cos(theta)
#                              y <- r * sin(theta)
#                              
#                              # Store the positions
#                              positions[i, ] <- c(x, y)
#                            }
#                            
#                            # Scale the node positions to the desired radius
#                            positions <- positions * radius / max(abs(positions))
#                            
#                            # Create the torus by wrapping the node positions around the edges
#                            positions[positions > radius] <- 2 * radius - positions[positions > radius]
#                            positions[positions < -radius] <- -2 * radius - positions[positions < -radius]
#                            
#                            # Create a matrix to store the edges
#                            edges <- matrix(NA, ncol = 2, nrow = 0)
#                            
#                            # Create the edges
#                            for (i in 1:(n_nodes - 1)) {
#                              # Get the neighbors of node i
#                              neighbors <- c(i + 1:n_nodes)
#                              
#                              # Calculate the distances between node i and all other nodes
#                              distances <- sqrt(rowSums((positions[neighbors, ] - positions[i, ])^2))
#                              
#                              # Connect node i to its nearest neighbors
#                              nearest_neighbors <- neighbors[order(distances)][1:floor(sqrt(n_nodes))]
#                              new_edges <- cbind(rep(i, length(nearest_neighbors)), nearest_neighbors)
#                              edges <- rbind(edges, new_edges)
#                              
#                              # Implement small world effects by rewiring some edges
#                              for (j in nearest_neighbors) {
#                                if (runif(1) < p_rewiring) {
#                                  rewired_neighbor <- sample(neighbors[neighbors != j], 1)
#                                  edges[edges[, 1] == i & edges[, 2] == j, 2] <- rewired_neighbor
#                                  edges[edges[, 2] == i & edges[, 1] == j, 1] <- rewired_neighbor
#                                }
#                              }
#                            }
#                            
#                            # Plot the generated network
#                            plot(positions[, 1], positions[, 2], xlab = "X", ylab = "Y", main = "Generative Model on a Torus",
#                                 
#                                 
#                                 ########%%%#########################%
#                                 #Define a function to generate a 2D Poisson process on a torus
#                                 generate_poisson_process <- function(lambda, L){
#                                   # Define the torus boundary
#                                   boundary <- c(0,L)
#                                   
#                                   # Generate Poisson process on the square
#                                   x <- runif(poisson(lambda * L^2), 0, L)
#                                   y <- runif(poisson(lambda * L^2), 0, L)
#                                   
#                                   # Implement torus topology by setting points outside the boundary back to the boundary
#                                   x[x<0] <- boundary[1]
#                                   x[x>boundary[2]] <- boundary[2]
#                                   y[y<0] <- boundary[1]
#                                   y[y>boundary[2]] <- boundary[2]
#                                   
#                                   return(cbind(x, y))
#                                 }
#                                 
#                                 # Define a function to generate connections between nodes on a torus
#                                 generate_connections <- function(points, k){
#                                   # Initialize the adjacency matrix
#                                   adjacency_matrix <- matrix(0, nrow = length(points), ncol = length(points))
#                                   
#                                   # Calculate the Euclidean distance between each pair of points
#                                   distances <- as.matrix(dist(points))
#                                   
#                                   # Define a function to handle the torus topology
#                                   handle_torus_topology <- function(distance, L){
#                                     if(distance > L/2){
#                                       distance <- L - distance
#                                     }
#                                     return(distance)
#                                   }
#                                   
#                                   # Iterate over each pair of points
#                                   for(i in 1:(length(points) - 1)){
#                                     for(j in (i + 1):length(points)){
#                                       # Handle the torus topology for the distance between the two points
#                                       distance <- handle_torus_topology(distances[i,j], L = max(points))
#                                       
#                                       # Connect the two points with probability proportional to their distance
#                                       if(runif(1) < (1 / (distance + 1))^k){
#                                         adjacency_matrix[i,j] <- 1
#                                         adjacency_matrix[j,i] <- 1
#                                       }
#                                     }
#                                   }
#                                   
#                                   return(adjacency_matrix)
#                                 }
#                                 
#                                 # Define a function to generate a scale-free network on a torus
#                                 generate_scale_free_network <- function(N, lambda, k, L){
#                                   # Generate the Poisson process on the torus
#                                   points <- generate_poisson_process(lambda, L)
#                                   
#                                   # Generate the initial connections on the torus
#                                   adjacency_matrix <- generate_connections(points, k)
#                                   
#                                   # Initialize the degree sequence
#                                   degree_sequence <- rowSums(adjacency_matrix)
#                                   
#                                   # Iterate until the number of nodes reaches N
#                                   while(length(degree_sequence) < N){
#                                     # Generate a new node with the Poisson process
#                                     new_point <- generate_poisson_process(lambda, L)
#                                     
#                                     # Generate connections for the new node
#                                     new_adjacency_matrix <- generate_connect
#                                     
#                                     ######%%%%%###################
#                                     # Generate Initial Nodes
#                                     generate_initial_nodes <- function(n, lambda) {
#                                       # Create a 2D Poisson Point Process with parameter lambda
#                                       x <- rpois(n, lambda)
#                                       y <- rpois(n, lambda)
#                                       data.frame(x, y)
#                                     }
#                                     
#                                     # Random Connection of Nodes
#                                     connect_random_nodes <- function(nodes, p) {
#                                       # Connect each node to p*n other nodes randomly
#                                       edges <- c()
#                                       for (i in 1:nrow(nodes)) {
#                                         connected_nodes <- sample(1:nrow(nodes), size = floor(p*nrow(nodes)), replace = TRUE)
#                                         for (j in connected_nodes) {
#                                           if (i != j) {
#                                             edges <- rbind(edges, c(i, j))
#                                           }
#                                         }
#                                       }
#                                       as.matrix(edges)
#                                     }
#                                     
#                                     # Scale Free Degree Distribution
#                                     create_scale_free_network <- function(n, m0, m, p) {
#                                       # Create the initial nodes
#                                       nodes <- generate_initial_nodes(n, p)
#                                       # Connect the initial nodes randomly
#                                       edges <- connect_random_nodes(nodes, m0/n)
#                                       # Add m-m0 edges using linear preferential attachment
#                                       for (i in 1:(m-m0)) {
#                                         attach_probs <- degree(edges)
#                                         attach_probs <- attach_probs/sum(attach_probs)
#                                         attach_node <- sample(1:nrow(nodes), size = 1, prob = attach_probs)
#                                         new_node <- nrow(nodes) + 1
#                                         nodes <- rbind(nodes, c(0,0))
#                                         edges <- rbind(edges, c(attach_node, new_node))
#                                       }
#                                       edges
#                                     }
#                                     
#                                     # Small World Effects
#                                     create_small_world_network <- function(n, m, p, k) {
#                                       # Create a regular lattice network
#                                       edges <- create_lattice_network(n, m, p)
#                                       # Rewire each edge with probability k
#                                       for (i in 1:nrow(edges)) {
#                                         rewire_prob <- runif(1)
#                                         if (rewire_prob < k) {
#                                           node1 <- edges[i,1]
#                                           node2 <- edges[i,2]
#                                           new_node2 <- sample(1:n, size = 1, replace = TRUE)
#                                           while (new_node2 == node1) {
#                                             new_node2 <- sample(1:n, size = 1, replace = TRUE)
#                                           }
#                                           edges[i,2] <- new_node2
#                                         }
#                                       }
#                                       edges
#                                     }
#                                     
#                                     # Example Usage
#                                     set.seed(123)
#                                     n <- 100
#                                     m0 <- 2
#                                     m <- n*m0
#                                     p <- 0.01
#                                     k <- 0.1
#                                     edges <- create_small_world_network(n, m, p, k)
#                                     #############%%%%%%%%%##########%%#
#                                     # Function to generate a 2D Poisson Point Process
#                                     generate_points <- function(lambda, size){
#                                       # Generate Poisson random variables for x and y coordinates
#                                       x <- rpois(size, lambda)
#                                       y <- rpois(size, lambda)
#                                       
#                                       # Return a data frame of points
#                                       return(data.frame(x = x, y = y))
#                                     }
#                                     
#                                     # Function to connect nodes based on the Euclidean distance
#                                     connect_nodes <- function(points, radius){
#                                       # Get the number of points
#                                       n <- nrow(points)
#                                       
#                                       # Initialize an empty adjacency matrix
#                                       adjacency_matrix <- matrix(0, n, n)
#                                       
#                                       # Iterate over all points and connect nodes if they are within the radius
#                                       for(i in 1:n){
#                                         for(j in (i+1):n){
#                                           distance <- sqrt((points[i, 1] - points[j, 1])^2 + (points[i, 2] - points[j, 2])^2)
#                                           if(distance <= radius){
#                                             adjacency_matrix[i, j] <- 1
#                                             adjacency_matrix[j, i] <- 1
#                                           }
#                                         }
#                                       }
#                                       
#                                       # Return the adjacency matrix
#                                       return(adjacency_matrix)
#                                     }
#                                     
#                                     # Function to generate a scale-free degree distribution
#                                     generate_scale_free_distribution <- function(n, m0, m){
#                                       # Initialize the degree sequence
#                                       degree_sequence <- rep(0, n)
#                                       
#                                       # Add the initial nodes with m0 edges
#                                       degree_sequence[1:m0] <- m
#                                       
#                                       # Iteratively add nodes and rewire edges
#                                       for(i in (m0+1):n){
#                                         for(j in 1:(m/2)){
#                                           k <- sample(1:(i-1), 1, prob = degree_sequence[1:(i-1)] / sum(degree_sequence[1:(i-1)]))
#                                           degree_sequence[k] <- degree_sequence[k] + 1
#                                           degree_sequence[i] <- degree_sequence[i] + 1
#                                         }
#                                       }
#                                       
#                                       # Return the degree sequence
#                                       return(degree_sequence)
#                                     }
#                                     
#                                     # Function to generate a small-world network
#                                     generate_small_world_network <- function(n, k, p){
#                                       # Generate a regular lattice network
#                                       adjacency_matrix <- matrix(0, n, n)
#                                       for(i in 1:n){
#                                         for(j in (i+1):(i+k/2)){
#                                           if(j <= n){
#                                             adjacency_matrix[i, j] <- 1
#                                             adjacency_matrix[j, i] <- 1
#                                           }
#                                         }
#                                       }
#                                       
#                                       # Re-wire edges with probability p
#                                       for(i in 1:n){
#                                         for(j in (i+1):n){
#                                           if(adjacency_matrix[i, j] == 1){
#                                             if(runif(1) < p){
#                                               k
#                                               ##########%%%%%###################%
#                                               # Generate node positions according to a 2D Poisson Point Process
#                                               generate_node_positions <- function(n, intensity){
#                                                 # Generate n points in 2D space with intensity parameter 'intensity'
#                                                 x <- rpois(n, intensity)
#                                                 y <- rpois(n, intensity)
#                                                 return(cbind(x, y))
#                                               }
#                                               
#                                               # Generate random connections between nodes
#                                               generate_random_connections <- function(n, p){
#                                                 # Generate a random adjacency matrix of size n x n with probability p
#                                                 # for each node to be connected to another node
#                                                 adjacency_matrix <- matrix(0, n, n)
#                                                 for (i in 1:(n-1)) {
#                                                   for (j in (i+1):n) {
#                                                     if (runif(1) < p) {
#                                                       adjacency_matrix[i, j] <- 1
#                                                       adjacency_matrix[j, i] <- 1
#                                                     }
#                                                   }
#                                                 }
#                                                 return(adjacency_matrix)
#                                               }
#                                               
#                                               # Generate a scale-free network by preferential attachment
#                                               generate_scale_free_network <- function(n, m0){
#                                                 # Generate a scale-free network with n nodes and initial m0 connections
#                                                 adjacency_matrix <- matrix(0, n, n)
#                                                 degree_sequence <- rep(0, n)
#                                                 for (i in 1:(m0-1)) {
#                                                   adjacency_matrix[i, i+1] <- 1
#                                                   adjacency_matrix[i+1, i] <- 1
#                                                   degree_sequence[i] <- 1
#                                                   degree_sequence[i+1] <- 1
#                                                 }
#                                                 for (i in (m0+1):n) {
#                                                   prob <- degree_sequence[1:(i-1)]/sum(degree_sequence[1:(i-1)])
#                                                   new_connections <- sample(1:(i-1), m0, replace = TRUE, prob = prob)
#                                                   for (j in new_connections) {
#                                                     adjacency_matrix[i, j] <- 1
#                                                     adjacency_matrix[j, i] <- 1
#                                                     degree_sequence[i] <- degree_sequence[i] + 1
#                                                     degree_sequence[j] <- degree_sequence[j] + 1
#                                                   }
#                                                 }
#                                                 return(adjacency_matrix)
#                                               }
#                                               
#                                               # Generate a small-world network by rewiring with probability p
#                                               generate_small_world_network <- function(adjacency_matrix, p){
#                                                 n <- dim(adjacency_matrix)[1]
#                                                 for (i in 1:(n-1)) {
#                                                   for (j in (i+1):n) {
#                                                     if (adjacency_matrix[i, j] == 1 && runif(1) < p) {
#                                                       adjacency_matrix[i, j] <- 0
#                                                       adjacency_matrix[j, i] <- 0
#                                                       new_connection_i <- sample(1:(n-2), 1) + (new_connection_i
#                                                                                                 
#                                                                                                 #########%%%%#####################
#                                                                                                 
#                                                                                                 # Generate a 2D Poisson Point Process
#                                                                                                 poisson_point_process <- function(lambda, n, min, max) {
#                                                                                                   # Generate a uniform distribution over the given range
#                                                                                                   unif_dist <- runif(n, min, max)
#                                                                                                   # Calculate the exponential distribution with rate lambda
#                                                                                                   exp_dist <- rexp(n, lambda)
#                                                                                                   # Combine the uniform and exponential distributions to get x and y coordinates
#                                                                                                   x <- unif_dist
#                                                                                                   y <- exp_dist / lambda
#                                                                                                   return(cbind(x, y))
#                                                                                                 }
#                                                                                                 
#                                                                                                 # Connect nodes with random edges
#                                                                                                 random_connection <- function(n, nodes) {
#                                                                                                   edges <- matrix(ncol = 2, nrow = n)
#                                                                                                   for (i in 1:(n-1)) {
#                                                                                                     for (j in (i+1):n) {
#                                                                                                       if (runif(1) < 0.5) {
#                                                                                                         edges[i, j] <- 1
#                                                                                                         edges[j, i] <- 1
#                                                                                                       }
#                                                                                                     }
#                                                                                                   }
#                                                                                                   return(edges)
#                                                                                                 }
#                                                                                                 
#                                                                                                 # Generate a scale-free degree distribution
#                                                                                                 scale_free_degree_distribution <- function(n, m0, m) {
#                                                                                                   # Initialize the degree distribution
#                                                                                                   degree_distribution <- rep(0, n)
#                                                                                                   # Calculate the probability distribution
#                                                                                                   p <- (1:n)^(-3/2)
#                                                                                                   p <- p / sum(p)
#                                                                                                   # Generate the degree sequence
#                                                                                                   for (i in 1:n) {
#                                                                                                     degree_distribution[i] <- sum(rpois(m0, p[i]*m))
#                                                                                                   }
#                                                                                                   return(degree_distribution)
#                                                                                                 }
#                                                                                                 
#                                                                                                 # Generate a small-world graph
#                                                                                                 small_world_graph <- function(n, nodes, edges, degree_distribution) {
#                                                                                                   # Re-wire the edges to introduce small-world effects
#                                                                                                   for (i in 1:n) {
#                                                                                                     for (j in (i+1):n) {
#                                                                                                       if (edges[i, j] == 1) {
#                                                                                                         # Calculate the distance between nodes i and j
#                                                                                                         d <- sqrt((nodes[i, 1] - nodes[j, 1])^2 + (nodes[i, 2] - nodes[j, 2])^2)
#                                                                                                         # Rewire the edge with probability p
#                                                                                                         p <- 1 / (1 + exp(-d))
#                                                                                                         if (runif(1) < p) {
#                                                                                                           # Choose a new node at random
#                                                                                                           k <- sample(1:n, 1)
#                                                                                                           while (k == i || k == j) {
#                                                                                                             k <- sample(1:n, 1)
#                                                                                                           }
#                                                                                                           # Connect nodes i and k
#                                                                                                           edges[i, k] <- 1
#                                                                                                           edges[k, i] <- 1
#                                                                                                           # Remove the edge between nodes i and j
#                                                                                                           edges[i, j] <- 0
#                                                                                                           edges[j, i] <- 0
#                                                                                                         }
#                                                                                                       }
#                                                                                                     }
#                                                                                                   }
#                                                                                                   return(edges)
#                                                                                                 }
#                                                                                                 
#                                                                                                 # Define the parameters for the generative model
#                                                                                                 n <- 100 # number of nodes
#                                                                                                 m0 <- 2 # initial number of edges per node
#                                                                                                 m <- 2 # number of edges to add in the scale-free step
#                                                                                                 lambda <- 0
#                                                                                                 
#                                                                                                 ###########%%%%%%##################
#                                                                                                 
#                                                                                                 # Define the parameters for the generative model
#                                                                                                 num_nodes <- 100 # number of nodes in the graph
#                                                                                                 scale_factor <- 0.1 # controls the spatial distribution of points
#                                                                                                 k_value <- 2 # number of connections per node
#                                                                                                 power_law_exponent <- 2 # degree distribution follows power-law with this exponent
#                                                                                                 small_world_factor <- 0.1 # controls the amount of rewiring for small-world effects
#                                                                                                 
#                                                                                                 # Generate node positions using a 2D Poisson Point Process
#                                                                                                 library(spatstat)
#                                                                                                 ppp <- rpoispp(num_nodes, lambda = scale_factor)
#                                                                                                 
#                                                                                                 # Define a function to calculate distances between nodes
#                                                                                                 node_distance <- function(i, j) {
#                                                                                                   sqrt((ppp$x[i] - ppp$x[j])^2 + (ppp$y[i] - ppp$y[j])^2)
#                                                                                                 }
#                                                                                                 
#                                                                                                 # Initialize the adjacency matrix with all zeros
#                                                                                                 adj_matrix <- matrix(0, nrow = num_nodes, ncol = num_nodes)
#                                                                                                 
#                                                                                                 # Define a function to make connections between nodes
#                                                                                                 make_connections <- function(node_index) {
#                                                                                                   # Create a list of nodes that haven't been connected to node_index
#                                                                                                   unconnected_nodes <- which(colSums(adj_matrix) == 0)
#                                                                                                   unconnected_nodes <- unconnected_nodes[unconnected_nodes != node_index]
#                                                                                                   
#                                                                                                   # Sort the unconnected nodes based on distance from node_index
#                                                                                                   distances <- sapply(unconnected_nodes, function(x) node_distance(node_index, x))
#                                                                                                   distances_sorted <- sort(distances)
#                                                                                                   closest_nodes <- unconnected_nodes[order(distances)]
#                                                                                                   
#                                                                                                   # Connect node_index to k closest nodes
#                                                                                                   for (i in 1:k_value) {
#                                                                                                     adj_matrix[node_index, closest_nodes[i]] <- 1
#                                                                                                     adj_matrix[closest_nodes[i], node_index] <- 1
#                                                                                                   }
#                                                                                                 }
#                                                                                                 
#                                                                                                 # Connect each node to k closest nodes
#                                                                                                 for (i in 1:num_nodes) {
#                                                                                                   make_connections(i)
#                                                                                                 }
#                                                                                                 
#                                                                                                 # Implement small-world rewiring
#                                                                                                 rewiring_prob <- small_world_factor / (num_nodes - 1)
#                                                                                                 for (i in 1:num_nodes) {
#                                                                                                   for (j in (i+1):num_nodes) {
#                                                                                                     if (adj_matrix[i, j] == 1) {
#                                                                                                       if (runif(1) < rewiring_prob) {
#                                                                                                         adj_matrix[i, j] <- 0
#                                                                                                         adj_matrix[j, i] <- 0
#                                                                                                         make_connections(i)
#                                                                                                         make_connections(j)
#                                                                                                       }
#                                                                                                     }
#                                                                                                   }
#                                                                                                 }
#                                                                                                 
#                                                                                                 # Calculate the degree distribution and plot it
#                                                                                                 degree_distribution <- sort(colSums(adj_matrix), decreasing = TRUE)
#                                                                                                 plot(degree_distribution, log = "xy", xlab = "Node Index", ylab = "Degree",
#                                                                                                      main = "Degree Distribution")
#                                                                                                 #######%%%#########################
#                                                                                                 
#                                                                                                 # Generate a 2D Poisson Point Process
#                                                                                                 generate_poisson_process <- function(lambda, size) {
#                                                                                                   points <- matrix(rpois(2 * size, lambda), ncol = 2)
#                                                                                                   return(points)
#                                                                                                 }
#                                                                                                 
#                                                                                                 # Generate random connections between nodes
#                                                                                                 generate_random_connections <- function(points) {
#                                                                                                   size <- nrow(points)
#                                                                                                   connections <- matrix(0, nrow = size, ncol = size)
#                                                                                                   
#                                                                                                   for (i in 1:(size - 1)) {
#                                                                                                     for (j in (i + 1):size) {
#                                                                                                       connections[i, j] <- connections[j, i] <- sample(c(0, 1), 1)
#                                                                                                     }
#                                                                                                   }
#                                                                                                   
#                                                                                                   return(connections)
#                                                                                                 }
#                                                                                                 
#                                                                                                 # Generate a scale-free degree distribution
#                                                                                                 generate_scale_free_degree_distribution <- function(connections) {
#                                                                                                   degree_distribution <- rowSums(connections)
#                                                                                                   size <- length(degree_distribution)
#                                                                                                   for (i in 1:size) {
#                                                                                                     for (j in 1:degree_distribution[i]) {
#                                                                                                       node <- sample(1:size, 1, prob = degree_distribution / sum(degree_distribution))
#                                                                                                       connections[i, node] <- connections[node, i] <- 1
#                                                                                                     }
#                                                                                                   }
#                                                                                                   
#                                                                                                   return(connections)
#                                                                                                 }
#                                                                                                 
#                                                                                                 # Generate small-world effects
#                                                                                                 generate_small_world_effects <- function(connections) {
#                                                                                                   size <- nrow(connections)
#                                                                                                   for (i in 1:(size - 1)) {
#                                                                                                     for (j in (i + 1):size) {
#                                                                                                       if (connections[i, j] == 1) {
#                                                                                                         rewire_prob <- runif(1)
#                                                                                                         if (rewire_prob <= 0.1) {
#                                                                                                           node <- sample(1:size, 1, prob = degree_distribution / sum(degree_distribution))
#                                                                                                           connections[i, j] <- connections[j, i] <- 0
#                                                                                                           connections[i, node] <- connections[node, i] <- 1
#                                                                                                         }
#                                                                                                       }
#                                                                                                     }
#                                                                                                   }
#                                                                                                   
#                                                                                                   return(connections)
#                                                                                                 }
#                                                                                                 
#                                                                                                 # Generate a graph with the desired properties
#                                                                                                 generate_graph <- function(lambda, size) {
#                                                                                                   points <- generate_poisson_process(lambda, size)
#                                                                                                   connections <- generate_random_connections(points)
#                                                                                                   connections <- generate_scale_free_degree_distribution(connections)
#                                                                                                   connections <- generate_small_world_effects(connections)
#                                                                                                   
#                                                                                                   return(list(points = points, connections = connections))
#                                                                                                 }
#                                                                                                 
#                                                                                                 # Example usage
#                                                                                                 set.seed(123)
#                                                                                                 graph <- generate_graph(lambda = 10, size = 100)
#                                                                                                 points <- graph$points
#                                                                                                 connections <- graph$connections
#                                                                                                 
#                                                                                                 ###################%%%#############
#                                                                                                 
#                                                                                                 # Function to generate a 2D Poisson point process
#                                                                                                 generate_poisson_points <- function(lambda, xmin, xmax, ymin, ymax) {
#                                                                                                   # Generate Poisson points
#                                                                                                   n <- rpois(1, lambda * (xmax - xmin) * (ymax - ymin))
#                                                                                                   points <- matrix(runif(2 * n, min = c(xmin, ymin), max = c(xmax, ymax)), ncol = 2)
#                                                                                                   colnames(points) <- c("x", "y")
#                                                                                                   points
#                                                                                                 }
#                                                                                                 
#                                                                                                 # Function to generate the network
#                                                                                                 generate_network <- function(points, n_nodes, p_conn, avg_degree) {
#                                                                                                   # Generate the edges between nodes
#                                                                                                   edges <- matrix(numeric(0), ncol = 2)
#                                                                                                   for (i in 1:(n_nodes - 1)) {
#                                                                                                     for (j in (i + 1):n_nodes) {
#                                                                                                       if (runif(1) < p_conn) {
#                                                                                                         edges <- rbind(edges, c(i, j))
#                                                                                                       }
#                                                                                                     }
#                                                                                                   }
#                                                                                                   
#                                                                                                   # Generate the node degrees
#                                                                                                   degree_sequence <- rt(n_nodes, df = avg_degree / 2, ncp = avg_degree / 2)
#                                                                                                   while (sum(degree_sequence) %% 2 != 0) {
#                                                                                                     degree_sequence[1] <- degree_sequence[1] + 1
#                                                                                                   }
#                                                                                                   
#                                                                                                   # Shuffle the edges
#                                                                                                   edges <- edges[sample(nrow(edges)), ]
#                                                                                                   
#                                                                                                   # Connect nodes with the desired degree
#                                                                                                   i <- 1
#                                                                                                   for (d in degree_sequence) {
#                                                                                                     for (j in 1:(d / 2)) {
#                                                                                                       edges <- rbind(edges, c(i, nrow(edges) + j))
#                                                                                                     }
#                                                                                                     i <- i + 1
#                                                                                                   }
#                                                                                                   
#                                                                                                   # Apply small world effect
#                                                                                                   n_rewires <- round(n_nodes * avg_degree * p_conn)
#                                                                                                   for (i in 1:n_rewires) {
#                                                                                                     edge1 <- sample(1:nrow(edges), 1)
#                                                                                                     edge2 <- sample(1:nrow(edges), 1)
#                                                                                                     while (edges[edge1, 1] == edges[edge2, 1] || edges[edge1, 1] == edges[edge2, 2] || 
#                                                                                                            edges[edge1, 2] == edges[edge2, 1] || edges[edge1, 2] == edges[edge2, 2]) {
#                                                                                                       edge2 <- sample(1:nrow(edges), 1)
#                                                                                                     }
#                                                                                                     node1 <- sample(c(edges[edge1, 1], edges[edge1, 2]), 1)
#                                                                                                     node2 <- sample(c(edges[edge2, 1], edges[edge2, 2]), 1)
#                                                                                                     edges[edge1, ] <- c(node1, node2)
#                                                                                                     edges[edge2, ] <- c(node2, node1)
#                                                                                                   }
#                                                                                                   
#                                                                                                   # Plot the network
#                                                                                                   plot(points[, 1], points[, 2], xlim = c(min(points[, 1]) - 1, max(points[, 1]) +
#                                                                                                                                             
#                                                                                                                                             ###########%%%%%#########
#                                                                                                                                           
#                                                                                                                                           # Generate node positions using a 2D PPP
#                                                                                                                                           generate_positions <- function(n, lambda) {
#                                                                                                                                             x <- rpois(n, lambda)
#                                                                                                                                             y <- rpois(n, lambda)
#                                                                                                                                             return(cbind(x, y))
#                                                                                                                                           }
#                                                                                                                                           
#                                                                                                                                           # Generate graph with random connections
#                                                                                                                                           generate_random_graph <- function(n, p) {
#                                                                                                                                             adj_matrix <- matrix(0, n, n)
#                                                                                                                                             for (i in 1:(n - 1)) {
#                                                                                                                                               for (j in (i + 1):n) {
#                                                                                                                                                 if (runif(1) < p) {
#                                                                                                                                                   adj_matrix[i, j] <- 1
#                                                                                                                                                   adj_matrix[j, i] <- 1
#                                                                                                                                                 }
#                                                                                                                                               }
#                                                                                                                                             }
#                                                                                                                                             return(adj_matrix)
#                                                                                                                                           }
#                                                                                                                                           
#                                                                                                                                           # Generate graph with scale-free degree distribution
#                                                                                                                                           generate_scale_free_graph <- function(n, m0, m) {
#                                                                                                                                             adj_matrix <- matrix(0, n, n)
#                                                                                                                                             degrees <- numeric(n)
#                                                                                                                                             for (i in 1:m0) {
#                                                                                                                                               adj_matrix[i, i + 1:(i + m)] <- 1
#                                                                                                                                               adj_matrix[i + 1:(i + m), i] <- 1
#                                                                                                                                               degrees[i] <- m0
#                                                                                                                                               degrees[i + 1:(i + m)] <- degrees[i + 1:(i + m)] + 1
#                                                                                                                                             }
#                                                                                                                                             for (i in (m0 + 1):n) {
#                                                                                                                                               prob <- degrees / sum(degrees)
#                                                                                                                                               neighbors <- sample(1:n, m, replace = TRUE, prob = prob)
#                                                                                                                                               for (j in neighbors) {
#                                                                                                                                                 adj_matrix[i, j] <- 1
#                                                                                                                                                 adj_matrix[j, i] <- 1
#                                                                                                                                                 degrees[j] <- degrees[j] + 1
#                                                                                                                                               }
#                                                                                                                                               degrees[i] <- m
#                                                                                                                                             }
#                                                                                                                                             return(adj_matrix)
#                                                                                                                                           }
#                                                                                                                                           
#                                                                                                                                           # Generate graph with small-world effect
#                                                                                                                                           generate_small_world_graph <- function(n, k, p, q) {
#                                                                                                                                             adj_matrix <- generate_random_graph(n, p)
#                                                                                                                                             for (i in 1:n) {
#                                                                                                                                               neighbors <- sort(c(which(adj_matrix[i, ] == 1), i))
#                                                                                                                                               if (length(neighbors) > k) {
#                                                                                                                                                 neighbors <- neighbors[-(1:(length(neighbors) - k))]
#                                                                                                                                               }
#                                                                                                                                               for (j in neighbors) {
#                                                                                                                                                 if (j > i) {
#                                                                                                                                                   if (runif(1) < q) {
#                                                                                                                                                     new_neighbor <- sample(1:n, 1, replace = TRUE, exclude = c(i, neighbors))
#                                                                                                                                                     adj_matrix[i, new_neighbor] <- 1
#                                                                                                                                                     adj_matrix[new_neighbor, i] <- 1
#                                                                                                                                                   }
#                                                                                                                                                 }
#                                                                                                                                               }
#                                                                                                                                             }
#                                                                                                                                             return(adj_matrix)
#                                                                                                                                           }
#                                                                                                                                           
#                                                                                                                                           # Generate the graph
#                                                                                                                                           n <- 1000 # number of nodes
#                                                                                                                                           lambda <- 10 # parameter for 2D PPP
#                                                                                                                                           positions <- generate_positions(n, lambda)
#                                                                                                                                           

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sample code 12
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# install necessary packages
install.packages("igraph")
install.packages("spatstat")

# load packages
library(igraph)
library(spatstat)

# Generative Model

# Generate 2D Poisson Point Process
poisson_process <- rpp(n=1000, win=owin(c(0,1),c(0,1)))

# Generate scale-free degree distribution
g <- sample_pa(1000, m=1, directed=FALSE)

# Add small-world effect
g <- rewire(g, probability=0.1, loops=FALSE)

# Function to calculate Euclidean distance between two points
calculate_distance <- function(x1, y1, x2, y2) {
  return (sqrt((x2 - x1)^2 + (y2 - y1)^2))
}

# Assign spatial positions to nodes in the network
V(g)$x <- coordinates(poisson_process)[,1]
V(g)$y <- coordinates(poisson_process)[,2]

# Calculate the spatial distance between each pair of nodes
spatial_distances <- as.matrix(distances(g, v=V(g), weights=NULL, directed=FALSE))

# Classification Model

# Function to calculate the likelihood of a network being generated by the generative model
calculate_likelihood <- function(network, spatial_distances, model) {
  # Calculate the likelihood of the network having the same spatial distribution as the model
  spatial_likelihood <- ...
  
  # Calculate the likelihood of the network having the same degree distribution as the model
  degree_likelihood <- ...
  
  # Calculate the likelihood of the network having the same small-world effect as the model
  small_world_likelihood <- ...
  
  # Combine the three likelihoods to get the overall likelihood of the network being generated by the model
  overall_likelihood <- ...
  
  return (overall_likelihood)
}

# Classify an empirical network based on the generative model
classify_network <- function(network) {
  # Calculate the likelihood of the network being generated by each model
  likelihoods <- c(calculate_likelihood(network, spatial_distances, "model1"), 
                   calculate_likelihood(network, spatial_distances, "model2"), 
                   calculate_likelihood(network, spatial_distances, "model3"))
  
  # Identify the model with the highest likelihood
  classified_model <- names(likelihoods)[which.max(likelihoods)]
  
  return (classified_model)
}

############%%%%################
library(igraph)
library(spatstat)

# Generate 2D Poisson Point Process
generate_ppp <- function(lambda, size){
  ppp(runifpoint(lambda * size^2), window = owin(c(0, size), c(0, size)))
}

# Generate Random Graph from 2D PPP
generate_random_graph <- function(lambda, size, p){
  ppp <- generate_ppp(lambda, size)
  g <- graph.ppp(ppp, radius = size * p, fun = mean, apply = mean)
  g
}

# Generate Scale-Free Network
generate_scale_free_network <- function(n, m0, m, p){
  g <- barabasi.game(n, m0, m, directed = F)
  g <- simplify(g)
  g <- rewire(g, p)
  g
}

# Generate Small-World Network
generate_small_world_network <- function(n, k, p){
  g <- watts.strogatz.game(1, n, k, p)
  g
}

# Tune parameters to generate different models
lambda <- 0.1
size <- 100
p_ppp <- 0.01

m0 <- 1
m <- 1
p_sf <- 0.01

k <- 10
p_sw <- 0.1

g_ppp <- generate_random_graph(lambda, size, p_ppp)
g_sf <- generate_scale_free_network(n = 1000, m0 = m0, m = m, p = p_sf)
g_sw <- generate_small_world_network(n = 1000, k = k, p = p_sw)

# Plot the generated networks
plot(g_ppp, vertex.size = 3, vertex.label = "")
plot(g_sf, vertex.size = 3, vertex.label = "")
plot(g_sw, vertex.size = 3, vertex.label = "")

########%%%%%%###############
library(caret)

# Prepare the data
data <- data.frame(type = c(rep("PPP", nrow(g_ppp)), 
                            rep("SF", nrow(g_sf)), 
                            rep("SW", nrow(g_sw))), 
                   features = c(degree(g_ppp), degree(g_sf), degree(g_sw)))

# Train and Test the Model
set.seed(123)
index <- createDataPartition(data$type, p = 0.7, list = F)
training <- data[index,]
testing <- data[-index,]
model <- train(type ~ features, data = training, method = "rf")
predictions <- predict(model, testing)
confusionMatrix(predictions, testing$type)

######%%%%%######################
# Generative Model for Network Formation

# Generate a spatial distribution of nodes
set.seed(123)
N <- 100 # number of nodes

# Generate x and y positions for nodes
x <- rnorm(N, mean = 0, sd = 10)
y <- rnorm(N, mean = 0, sd = 10)

# Calculate the Euclidean distances between nodes
distances <- matrix(0, nrow = N, ncol = N)
for (i in 1:N) {
  for (j in 1:N) {
    distances[i,j] <- sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2)
  }
}

# Generate links between nodes based on a 2D Poisson Point Process
p <- 0.1 # probability of forming a link
links <- matrix(0, nrow = N, ncol = N)
for (i in 1:N) {
  for (j in (i+1):N) {
    links[i,j] <- rbinom(1, 1, p = p * exp(-distances[i,j]))
    links[j,i] <- links[i,j]
  }
}

# Generate the degree distribution of nodes
degrees <- rowSums(links)

# Generate links between high degree nodes with a probability proportional to their degrees
p_high_degree <- 0.5 # probability of forming a link between high degree nodes
for (i in 1:N) {
  for (j in (i+1):N) {
    if (degrees[i] > mean(degrees) & degrees[j] > mean(degrees)) {
      links[i,j] <- links[i,j] + rbinom(1, 1, p = p_high_degree * degrees[i] * degrees[j])
      links[j,i] <- links[i,j]
    }
  }
}

# Plot the network
library(ggplot2)
network_data <- data.frame(x, y, links = as.vector(t(links)))
ggplot(network_data, aes(x, y)) +
  geom_point(aes(color = links)) +
  scale_color_gradient2(low = "red", high = "blue", mid = "white", midpoint = 0.5) +
  theme_classic() +
  ggtitle("Generated Network")

########%%%%%%%%%################

# Generate random points in 2d space
generate_random_points <- function(n, mean = c(0,0), sd = c(1,1)) {
  x <- rnorm(n, mean = mean[1], sd = sd[1])
  y <- rnorm(n, mean = mean[2], sd = sd[2])
  return(cbind(x, y))
}

# Generate adjacency matrix
generate_adjacency_matrix <- function(points, k, p) {
  n <- nrow(points)
  adjacency_matrix <- matrix(0, nrow = n, ncol = n)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      distance <- sqrt((points[i, 1] - points[j, 1])^2 + (points[i, 2] - points[j, 2])^2)
      if (distance <= k) {
        adjacency_matrix[i, j] <- 1
        adjacency_matrix[j, i] <- 1
      }
    }
  }
  return(adjacency_matrix)
}

# Generate scale-free network
generate_scale_free_network <- function(adjacency_matrix, m0, m) {
  n <- nrow(adjacency_matrix)
  for (i in (m0 + 1):n) {
    connections <- sample(1:(i - 1), size = m, replace = TRUE, prob = colSums(adjacency_matrix[1:(i - 1), 1:(i - 1)]))
    for (j in connections) {
      adjacency_matrix[i, j] <- 1
      adjacency_matrix[j, i] <- 1
    }
  }
  return(adjacency_matrix)
}

# Main function
generate_network <- function(n, mean = c(0,0), sd = c(1,1), k, p, m0, m) {
  points <- generate_random_points(n, mean, sd)
  adjacency_matrix <- generate_adjacency_matrix(points, k, p)
  scale_free_network <- generate_scale_free_network(adjacency_matrix, m0, m)
  return(scale_free_network)
}

#######%%%%%%################
library(MASS) # for mvrnorm function
library(matrixStats) # for rowMaxs function

# Define parameters
n <- 500 # number of nodes to generate
k <- 5 # number of initial connections for each new node
beta <- 0.5 # controls the preference for short spatial distances
lambda <- 100 # controls the density of nodes in the 2D space
gamma <- 2 # power law exponent for the degree distribution

# Generate 2D positions using Poisson point process
positions <- mvrnorm(n, mu=c(0,0), Sigma=diag(2)/lambda)
x <- positions[,1]
y <- positions[,2]

# Initialize adjacency matrix and degrees vector
adj_matrix <- matrix(0, n, n)
degrees <- rep(0, n)

# Add nodes and connections
for (i in 1:n) {
  # Connect to k nearest neighbors based on distance and degree
  dists <- sqrt((x-x[i])^2 + (y-y[i])^2)
  neighbors <- sort(order(dists, degrees), decreasing=FALSE)[2:(k+1)]
  adj_matrix[i,neighbors] <- 1
  adj_matrix[neighbors,i] <- 1
  degrees[c(i,neighbors)] <- degrees[c(i,neighbors)] + 1
  
  # Add new node with position based on spatial preference
  dist_weights <- exp(-beta*dists)
  new_x <- sum(x*dist_weights) / sum(dist_weights)
  new_y <- sum(y*dist_weights) / sum(dist_weights)
  positions <- rbind(positions, c(new_x, new_y))
  x <- c(x, new_x)
  y <- c(y, new_y)
  
  # Update adjacency matrix and degrees vector for new node
  adj_matrix <- rbind(adj_matrix, rep(0, n))
  adj_matrix <- cbind(adj_matrix, rep(0, n+1))
  degrees <- c(degrees, 0)
}

# Generate degree distribution following power law
degree_probs <- seq(1, n)^(-gamma)
degree_probs <- degree_probs / sum(degree_probs)
degrees <- sample(seq(1, n), n, replace=TRUE, prob=degree_probs)

# Add small world effect by rewiring edges
rewire_prob <- 0.2
for (i in 1:n) {
  for (j in (i+1):n) {
    if (adj_matrix[i,j] == 1 && runif(1) < rewire_prob) {
      # Find potential rewiring partners
      non_neighbors <- which(adj_matrix[i,] == 0 & degrees > 0)
      potential_partners <- non_neighbors[non_neighbors != j]
      if (length(potential_partners) > 0) {
        # Choose random rewiring partner and rewire
        partner <- sample(potential_partners, 1)
        adj_matrix[i,j] <- 0
        adj_matrix[j,i] <- 0
        adj_matrix[i,partner] <- 1
        adj_matrix[partner,i] <- 1
        degrees[c(i,j)] <- degrees[c(i,j)] - 1
        degrees[c(i,partner)] <- degrees[c(i,partner)] + 1
        degrees[c(j,partner)] <- degrees[c(j,partner)] + 1
      }
    }
  }
}

# Generate random network and classify using machine learning
empirical_network <- matrix(rbinom(n*n, 1, 0.5), n, n
                            
######%%%%######################
# library(pracma) # for distance calculation
#                             
#                             # Set model parameters
#                             n <- 500 # number of nodes
#                             m <- 5 # number of initial connections per node
#                             p <- 0.05 # probability of adding new connections
#                             alpha <- 1 # degree distribution exponent
#                             beta <- 0.5 # spatial distance exponent
#                             delta <- 0.1 # small-world rewiring probability
#                             
#                             # Initialize graph and spatial positions
#                             g <- vector(mode = "list", length = n)
#                             x <- runif(n)
#                             y <- runif(n)
#                             
#                             # Initialize first m nodes
#                             for (i in 1:m) {
#                               g[[i]] <- sample(setdiff(1:n, i), m-1)
#                             }
#                             
#                             # Add new nodes and connections
#                             for (i in (m+1):n) {
#                               # Calculate distances to existing nodes
#                               d <- dist2(cbind(x, y), cbind(x[g[[1:i-1]]], y[g[[1:i-1]]]))
#                               # Calculate probabilities based on distance and degree
#                               prob <- (1/d)^beta * (degree(g[1:i-1])^alpha)/sum(degree(g[1:i-1])^alpha)
#                               # Choose m nodes to connect to with highest probability
#                               new_connections <- sample(1:i-1, m, prob = prob, replace = FALSE)
#                               # Add new node and connections
#                               g[[i]] <- new_connections
#                               g[[new_connections]] <- c(g[[new_connections]], i)
#                               # Add small-world rewiring
#                               if (runif(1) < delta) {
#                                 old_connections <- setdiff(1:i-1, new_connections)
#                                 for (j in new_connections) {
#                                   if (runif(1) < delta) {
#                                     old_node <- sample(old_connections, 1)
#                                     g[[j]] <- setdiff(g[[j]], old_node)
#                                     g[[old_node]] <- setdiff(g[[old_node]], j)
#                                     g[[j]] <- c(g[[j]], old_node)
#                                     g[[old_node]] <- c(g[[old_node]], j)
#                                   }
#                                 }
#                               }
#                               # Add new spatial position
#                               x[i] <- runif(1, min(x), max(x))
#                               y[i] <- runif(1, min(y), max(y))
#                             }
#                             
#                             # Plot the resulting graph and spatial positions
#                             plot(x, y, pch = 20, col = "blue", cex = 0.5)
#                             for (i in 1:n) {
#                               lines(x[c(i, g[[i]])], y[c(i, g[[i]])], col = "gray", lwd = 0.5)
#                             }
#                             
#                             #####%%%%%%%%###################
#                             # Parameters
#                             n <- 100 # number of nodes
#                             alpha <- 2 # power law exponent for degree distribution
#                             p <- 0.05 # probability of connecting to an existing node
#                             r <- 0.1 # spatial decay parameter
#                             sigma <- 0.1 # standard deviation of spatial noise
#                             
#                             # Generate spatial positions for nodes
#                             x <- runif(n)
#                             y <- runif(n)
#                             for (i in 2:n) {
#                               # Generate new position by adding spatial noise to previous position
#                               x[i] <- x[i-1] + rnorm(1, mean=0, sd=sigma)
#                               y[i] <- y[i-1] + rnorm(1, mean=0, sd=sigma)
#                             }
#                             
#                             # Initialize adjacency matrix
#                             adj <- matrix(0, nrow=n, ncol=n)
#                             
#                             # Add edges iteratively
#                             for (i in 1:n) {
#                               # Calculate distances and degrees to other nodes
#                               dist <- sqrt((x[i]-x)^2 + (y[i]-y)^2)
#                               deg <- colSums(adj) + rowSums(adj)
#                               
#                               # Calculate probability of connecting to each other node
#                               prob <- p * (deg + 1)^(-alpha) * exp(-dist/r)
#                               
#                               # Randomly connect to some nodes
#                               for (j in 1:n) {
#                                 if (i != j & runif(1) < prob[j]) {
#                                   adj[i,j] <- 1
#                                   adj[j,i] <- 1
#                                 }
#                               }
#                             }
#                             
#                             # Calculate degree distribution
#                             deg <- colSums(adj)
#                             deg_dist <- table(deg)
#                             
#                             # Plot network
#                             plot(x, y, pch=19, cex=0.5)
#                             for (i in 1:n) {
#                               for (j in (i+1):n) {
#                                 if (adj[i,j] == 1) {
#                                   lines(c(x[i], x[j]), c(y[i], y[j]))
#                                 }
#                               }
#                             }
#                             
#                             ########%%%%%###################
#                             # Load required libraries
#                             library(igraph)
#                             library(tidyverse)
#                             library(randomcoloR)
#                             
#                             # Set up the generative models
#                             g1 <- sample_gnm(20, 30) # Example of a random graph model
#                             g2 <- sample_pa(20, directed = TRUE) # Example of a preferential attachment model
#                             g3 <- sample_degseq(degree(g1)) # Example of a configuration model
#                             
#                             # Create a list of the generative models and their names
#                             generative_models <- list(g1, g2, g3)
#                             model_names <- c("Random graph", "Preferential attachment", "Configuration model")
#                             
#                             # Create a function to calculate the similarity between the input graph and the generative models
#                             similarity <- function(input_graph, generative_models) {
#                               sim <- sapply(generative_models, function(g) {
#                                 vf2 <- graph.isomorphism.vf2(input_graph, g)
#                                 if (is.na(vf2)) {
#                                   return(0)
#                                 } else {
#                                   return(1)
#                                 }
#                               })
#                               return(sim)
#                             }
#                             
#                             # Load the empirical network
#                             empirical_network <- read.graph("path/to/empirical_network.txt", format = "edgelist")
#                             
#                             # Calculate the similarity between the empirical network and the generative models
#                             sim <- similarity(empirical_network, generative_models)
#                             
#                             # Identify the most similar model
#                             most_similar_model <- model_names[which.max(sim)]
#                             
#                             # Print the result
#                             cat("The empirical network is most similar to the", most_similar_model, "model.")
#                             
#                             ########%%%%%%######################
#                             # Load necessary libraries
#                             library(igraph)
#                             library(dplyr)
#                             library(caret)
#                             
#                             # Load data
#                             network_data <- read.graph("network_data.graphml", format = "graphml")
#                             
#                             # Calculate network statistics
#                             stats <- c(degree(network_data), closeness(network_data), betweenness(network_data))
#                             
#                             # Create labels for each network type
#                             labels <- rep(c("model1", "model2", "model3"), each = 100)
#                             
#                             # Combine statistics and labels into a data frame
#                             df <- data.frame(stats, labels)
#                             
#                             # Split data into training and testing sets
#                             train_idx <- createDataPartition(df$labels, p = 0.8, list = FALSE)
#                             train <- df[train_idx,]
#                             test <- df[-train_idx,]
#                             
#                             # Train a machine learning model
#                             model <- train(labels ~ ., data = train, method = "svmRadial")
#                             
#                             # Predict labels for test set
#                             predictions <- predict(model, newdata = test)
#                             
#                             # Evaluate model performance
#                             confusionMatrix(predictions, test$labels)
#                             
#                             ######%%%%%%##################
#                             # Load necessary packages
#                             library(igraph)
#                             library(randomForest)
#                             
#                             # Load the empirical network
#                             empirical_network <- read.graph("empirical_network.csv", format="edgelist", directed=FALSE)
#                             
#                             # Calculate network features
#                             network_features <- cbind(degree(empirical_network), betweenness(empirical_network), closeness(empirical_network))
#                             
#                             # Load the synthetic networks
#                             synthetic_networks_1 <- read.graph("synthetic_network_1.csv", format="edgelist", directed=FALSE)
#                             synthetic_networks_2 <- read.graph("synthetic_network_2.csv", format="edgelist", directed=FALSE)
#                             synthetic_networks_3 <- read.graph("synthetic_network_3.csv", format="edgelist", directed=FALSE)
#                             
#                             # Calculate network features for synthetic networks
#                             synthetic_features_1 <- cbind(degree(synthetic_networks_1), betweenness(synthetic_networks_1), closeness(synthetic_networks_1))
#                             synthetic_features_2 <- cbind(degree(synthetic_networks_2), betweenness(synthetic_networks_2), closeness(synthetic_networks_2))
#                             synthetic_features_3 <- cbind(degree(synthetic_networks_3), betweenness(synthetic_networks_3), closeness(synthetic_networks_3))
#                             
#                             # Combine the synthetic features into a single dataset
#                             synthetic_features <- rbind(synthetic_features_1, synthetic_features_2, synthetic_features_3)
#                             
#                             # Add labels to indicate which synthetic network each set of features came from
#                             labels <- c(rep(1, nrow(synthetic_features_1)), rep(2, nrow(synthetic_features_2)), rep(3, nrow(synthetic_features_3)))
#                             
#                             # Combine the synthetic features and labels into a single data frame
#                             synthetic_data <- data.frame(synthetic_features, labels)
#                             
#                             # Train a random forest classifier using the synthetic data
#                             model <- randomForest(labels ~ ., data=synthetic_data)
#                             
#                             # Use the trained model to predict which synthetic network the empirical network was generated from
#                             prediction <- predict(model, newdata=network_features)
#                             
#                             ####%%%%%%########%%%%###########
#                             # Load required libraries
#                             library(matrixStats)
#                             library(pracma)
#                             
#                             # Function to calculate network metrics
#                             network_metrics <- function(mat) {
#                               # Calculate degree distribution
#                               degree_distribution <- rowSums(mat)
#                               
#                               # Calculate average degree
#                               avg_degree <- mean(degree_distribution)
#                               
#                               # Calculate average clustering coefficient
#                               clustering_coefs <- sapply(1:nrow(mat), function(i) {
#                                 neighbors <- which(mat[i, ] == 1)
#                                 if (length(neighbors) < 2) return(0)
#                                 numerator <- sum(mat[neighbors, neighbors]) / 2
#                                 denominator <- choose(length(neighbors), 2)
#                                 numerator / denominator
#                               })
#                               avg_clustering_coef <- mean(clustering_coefs, na.rm = TRUE)
#                               
#                               # Return results
#                               return(list(degree_distribution = degree_distribution,
#                                           avg_degree = avg_degree,
#                                           avg_clustering_coef = avg_clustering_coef))
#                             }
#                             
#                             # Function to classify network as generated by one of the synthetic generative models
#                             classify_network <- function(mat) {
#                               # Calculate network metrics
#                               network_metrics <- network_metrics(mat)
#                               
#                               # Classify network based on network metrics
#                               if (network_metrics$avg_degree >= 8 &&
#                                   network_metrics$avg_clustering_coef >= 0.3) {
#                                 return("Small-World")
#                               } else if (network_metrics$avg_degree >= 8 &&
#                                          network_metrics$avg_clustering_coef < 0.3) {
#                                 return("Random")
#                               } else if (network_metrics$avg_degree < 8 &&
#                                          network_metrics$avg_clustering_coef >= 0.3) {
#                                 return("Regular")
#                               } else {
#                                 return("Unknown")
#                               }
#                             }
#                             
#                             ###########%%%%%%%%%%%#############
#                             # Load required libraries
#                             library(caret)
#                             library(randomForest)
#                             
#                             # Load the empirical network data and the corresponding labels
#                             empirical_data <- read.csv("empirical_data.csv")
#                             empirical_labels <- read.csv("empirical_labels.csv")
#                             
#                             # Convert the labels to a factor variable
#                             empirical_labels$label <- as.factor(empirical_labels$label)
#                             
#                             # Split the data into training and testing sets
#                             set.seed(123)
#                             train_index <- createDataPartition(empirical_labels$label, p = 0.7, list = FALSE)
#                             train_data <- empirical_data[train_index, ]
#                             train_labels <- empirical_labels[train_index, ]
#                             test_data <- empirical_data[-train_index, ]
#                             test_labels <- empirical_labels[-train_index, ]
#                             
#                             # Train a random forest classifier on the training data
#                             model <- randomForest(label ~ ., data = train_data, ntree = 100, importance = TRUE)
#                             
#                             # Make predictions on the test data
#                             predictions <- predict(model, newdata = test_data)
#                             
#                             # Evaluate the accuracy of the model
#                             accuracy <- mean(predictions == test_labels$label)
#                             print(accuracy)
#                             
#                             ##########%%%%%%%###############
#                             # Load necessary libraries
#                             library(caret)
#                             library(tidyverse)
#                             
#                             # Generate sample data
#                             set.seed(123)
#                             generative_model1 <- erdos.renyi.game(100, 0.1)
#                             generative_model2 <- barabasi.game(100, m=1)
#                             
#                             # Convert adjacency matrices to data frames
#                             df1 <- as.data.frame(get.adjacency(generative_model1))
#                             df2 <- as.data.frame(get.adjacency(generative_model2))
#                             
#                             # Bind the two data frames into one and add a response variable column
#                             df_all <- bind_rows(list(df1=df1, df2=df2)) %>%
#                               mutate(response = rep(c("model1", "model2"), each = nrow(df1)))
#                             
#                             # Split the data into training and testing sets
#                             trainIndex <- createDataPartition(df_all$response, p = 0.8, list = FALSE, times = 1)
#                             training <- df_all[ trainIndex,]
#                             testing  <- df_all[-trainIndex,]
#                             
#                             # Train a random forest classifier
#                             model_rf <- train(response ~ ., data = training, method = "rf")
#                             
#                             # Make predictions on the testing data
#                             predictions <- predict(model_rf, newdata = testing)
#                             
#                             # Evaluate the model's performance
#                             confusionMatrix(predictions, testing$response)
#                             
#                             #######%%%%%%##################
#                             
#                             # Load required libraries
#                             library(randomForest)
#                             library(caret)
#                             
#                             # Load the real network data and synthetic network data
#                             real_data <- read.csv("real_network_data.csv")
#                             synthetic_data <- read.csv("synthetic_network_data.csv")
#                             
#                             # Combine the data into one data frame
#                             data <- rbind(real_data, synthetic_data)
#                             
#                             # Create a target variable with 0 for real data and 1 for synthetic data
#                             target <- c(rep(0, nrow(real_data)), rep(1, nrow(synthetic_data)))
#                             
#                             # Split the data into training and testing sets
#                             set.seed(123)
#                             indices <- createDataPartition(target, p = 0.7, list = FALSE)
#                             training_data <- data[indices, ]
#                             testing_data <- data[-indices, ]
#                             training_target <- target[indices]
#                             testing_target <- target[-indices]
#                             
#                             # Train a Random Forest classifier on the training data
#                             model <- randomForest(training_data, training_target, ntree = 500)
#                             
#                             # Predict the class of the testing data using the trained model
#                             predictions <- predict(model, testing_data)
#                             
#                             # Evaluate the model using confusion matrix
#                             confusion_matrix <- table(testing_target, predictions)
#                             print(confusion_matrix)
#                             
#                             # Calculate the accuracy of the model
#                             accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
#                             print(accuracy)
#                             
#                             #####%%%%###################
#                             # Import required libraries
#                             import networkx as nx
#                             import numpy as np
#                             import pandas as pd
#                             from sklearn.model_selection import train_test_split
#                             from sklearn.preprocessing import StandardScaler
#                             from sklearn.neural_network import MLPClassifier
#                             from sklearn.metrics import accuracy_score
#                             
#                             # Load the real network data into a NetworkX graph object
#                             G = nx.read_edgelist('real_network.txt', create_using=nx.Graph(), nodetype=int)
#                             
#                             # Extract graph features from the network
#                             degree_sequence = [d for n, d in G.degree()]
#                             clustering_coefficients = nx.clustering(G).values()
#                             
#                             # Combine the extracted features into a feature matrix
#                             features = np.column_stack([degree_sequence, clustering_coefficients])
#                             
#                             # Load the synthetic network data and extract graph features
#                             synthetic_networks = []
#                             for i in range(1, 11):
#                               G = nx.read_edgelist('synthetic_network_' + str(i) + '.txt', create_using=nx.Graph(), nodetype=int)
#                             degree_sequence = [d for n, d in G.degree()]
#                             clustering_coefficients = nx.clustering(G).values()
#                             synthetic_networks.append(np.column_stack([degree_sequence, clustering_coefficients]))
#                             
#                             # Combine the extracted features of all synthetic networks into a feature matrix
#                             synthetic_features = np.concatenate(synthetic_networks)
#                             
#                             # Combine the real network features and synthetic network features into a single feature matrix
#                             all_features = np.concatenate([features, synthetic_features])
#                             
#                             # Scale the features using StandardScaler
#                             scaler = StandardScaler()
#                             scaled_features = scaler.fit_transform(all_features)
#                             
#                             # Create a target vector for the real and synthetic networks
#                             real_network_target = np.zeros(features.shape[0])
#                             synthetic_network_target = np.ones(synthetic_features.shape[0])
#                             all_targets = np.concatenate([real_network_target, synthetic_network_target])
#                             
#                             # Split the features and targets into training and testing sets
#                             X_train, X_test, y_train, y_test = train_test_split(scaled_features, all_targets, test_size=0.2)
#                             
#                             # Train a Multi-Layer Perceptron classifier on the training data
#                             mlp = MLPClassifier(hidden_layer_sizes=(100,100,100), max_iter=500, alpha=0.0001,
#                                                 solver='sgd', verbose=10,  random_state=21, tol=0.000000001)
#                             mlp.fit(X_train, y_train)
#                             
#                             # Predict the classes of the test data using the trained classifier
#                             y_pred = mlp.predict(X_test)
#                             
#                             # Calculate the accuracy of the classifier
#                             accuracy = accuracy_score(y_test, y_pred)
#                             
#                             ########%%%%%%#################
#                             # Load necessary libraries
#                             library(caret)
#                             library(ggplot2)
#                             
#                             # Load real network data
#                             data("karate")
#                             real_network <- as.matrix(get.adjacency(karate))
#                             
#                             # Generate synthetic networks using Barabasi-Albert model
#                             library(igraph)
#                             set.seed(123)
#                             synthetic_network <- generate_synthetic(real_network, model = "barabasi.game")
#                             
#                             # Extract graph features
#                             real_features <- graph_features(real_network)
#                             synthetic_features <- graph_features(synthetic_network)
#                             
#                             # Combine real and synthetic features into a single data frame
#                             all_features <- rbind(real_features, synthetic_features)
#                             
#                             # Create target variable indicating real or synthetic network
#                             real_target <- rep(1, nrow(real_features))
#                             synthetic_target <- rep(0, nrow(synthetic_features))
#                             all_target <- c(real_target, synthetic_target)
#                             
#                             # Split data into training and testing sets
#                             set.seed(456)
#                             split <- createDataPartition(all_target, p = 0.7, list = FALSE)
#                             train_features <- all_features[split, ]
#                             train_target <- all_target[split]
#                             test_features <- all_features[-split, ]
#                             test_target <- all_target[-split]
#                             
#                             # Train a random forest classifier
#                             model <- train(x = train_features, y = train_target, method = "rf")
#                             
#                             # Predict on the test set
#                             prediction <- predict(model, newdata = test_features)
#                             
#                             # Calculate accuracy
#                             accuracy <- mean(prediction == test_target)
#                             print(paste("Accuracy:", accuracy))
#                             
#                             ###################################
#                             # Load required packages
#                             library(data.table)
#                             library(caret)
#                             library(randomForest)
#                             
#                             # Generate Example Graphs
#                             set.seed(123)
#                             
#                             # Erdos-Renyi
#                             erdos_renyi_graph <- erdos.renyi.game(100, 0.05, type = "gnp", directed = FALSE)
#                             erdos_renyi_adj <- get.adjacency(erdos_renyi_graph, sparse = FALSE)
#                             
#                             # Small-World
#                             small_world_graph <- watts.strogatz.game(1, 100, 5, 0.05, circular = FALSE)
#                             small_world_adj <- get.adjacency(small_world_graph, sparse = FALSE)
#                             
#                             # Scale-Free
#                             scale_free_graph <- barabasi.game(100, m = 3, directed = FALSE)
#                             scale_free_adj <- get.adjacency(scale_free_graph, sparse = FALSE)
#                             
#                             # Spatial
#                             x <- rnorm(100)
#                             y <- rnorm(100)
#                             d <- as.dist(squareform(dist(cbind(x, y))))
#                             spatial_graph <- graph.adjacency(as.matrix(d) <= 1, mode = "undirected", weighted = NULL)
#                             spatial_adj <- get.adjacency(spatial_graph, sparse = FALSE)
#                             
#                             # Generate Graph Features
#                             graph_features <- data.table(
#                               type = c(rep("Erdos-Renyi", 100), rep("Small-World", 100), rep("Scale-Free", 100), rep("Spatial", 100)),
#                               degree = c(degree(erdos_renyi_graph), degree(small_world_graph), degree(scale_free_graph), degree(spatial_graph)),
#                               closeness = c(closeness(erdos_renyi_graph), closeness(small_world_graph), closeness(scale_free_graph), closeness(spatial_graph)),
#                               betweenness = c(betweenness(erdos_renyi_graph), betweenness(small_world_graph), betweenness(scale_free_graph), betweenness(spatial_graph)),
#                               eigen_centrality = c(evcent(erdos_renyi_graph)$vector, evcent(small_world_graph)$vector, evcent(scale_free_graph)$vector, evcent(spatial_graph)$vector),
#                               clustering_coefficient = c(transitivity(erdos_renyi_graph), transitivity(small_world_graph), transitivity(scale_free_graph), transitivity(spatial_graph)),
#                               graph_density = c(graph.density(erdos_renyi_adj), graph.density(small_world_adj), graph.density(scale_free_adj), graph.density(spatial_adj))
#                             )
#                             
#                             # Split Data into Train and Test Sets
#                             train_index <- createDataPartition(graph_features$type, p = 0.7, list = FALSE)
#                             train_data <- graph_features[train_index,]
#                             test_data
#                             
#                             ###############%%%%%%###############
#                             # Load necessary libraries
#                             library(dplyr)
#                             library(caret)
#                             library(e1071)
#                             
#                             # Generate sample data
#                             set.seed(123)
#                             erdos_renyi <- erdos.renyi.game(1000, 0.1)
#                             small_world <- watts.strogatz.game(1, 1000, 6, 0.2)
#                             scale_free <- barabasi.game(1000, m = 5)
#                             spatial_network <- randomGeometricGraph(1000, 0.1)
#                             
#                             # Extract graph features
#                             erdos_renyi_features <- as.data.frame(t(get.adjacency(erdos_renyi, sparse = F)))
#                             small_world_features <- as.data.frame(t(get.adjacency(small_world, sparse = F)))
#                             scale_free_features <- as.data.frame(t(get.adjacency(scale_free, sparse = F)))
#                             spatial_network_features <- as.data.frame(t(get.adjacency(spatial_network, sparse = F)))
#                             
#                             # Create response variable
#                             erdos_renyi_response <- rep("Erdos-Renyi", 1000)
#                             small_world_response <- rep("Small-World", 1000)
#                             scale_free_response <- rep("Scale-Free", 1000)
#                             spatial_network_response <- rep("Spatial", 1000)
#                             
#                             # Combine features and response into one dataset
#                             data <- rbind(cbind(erdos_renyi_features, response = erdos_renyi_response),
#                                           cbind(small_world_features, response = small_world_response),
#                                           cbind(scale_free_features, response = scale_free_response),
#                                           cbind(spatial_network_features, response = spatial_network_response))
#                             
#                             # Split data into training and testing sets
#                             set.seed(456)
#                             train_index <- createDataPartition(data$response, p = 0.7, list = F)
#                             train_data <- data[train_index, ]
#                             test_data <- data[-train_index, ]
#                             
#                             # Train a SVM model
#                             svm_model <- svm(response ~ ., data = train_data, kernel = "linear")
#                             
#                             # Make predictions on the test data
#                             predictions <- predict(svm_model, test_data[, -ncol(test_data)])
#                             
#                             # Evaluate the model's performance
#                             confusionMatrix(predictions, test_data$response)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sample code 13
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# library(igraph)
# 
# # Function to generate a graph with similar statistical and structural properties as the empirical graph
# generate_graph <- function(emp_graph) {
#   # Generate a random graph with the same number of nodes as the empirical graph
#   n <- vcount(emp_graph)
#   rand_graph <- erdos.renyi.game(n, 0.5, directed = FALSE)
#   
#   # Generate the degree sequence for the random graph
#   degree_seq <- degree(emp_graph)
#   
#   # Use the degree sequence to generate a configuration model graph
#   rand_graph_cm <- sample_degseq(degree_seq, method = "simple.no.multiple")
#   
#   # Apply a rewiring algorithm to the configuration model graph to generate a graph with similar structural properties
#   rand_graph_rewire <- rewiring_game(rand_graph_cm, 100)
#   
#   # Apply a community detection algorithm to the generated graph and compare the community structure with that of the empirical graph
#   rand_graph_comm <- cluster_louvain(rand_graph_rewire)
#   emp_graph_comm <- cluster_louvain(emp_graph)
#   
#   # If the community structure is different, repeat the process
#   while (length(unique(rand_graph_comm$membership)) != length(unique(emp_graph_comm$membership))) {
#     rand_graph_cm <- sample_degseq(degree_seq, method = "simple.no.multiple")
#     rand_graph_rewire <- rewiring_game(rand_graph_cm, 100)
#     rand_graph_comm <- cluster_louvain(rand_graph_rewire)
#   }
#   
#   return(rand_graph_rewire)
# }
# 
# # Generate a set of 10 graphs with similar statistical and structural properties as the empirical graph
# emp_graph <- barabasi.game(100, m = 3)
# set.seed(123)
# gen_graphs <- replicate(10, generate_graph(emp_graph))
# 
# # Plot the empirical graph and the generated graphs
# par(mfrow = c(3, 4))
# plot(emp_graph, main = "Empirical Graph")
# for (i in 1:10) {
#   plot(gen_graphs[[i]], main = paste0("Generated Graph ", i))
# }
# 
# ##############%%%%%####################
# What are the properties that a best generative model should have in order to be considered as a generative model for empirical graphs. Code in R built from scratch to implement a single generative model that has all these properties without using stochastic block matrix or igraph
# 
# generate_graph <- function(n, p) {
#   A <- matrix(0, nrow=n, ncol=n)
#   for (i in 1:(n-1)) {
#     for (j in (i+1):n) {
#       if (runif(1) < p) {
#         A[i,j] <- 1
#         A[j,i] <- 1
#       }
#     }
#   }
#   return(A)
# }
# 
# ####%%%%%%%%########################
# generate_graph <- function(n, p, q, alpha, beta) {
#   # n: the number of vertices
#   # p: the probability of an edge between vertices in the same community
#   # q: the probability of an edge between vertices in different communities
#   # alpha: a parameter that controls the clustering coefficient
#   # beta: a parameter that controls the degree distribution
#   
#   # Generate the community structure
#   com <- sample(1:3, n, replace=TRUE, prob=c(0.4, 0.3, 0.3))
#   
#   # Generate the adjacency matrix
#   A <- matrix(0, n, n)
#   for (i in 1:(n-1)) {
#     for (j in (i+1):n) {
#       if (com[i] == com[j]) {
#         A[i,j] <- rbinom(1, 1, p)
#         A[j,i] <- A[i,j]
#       } else {
#         A[i,j] <- rbinom(1, 1, q)
#         A[j,i] <- A[i,j]
#       }
#     }
#   }
#   
#   # Apply the configuration model to match the degree distribution
#   d <- degree(A)
#   D <- matrix(0, n, n)
#   for (i in 1:n) {
#     for (j in (i+1):n) {
#       D[i,j] <- rbinom(1, d[i]*d[j]/sum(d), beta)
#       D[j,i] <- D[i,j]
#     }
#   }
#   A <- A*D
#   
#   # Apply the rewiring algorithm to increase the clustering coefficient
#   for (i in 1:n) {
#     neighbors <- which(A[i,] == 1)
#     for (j in 1:(length(neighbors)-1)) {
#       for (k in (j+1):length(neighbors)) {
#         if (rbinom(1, 1, alpha) == 1) {
#           A[neighbors[j], neighbors[k]] <- 1
#           A[neighbors[k], neighbors[j]] <- 1
#         }
#       }
#     }
#   }
#   
#   # Return the adjacency matrix
#   return(A)
# }
# 
# ##########%%%%%%%%##################
# # Generate an EP graph with n nodes and degree k
# generate_ep_graph <- function(n, k) {
#   # Initialize the adjacency matrix
#   A <- matrix(0, nrow = n, ncol = n)
#   # Set the first k nodes to be fully connected
#   A[1:k, 1:k] <- 1
#   # Initialize the list of frontier nodes
#   frontier <- 1:k
#   # Iterate until all nodes have degree k
#   while (length(frontier) > 0) {
#     # Choose a random frontier node
#     i <- sample(frontier, 1)
#     # Find its neighbors
#     neighbors <- which(A[i,] == 1)
#     # Choose a random neighbor
#     j <- sample(neighbors, 1)
#     # Find the neighbors of j that are not already connected to i
#     available_neighbors <- setdiff(which(A[j,] == 1), neighbors)
#     # If there are no available neighbors, remove j from the frontier
#     if (length(available_neighbors) == 0) {
#       frontier <- setdiff(frontier, j)
#     } else {
#       # Choose a random available neighbor
#       l <- sample(available_neighbors, 1)
#       # Connect i to l
#       A[i,l] <- 1
#       A[l,i] <- 1
#       # Add l to the frontier if it now has degree k
#       if (sum(A[l,]) == k) {
#         frontier <- c(frontier, l)
#       }
#     }
#   }
#   # Return the adjacency matrix
#   return(A)
# }
# 
# ########%%%%%%###################%
# library(Matrix)
# 
# # Initialize the adjacency matrix
# n <- 100
# A <- Matrix(0, n, n)
# 
# # Set the expansion parameter
# alpha <- 0.8
# 
# # Generate the graph using the EP algorithm
# for (i in 2:n) {
#   # Choose a random set of neighbors for the new node
#   neighbor_set <- sample(1:(i-1), floor(alpha*(i-1)))
#   
#   # Add the new node and its neighbors to the graph
#   A[i, neighbor_set] <- 1
#   A[neighbor_set, i] <- 1
# }
# 
# # Print the adjacency matrix
# A
# 
# 
# #########%%%%%%%%%#################
# # Function to generate expander propagation graph
# expander_propagation_graph <- function(n, p, q) {
#   # Initialize graph
#   g <- matrix(0, n, n)
#   # Create expander subgraph
#   k <- round(log(n))
#   expandersubgraph <- matrix(0, k, k)
#   for (i in 1:k) {
#     for (j in (i+1):k) {
#       if (runif(1) < p) {
#         expandersubgraph[i,j] <- 1
#         expandersubgraph[j,i] <- 1
#       }
#     }
#   }
#   # Connect expander subgraph to main graph
#   for (i in (k+1):n) {
#     for (j in 1:k) {
#       if (runif(1) < q) {
#         g[i,j] <- 1
#         g[j,i] <- 1
#       }
#     }
#   }
#   # Propagate edges from expander subgraph to main graph
#   for (i in (k+1):n) {
#     candidates <- which(g[i,] == 0 & rowSums(expandersubgraph[g[1:k,] == 1,] == 1) >= 2)
#     if (length(candidates) > 0) {
#       j <- sample(candidates, 1)
#       g[i,j] <- 1
#       g[j,i] <- 1
#     }
#   }
#   return(g)
# }
# 
# # Function to compute degree distribution of graph
# degree_distribution <- function(g) {
#   degrees <- rowSums(g)
#   hist(degrees, breaks = seq(0, max(degrees)+1, by = 1), freq = FALSE, main = "Degree distribution", xlab = "Degree")
# }
# 
# # Generate graph with 1000 nodes, p=0.5 and q=0.1
# set.seed(123)
# g <- expander_propagation_graph(1000, 0.5, 0.1)
# 
# # Compute degree distribution of graph
# degree_distribution(g)
# 

