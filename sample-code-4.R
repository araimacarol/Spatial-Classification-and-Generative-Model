# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Code in R to grow a spatial scale-free expander graph with the following properties. This code should be built from scratch and not contain any igraph packages etc.
# 
# 1) Nodes should be generated on a unit square where the initial set of nodes are distributed spatially with an underlying spatial structure, e.g. by not permitting nodes to land too near each other, or not too far from each other using an in built 'quadtree' algorithm built from scratch
# 
# 2) Use the distance matrix to create an adjacency matrix for the graph.o create an adjacency matrix based on the distance matrix, you can define a cutoff distance 'r' and set the adjacency matrix element A_{ij} to 1 if dist_{ij} <= 'r' and 0 otherwise . This parameter 'r' controls the strength of the spatial distance effect 
# 
# 3) This network should be grown  and for each of these nodes added to the network, add edges based on a 'probability of attachment to other nodes' favouring (i)short spatial distances between the new node (s) and existing nodes  * (ii)  community structure denoted by probability 'mu' that allows you to tune the number of communities in the graph and the degree of connectivity within and between communities. *(iii) a scale-free preferential attachment with power law parameter 'beta' and the parameter 'm' that controls the number of edges at each instance that define the attachment types that we can tune to get 'linear preferential attachment' (one edge added each instance) and any other form of attachment (multiple nodes). Include 'alpha': the parameter that controls the strength of the degree effect 
# 
# 4)Add a rewiring probability 'prewire' for small world effect to the resulting graph to enhance the connectivity to long range nodes. Thus  you can randomly rewire each edge with probability 'prewire' to a random node within a certain distance
# 
# 5) Add arguments to create edge weights, node attributes etc
# 
# 6)Check the expansion properties of the graph by computing the Laplacian matrix. The Laplacian matrix represents the graph's connectivity and its ability to spread information. A good spatial graph should have a small eigenvalue gap, indicating strong expansion properties. 
# 7)This graph should be a single function where we can tune, 'beta', N, 'alpha', 'mu' , 'prewire', 'm', to generate different graph models+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+Code in R to grow a spatial scale-free expander graph with the following properties. This code should be built from scratch and not contain any igraph packages etc.

# 1) Nodes should be generated on a unit square where the initial set of nodes are distributed spatially with an underlying spatial structure, e.g. by not permitting nodes to land too near each other, or not too far from each other
# 
# 2) Use the distance matrix to create an adjacency matrix for the graph. Cr eate an adjacency matrix based on the distance matrix, you can define a cutoff distance 'r' and set the adjacency matrix element A_{ij} to 1 if dist_{ij} <= 'r' and 0 otherwise. This parameter 'r' controls the strength of the spatial distance effect. 
# 
# 3) This network should be grown from scratch where each node added at an instance to the network be categorize into non overlapping communities,  where nodes within the same community and between different communities connect with each other with an attachment probability matrix 'P_ij' favoring   (i)short spatial distances (constrained by 'r') and a scale-free degree preferential attachment with power law parameter 'beta'.  Let 'alpha': the parameter that controls the strength of the degree effect.  The attachment probability matrix 'P_ij' should be symmetric if we want erdos renyi random graph but otherwise non-symmetric.
# 
# 4)Add a rewiring probability 'prewire' for small world effect to the resulting graph to enhance the connectivity to long range nodes. Thus, you can randomly rewire each edge with probability 'prewire' to a random node within a certain distance.
# 
# 5)This graph should be a single function where we can tune, 'beta', N, 'alpha', 'mu' , 'prewire', 'm', to generate different graph models
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Code in R to grow a spatial scale-free graph with the following properties. This code should be built from scratch and not contain any igraph packages etc.
# 
# 1) Nodes should be generated on a unit square where the initial set of nodes are distributed spatially with an underlying spatial structure.
# 
# 2) Initialize the network with a single node and place it in an initial community /block
# 
# 3) For each new node now added to the network, do the following:
#   
#   i) Step 1: Create an arbitrary number of communities/blocks to partition or divide the network into
# 
# ii) Step 2: After creating the arbitrary number of communities, define a community (or block) probability attachment matrix, which specifies the probability of an edge between any two nodes in the same community/block and nodes in different communities. This block/community probability matrix should have size equal to the number of communities/blocks and must favour (i) short 'spatial distances' constrained by parameter 'r' (the parameter that controls the degree of sparsity) and (ii) a scale-free degree preferential attachment with power law parameter 'beta' and m the parameter that controls the number of edges added at each instance. This probability attachment matrix should depend on all of the degree, distance and community affiliation
# 
# iii) Step 3: After defining the probability matrix, generate community/ block assignment for each new node. Thus, assign each of the new node added to the network to a community or block using a 2D poison point process
# 
# iv) Step 4: Add edges between the new node and existing nodes based on the community/block assignments and community/block probability matrix.
# 
# 4)Add a rewiring probability 'prewire' for small world effect to the resulting graph to enhance the connectivity to long range nodes. Thus, you can randomly rewire each edge with probability 'prewire' to a random node within a certain distance.
# 
# 4)This graph should be a single function where we can tune, 'beta', N,  'prewire', 'm', to generate different graph models
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Code in R to grow a spatial scale-free graph with the following properties. This code should be built from scratch and not contain any igraph packages etc.
# 
# 1) Nodes should be generated on a unit square where the initial set of nodes are distributed spatially with an underlying spatial structure.
# 
# 2) Initialize the network with a single node and place it in an initial community /block
# 
# 3) For each new node now added to the network, do the following:
#   
#   i) Step 1: Create an arbitrary number of communities/blocks to partition or divide the network into using the parameter "numofComm"
# 
# ii) Step 2: After creating the arbitrary number of communities, define a community (or block) probability attachment matrix denoted as 'P_ij' , which specifies the probability of an edge between any two nodes in the same community/block and nodes in different communities. The block/community probability matrix should have size equal to the number of communities/blocks.
# 
# iii) Step 3: We assume that the probability of attachment matrix for connections between nodes within the same communities for the different communities/blocks is the same. Therefore, the within community attachment probability matrix should depend solely on the following: (A) short 'spatial distances' between nodes constrained by parameter 'r' (the parameter that controls the degree of sparsity) and (B) a scale-free degree preferential attachment with power law parameter 'beta' and m the parameter that controls the number of edges added at each instance.
# 
# iv) Step 4: The between community probability of attachment is different from within community attachment. Therefore, the between probability of attachment should depend on the following: (C) rewiring probability 'prewire' constrained by 'r', (D) short 'spatial distances' between nodes constrained by parameter 'r' (the parameter that controls the degree of sparsity) and (B) a scale-free degree preferential attachment with power law parameter 'beta' and m the parameter that controls the number of edges added at each instance. 
# 
# v) Step 5: After defining the probability matrix, generate community/ block assignment for each new node. Thus, assign each of the new node added to the network to a community or block using a 2D poison point process.
# 
# vi) Step 6: Add edges between the new node and existing nodes based on the community/block assignments and community/block probability matrix for within and between communities.
# 
# 4)Add a rewiring probability 'prewire' for small world effect to the entire resulting graph to enhance the connectivity to long range nodes. Thus, you can randomly rewire each edge with probability 'prewire' to a random node within a certain distance.
# 
# 5)This graph should be a single function where we can tune, 'beta', N,  "numofComm", 'prewire', 'm', to generate different graph models
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Code in R to grow a spatial scale-free graph with the following properties. This code should be built from scratch and not contain any igraph packages etc.
# 
# 1) Nodes should be generated on a unit square where the initial set of nodes are distributed spatially with an underlying spatial structure. Define a distance matrix between nodes as 'dist_mat'
# 
# 2) Initialize the network with a single node. eg adj_mat <- matrix(0, nrow = N, ncol = N), therefore, adj_mat[1, 1] <- 1 and place it in an initial community /block.
# 
# 3) For each new node now added to the network, do the following:
#   
#   i) Step 1: Create an arbitrary number of communities/blocks to partition or divide the network into using the parameter "numofComm"
# 
# ii) Step 2: After creating the arbitrary number of communities, define a community (or block) probability attachment matrix denoted as 'P_ij' , which specifies the probability of an edge between any two nodes in the same community/block and nodes in different communities. The block/community probability matrix should have size equal to the number of communities/blocks.
# 
# iii) Step 3: We assume that the probability of attachment matrix for connections between nodes within the same communities for the different communities/blocks is the same. Therefore, the within community attachment probability matrix should depend solely on the following: (A) short 'spatial distances' between nodes constrained by parameter 'r' (the parameter that controls the degree of sparsity) and (B) a scale-free degree preferential attachment with power law parameter 'beta' and m the parameter that controls the number of edges added at each instance. The distance probability for within community attachment should be defined as spat_distprobs=1/(1+(dist_mat[1:(i-1), i])^beta) if and only if dist_mat[1:(i-1), i] <= r. The degree probability for within community attachment should be define as
# deg_probs=(degrees^alpha)/(sum(degrees^alpha)), where degrees=rowSums(adj_mat[1:(i-1), ])
# 
# iv) Step 4: The between community probability of attachment is different from within community attachment. Therefore, the between probability of attachment should depend on the following: (C) rewiring probability 'prewire' constrained by 'r', (D) short 'spatial distances' between nodes constrained by parameter 'r' (the parameter that controls the degree of sparsity) and (B) a scale-free degree preferential attachment with power law parameter 'beta' and m the parameter that controls the number of edges added at each instance. Similarly, the distance probability for between community attachment should be defined as spat_distprobs=1/(1+(dist_mat[1:(i-1), i])^beta) if and only if dist_mat[1:(i-1), i] <= r. The degree probability for between community attachment should be define as
# deg_probs=(degrees^alpha)/(sum(degrees^alpha)), where degrees=rowSums(adj_mat[1:(i-1), ])
# 
# v) Step 5: After defining the probability matrix, generate community/ block assignment for each new node. Thus, assign each of the new node added to the network to a community or block using a 2D poison point process.
# 
# vi) Step 6: Add edges between the new node and existing nodes based on the community/block assignments and community/block probability matrix for within and between communities.
# 
# 4)Add a rewiring probability 'prewire' for small world effect to the entire resulting graph to enhance the connectivity to long range nodes. Thus, you can randomly rewire each edge with probability 'prewire' to a random node within a certain distance.
# 
# 5)This graph should be a single function where we can tune, 'beta', N,  "numofComm", 'prewire', 'm', to generate different graph models
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to generate a spatial scale-free expander graph
# Arguments:
#   N: number of nodes
#   lambda: intensity parameter of Poisson point process
#   r: cutoff distance for spatial adjacency
#   mu: probability of attachment to nodes within community
#   beta: power law parameter for distance-dependent probability
#   m: number of edges to attach to new node at each step
#   alpha: strength of degree effect in attachment probability
#   sigma: decay of degree effect with distance
#   L: number of rewiring steps for small world effect
#   p: rewiring probability
#   node_attrs: list of node attributes
#   edge_weights: matrix of edge weights
# Output:
#   A list containing the adjacency matrix, node attributes, and edge weights
spatial_scale_free_expander <- function(N, lambda, r, mu, beta, m, alpha, sigma, L, p,
                                        node_attrs = NULL, edge_weights = NULL) {
  # Generate nodes using 2D Poisson point process on d-dimensional torus
  d <- 2 # dimension
  grid_size <- sqrt(N/lambda) # size of torus
  points <- matrix(runif(N*d), ncol = d) * grid_size
  torus_dist <- function(x, y) {
    dx <- abs(x - y)
    dx <- dx - grid_size * round(dx/grid_size)
    sqrt(sum(dx^2))
  }
  # Initialize adjacency matrix and degree vector
  A <- matrix(0, N, N)
  degree <- rep(0, N)
  # Compute distance matrix and minimum distance between nodes
  dist_mat <- as.matrix(dist(points, diag = TRUE, upper = TRUE))
  min_dist <- min(dist_mat[dist_mat > 0])
  # Add nodes one by one
  for (i in 1:N) {
    # Compute attachment probabilities for existing nodes
    if (i <= m) {
      p_attach <- rep(1/N, N) # first m nodes attach to all nodes uniformly
    } else {
      dist_i <- dist_mat[i, ]
      degree_i <- degree
      degree_i[i] <- 0 # exclude self from degree calculation
      # Compute distance-dependent probability function
      p_dist <- 1/(dist_i^beta + 1e-16)
      # Compute degree-dependent probability function
      p_degree <- (degree_i^alpha + 1e-16)/(sum(degree_i^alpha) + N*1e-16)
      # Combine distance and degree probabilities with spatial community structure
      p_attach <- mu * p_dist * p_degree
      p_attach <- p_attach + (1-mu) * rep(1/N, N)
    }
    # Sample m nodes to attach to with probability proportional to p_attach
    attach_nodes <- sample.int(N, m, replace = TRUE, prob = p_attach)
    # Add edges between new node and selected nodes
    for (j in attach_nodes) {
      if (i != j & A[i, j] == 0) { # check for self-loops and existing edges
        # Compute edge weight as inverse of distance
        edge_weight <- 1/(torus_dist(points[i, ], points[j, ])^sigma + 1e-16)
        A[i, j] <- edge_weight
        A[j, i] <- edge_weight
        degree[i] <- degree[i] + 1
        degree[j] <- degree[j] + 1

        
                
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Function to generate a spatial scale-free expander graph
        spatial_sf_expander_graph <- function(N, lambda, r, beta, alpha, mu, L, m, p, dim = 2, weighted = FALSE, node_atts = FALSE) {
          
          # Generate points using Poisson point process
          points <- matrix(rpois(N*dim, lambda), ncol = dim)
          
          # Define quadtree algorithm to distribute points spatially
          quadtree <- function(pts, r, n) {
            if(n <= 1) {
              return(pts)
            } else {
              # Divide points into quadrants
              xmid <- median(pts[, 1])
              ymid <- median(pts[, 2])
              quad1 <- pts[pts[, 1] <= xmid & pts[, 2] <= ymid, ]
              quad2 <- pts[pts[, 1] <= xmid & pts[, 2] > ymid, ]
              quad3 <- pts[pts[, 1] > xmid & pts[, 2] <= ymid, ]
              quad4 <- pts[pts[, 1] > xmid & pts[, 2] > ymid, ]
              # Recursively apply quadtree to each quadrant
              q1 <- quadtree(quad1, r/2, n/4)
              q2 <- quadtree(quad2, r/2, n/4)
              q3 <- quadtree(quad3, r/2, n/4)
              q4 <- quadtree(quad4, r/2, n/4)
              # Combine points from all quadrants
              return(rbind(q1, q2, q3, q4))
            }
          }
          
          # Use quadtree algorithm to distribute points
          points <- quadtree(points, r, N)
          
          # Calculate distance matrix between points
          dist_mat <- as.matrix(dist(points))
          
          # Create adjacency matrix based on distance matrix and cutoff distance r
          adj_mat <- ifelse(dist_mat <= r, 1, 0)
          diag(adj_mat) <- 0
          
          # Add nodes and edges based on preferential attachment
          for (i in 1:N) {
            degrees <- colSums(adj_mat)
            distances <- dist_mat[i, ]
            prob <- (degrees^alpha) * (1/(distances + 1e-6)^beta)
            prob[i] <- 0
            prob <- prob/sum(prob)
            new_edges <- sample(1:N, m, prob = prob)
            adj_mat[i, new_edges] <- 1
            adj_mat[new_edges, i] <- 1
          }
          
          # Add community structure
          if(mu > 0 & mu < 1) {
            # Divide nodes into communities
            num_communities <- round(mu*N)
            communities <- rep(1:num_communities, each = ceiling(N/num_communities))[1:N]
            # Add edges within communities
            for(c in 1:num_communities) {
              comm_nodes <- which(communities == c)
              comm_degrees <- colSums(adj_mat[comm_nodes, ])
              for(n in comm_nodes) {
                comm_prob <- comm_degrees/n
                comm_prob[n] <- 0
                comm_prob <- comm_prob/sum(comm_prob)
                new_edges <- sample(comm_nodes, m, prob = comm_prob)
                adj_mat[n, new_edges] <- 1
                adj_mat[new_edges, n] <- 1
                
                
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to grow a spatial scale-free expander graph
spatial_scale_free_expander <- function(N, lambda, beta, sigma, alpha, mu, L, p, m) {
  
  # Step 1: Generate points on a torus using 2D Poisson point process
  coords <- matrix(runif(2*N), ncol = 2)
  coords <- coords + rep(c(0,1), each = N/2) # periodic boundary condition
  
  # Step 2: Calculate distance matrix and minimum distance
  dist_mat <- as.matrix(dist(coords, diag = TRUE, upper = TRUE))
  min_dist <- min(dist_mat[dist_mat > 0])
  
  # Step 3: Create adjacency matrix based on distance matrix
  cutoff <- L * min_dist
  adj_mat <- ifelse(dist_mat <= cutoff, 1, 0)
  diag(adj_mat) <- 0
  
  # Step 4: Add nodes and edges based on preferential attachment and spatial constraints
  for (i in 1:(N - 1)) {
    prob_attach <- rep(0, N)
    dist_vec <- dist_mat[i, ]
    deg_vec <- colSums(adj_mat)
    for (j in 1:(i-1)) {
      prob_attach[j] <- (dist_vec[j]^(-beta)) * (deg_vec[j]^alpha)
    }
    for (j in (i+1):N) {
      prob_attach[j] <- (dist_vec[j]^(-beta)) * (deg_vec[j]^alpha)
    }
    # Apply community structure
    if (runif(1) < mu) {
      comm_size <- rpois(1, N/2)
      comm_nodes <- sample(1:N, comm_size)
      prob_attach[-comm_nodes] <- 0
    }
    # Apply preferential attachment
    if (m == 1) {
      new_edge <- sample(1:(i-1), 1, prob = prob_attach[1:(i-1)]/sum(prob_attach[1:(i-1)]))
      adj_mat[i, new_edge] <- 1
      adj_mat[new_edge, i] <- 1
    } else {
      for (j in 1:m) {
        new_edge <- sample(1:(i-1), 1, prob = prob_attach[1:(i-1)]/sum(prob_attach[1:(i-1)]))
        adj_mat[i, new_edge] <- 1
        adj_mat[new_edge, i] <- 1
      }
    }
  }
  
  # Step 5: Randomly rewire edges for small world effect
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (adj_mat[i,j] == 1) {
        if (runif(1) < p) {
          rew_dist <- runif(1) * cutoff
          rew_nodes <- which(dist_mat[i,] >= rew_dist & adj_mat[i,] == 0)
          if (length(rew_nodes) > 0) {
            rew_node <- sample(rew_nodes, 1)
            adj_mat[i,j] <- 0
            adj_mat[j,i] <- 0
            adj_mat[i,rew_node] <- 1
            adj_mat[rew_node,i] <- 1
          }
        }
      }
    }
  }
  
  # Step 6: Add edge weights and node attributes if desired
  
  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+samplecode
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
generate_spatial_scale_free_expander_graph <- function(N, lambda, L, beta, sigma, alpha, mu, m, p) {
  # Generate the initial set of nodes
  points <- runif(N*2, min=0, max=L)
  x <- points[1:N]
  y <- points[(N+1):(2*N)]
  
  # Generate a distance matrix between all pairs of nodes
  dist_mat <- as.matrix(dist(cbind(x, y), method="euclidean"))
  
  # Create an adjacency matrix based on the distance matrix
  r <- sigma * log(N) / sqrt(N)
  A <- as.matrix(dist_mat <= r)
  diag(A) <- 0
  
  # Initialize the degree sequence and the list of edges
  degrees <- colSums(A)
  edges <- matrix(0, nrow=0, ncol=2)
  #exp(-sigma*dist_mat[i,])
  # Add nodes and edges using the Barabasi-Albert model
  for (i in seq(N-m)) {
    # Calculate the attachment probability for each existing node
    dist_probs <- 1 / (dist_mat[i,]^beta)
    deg_probs <- (degrees[1:(i-1)] / sum(degrees[1:(i-1)]))^alpha
    attach_probs <- (1 - mu) * deg_probs + (dist_probs *mu)}
    
    # Choose m nodes to attach to
    attach_nodes <- sample.int(i, size=m, replace=FALSE, prob=attach_probs)
    
    # Add edges to the selected nodes
    edges_to_add <- cbind(rep(i, m), attach_nodes)
    edges <- rbind(edges, edges_to_add)
    degrees[i] <- degrees[i] + m
    degrees[attach_nodes] <- degrees[attach_nodes] + 1
  }
  
  # Rewire edges for small world effect
  for (i in seq(nrow(edges))) {
    if (runif(1) < p) {
      dists_to_other_nodes <- dist_mat[edges[i,1],] + dist_mat[edges[i,2],]
      dists_to_other_nodes[edges[i,1]] <- Inf
      dists_to_other_nodes[edges[i,2]] <- Inf
      nodes_to_rewire <- which(dists_to_other_nodes < r)
      if (length(nodes_to_rewire) > 0) {
        edges[i,2] <- sample(nodes_to_rewire, size=1)
      }
    }
  }
  
  # Convert the edge list to an adjacency matrix
  A <- matrix(0, nrow=N, ncol=N)
  for (i in seq(nrow(edges))) {
    A[edges[i,1], edges[i,2]] <- 1
    A[edges[i,2], edges[i,1]] <- 1
  }
  
  # Add edge weights and node attributes
  # ...
  
  # Compute the Laplacian matrix and check the expansion properties
  D <- diag(colSums(A))
  L <- D - A
  eigenvalues <- eigen(L)$values
  eigen_gap <- eigenvalues[2] - eigenvalues[1]
  
  return(list(adj_matrix=A, laplacian_matrix=L, eigen_gap=eigen_gap))
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+samplecode
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
#   sigma: parameter controlling decay of degree effect
#   mu: parameter controlling number of communities
#   L: number of iterations for Barabasi-Albert model
#   p: rewiring probability for small world effect
#   node_attrs: list of node attributes to add to graph
#   edge_weights: logical indicating whether to add edge weights
#
# Returns:
#   A graph object of class igraph

spatial_scale_free_expander <- function(N, lambda, d, r, beta, m, alpha, sigma, mu, L, p, node_attrs = NULL, edge_weights = FALSE) {
  # Generate points using 2D Poisson point process on torus
  points <- matrix(runif(N*d), ncol = d)
  points <- points - floor(points) # Periodic boundary condition
  points <- points + jitter(matrix(0, nrow = N, ncol = d), amount = lambda, factor = 1) # Add spatial structure
  # Calculate distance matrix
  dist_mat <- as.matrix(dist(points, diag = TRUE, upper = TRUE))
  # Create adjacency matrix based on distance matrix
  adj_mat <- ifelse(dist_mat <= r, 1, 0)
  diag(adj_mat) <- 0 # No self-loops
  # Initialize graph with all nodes and edges
  g <- graph.adjacency(adj_mat, mode = "undirected")
  # Add node attributes if provided
  if (!is.null(node_attrs)) {
    for (name in names(node_attrs)) {
      g[[name]] <- node_attrs[[name]]
    }
  }
  # Add edges based on Barabasi-Albert model with distance-dependent probability
  for (i in 1:(N - m)) {
    # Probability of attachment to existing nodes
    prob <- degree(g)^alpha * (1/dist_mat[i,]^beta) * exp(-sigma*dist_mat[i,])
    # Exclude existing neighbors and nodes already connected
    prob[neighbors(g, i) | seq_len(N) <= i] <- 0 
    # Choose m neighbors based on probability of attachment
    for (j in 1:m) {
      nbr <- sample.int(N, size = 1, prob = prob/sum(prob))
      g <- add_edges(g, c(i, nbr))
    }
  }
  # Add edges to remaining nodes with linear preferential attachment
  for (i in (N - m + 1):N) {
    # Probability of attachment to existing nodes
    prob <- degree(g)^alpha * (1/dist_mat[i,]^beta) * exp(-sigma*dist_mat[i,])
    # Exclude existing neighbors and nodes already connected
    prob[neighbors(g, i) | seq_len(N) <= i] <- 0 
    # Choose one neighbor based on probability of attachment
    nbr <- sample.int(N, size = 1, prob = prob/sum(prob))
    g <- add_edges(g, c(i, nbr))
  }
  # Rewire edges with small world probability
  for (i in 1:E(g
                
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+samplecode
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to generate a spatial scale-free expander graph
# with the given parameters
spatial_scale_free_graph <- function(N, lambda, L, beta, sigma, alpha, mu, m, p, weighted = FALSE) {
  
  # Generate points using a Poisson point process on a torus
  points <- matrix(runif(N*2), ncol = 2) * L
  if (weighted) {
    weights <- runif(N)
  }
  
  # Calculate the distance matrix between all pairs of points
  dist_matrix <- as.matrix(dist(points, method = "euclidean", diag = TRUE, upper = TRUE))
  
  # Create an adjacency matrix based on the distance matrix
  r <- log(N)/N
  adj_matrix <- ifelse(dist_matrix <= r, 1, 0)
  
  # Add nodes and edges based on preferential attachment
  for (i in 1:(N-m)) {
    # Calculate the probability of attachment to each existing node
    dist_vec <- dist_matrix[i,]
    deg_vec <- rowSums(adj_matrix)
    deg_vec[i] <- 0
    p_vec <- (dist_vec^(-beta)) * (deg_vec^alpha)
    if (mu < 1) {
      # Add community structure to the probability of attachment
      community_vec <- rep(1, N)
      community_vec[i:(i+floor((1-mu)*N))] <- 2
      community_prob <- matrix(1, nrow = 2, ncol = 2)
      community_prob[1,1] <- (1-mu)^2
      community_prob[2,2] <- mu^2
      community_prob[1,2] <- community_prob[2,1] <- mu*(1-mu)
      community_vec_prob <- community_prob[community_vec, community_vec]
      p_vec <- p_vec * community_vec_prob
    }
    # Choose the nodes to attach to based on the probability of attachment
    m_idx <- sample.int(N-i, m, replace = TRUE, prob = p_vec/(sum(p_vec)))
    new_edges <- cbind(rep(i, m), i+m_idx)
    adj_matrix[new_edges] <- 1
    adj_matrix[new_edges[,2],new_edges[,1]] <- 1
  }
  
  # Add small world effect by rewiring edges
  for (i in 1:N) {
    for (j in 1:N) {
      if (i < j & adj_matrix[i,j] == 1) {
        if (runif(1) < p) {
          # Find a random node within a certain distance
          max_dist <- L/2
          idx <- which(dist_matrix[i,] <= max_dist)
          idx <- idx[idx != i & adj_matrix[i,idx] == 0]
          if (length(idx) > 0) {
            
            

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+samplecode
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+samplecode
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+samplecode
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+samplecode
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+samplecode
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+samplecode
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# Function to grow a spatial scale-free expander graph
spatial_scale_free_expander <- function(N, lambda, r, beta, sigma, alpha, mu, L, p, m, seed = 1) {
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Generate nodes using a 2D Poisson point process on a torus with intensity parameter lambda
  points <- matrix(runif(N*2), ncol=2)
  points <- points + rep(c(0, 1), each = N/2)
  points <- points * matrix(c(L, L/2, 0, sqrt(3)/2), ncol = 2, byrow = TRUE)
  
  # Add spatial structure to the distribution by jittering the points
  jitter_x <- runif(N, min = -r/2, max = r/2)
  jitter_y <- runif(N, min = -r/2, max = r/2)
  points[,1] <- points[,1] + jitter_x
  points[,2] <- points[,2] + jitter_y
  
  # Calculate the distance matrix between all pairs of points and the minimum distance between pairs of nodes
  dist_mat <- as.matrix(dist(points))
  min_dist <- min(dist_mat[dist_mat > 0])
  
  # Create an adjacency matrix based on the distance matrix
  adj_mat <- ifelse(dist_mat <= r, 1, 0)
  diag(adj_mat) <- 0
  
  # Add nodes and edges based on preferential attachment with spatial constraints
  for (i in 1:(N-m)) {
    # Create a probability vector for the new node to attach to existing nodes
    prob_vec <- colSums(adj_mat) / sum(adj_mat)
    dist_vec <- ifelse(dist_mat[i,] > 0, 1 / (dist_mat[i,]^beta), 0)
    prob_vec <- prob_vec * dist_vec
    # Add edges to existing nodes
    for (j in 1:m) {
      # Choose an existing node to attach to based on the probability vector
      attach_node <- sample(1:N, size = 1, prob = prob_vec)
      # Add edge between the new node and the chosen node
      adj_mat[i, attach_node] <- 1
      adj_mat[attach_node, i] <- 1
    }
  }
  
  # Add rewiring probability for small world effect
  for (i in 1:N) {
    for (j in (i+1):N) {
      # Check if there is an edge between nodes i and j
      if (adj_mat[i, j] == 1) {
        # Generate a random number to decide whether to rewire the edge
        if (runif(1) < p) {
          # Find nodes within a certain distance to rewire the edge to
          dist_to_node <- dist_mat[i,] + dist_mat[j,] - 2*dist_mat[i,j]
          dist_to_node[which(dist_to_node <= min_dist)] <- 0
          # Choose a random node to rewire the edge to
          rewire_node <- sample(which(dist_to_node > 0), size = 1)
          # Rewire the edge
          adj_mat[i, j] <- 0
          adj_mat[j, i] <- 0
          adj_mat[i, rewire_node] <- 1
          adj_mat[rewire_node, i] <- 1
        }
        
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+samplecode
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


spatial_scale_free_expander <- function(N, lambda, L, beta, sigma, alpha, m, p, mu = 0, communities = NULL, edge_weights = FALSE, node_attributes = FALSE) {
  
  # Generate points on a torus using a Poisson point process
  points <- matrix(runif(N*L), ncol = L) * c(2*pi, rep(1, L-1))
  points <- matrix(apply(points, 2, function(x) jitter(x, factor = 0.1)), ncol = L)
  points <- apply(points, 2, function(x) x - min(x))
  points <- points / apply(points, 2, max)
  
  # Calculate the distance matrix between all pairs of points
  dist_matrix <- as.matrix(dist(points, diag = TRUE, upper = TRUE))
  min_dist <- min(dist_matrix[dist_matrix > 0])
  
  # Create adjacency matrix based on distance matrix and cutoff distance 'r'
  r <- min_dist / 2
  adj_matrix <- as.matrix(dist_matrix <= r)
  
  # Initialize graph with m nodes and a fully connected subgraph
  graph <- matrix(0, nrow = N, ncol = N)
  for (i in 1:m) {
    for (j in (i+1):m) {
      graph[i,j] <- 1
      graph[j,i] <- 1
    }
  }
  
  # Sequentially add nodes with preferential attachment
  for (i in (m+1):N) {
    # Calculate probability of attachment to existing nodes based on degree and distance
    dist_probs <- 1 / (dist_matrix[i,] ^ beta)
    deg_probs <- (colSums(graph) + alpha) / (sum(colSums(graph)) + alpha*N)
    probs <- dist_probs * deg_probs
    probs[communities != communities[i]] <- probs[communities != communities[i]] * (1-mu)
    probs[communities == communities[i]] <- probs[communities == communities[i]] * m / (sum(probs[communities == communities[i]]) + 1e-9)
    probs[i] <- 0 #
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    generate_graph <- function(beta, N, lambda, sigma, alpha, mu, L, p, m, add_weights = FALSE, add_attributes = FALSE) {
      # Generate N points using a 2D Poisson point process with intensity parameter lambda
      points <- matrix(runif(N*2), ncol = 2) * L
      points <- points + matrix(runif(N*2), ncol = 2) * (L / sqrt(lambda))
      points <- points %% L
      
      # Calculate distance matrix
      dist_matrix <- as.matrix(dist(points))
      
      # Create adjacency matrix based on distance matrix and cutoff distance r
      r <- sigma * mean(dist_matrix)
      adj_matrix <- ifelse(dist_matrix <= r, 1, 0)
      diag(adj_matrix) <- 0
      
      # Grow graph using Barabasi-Albert model with preferential attachment based on degree and distance
      for (i in 2:N) {
        # Calculate distance-dependent probability function
        p_ij <- 1 / (dist_matrix[i,] ^ beta)
        
        # Calculate attachment probabilities based on degree, distance, and community structure (if specified)
        if (mu > 0) {
          comm_size <- floor(N / mu)
          comm_ids <- rep(1:mu, each = comm_size)
          comm_ids <- c(comm_ids, rep(comm_ids[1], N %% mu))
          same_comm <- (comm_ids == comm_ids[i])
          p_attach <- ifelse(same_comm, m, 1 - mu)
          p_attach[which(adj_matrix[i,] == 1)] <- 0
          p_attach <- p_attach * p_ij
        } else {
          p_attach <- adj_matrix[i,] * p_ij
          p_attach[which(adj_matrix[i,] == 1)] <- 0
        }
        
        # Add m edges to existing nodes based on attachment probabilities
        for (j in 1:m) {
          if (sum(p_attach) == 0) {
            break
          }
          j_attach <- sample(1:N, size = 1, prob = p_attach / sum(p_attach))
          adj_matrix[i, j_attach] <- 1
          adj_matrix[j_attach, i] <- 1
          p_attach[j_attach] <- 0
          if (p > 0 && runif(1) < p) {
            k <- sample(which(dist_matrix[i,] > r), size = 1)
            adj_matrix[i, k] <- 1
            adj_matrix[k, i] <- 1
          }
        }
      }
      
      # Add edge weights (if specified)
      if (add_weights) {
        edge_weights <- matrix(runif(sum(adj_matrix)), ncol = 1)
        adj_matrix <- adj_matrix * edge_weights
      }
      
      # Add node attributes (if specified)
      if (add_attributes) {
        node_attributes <- matrix(runif(N), ncol = 1)
      }
      
      # Compute Laplacian matrix and check eigenvalue gap
      degree_matrix <- diag(rowSums(adj_matrix))
      laplacian_matrix <- degree_matrix - adj_matrix
      eigenvalues <- eigen(laplacian_matrix)$values
      eig_gap <- eigenvalues[2] - eigenvalues[1]
      
      # Return graph and relevant information
      graph
      
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      # Function to grow a spatial scale-free expander graph
      # with the specified properties
      spatial_scale_free_expander <- function(N, lambda, beta, sigma, alpha, mu, L, p, m, weight=FALSE, attribute=FALSE){
        
        # Define the d-dimensional torus
        d <- 2
        torus <- function(x) pmax(pmin(x, L), -L)
        
        # Generate N nodes using a Poisson point process on the torus
        points <- torus(matrix(rnorm(N*d, mean=0, sd=1/lambda), ncol=d))
        
        # Perturb the nodes using jittering to add spatial structure
        points <- points + matrix(runif(N*d, min=-0.5, max=0.5), ncol=d)
        
        # Calculate the distance matrix
        dist_mat <- as.matrix(dist(points))
        
        # Define the cutoff distance
        r <- mean(dist_mat) / 2
        
        # Create the adjacency matrix based on the distance matrix
        adj_mat <- ifelse(dist_mat <= r, 1, 0)
        
        # Grow the graph using Barabasi-Albert model
        for (i in 1:(N-m)){
          
          # Calculate the probability of attachment to other nodes
          dist_i <- dist(points[i,], points)
          p_i <- 1/(dist_i^beta + sigma*rowSums(adj_mat[i,]*points[i,]^2)^alpha)
          p_i[i] <- 0
          
          # Choose m nodes to attach to based on the probability of attachment
          attach_to <- sample(1:(i-1), m, replace=TRUE, prob=p_i/(sum(p_i)))
          
          # Add edges to the selected nodes
          adj_mat[i,attach_to] <- 1
          adj_mat[attach_to,i] <- 1
          
          # Add small world effect by rewiring edges with probability p
          if (p > 0){
            for (j in attach_to){
              if (runif(1) < p){
                k <- sample(setdiff(1:N, c(i,j)), 1)
                adj_mat[i,j] <- 0
                adj_mat[j,i] <- 0
                adj_mat[i,k] <- 1
                adj_mat[k,i] <- 1
              }
            }
          }
        }
        
        # Create graph object
        g <- graph.adjacency(adj_mat, mode="undirected")
        
        # Add weights and attributes if specified
        if (weight){
          E(g)$weight <- runif(ecount(g), min=0, max=1)
        }
        if (attribute){
          V(g)$attribute <- sample(1:mu, N, replace=TRUE)
        }
        
        # Compute Laplacian matrix
        Lap <- laplacian_matrix(g)
        
        # Return graph object and Laplacian matrix
        return(list(graph=g, Laplacian=Lap))
      }
      
      # Function to compute Laplacian matrix
      laplacian_matrix <- function(g){
        
        # Get adjacency matrix
        adj_mat <- get.adjacency(g)
        
        # Get degree matrix
        deg <- degree(g)
        D <- diag(deg)
        
        # Compute Laplacian matrix
        Lap <- D - adj_mat
        
        return(Lap)
      }
      
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Function to grow a spatial scale-free expander graph
      # Arguments:
      # N : Number of nodes to be generated
      # lambda : Intensity parameter for Poisson point process
      # r : Cutoff distance for adjacency matrix
      # beta : Power law parameter for preferential attachment
      # mu : Parameter for community structure
      # sigma : Parameter controlling the decay of degree effect with distance
      # alpha : Parameter controlling the strength of degree effect
      # L : Number of rewiring iterations for small world effect
      # p : Probability of rewiring
      # m : Number of edges to attach for each new node
      # Attributes: 
      # weights : Edge weights based on distance
      # nodes : Node attributes with spatial coordinates
      # laplacian : Laplacian matrix of the graph
      grow_spatial_scale_free_graph <- function(N, lambda, r, beta, mu, sigma, alpha, L, p, m){
        # Initialize the graph
        adj_matrix <- matrix(0, N, N)
        nodes <- matrix(0, N, 2)
        weights <- matrix(Inf, N, N)
        
        # Generate nodes using a 2D Poisson point process
        points <- matrix(runif(2*N), N, 2)
        points <- points * c(1/lambda, 1/lambda)
        points <- points %% matrix(rep(c(1, 1), N), N, 2) # periodic boundary condition
        
        # Add nodes to the graph one at a time
        for (i in 1:N){
          nodes[i,] <- points[i,] # save the coordinates of the node
          
          # Calculate the distance matrix between nodes
          dist_matrix <- as.matrix(dist(nodes[1:i,]))
          
          # Create the adjacency matrix based on the distance matrix
          adj_matrix[i,1:(i-1)] <- as.numeric(dist_matrix[i,1:(i-1)] <= r)
          adj_matrix[1:(i-1),i] <- adj_matrix[i,1:(i-1)]
          
          # Calculate the edge weights based on distance
          weights[i,1:(i-1)] <- 1/(dist_matrix[i,1:(i-1)]^sigma)
          weights[1:(i-1),i] <- weights[i,1:(i-1)]
          
          # Perform preferential attachment
          if (i > 1){
            # Calculate the attachment probability based on distance and degree
            prob <- rep(0, i-1)
            for (j in 1:(i-1)){
              dist_factor <- 1/(dist_matrix[i,j]^beta)
              degree_factor <- (adj_matrix[i,j] + adj_matrix[j,i] + 1)^alpha
              community_factor <- (mu + (1-mu)*(1 - (sum(adj_matrix[i,]) + sum(adj_matrix[j,]))/(2*(i-1))))
              prob[j] <- dist_factor * degree_factor * community_factor
            }
            prob <- prob/sum(prob)
            
            # Perform preferential attachment with m edges
            for (k in 1:m){
              attach_node <- sample(1:(i-1), 1, prob = prob)
              adj_matrix[i,attach_node] <- 1
              adj_matrix[attach_node,i] <- 1
              weights[i,attach_node] <- 1/(dist_matrix[i,attach_node]^sigma)
              weights[attach_node,i] <- weights[i,attach_node]
            }
          }
        }
      }     
        # Add small world effect by rewiring edges
        
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Function to generate a spatial scale-free expander graph
        # with specified parameters
        # Arguments:
        #   - N: number of nodes
        #   - lambda: intensity parameter of Poisson point process
        #   - L: size of torus in each dimension
        #   - r: cutoff distance for adjacency matrix
        #   - beta: power law parameter for preferential attachment
        #   - m: number of edges to attach for each new node
        #   - alpha: degree effect parameter for preferential attachment
        #   - sigma: decay parameter for degree effect
        #   - mu: community structure parameter
        #   - p: rewiring probability for small world effect
        #   - weight: function to compute edge weights
        #   - node_attr: function to compute node attributes
        # Output:
        #   - adjacency matrix of the graph
        expand_graph <- function(N, lambda, L, r, beta, m, alpha, sigma, mu, p, weight=NULL, node_attr=NULL) {
          
          # Generate initial points with Poisson point process on torus
          set.seed(123) # for reproducibility
          points <- matrix(rpois(N*L^2/lambda), ncol=2) %% L
          if (!is.null(node_attr)) {
            node_attr <- apply(points, 1, node_attr)
          }
          
          # Add edges one at a time using Barabasi-Albert model with spatial constraints
          A <- matrix(0, nrow=N, ncol=N)
          degrees <- rep(0, N)
          for (i in seq_len(N)[-1]) {
            # Compute distance-dependent attachment probabilities
            dists <- apply(points[seq_len(i-1),], 1, function(x) sqrt(rowSums((x-points[i,])^2)))
            attach_probs <- 1/dists^beta}
            attach_probs[is.infinite(attach_probs)] <- 0 # handle divide-by-zero errors
            # Compute degree-dependent attachment probabilities
            degree_probs <- (degrees[seq_len(i-1)] + alpha) / sum(degrees[seq_len(i-1)] + alpha)
            # Combine spatial and degree probabilities
            probs <- mu * attach_probs * degree_probs + (1-mu) * degree_probs
            # Choose m nodes to attach to
            if (m == 1) {
              attach_nodes <- sample(seq_len(i-1), size=m, prob=probs)
            } else {
              attach_nodes <- numeric(m)
              for (j in seq_len(m)) {
                attach_nodes[j] <- sample(seq_len(i-1), size=1, prob=probs)
                probs[attach_nodes[j]] <- 0 # avoid selecting same node twice
              }
            }
            # Add edges to the selected nodes
            for (j in attach_nodes) {
              A[i,j] <- 1
              A[j,i] <- 1
              degrees[i] <- degrees[i] + 1
              degrees[j] <- degrees[j] + 1
            }
            # Add new point with jittering
            new_point <- rnorm(2, sd=0.1) + points[i,]
            new_point <- new_point %% L
            points <- rbind(points, new_point)
            if (!is.null(node_attr)) {
              node_attr <- rbind(node_attr, node_attr[i,])
            }
          }
          
          # Compute edge weights if specified
          if (!is.null(weight)) {
            weights <- apply(A, c(1,
                                  
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Function to generate a spatial scale-free expander graph
spatial_scalefree_expander <- function(N, lambda, beta, L, m, p = 0, mu = 0, probWithin = 0.5, probBetween = 0.1, create_weights = FALSE, create_node_attributes = FALSE) {
  
  # Generate points using a 2D Poisson point process with intensity parameter lambda on a torus
  points <- matrix(runif(N * 2, 0, L), N, 2) # generate N points in a 2D torus
  d <- as.matrix(dist(points, method = "euclidean")) # calculate distance matrix
  
  # Apply periodic boundary condition
  d[d > L / 2] <- L - d[d > L / 2]
  
  # Add spatial structure to the distribution by jittering
  points <- points + matrix(rnorm(N * 2, 0, lambda), N, 2)
  
  # Initialize adjacency matrix
  adj_matrix <- matrix(0, N, N)
  
  # Add edges based on preferential attachment
  # First, connect the first two nodes
  adj_matrix[1, 2] <- 1
  adj_matrix[2, 1] <- 1
  
  # Initialize degree and distance vectors
  degree <- rep(1, N)
  distance <- d[1, ]
  
  # Loop over remaining nodes
  for (i in 3:N) {
    # Add m edges to existing nodes with probability proportional to their degree and distance
    probs <- (degree^mu / distance^beta) * (1 - mu) # probability of attachment to existing nodes
    within_idx <- which(i %% (N/probWithin) == seq(0, N/probWithin - 1) * (N/probWithin)) # indices of nodes within the same community
    between_idx <- setdiff(1:N, within_idx) # indices of nodes in other communities
    probs[between_idx] <- probs[between_idx] * (1 - probBetween) # reduce probability of attachment to nodes in other communities
    probs[within_idx] <- probs[within_idx] * probBetween # increase probability of attachment to nodes in the same community
    
    if (m > 1) { # add multiple edges
      for (j in 1:m) {
        new_edge <- sample(1:N, size = 1, prob = probs)
        adj_matrix[i, new_edge] <- 1
        adj_matrix[new_edge, i] <- 1
        degree[c(i, new_edge)] <- degree[c(i, new_edge)] + 1
        distance[c(i, new_edge)] <- d[i, new_edge]
      }
    } else { # add one edge
      new_edge <- sample(1:N, size = 1, prob = probs)
      adj_matrix[i, new_edge] <- 1
      adj_matrix[new_edge, i] <- 1
      degree[c(i, new_edge)] <- degree[c(i, new_edge)] + 1
      distance[c(i, new_edge)] <- d[i, new_edge]
    }
    
    # Rewire edges with probability p
    # for (j in 1:(i-1)) {
    #   if (adj_matrix[i, j] == 1 & runif(1) < p) {
    #     new_edge <- sample(setdiff(1:N, c(i, j)), size = 1)
    #     adj_matrix[i, j] <- 0
    #     
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# dist_probs=exp(-dist_mat[i,] / (min_dist))
# pref_probs=(deg[1:(i+1)]^(-beta))
# # Compute the attachment probabilities for the new node
# #### (ii) p_i <- (deg/(2*sum(deg)))*(p_ij^prewire)#change p to prewire
# p_i=mu*dist_probs[1:(i+1)]*pref_probs + (1-mu)*rep(sum(dist_probs[1:(i+1)]*pref_probs)/sum(pref_probs), (i+1)) 
# p_i[which(!is.finite(p_i))] <- 0
# #### (iii) p_i[which(!is.finite(p_i))] <- 0
# #   print(p_i)
# 
# #Add m edges to existing nodes with probability proportional to their degree and distance
# for (j in 1:m){
#   # Choose a node to attach to
#   node_to_attach <- sample(1:(i+1), size=m, replace = TRUE, prob=p_i)



# Function to generate a Spatial Scale-Free Graph
# Arguments:
#  N: number of nodes
#  d: number of dimensions for the torus (default=2)
#  lambda: intensity parameter for the Poisson point process (default=1)
#  L: size of the torus (default=1)
#  beta: power-law exponent for preferential attachment (default=1)
#  m: number of edges added for each new node (default=1)
#  p: rewiring probability (default=0)
#  probWithin: probability of connecting nodes within communities (default=0)
#  probBetween: probability of connecting nodes between communities (default=0)
#  attributes: list of node attributes (default=NULL)
#  edge_weights: logical indicating whether to include edge weights (default=FALSE)
# Returns:
#  a graph object of class "igraph"
spatial_sf_graph <- function(N, d=2, lambda=1, L=1, beta=1, m=1, p=0, 
                probWithin=0, probBetween=0, attributes=NULL, edge_weights=FALSE) {
  # Generate N points on a d-dimensional torus using Poisson point process
  points <- matrix(runif(N*d), ncol=d) * L
  if (d == 2) {
    # Apply spatial jittering for 2D torus
    points <- points + matrix(runif(N*d), ncol=d) * lambda/L
  } else {
    # Apply spatial jittering for higher dimensions
    points <- points + matrix(rnorm(N*d), ncol=d) * lambda/L
  }
  # Create an adjacency matrix for the graph
  adj_matrix <- matrix(0, nrow=N, ncol=N)
  # Compute the distance matrix between all pairs of points
  dist_matrix <- as.matrix(dist(points, diag=TRUE, upper=TRUE))
  min_dist <- min(dist_matrix[dist_matrix > 0]) # Minimum non-zero distance
  # Initialize the graph with the first node
  #Create adjacency matrix with preferential attachment and small-world rewiring
  adj_matrix <- matrix(0, nrow = N, ncol = N)
  #graph <- numeric(N) # Track degree of each node
  #graph[1] <- 1 # Start with a single node
   Graph <- graph.empty(N, directed=FALSE)
   Graph <- add.edges(Graph,c(1,2),attr=list(width=10))
  # Loop through the rest of the nodes
  for (i in 2:N) {
    # Compute the probability of attachment to existing nodes
    prob_attach <- rep(0, i-1)
    for (j in 1:(i-1)) {
      # Compute the distance between the new node and existing nodes
      d_ij <- dist_matrix[i, j]
      # Compute the degree of the existing node using its index
      k_j <- degree(Graph, j)
      # Compute the product of distance-based and scale-free attachment probabilities
      prob_attach[j] <- (d_ij^-beta) * (k_j^beta)
    }
    # Normalize the attachment probabilities
    prob_attach <- prob_attach/sum(prob_attach)
    # Choose m existing nodes to attach to based on their probabilities
    attach_to <- sample(1:(i-1), m, replace=TRUE, prob=prob_attach)}
    # Add edges between the new node and the chosen nodes
    for (j in attach_to) {
      # Check if an edge already exists
      # Small-world rewiring with probability p
      if (p > 0 && runif(1) < p) {
        neighbors <- which(adj_matrix[i,j] == 1)
        if (length(neighbors) > 0) {
          new_neighbor <- sample(neighbors, 1)
          adj_matrix[j, new_neighbor] <- 0
          adj_matrix[new_neighbor, j] <- 0
          available_nodes <- setdiff(1:N, c(j, neighbors))
          new_neighbor <- sample(available_nodes, 1)
          adj_matrix[j, new_neighbor] <- 1
          adj_matrix[new_neighbor, j] <- 1
        }
      }
    }
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
grow_spatial_scale_free_graph <- function(beta, N, lambda, L, p, m, probWithin = 0.5, probBetween = 0.1, add_edge_weights = FALSE, add_node_attributes = FALSE) {
     # Generate points using a 2D Poisson point process with intensity parameter lambda on a d-dimensional torus
     coords <- matrix(runif(N*L, 0, L), ncol = L)
     coords_jittered <- coords + runif(N*L, -0.5, 0.5)/lambda
     coords_jittered <- coords_jittered %% L
     
     # Calculate the distance matrix between all pairs of points
     dist_mat <- as.matrix(dist(coords_jittered))
     
     # Create an adjacency matrix for the graph using preferential attachment
     adj_mat <- matrix(0, nrow = N, ncol = N)
     degree <- numeric(N)
     degree[1] <- 1
     for (i in 2:N) {
       # Compute probability of attachment based on spatial distance and preferential attachment
       prob_attach <- (dist_mat[i,1:(i-1)] + 1e-6) ^ (-beta) * degree[1:(i-1)]
       prob_attach[i] <- 0 # exclude self-attachment
       probsAttach <- prob_attach/sum(prob_attach)
       x2=probsAttach
       # Sample m nodes to attach to
       attach_to <- sample(1:N, size=1, replace = TRUE, prob = x2)
       
       # Connect to selected nodes with probability probWithin or probBetween
       for (j in attach_to) {
         if (runif(1) < probWithin) {
           adj_mat[i,j] <- 1
           adj_mat[j,i] <- 1
         } else if (runif(1) < probBetween) {
           available_nodes <- setdiff(1:N, c(i, j, which(adj_mat[i,] == 1), which(adj_mat[j,] == 1)))
           if (length(available_nodes) > 0) {
             to_node <- sample(available_nodes, 1)
             adj_mat[i,to_node] <- 1
             adj_mat[to_node,i] <- 1
           }
         }
       }
       
       # Update degree vector
       degree[i] <- sum(adj_mat[i,])
       degree[attach_to] <- degree[attach_to] + 1
     }
     
     # Rewire edges with probability p for small-world effect
     for (i in 1:(N-1)) {
       for (j in (i+1):N) {
         if (adj_mat[i,j] == 1 && runif(1) < p) {
           available_nodes <- setdiff(1:N, c(i, j, which(adj_mat[i,] == 1), which(adj_mat[j,] == 1)))
           if (length(available_nodes) > 0) {
             to_node <- sample(available_nodes, 1)
             adj_mat[i,j] <- 0
             adj_mat[j,i] <- 0
             adj_mat[i,to_node] <- 1
             adj_mat[to_node,i] <- 1
           }
         }
       }
     }
     
     # Add edge weights if requested
     if (add_edge_weights) {
       edge_weights <- runif(sum(adj_mat), 0, 1)
       adj_mat[adj_mat == 1] <- edge_weights
     }
     
     # Add node attributes if requested
}    
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   # Function to generate a spatial scale-free expander graph
   spatial_scale_free_expander <- function(N, lambda, L, beta, m, p = 0, 
                                           probWithin = 0.5, probBetween = 0.1,
                                           add_edge_weight = FALSE, 
                                           add_node_attr = FALSE) {
     
     # Set seed for reproducibility
     #set.seed(123)
     
     # Generate N points using a 2D Poisson point process with intensity parameter lambda
     points <- matrix(runif(2*N, min=0, max=L), ncol=2)
     
     # Apply periodic boundary condition to ensure points are distributed on a torus
     points <- points %% L
     
     # Add spatial structure to the point distribution using jittering
     jitter_factor <- 0.1 # adjust this to control the degree of perturbation
     points <- points + matrix(runif(2*N, min=-jitter_factor, max=jitter_factor), ncol=2)
     
     
     # Calculate distance matrix between all pairs of points
     dist_mat <- as.matrix(dist(points))
     
     # Set the adjacency matrix element A_{ij} to 1 if dist_{ij} <= d_cutoff and 0 otherwise
     d_cutoff <- 0.2 # adjust this to control the degree of sparsity
     adj_mat <- as.matrix(dist_mat <= d_cutoff)
     
     # Initialize the degree vector and the edge list
     deg <- colSums(adj_mat)
     edge_list <- matrix(ncol=2, nrow=0)
     
     # Grow the network with nodes added one at a time
     for (i in 1:(N-2)) {
       # Compute the distance-dependent probability function
       p_ij <- 1/dist_mat[,(i+1)]^beta
       # Compute the attachment probabilities for the new node
       p_i <- (deg/(2*sum(deg)))*(p_ij^p)
      # p_i[i]=0}
       p_i[which(adj_mat[,i+1])] <- 0
       #node_to_attach <- sample(1:N, size=1, replace = TRUE, prob=p_i)}
       # Add m edges to existing nodes with probability proportional to their degree and distance
       for (j in 1:m){
         #x1=1:(i+2)
         # Choose a node to attach to
         node_to_attach <- sample(1:N, size=m, replace = TRUE, prob=p_i)
         # Add the edge if it doesn't already exist
         if (adj_mat[(i+1), node_to_attach] == 0) {
           edge_list <- rbind(edge_list, c(i+1, node_to_attach))
           degree[i+1] <- degree[i+1] + 1
           degree[node_to_attach] <- degree[node_to_attach] + 1
           adj_mat[i+1, node_to_attach] <- 1
           adj_mat[node_to_attach, i+1] <- 1
         }
     }}
         
         # # Rewire the edge with probability p
         # if (runif(1) < p) {
         #   # Choose a random node within a certain distance
         #   nodes_within_dist <- which(dist_mat[,i+1] <= d_cutoff & degree > 0)
         #   node_to_rewire <- sample(nodes_within_dist, size=1)
         #   # Rewire the edge if it doesn't already exist
         #   if (adj_mat[i+1, node_to_rewire] == 0) {
         #     adj_mat[i+1, node_to_attach] <- 0
         #     adj_mat[node_to_attach, i+1] <- 0
         #     adj_mat[i+
                       
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     library(Matrix)
     
     generate_points <- function(N, L, lambda, d) {
       # Generate N points using a 2D Poisson point process with intensity parameter lambda
       # on a d-dimensional torus with side length L
       
       # Calculate expected number of points
       expected_points <- floor(lambda * L^d)
       
       # Generate Poisson process
       x <- matrix(runif(N*d), ncol=d) * L
       points <- x[runif(expected_points) < lambda, ]
       
       # Add spatial structure to the points by jittering
       jitter_amount <- L * 0.01
       jitter <- matrix(runif(N*d, -jitter_amount, jitter_amount), ncol=d)
       points <- points + jitter
       
       # Apply periodic boundary condition to ensure torus topology
       points <- points %% L
       
       return(points)
     }
     
     dist_matrix <- function(points) {
       # Calculate distance matrix between all pairs of points
       dist_mat <- as.matrix(dist(points))
       return(dist_mat)
     }
     
     adj_matrix <- function(dist_mat, r) {
       # Create adjacency matrix based on distance matrix and cutoff distance
       adj_mat <- ifelse(dist_mat <= r, 1, 0)
       return(adj_mat)
     }
     
     preferential_attachment <- function(adj_mat, m, beta) {
       # Perform preferential attachment to grow network
       
       N <- nrow(adj_mat)
       
       # Initialize with fully connected graph of m+1 nodes
       G <- matrix(1, nrow=m+1, ncol=m+1) - diag(1, m+1)
       degree <- rep(m, m+1)
       
       # Loop over remaining nodes to add
       for (i in (m+2):N) {
         # Calculate probability of attaching to existing nodes
         p <- (degree^beta) / sum(degree^beta)
         
         # Select m nodes to attach to
         attach_nodes <- sample(1:(m+1), m, replace=T, prob=p)
         
         # Add edges to selected nodes
         G[i, attach_nodes] <- 1
         G[attach_nodes, i] <- 1
         
         # Update degree of attached nodes
         degree[attach_nodes] <- degree[attach_nodes] + 1
         degree[i] <- m}
         
         # Rewire edges with probability p_rewire
         p_rewire <- 0.1
         for (j in attach_nodes) {
           if (runif(1) < p_rewire) {
             # Select random node within cutoff distance
             dist_j <- which(dist_mat[j,] <= d_cutoff)
             rewire_node <- sample(dist_j, 1)
             while (rewire_node == i || adj_mat[i, rewire_node] == 1) {
               rewire_node <- sample(dist_j, 1)
             }
             
             # Rewire edge
             G[i, j] <- 0
             G[j, i] <- 0
             G[i, rewire_node] <- 1
             G[rewire_node, i] <- 1
           }
         }
       }
       
       return(G)
     }
     
     stochastic_block_model <- function(N, K, probWithin, probBetween) {
       # Perform stochastic block model to add community structure
       
       # Assign nodes to communities
       community_size <- rep(floor(N/K), K)
       community_size[1:(N %% K)] <- community_size[1:(N %% K)] + 1
       community_assignments <- rep(1
                                    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# function to generate spatial scale-free graph
spatial_sf_graph <- function(beta, N, lambda, L, p, probWithin, probBetween, weights = FALSE, nodeAttrs = FALSE) {
  
  # generate nodes using 2D Poisson point process on a torus
  coords <- matrix(runif(N*2, 0, L), ncol=2)
  
  # calculate distance matrix between all pairs of nodes
  dist_mat <- as.matrix(dist(coords, method = "euclidean"))
  dist.min <- pmin(dist_mat, L - dist_mat)#minimum distance
  # initialize adjacency matrix
  adj_mat <- matrix(0, nrow = N, ncol = N)
  
  # calculate probability of attachment for each node
  pa_prob <- rep(1, N)
  
  for (i in 2:N) {
    # preferential attachment
    pa_prob[i] <- sum(pa_prob[1:(i-1)]^beta)
    # spatial attachment
    spatial_dist <- dist.min[1:(i-1), i]
    pa_prob[i] <- pa_prob[i] * exp(-spatial_dist / (lambda * L^(1/2)))}
    # normalize probabilities
    pa_prob[i] <- pa_prob[i] / sum(pa_prob[1:i])
  }
  
  # add nodes to the graph sequentially
  for (i in 2:N) {
    # preferential attachment
    pa_sample <- sample(1:(i-1), size = 1, prob = pa_prob[1:(i-1)])
    adj_mat[pa_sample, i] <- 1
    adj_mat[i, pa_sample] <- 1
    # spatial attachment
    spatial_sample <- sample(1:(i-1), size = 1, prob = exp(-dist_mat[1:(i-1), i] / (lambda * L^(1/2))))
    adj_mat[spatial_sample, i] <- 1
    adj_mat[i, spatial_sample] <- 1
  }
  
  # rewire edges with small probability
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (adj_mat[i,j] == 1 & runif(1) < p) {
        adj_mat[i,j] <- 0
        adj_mat[j,i] <- 0
        rewired <- FALSE
        while (!rewired) {
          k <- sample(1:N, size = 1)
          if (k != i & k != j & adj_mat[i,k] == 0 & adj_mat[k,j] == 0) {
            adj_mat[i,k] <- 1
            adj_mat[k,i] <- 1
            adj_mat[k,j] <- 1
            adj_mat[j,k] <- 1
            rewired <- TRUE
          }
        }
      }
    }
  }
  
  # add community structure
#  if (probWithin > 0 & probBetween > 0) {

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    library(Matrix)
    library(mvtnorm)
    
    growSpatialScaleFreeGraph <- function(N, lambda, L, beta, p = 0, probWithin = 0, probBetween = 0, edgeWeights = FALSE, nodeAttributes = FALSE) {
      
      # Generate points using a Poisson point process on a d-dimensional torus
      points <- rmvnorm(N, rep(0, L), diag(rep(1, L)))
      points <- points %% 1 # wrap around to torus
      
      # Calculate distance matrix
      distMat <- as.matrix(dist(points, diag = TRUE, upper = TRUE))
      
      # Initialize adjacency matrix
      adjMat <- matrix(0, nrow = N, ncol = N)
      
      # Attach first node
      adjMat[1, 1] <- 1
      
      # Create function for calculating attachment probabilities
      calcProb <- function(i, j) {
        dist_ij <- distMat[i, j]
        degree_j <- sum(adjMat[j, ])
        return (probWithin * adjMat[i, j] + probBetween * (1 - adjMat[i, j])) * (dist_ij ^ (-beta)) * (degree_j ^ beta)
      }
      
      # Attach remaining nodes
      for (i in 2:N) {
        # Calculate attachment probabilities
        probs <- sapply(1:(i-1), calcProb, i = i)
        probs[is.na(probs)] <- 0 # set probabilities to 0 for disconnected nodes
        
        # Apply rewiring probability for small world effect
        if (p > 0) {
          for (j in 1:(i-1)) {
            if (runif(1) < p) {
              probs[j] <- calcProb(i, sample(1:(i-1), 1))
            }
          }
        }
        
        # Normalize probabilities and randomly attach node
        probs <- probs / sum(probs)
        j <- sample(1:(i-1), 1, prob = probs)
        adjMat[i, j] <- 1
        adjMat[j, i] <- 1
      }
      
      # Create edge weights
      if (edgeWeights) {
        w <- runif(sum(adjMat), 0, 1)
        adjMat[adjMat == 1] <- w
      }
      
      # Create node attributes
      if (nodeAttributes) {
        attributes <- matrix(runif(N * 2), nrow = N, ncol = 2)
        return (list(adjMat = adjMat, attributes = attributes))
      }
      
      return (adjMat)
    }
    
    # Compute Laplacian matrix
    computeLaplacian <- function(adjMat) {
      degVec <- rowSums(adjMat)
      degMat <- diag(degVec)
      lapMat <- degMat - adjMat
      return (lapMat)
    }
    
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
growSpatialScaleFreeGraph <- function(N, lambda, L, beta, p = 0, probWithin = 0.5, probBetween = 0.1, edgeWeights = FALSE, nodeAttributes = FALSE) {
      
# Generate N points using 2D Poisson point process with intensity lambda on a torus of size L^d
  d <- 2 # Dimension
  xy <- matrix(runif(N*d, 0, L), ncol = d) # Generate random coordinates in [0,L]^d
  xy_jittered <- xy + runif(N*d, -lambda/2, lambda/2) # Jitter points with uniform noise of amplitude lambda
      
  # Calculate distance matrix between all pairs of points and minimum distance between pairs of nodes
  dist_matrix <- as.matrix(dist(xy_jittered, method = "euclidean", diag = TRUE, upper = TRUE))
  min_dist <- min(dist_matrix[dist_matrix > 0]) # Minimum non-zero distance
      
  # Create adjacency matrix with preferential attachment and small-world rewiring
  adj_matrix <- matrix(0, nrow = N, ncol = N)
  degrees <- numeric(N) # Track degree of each node
  degrees[1] <- 1 # Start with a single node
    for (i in 2:N) {
        # Calculate probability of attachment to each existing node
      dist_probs <- exp(-dist_matrix[i,]) # Favor short spatial distances
      pref_probs <- (degrees[1:(i-1)]^(-beta)) # Scale-free preferential attachment
      # Combine probabilities
      # attach_probs <- probWithin*dist_probs[1:i-1]*pref_probs + 
      attach_probs= rep(sum(dist_probs[1:(i-1)]*pref_probs)/sum(pref_probs), i-1)}
        
      # Choose nodes to attach to based on the attachment probabilities
      attach_nodes <- sample(1:(i-1), 1, replace=F, prob = attach_probs)
      adj_matrix[i, attach_nodes] <- 1
      adj_matrix[attach_nodes, i] <- 1
      degrees[c(i, attach_nodes)] <- degrees[c(i, attach_nodes)] + 1
        
        # Small-world rewiring with probability p
      if (p > 0 && runif(1) < p) {
        neighbors <- which(adj_matrix[i,] == 1)
        if (length(neighbors) > 0) {
          new_neighbor <- sample(neighbors, 1)
          adj_matrix[i, new_neighbor] <- 0
          adj_matrix[new_neighbor, i] <- 0
          available_nodes <- setdiff(1:N, c(i, neighbors))
          new_neighbor <- sample(available_nodes, 1)
          adj_matrix[i, new_neighbor] <- 1
          adj_matrix[new_neighbor, i] <- 1
        }
        }
      }
      
      # Create edge weights if requested
      if (edgeWeights) {
        adj_matrix[adj_matrix == 1] <- runif(sum(adj_matrix == 1))
      }
      
      # Create node attributes if requested
      if (nodeAttributes) {
        node_attrs <- data.frame(x = xy_jittered[,1], y = xy_jittered[,2], degree = degrees)
        return(list(adj_matrix = adj_matrix, node_attrs = node_attrs))
      } 
      else {
        return(adj_matrix)
      }
  ### Graph object
  # GraphModel <- graph.adjacency(as.matrix(adj_matrix), mode="undirected")
  # GraphModel=simplify(GraphModel,remove.loops = TRUE,remove.multiple = TRUE)
  # graph <- list()
  # graph$adj_mat=adj_mat
  # graph$GraphObject <- GraphModel
  # graph$ExpansionRatio=expansion_ratio
 # return(GraphModel)
}      

g1=growSpatialScaleFreeGraph(N=100, lambda=30, L=2, beta=8, p = 0.9, probWithin = 0.1, probBetween = 0.1, edgeWeights = FALSE, nodeAttributes = FALSE) 
GraphModel <- graph.adjacency(as.matrix(g1), mode="undirected")
GraphModel=simplify(GraphModel,remove.loops = TRUE,remove.multiple = TRUE)
plot(GraphModel,vertex.label=NA,vertex.size=2)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
grow_spatial_scale_free_graph <- function(N, beta, lambda, L, p, probWithin = NULL, probBetween = NULL, create_edge_weights = FALSE, create_node_attributes = FALSE) {
  # Generate points on torus with jittering
  d <- 2
  grid <- seq(0, 1, length.out = floor(sqrt(lambda*N)))
  points <- expand.grid(replicate(d, grid, simplify = FALSE))
  jittered_points <- points + runif(d*N, -L/2, L/2)
  jittered_points <- jittered_points %% L
  
  # Calculate distance matrix
  dist_mat <- as.matrix(dist(jittered_points))
  diag(dist_mat) <- Inf
  min_dists <- apply(dist_mat, 1, min)
  
  # Initialize adjacency matrix
  adj_mat <- matrix(0, nrow = N, ncol = N)
  
  # Add nodes sequentially with preferential attachment
  for (i in 1:N) {
    # Probability of attachment based on distance and degree
    dist_probs <- exp(-min_dists^beta) * rowSums(adj_mat)^(-1)
    dist_probs[i] <- 0
    dist_probs <- dist_probs / sum(dist_probs)
    
    # Preferential attachment
    new_edges <- sample(1:N, size = 5, replace = TRUE, prob = dist_probs)
    new_edges <- new_edges[new_edges != i]
    
    # Add edges with rewiring probability
    for (j in new_edges) {
      if (runif(1) < p) {
        possible_edges <- setdiff(1:N, c(i, j, which(adj_mat[i,] == 1), which(adj_mat[j,] == 1)))
        if (length(possible_edges) > 0) {
          rewired_node <- sample(possible_edges, size = 1)
          adj_mat[i, j] <- 0
          adj_mat[j, i] <- 0
          adj_mat[i, rewired_node] <- 1
          adj_mat[rewired_node, i] <- 1
        }
      } else {
        adj_mat[i, j] <- 1
        adj_mat[j, i] <- 1
      }
    }
  }
  
  # Create community structure if specified
  if (!is.null(probWithin) & !is.null(probBetween)) {
    n_communities <- ceiling(sqrt(N))
    community_labels <- rep(1:n_communities, each = N/n_communities)
    within_community_probs <- matrix(probWithin, nrow = n_communities, ncol = n_communities)
    diag(within_community_probs) <- 1
    between_community_probs <- matrix(probBetween, nrow = n_communities, ncol = n_communities)
    diag(between_community_probs) <- 0
    
    for (i in 1:N) {
      for (j in 1:N) {
        if (i != j) {
          if (community_labels[i] == community_labels[j]) {
            if (runif(1) < within_community_probs[community_labels[i], community_labels[j]]) {
              adj_mat[i, j] <- 1
              adj_mat[j, i] <- 1
            }
          } else {
            if (runif(1) < between_community_probs[community_labels[i], community_labels[j]]) {
              adj_mat[i, j] <- 1
              adj_mat[j, i] <- 1
            }
          }
        }
      }
    }
  }
  
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
        