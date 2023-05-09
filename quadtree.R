#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
setwd("C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi-New")
source("SPATIAL-PIPELINE-NEW.R")
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Spatial Expander Propagation Graph with weighted edges, 
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

spatial_scale_free_expander <- function(N, beta, m,delta, prewire = 0.1, 
                                        mu=0.2,add_edge_weight = FALSE, 
                                        add_node_attr = FALSE,r=0.1,alpha=0.2) {
  #Notes:
  ##--(1) use of of N or lambda
  ##--(2) fix L or ignore
  ###--(3) mu > 0 for preferential attachment: mu captures what p does in the doi:10.1088/1367-2630/9/6/190
  ###--(4) use unit square instead of torus
  
  ##--Generate points on a unit square
  points <- matrix(runif(N*2), ncol=2)
  
  #--Calculate distance matrix between all pairs of points
  #dist_mat <- as.matrix(dist(points))
  
  
  #--Set the adjacency matrix element A_{ij} to 1 if dist_{ij} <= r and 0 otherwise
  #adj_mat <- ifelse(dist_mat <=r,1,0)
  #excludes self from degree 
  #diag(adj_mat) <- 0 
  
  #--Initialize the degree vector and the edge list
  #degree <- colSums(adj_mat)
  edge_list <- matrix(ncol=2, nrow=0)
  
  # create initial random graph with 2 nodes
  # calculate distance matrix between all pairs of points
  dist_mat <- as.matrix(dist(points))
  
  
  # create adjacency matrix based on distance matrix
  adj_mat <- ifelse(dist_mat <= r, 1, 0)
  
  #excludes self loops 
  diag(adj_mat) <- 0 
  
  
  # initialize graph with a single node and no edges
  graph <- list()
  graph$adj_mat <- matrix(0, nrow = N, ncol = N)
  graph$adj_mat[1, 1] <- 1
  graph$degrees <- rep(0, N)
  
  # Grow the network with new nodes added one at a time
  for (i in 2:N) {
    # calculate probability of attaching to each existing node
    probs <- rep(0, (i-1))
    for (j in 1:(i-1)) {
      # spatial distance effect
      d_ij <- dist(points[c(i, j), ])
      spatial_prob <- ifelse(d_ij <= r, 1/1+(d_ij^alpha), 0)
      # community structure effect
      # community_prob <- ifelse(runif(1) < delta, adj_mat[i, j], 0)
      # preferential attachment effect
      degree_effect <- graph$degrees[j]^beta
      # total probability of attaching to node j
      probs[j] <- (spatial_prob + 0.3 +
                     degree_effect) / sum(spatial_prob + 0.3 + degree_effect)
      # randomly select an existing node to attach to based on the probabilities
      attach_node <- sample(1:(i-1), size = 1, prob = probs)
      # add edge between new node and attach_node
      graph$adj_mat[i, attach_node] <- 1
      graph$adj_mat[attach_node, i] <- 1
      graph$degrees[i] <- graph$degrees[i] + 1
      graph$degrees[attach_node] <- graph$degrees[attach_node] + 1}}
  
  # small world effect: randomly rewire edge with probability p
  if (runif(1) < p) {
    possible_rewire_nodes <- setdiff(1:N, c(i, attach_node, which(graph$adj_matrix[i, ] == 1)))
    if (length(possible_rewire_nodes) > 0) {
      rewire_node <- sample(possible_rewire_nodes, size = 1)
      graph$adj_matrix[i, attach_node] <- 0
      graph$adj_matrix[attach_node, i] <- 0
      graph$adj_matrix[i, rewire_node] <- 1
      graph$adj_matrix[rewire_node, i] <- 1
      graph$degrees[i] <- graph$degrees[i] - 1
      graph$degrees[attach_node] <- graph$degrees[attach_node] - 1
      graph$degrees[rewire_node] <- graph}}}}

#--Compute the distance-dependent probability function
#dist_probs <- 1/(1+(dist_mat[i,])^beta+1e-10) 
#--Compute the degree-dependent probability function
#deg_probs<- (degree^alpha)/(sum(degree^alpha))
#--compute community dependent probability
#comm_prob <- delta*comm + (1-delta)*rep(1,N)

# Attachment probabilities for new nodes
f_ij=mu*dist_probs + (1-mu)*deg_probs
f_ij[i]=0 # probability of attachment to itself is zero 

probattach=f_ij/sum(f_ij)#normalize attachment probabilities

probattach[which(probattach == Inf)] <- 0 # remove infinite probabilities


#Add m edges to existing nodes with probability proportional to their degree and distance
for (j in 1:m){
  # Choose a node to attach to
  node_to_attach <- sample(1:N, size=m, replace = TRUE, prob=probattach)
  
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

g=spatial_scale_free_expander(N=100, beta=4, m=1, prewire = 0.1, 
                              add_edge_weight = FALSE,mu=0, 
                              add_node_attr = T,r=0.2,alpha=2,delta = 0.3)

lo=layout.norm(as.matrix(cbind(g$NodeAttributes$x,g$NodeAttributes$y)))

plot(g$GraphObject,vertex.label=NA,vertex.size=2,layout=lo)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
spatial_scale_free_graph <- function(N, beta, lambda, alpha, delta, L, p, m) {
  # Step 1: Generate N points uniformly at random in a unit square
  points <- matrix(runif(N*2), ncol=2)
  
  # Step 2: Compute the distance matrix between all pairs of points
  dist_mat <- as.matrix(dist(points))
  
  # Step 3: Create an adjacency matrix based on the distance matrix, with a cutoff distance r controlled by lambda
  adj_mat <- ifelse(dist_mat <= lambda, 1, 0)
  
  # Step 4: Grow the graph with preferential attachment
  # Initialize the graph with m nodes and a fully connected subgraph
  graph <- matrix(0, nrow=N, ncol=N)
  for (i in 1:m) {
    for (j in (i+1):m) {
      graph[i, j] <- 1
      graph[j, i] <- 1
    }
  }
  # Keep track of node degrees
  degrees <- rep(m, m)
  # Loop over the remaining nodes and add edges using preferential attachment
  for (i in (m+1):N) {
    # Compute the distance-dependent probability function
    prob <- 1/(dist_mat[i,]^beta)
    # Compute the attachment probabilities for each existing node
    attach_probs <- (degrees * prob) + (delta * adj_mat[i,])
    # Normalize the probabilities
    attach_probs <- attach_probs / sum(attach_probs)
    # Sample m existing nodes to attach to, with replacement and using the computed probabilities
    attach_nodes <- sample.int(N, size=m, replace=TRUE, prob=attach_probs)
    # Add edges to the selected nodes and update their degrees
    for (j in attach_nodes) {
      graph[i, j] <- 1
      graph[j, i] <- 1
      degrees[i] <- degrees[i] + 1
      degrees[j] <- degrees[j] + 1
    }
  }
  
  # Step 5: Rewire edges with probability p
  for (i in 1:N) {
    for (j in (i+1):N) {
      if (graph[i, j] == 1 && runif(1) < p) {
        # Find a new node within a certain distance to rewire to
        dist_to_j <- dist_mat[j,]
        new_node <- sample.int(N, size=1, prob=(1/dist_to_j^L) / sum(1/dist_to_j^L))
        # Rewire the edge
        graph[i, j] <- 0
        graph[j, i] <- 0
        graph[i, new_node] <- 1
        graph[new_node, i] <- 1
      }
    }
  }
  
  # Step 6: Add edge weights and node attributes if desired
  
  # Step 7: Check the expansion properties of the graph by computing the Laplacian matrix
  laplacian <- diag(colSums(graph)) - graph
  eigenvalues <- eigen(laplacian, symmetric=TRUE)$values
  expansion <- eigenvalues[2]
  print(paste0("Expansion: ", expansion))
  
  # Step 8: Return the adjacency matrix
  return(adj_mat)
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
generate_spatial_scale_free_graph <- function(N, r, m, alpha, delta, p, beta) {
  # Generate initial nodes in a unit square
  nodes <- matrix(runif(N*2), ncol=2)
  
  # Compute the distance matrix
  dist_mat <- as.matrix(dist(nodes))
  
  # Initialize adjacency matrix with nodes within distance r connected
  adj_mat <- as.matrix(dist_mat <= r)
  diag(adj_mat) <- 0
  
  # Add nodes and edges sequentially
  for (i in seq_len(N-1)) {
    # Compute the probability of attachment for each existing node
    degree_dist <- colSums(adj_mat)
    spatial_dist <- dist_mat[i,]
    prob_attach <- (spatial_dist^(-beta))*(degree_dist^alpha)
    prob_attach[which(adj_mat[i,]==1)] <- 0 # already connected nodes have zero probability
    prob_attach <- prob_attach/sum(prob_attach) # normalize probabilities
    
    # Choose m nodes to connect to with probability proportional to prob_attach
    to_connect <- sample(seq_len(N), size=m, prob=prob_attach, replace=TRUE)
    
    # Connect the new node to the chosen nodes
    adj_mat[i,to_connect] <- 1
    adj_mat[to_connect,i] <- 1
  }
  
  # Add community structure with probability delta
  comm <- rep(1, N)
  num_comm <- round(N*delta)
  comm[1:num_comm] <- sample(seq_len(num_comm))
  comm <- sample(comm)
  for (i in seq_len(N)) {
    for (j in seq_len(N)) {
      if (comm[i] == comm[j]) {
        if (runif(1) < delta) {
          adj_mat[i,j] <- 1
          adj_mat[j,i] <- 1
        }
      } else {
        if (runif(1) < (1-delta)) {
          adj_mat[i,j] <- 1
          adj_mat[j,i] <- 1
        }
      }
    }
  }
  
  # Rewire edges with probability p
  for (i in seq_len(N)) {
    for (j in seq_len(N)) {
      if (adj_mat[i,j] == 1 && runif(1) < p) {
        # Find a new node within distance r to rewire the edge to
        candidates <- which(dist_mat[i,] <= r & adj_mat[i,]==0)
        if (length(candidates) > 0) {
          new_j <- sample(candidates, size=1)
          adj_mat[i,j] <- 0
          adj_mat[j,i] <- 0
          adj_mat[i,new_j] <- 1
          adj_mat[new_j,i] <- 1
        }
      }
    }
  }
  
  # Create edge weights
  edge_weights <- adj_mat*runif(N*N, min=0.5, max=1)
  
  # Create node attributes
  node_attrs <- data.frame(x=nodes[,1], y=nodes[,2])
  
  # Compute Laplacian matrix
  degrees <- rowSums(adj_mat)
  D <- diag(degrees)
  L <- D - adj_mat
  evals <- eigen(L, only.values = TRUE)$values
  gap <- evals[2] - evals[1]
  
  # Return graph object
  graph <- list(adj_mat=adj_mat, edge_weights=edge_weights, node_attrs=node_attrs, laplacian=L
                

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
grow_spatial_graph <- function(beta, N, r, alpha, p, m, weight_func=NULL, node_attr_func=NULL) {
  # Generate random spatial points on a unit square
  points <- matrix(runif(N*2), ncol=2)
  
  # Compute distance matrix
  dist_mat <- as.matrix(dist(points))
  
  # Initialize the network with one node
  adjacency_matrix <- matrix(0, nrow=N, ncol=N)
  adjacency_matrix[1,1] <- 1
  
  # Barabasi-Albert model to add nodes and edges
  for (i in 2:N) {
    # Compute probabilities of attachment
    d <- dist_mat[1:i-1, i]
    deg <- rowSums(adjacency_matrix[1:i-1, 1:i-1])}
    prob <- (d^(-beta))*(deg^alpha)
    prob <- prob/sum(prob)
    
    # Add m edges
    for (j in 1:m) {
      # Choose a node to attach to
      attach_node <- sample(1:i-1, size=1, prob=prob)
      
      # Rewire with probability p
      if (runif(1) < p) {
        # Choose a random node within a certain distance to rewire to
        candidates <- which(dist_mat[i,] <= r & adjacency_matrix[i,] == 0)
        rewire_node <- sample(candidates, size=1)
        adjacency_matrix[i, rewire_node] <- 1
        adjacency_matrix[rewire_node, i] <- 1
        adjacency_matrix[i, attach_node] <- 0
        adjacency_matrix[attach_node, i] <- 0
      } else {
        adjacency_matrix[i, attach_node] <- 1
        adjacency_matrix[attach_node, i] <- 1
      }
    }
    adjacency_matrix[i,i] <- 1
  }
  
  # Add edge weights if specified
  if (!is.null(weight_func)) {
    weights <- matrix(0, nrow=N, ncol=N)
    for (i in 1:N) {
      for (j in 1:N) {
        if (adjacency_matrix[i,j] == 1) {
          weights[i,j] <- weight_func(dist_mat[i,j])
        }
      }
    }
    return(graph_from_adjacency_matrix(adjacency_matrix, weighted=TRUE, diag=FALSE, edge.attr=list(weight=weights)))
  } else {
    g <- graph_from_adjacency_matrix(adjacency_matrix, weighted=FALSE, diag=FALSE)
  }
  
  # Add node attributes if specified
  if (!is.null(node_attr_func)) {
    node_attrs <- data.frame(matrix(0, nrow=N, ncol=length(node_attr_func(1))))
    for (i in 1:N) {
      node_attrs[i,] <- node_attr_func(points[i,])
    }
    set_vertex_attr(g, name="attributes", value=node_attrs)
  }
  
  # Compute Laplacian matrix and check expansion properties
  L <- laplacian_matrix(g, normalized=FALSE)
  eigenvals <- eigen(L, symmetric=TRUE)$values
  gap <- eigenvals[2]
  if (gap > 0.1) {
    warning("Large eigenvalue gap: ", gap)
  }
  
  return(g)
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Define the quadtree data structure
quadtree <- function(xmin, ymin, xmax, ymax, maxdepth = 5) {
  structure(list(
    xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax,
    depth = 0, count = 0, nodes = NULL, leaf = TRUE,
    maxdepth = maxdepth
  ), class = "quadtree")
}

# Define the function to insert a node into the quadtree
insert <- function(qtree, node) {
  if (node$x < qtree$xmin || node$x > qtree$xmax ||
      node$y < qtree$ymin || node$y > qtree$ymax) {
    # Node is outside the bounds of this quadrant
    return(qtree)
  }
  if (qtree$leaf && qtree$count < 4) {
    # This quadrant is a leaf and has room for another node
    qtree$count <- qtree$count + 1
    qtree$nodes <- c(qtree$nodes, list(node))
    return(qtree)
  }
  if (qtree$leaf && qtree$count == 4) {
    # This quadrant is a leaf and is full, so split it
    qtree$leaf <- FALSE
    for (i in 1:4) {
      # Create subquadrants
      xmid <- (qtree$xmin + qtree$xmax) / 2
      ymid <- (qtree$ymin + qtree$ymax) / 2
      if (i == 1) {
        subqtree <- quadtree(qtree$xmin, qtree$ymin, xmid, ymid, maxdepth = qtree$maxdepth)
      } else if (i == 2) {
        subqtree <- quadtree(xmid, qtree$ymin, qtree$xmax, ymid, maxdepth = qtree$maxdepth)
      } else if (i == 3) {
        subqtree <- quadtree(qtree$xmin, ymid, xmid, qtree$ymax, maxdepth = qtree$maxdepth)
      } else {
        subqtree <- quadtree(xmid, ymid, qtree$xmax, qtree$ymax, maxdepth = qtree$maxdepth)
      }
      qtree$nodes <- c(qtree$nodes, list(subqtree))
    }
    # Move existing nodes into subquadrants
    for (i in 1:4) {
      qtree$nodes[[i]] <- insert(qtree$nodes[[i]], qtree$nodes[[qtree$count-3+i]])
    }
    qtree$count <- 1
    qtree$nodes <- c(qtree$nodes, list(node))
    return(qtree)
  }
  
  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Define a class for a quad tree node
  Node <- setRefClass("Node", fields = list(
    x = "numeric",
    y = "numeric",
    value = "any",
    nw = "Node",
    ne = "Node",
    sw = "Node",
    se = "Node"
  ))
  
  # Define a function to create a new node
  newNode <- function(x, y, value) {
    return(Node(x = x, y = y, value = value, nw = NULL, ne = NULL, sw = NULL, se = NULL))
  }
  
  # Define a class for a quad tree
  QuadTree <- setRefClass("QuadTree", fields = list(
    root = "Node",
    xMin = "numeric",
    yMin = "numeric",
    xMax = "numeric",
    yMax = "numeric",
    maxDepth = "numeric",
    capacity = "numeric"
  ))
  
  # Define a function to create a new quad tree
  newQuadTree <- function(xMin, yMin, xMax, yMax, maxDepth, capacity) {
    return(QuadTree(root = NULL, xMin = xMin, yMin = yMin, xMax = xMax, yMax = yMax, maxDepth = maxDepth, capacity = capacity))
  }
  
  # Define a function to insert a point into the quad tree
  QuadTree$insert <- function(this, x, y, value) {
    if (is.null(this$root)) {
      this$root <- newNode(x, y, value)
      return
    }
    
    node <- this$root
    depth <- 0
    
    while (!is.null(node)) {
      if (is.null(node$nw)) {
        if (is.null(node$value)) {
          node$x <- x
          node$y <- y
          node$value <- value
          return
        } else {
          nw <- newNode(node$x, node$y, node$value)
          ne <- newNode(this$xMin + (this$xMax - this$xMin) / 2, node$y, NULL)
          sw <- newNode(node$x, this$yMin + (this$yMax - this$yMin) / 2, NULL)
          se <- newNode(this$xMin + (this$xMax - this$xMin) / 2, this$yMin + (this$yMax - this$yMin) / 2, NULL)
          node$nw <- nw
          node$ne <- ne
          node$sw <- sw
          node$se <- se
          node$value <- NULL
        }
      }
      
      if (x < node$x) {
        if (y < node$y) {
          node <- node$sw
        } else {
          node <- node$nw
        }
      } else {
        if (y < node$y) {
          node <- node$se
        } else {
          node <- node$ne
        }
      }
      
      depth <- depth + 1
      
      if (depth > this$maxDepth) {
        return
      }
    }
  }
  
  # Define a function to query the quad tree
  QuadTree$query <- function(this, xMin, yMin, xMax, yMax) {
    results <- list()
    stack <- list(this$root)
    
    while (length(stack) > 0) {
      node <- stack[[length(stack)]]
      stack <- stack[1:length(stack) - 1]
      
      if (!is.null(node)) {
        if (node$x >= xMin && node$x <= xMax && node$y >= yMin && node$y
            

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Define the quadtree data structure
quadtree <- function(xmin, ymin, xmax, ymax, maxdepth = 5) {
  structure(list(
    xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax,
    depth = 0, count = 0, nodes = NULL, leaf = TRUE,
    maxdepth = maxdepth
  ), class = "quadtree")
}

# Define the function to insert a node into the quadtree
insert <- function(qtree, node) {
  if (node$x < qtree$xmin || node$x > qtree$xmax ||
      node$y < qtree$ymin || node$y > qtree$ymax) {
    # Node is outside the bounds of this quadrant
    return(qtree)
  }
  if (qtree$leaf && qtree$count < 4) {
    # This quadrant is a leaf and has room for another node
    qtree$count <- qtree$count + 1
    qtree$nodes <- c(qtree$nodes, list(node))
    return(qtree)
  }
  if (qtree$leaf && qtree$count == 4) {
    # This quadrant is a leaf and is full, so split it
    qtree$leaf <- FALSE
    for (i in 1:4) {
      # Create subquadrants
      xmid <- (qtree$xmin + qtree$xmax) / 2
      ymid <- (qtree$ymin + qtree$ymax) / 2
      if (i == 1) {
        subqtree <- quadtree(qtree$xmin, qtree$ymin, xmid, ymid, maxdepth = qtree$maxdepth)
      } else if (i == 2) {
        subqtree <- quadtree(xmid, qtree$ymin, qtree$xmax, ymid, maxdepth = qtree$maxdepth)
      } else if (i == 3) {
        subqtree <- quadtree(qtree$xmin, ymid, xmid, qtree$ymax, maxdepth = qtree$maxdepth)
      } else {
        subqtree <- quadtree(xmid, ymid, qtree$xmax, qtree$ymax, maxdepth = qtree$maxdepth)
      }
      qtree$nodes <- c(qtree$nodes, list(subqtree))
    }
    # Move existing nodes into subquadrants
    for (i in 1:4) {
      qtree$nodes[[i]] <- insert(qtree$nodes[[i]], qtree$nodes[[qtree$count-3+i]])
    }
    qtree$count <- 1
    qtree$nodes <- c(qtree$nodes, list(node))
    return(qtree)
  }

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Function to generate a spatial scale-free expander graph
  # Arguments:
  #   - n: number of nodes in the graph
  #   - r: cutoff distance for adjacency matrix based on distance matrix
  #   - beta: power law parameter for scale-free preferential attachment
  #   - m: number of edges added at each step for preferential attachment
  #   - alpha: strength of degree effect in preferential attachment
  #   - mu: probability of attachment within or between communities
  #   - prewire: probability of rewiring edges for small-world effect
  #   - weight_fun: function to generate edge weights based on distance or other factors
  #   - attr_fun: function to generate node attributes
  # Output: 
  #   - igraph object representing the spatial scale-free expander graph
  generate_spatial_scalefree_graph <- function(n, r, beta, m, alpha, mu, prewire, weight_fun = NULL, attr_fun = NULL) {
    
    # Generate initial placement of nodes using quadtree algorithm
    quadtree <- make_quadtree(cbind(runif(n), runif(n)))
    
    # Initialize adjacency matrix with all zeros
    adj_mat <- matrix(0, nrow = n, ncol = n)
    
    # Function to compute probability of attachment for a new node
    # based on distance, degree, and community structure
    prob_attach <- function(node, neighbors, degree_vec) {
      dist_vec <- sqrt(rowSums((node - neighbors)^2))
      dist_factor <- exp(-dist_vec/r)
      degree_factor <- degree_vec[neighbors]^alpha
      comm_factor <- ifelse(degree_vec[node] < m, 1, mu)
      return(dist_factor * degree_factor * comm_factor)
    }
    
    # Initialize graph with first node and attributes
    graph <- make_empty_graph(1)
    if (!is.null(attr_fun)) {
      graph <- set_vertex_attr(graph, "attr", value = attr_fun(1))
    }
    
    # Add nodes and edges one at a time
    for (i in 2:n) {
      
      # Compute probability of attachment for each existing node
      prob_vec <- prob_attach(cbind(runif(1), runif(1)), get_quadtree_points(quadtree), degree(graph))
      
      # Sample m nodes with probability proportional to attachment probability
      new_edges <- sample(1:(i-1), size = m, replace = TRUE, prob = prob_vec)
      
      # Add edges to new node and existing nodes
      edges <- cbind(rep(i, m), new_edges)
      if (!is.null(weight_fun)) {
        edge_weights <- weight_fun(get_vertex_attr(graph, "attr")[edges[,1],], get_vertex_attr(graph, "attr")[edges[,2],])
        graph <- add_edges(graph, edges, weight = edge_weights)
      } else {
        graph <- add_edges(graph, edges)
      }
      
      # Add node to quadtree
      quadtree <- insert_into_quadtree(quadtree, cbind(runif(1), runif(1)))
      
      # Add node and attributes to graph
      graph <- add_vertices(graph, 1)
      if (!is.null(attr_fun)) {
        graph <- set_vertex_attr(graph, "attr", value = rbind(get_vertex_attr(graph, "attr"), attr_fun(i)))
      }
      
      # Rewire edges with small-world probability
      if (prewire > 0) {
        edges <- get_edgelist(graph)
        rewire_mask <- runif(nrow(edges)) < prewire
        
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # function to generate a spatial scale-free expander graph
        spatial_sf_expander_graph <- function(N, beta, alpha, mu, prewire, m, r, plot_graph = TRUE) {
          
          # generate initial nodes
          x <- runif(N)
          y <- runif(N)
          nodes <- data.frame(x = x, y = y)
          
          # generate adjacency matrix based on distance matrix
          dist_mat <- as.matrix(dist(nodes))
          adj_mat <- as.matrix(dist_mat <= r)
          
          # initialize the degree vector
          degrees <- rep(0, N)
          
          # loop over new nodes
          for (i in 2:N) {
            
            # calculate distances to existing nodes
            distances <- dist(nodes[i,], nodes) 
            
            # calculate probabilities of attachment
            prob_attach <- (distances^(-beta) + alpha*degrees + mu*adj_mat[i,])/sum(distances^(-beta) + alpha*degrees + mu*adj_mat[i,])
            
            # add new edges based on preferential attachment
            for (j in 1:m) {
              new_edge <- sample(1:N, size = 1, prob = prob_attach)
              adj_mat[i, new_edge] <- 1
              adj_mat[new_edge, i] <- 1
              degrees[i] <- degrees[i] + 1
              degrees[new_edge] <- degrees[new_edge] + 1
            }
            
            # rewire edges with small-world probability
            for (j in 1:N) {
              if (runif(1) < prewire) {
                if (adj_mat[i,j] == 1) {
                  new_edge <- sample(1:N, size = 1)
                  adj_mat[i,j] <- 0
                  adj_mat[j,i] <- 0
                  adj_mat[i,new_edge] <- 1
                  adj_mat[new_edge,i] <- 1
                }
              }
            }
          }
          
          # create the graph object
          g <- graph_from_adjacency_matrix(adj_mat)
          
          # calculate Laplacian matrix
          laplacian_mat <- laplacian_matrix(g)
          
          # plot the graph
          if (plot_graph) {
            plot(g, layout = layout.fruchterman.reingold)
          }
          
          # return the graph and Laplacian matrix
          list(graph = g, laplacian_matrix = laplacian_mat)
        }
        
        # example usage
        set.seed(123)
        spatial_sf_expander_graph(N = 100, beta = 1, alpha = 0.1, mu = 0.1, prewire = 0.05, m = 2, r = 0.1)
        

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


d <- 2 # Set dimensionality of torus
quadtree <- function(xlim, ylim, depth, min_size) {
  if (depth == 0) return(NULL)
  if (xlim[2] - xlim[1] <= min_size && ylim[2] - ylim[1] <= min_size) {
    # Create new node at center of quadtree
    return(list(x = (xlim[1] + xlim[2])/2, y = (ylim[1] + ylim[2])/2))
  }
  # Recursively subdivide quadtree
  xmid <- (xlim[1] + xlim[2])/2
  ymid <- (ylim[1] + ylim[2])/2
  return(c(quadtree(c(xlim[1], xmid), c(ymid, ylim[2]), depth-1, min_size),
           quadtree(c(xlim[1], xmid), c(ylim[1], ymid), depth-1, min_size),
           quadtree(c(xmid, xlim[2]), c(ymid, ylim[2]), depth-1, min_size),
           quadtree(c(xmid, xlim[2]), c(ylim[1], ymid), depth-1, min_size)))
}
min_size <- 1/(2*N) # Set minimum separation distance between nodes
nodes <- quadtree(c(0,1), c(0,1), L, min_size)

# Step 2: Calculate distance matrix and create adjacency matrix
dist_mat <- as.matrix(dist(t(as.data.frame(nodes))))
adj_mat <- as.matrix(dist_mat <= r)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Define the quadtree generation function
generate_quadtree <- function(min_coords, max_coords, min_dist, max_dist, depth) {
  # Initialize the quadtree node
  node <- list()
  node$min_coords <- min_coords
  node$max_coords <- max_coords
  
  # Determine if the node is a leaf or not
  if (depth == 0) {
    # Generate a random point within the node's boundaries
    point <- rep(0, length(min_coords))
    for (i in 1:length(min_coords)) {
      point[i] <- runif(1, min_coords[i], max_coords[i])
    }
    node$point <- point
  } else {
    # Generate four child nodes and recursively generate their children
    node$children <- list()
    mid_coords <- (min_coords + max_coords) / 2
    for (i in 0:1) {
      for (j in 0:1) {
        child_min_coords <- mid_coords + c(i, j) * (min_dist / 2)
        child_max_coords <- mid_coords + c(i + 1, j + 1) * (min_dist / 2)
        if (max_dist > min_dist && depth > 1) {
          child_max_dist <- max_dist / 2
        } else {
          child_max_dist <- max_dist
        }
        child_node <- generate_quadtree(child_min_coords, child_max_coords, min_dist, child_max_dist, depth - 1)
        node$children[[i + 2 * j + 1]] <- child_node
      }
    }
  }
  
  return(node)
}

# Define a function to generate points on a d-dimensional torus
generate_torus_points <- function(n, d, min_dist, max_dist) {
  # Generate the quadtree
  min_coords <- rep(0, d)
  max_coords <- rep(1, d)
  quadtree <- generate_quadtree(min_coords, max_coords, min_dist, max_dist, log2(n))
  
  # Traverse the quadtree and extract the leaf nodes
  leaf_nodes <- list()
  stack <- list(quadtree)
  while (length(stack) > 0) {
    node <- stack[[length(stack)]]
    stack <- stack[-length(stack)]
    if (length(node$point) > 0) {
      leaf_nodes[[length(leaf_nodes) + 1]] <- node
    } else {
      stack <- c(stack, node$children)
    }
  }
  
  # Extract the points from the leaf nodes and apply torus periodicity
  points <- matrix(0, nrow=n, ncol=d)
  for (i in 1:n) {
    node <- leaf_nodes[[i]]
    point <- node$point
    for (j in 1:d) {
      if (point[j] < 0) {
        point[j] <- point[j] + 1
      } else if (point[j] > 1) {
        point[j] <- point[j] - 1
      }
    }
    points[i,] <- point
  }
  
  return(points)
}

# Generate 100 points on a 2-dimensional torus with minimum distance 0.1 and maximum distance 0.5
points <- generate_torus_points(100, 2, 0.1, 0.5)
plot(points, pch=16)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# To generate nodes on a d-dimensional torus with an underlying spatial structure/information using a quadtree generation algorithm that has a constraint not to allow nodes to be too near or too far from each other in R, you can write your own code using the following steps:
#   
#   Define the dimensions of your torus, the number of nodes you want to generate, and the minimum and maximum distance between nodes.
# Create a function that generates a random point in the torus based on the dimensions provided. This can be done by generating a random coordinate for each dimension, and then wrapping the coordinates around the torus if they go beyond the boundaries.
# Implement a quadtree algorithm to partition the torus into smaller regions. This involves creating a root node that covers the entire torus, and recursively splitting it into four child nodes until each node contains a maximum number of points or is smaller than a minimum size. You can use the rpart package in R to implement the quadtree algorithm.
# Add a constraint to the quadtree algorithm to check the minimum and maximum distance between nodes in each node. If the distance is too small, generate a new point in that node. If the distance is too large, split the node into four child nodes and repeat the process.
# Repeat steps 3 and 4 until the desired number of nodes is generated.
# Here's some example code that implements these steps:
# 
# R
# Copy code
library(rpart)

# Define the dimensions of the torus
dims <- c(10, 10, 10)  # 3-dimensional torus
volume <- prod(dims)   # Total volume of the torus

# Define the number of nodes and minimum/maximum distance between nodes
n_nodes <- 100
min_dist <- 1
max_dist <- 3

# Generate a random point in the torus
generate_point <- function() {
  coords <- runif(length(dims), min = 0, max = 1)
  coords <- coords * dims
  coords <- coords %% dims  # Wrap coordinates around torus
  return(coords)
}

# Define a function to check the distance between two points
distance <- function(p1, p2) {
  diff <- abs(p1 - p2)
  diff <- pmin(diff, dims - diff)  # Handle torus wrapping
  return(sqrt(sum(diff^2)))
}

# Define a function to generate nodes using quadtree algorithm with distance constraint
generate_nodes <- function(n, min_dist, max_dist) {
  tree <- rpart(target ~ ., data = data.frame(coords = I(list()), target = logical()), control = rpart.control(minsplit = 1))
  while (nrow(tree$frame) < n) {
    # Split nodes that contain too many points
    too_many_points <- which(tree$frame$var == "")
    for (i in too_many_points) {
      coords <- tree$frame$coords[i]
      indices <- tree$where(tree$frame$data[["coords"]] == coords, use.nas = FALSE)
      if (nrow(indices) > 1) {
        tree <- prune(tree, prune.tree(tree, i))
        tree <- rpart(target ~ ., data = data.frame(coords = I(list()), target = logical()), control = rpart.control(minsplit = 1))
        break
      }
    }
    # Generate new points in nodes with insufficient distance between points
    for (i in seq_along(tree$frame$var)) {
      coords <- tree$frame$coords[i]
      points <- tree$frame$data[["coords"]][tree$where(tree$frame$data[["coords"]] == coords, use.nas = FALSE),]
      if (

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        generate_quadtree_torus <- function(d, n, min_dist, max_dist) {
          # Generate a set of random points in the d-dimensional torus
          points <- matrix(runif(n*d, min = 0, max = 1), ncol = d)
          for (i in 1:d) {
            points[,i] <- points[,i] * 2 * pi
          }
          cos_points <- cos(points)
          sin_points <- sin(points)
          radius <- 1 / sqrt(d)
          for (i in 1:n) {
            r <- radius * sqrt(runif())
            theta <- runif() * 2 * pi
            for (j in 1:d) {
              points[i,j] <- points[i,j] + r * cos(theta)
              if (points[i,j] < 0) {
                points[i,j] <- points[i,j] + 2 * pi
              }
              if (points[i,j] > 2 * pi) {
                points[i,j] <- points[i,j] - 2 * pi
              }
            }
          }
          
          # Build the quadtree using a recursive function
          quadtree <- build_quadtree(points, min_dist, max_dist)
          return(quadtree)
        }
        
        build_quadtree <- function(points, min_dist, max_dist) {
          n <- nrow(points)
          if (n == 0) {
            return(NULL)
          }
          if (n == 1) {
            return(list(point = points[1,], children = NULL))
          }
          center <- apply(points, 2, mean)
          ranges <- apply(points, 2, range)
          if (max(abs(ranges[2,] - ranges[1,])) <= min_dist) {
            return(list(point = center, children = NULL))
          }
          children <- list(NULL, NULL, NULL, NULL)
          for (i in 1:n) {
            child_points <- points[-i,]
            dists <- apply(child_points, 1, function(x) sqrt(sum((x - points[i,])^2)))
            close_points <- child_points[dists <= max_dist & dists >= min_dist,]
            for (j in 1:4) {
              quadrant <- get_quadrant(center, points[i,], j)
              quad_points <- close_points[rowSums(get_quadrant_ranges(quadrant, center, ranges)) == d,]
              if (nrow(quad_points) > 0) {
                children[[j]] <- c(children[[j]], build_quadtree(quad_points, min_dist, max_dist))
              }
            }
          }
          return(list(point = center, children = children))
        }
        
        get_quadrant <- function(center, point, i) {
          d <- length(center)
          center_mod <- center + pi
          point_mod <- point + pi
          center_mod[(i-1) %% d + 1] <- center[(i-1) %% d + 1]
          point_mod[(i-1) %% d + 1] <- center[(i-1) %% d + 1]
          return(as.logical(point_mod > center_mod))
        }
        
        get_quadrant_ranges <- function(quadrant, center, ranges) {
          d <- length(center)
          ranges_mod <- matrix(0, nrow = d, ncol = 2)
          ranges_mod[quadrant] <- ranges[quadrant,]
          ranges
          

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          # Function to generate a spatial scale-free expander graph
          # Inputs:
          # - N: number of nodes
          # - L: length of the unit square
          # - r: cutoff distance for adjacency matrix
          # - beta: power law parameter for preferential attachment
          # - alpha: parameter for the strength of degree effect
          # - sigma: parameter for the decay of degree effect with distance
          # - delta: probability of attachment within communities
          # - m: number of edges added at each attachment
          # - p: rewiring probability for small world effect
          # - node_attrs: a list of node attributes (optional)
          # - edge_weights: a matrix of edge weights (optional)
          # Output:
          # - An adjacency matrix of the generated graph
          spatial_scale_free_expander <- function(N, L, r, beta, alpha, sigma, delta, m, p, node_attrs = NULL, edge_weights = NULL) {
            # Generate the initial nodes using a quadtree algorithm
            quadtree <- quadtree_gen(N, L)
            coords <- quadtree_to_coords(quadtree)
            
            # Initialize the adjacency matrix
            adj_matrix <- matrix(0, N, N)
            
            # Add edges for each new node
            for (i in 2:N) {
              # Compute the probability of attachment to existing nodes
              dists <- dist(coords[1:i-1,])
              degrees <- rowSums(adj_matrix[1:i-1,1:i-1])
              degrees[degrees == 0] <- 1
              degree_weights <- (degrees^(-beta)) * (dists^(-alpha))
              comm_weights <- matrix(0, i-1, i-1)
              for (j in 1:(i-1)) {
                for (k in 1:(i-1)) {
                  if (quadtree[j] == quadtree[k]) {
                    comm_weights[j,k] <- delta
                  }
                }
              }
              attachment_weights <- degree_weights * comm_weights
              attachment_probs <- attachment_weights / sum(attachment_weights)
              
              # Add m edges to existing nodes with probability proportional to their attachment probabilities
              for (j in 1:m) {
                if (j == 1) {
                  # Add one edge using linear preferential attachment
                  idx <- sample(1:(i-1), size = 1, prob = attachment_probs)
                } else {
                  # Add multiple edges using regular preferential attachment
                  idx <- sample(1:(i-1), size = 1, prob = degrees[1:i-1])
                }
                adj_matrix[i, idx] <- 1
                adj_matrix[idx, i] <- 1
              }
              
              # Rewire edges with small world probability p
              for (j in 1:i) {
                for (k in (j+1):N) {
                  if (adj_matrix[j,k] == 1 && runif(1) < p) {
                    # Rewire edge (j,k) to a random node within distance r
                    dist_jk <- dists[j,k]
                    candidate_idxs <- which(dists[j,] <= r)
                    candidate_idxs <- candidate_idxs[candidate_idxs != j & candidate_idxs != k]
                    if (length(candidate_idxs) > 0) {
                      candidate_dists <- dists[k,candidate_idxs]
                      candidate_weights <- (candidate_dists^(-beta)) * (degrees[candidate_idxs]^(-alpha))
                      candidate_probs <- candidate_weights / sum(candidate_weights)
                      new_k <- sample(candidate_idxs, size = 1
                                      

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Function to grow a spatial scale-free expander graph
# Inputs:
#   - N: number of nodes to generate
#   - beta: power law parameter for degree distribution
#   - sigma: decay parameter for degree effect in attachment probability
#   - alpha: strength parameter for degree effect in attachment probability
#   - mu: probability of attachment to nodes in same community
#   - L: number of communities in the graph
#   - r: cutoff distance for spatial distance effect
#   - m: number of edges to add at each instance
#   - p: rewiring probability for small world effect
#   - lambda: distance parameter for distance-dependent attachment probability
# Output:
#   - adjacency matrix of the resulting graph
grow_spatial_scale_free_graph <- function(N, beta, sigma, alpha, mu, L, r, m, p, lambda) {
  # Generate nodes spatially using quad tree
  library(SpatialPosition)
  quad_tree <- spindex(createSP(n = N, type = "random", extent = c(0, 1, 0, 1)))
  nodes <- quad_tree$pts
  # Calculate distance matrix
  dist_matrix <- as.matrix(dist(nodes))
  # Create adjacency matrix based on distance matrix and cutoff distance r
  adj_matrix <- ifelse(dist_matrix <= r, 1, 0)
  diag(adj_matrix) <- 0 # remove self-loops
  # Initialize the graph with L communities
  community <- rep(1:L, each = N/L)
  # Grow the graph one node at a time
  for (i in 2:N) {
    # Compute attachment probabilities based on distance, degree, and community
    prob <- numeric(i-1)
    for (j in 1:(i-1)) {
      d_ij <- dist(nodes[i,], nodes[j,])
      if (d_ij <= r) {
        prob[j] <- (1 / (d_ij^beta)) * (adj_matrix[i,j] + 1)^(alpha-sigma*d_ij) * (mu^(community[i] == community[j]))
      } else {
        prob[j] <- 0
      }
    }
    # Add m edges based on preferential attachment
    for (k in 1:m) {
      # Choose a node to attach to based on the attachment probabilities
      j <- sample(1:(i-1), size = 1, prob = prob/sum(prob))
      # Rewire edge with probability p
      if (runif(1) < p) {
        # Choose a random node within a certain distance to rewire to
        candidates <- which(dist_matrix[i,] > r & dist_matrix[i,] <= lambda)
        if (length(candidates) > 0) {
          j <- sample(candidates, size = 1)
        }
      }
      # Add edge to the chosen node
      adj_matrix[i,j] <- 1
      adj_matrix[j,i] <- 1
    }
    # Assign node to a community based on the community of its neighbors
    community[i] <- sample(community[adj_matrix[i,]==1], size=1)
  }
  # Return adjacency matrix
  return(adj_matrix)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Function to generate a spatial scale-free expander graph
# with the given parameters
# Arguments:
#   N: number of nodes in the graph
#   beta: power-law exponent for preferential attachment
#   alpha: parameter that controls the strength of degree effect
#   delta: probability of attachment to nodes within/outside communities
#   r: cutoff distance for spatial effect
#   m: number of edges to add for each new node
#   p: probability of rewiring an edge to a random node
# Return:
#   igraph graph object

generate_spatial_scalefree_graph <- function(N, beta, alpha, delta, r, m, p) {
  # Generate initial network with two nodes connected by an edge
  g <- graph.empty(n = 2, directed = FALSE)
  E(g) <- rbind(c(1, 2), c(2, 1))
  V(g)$community <- c(1, 2)
  V(g)$degree <- c(1, 1)
  pos <- matrix(c(0, 0, 1, 1), ncol = 2)
  
  # Generate remaining nodes
  for (i in 3:N) {
    # Add new node
    g <- add.vertices(g, 1)
    V(g)$community[i] <- ifelse(runif(1) < delta, sample(unique(V(g)$community), 1), max(V(g)$community) + 1)
    V(g)$degree[i] <- 0
    
    # Calculate probabilities for attaching to existing nodes
    d <- as.matrix(dist(pos[1:i-1,]))
    p_dist <- 1 / (d^beta + 1e-6) # add small number to avoid division by zero
    p_deg <- V(g)$degree[1:i-1]^alpha
    p_comm <- ifelse(V(g)$community[1:i-1] == V(g)$community[i], delta, 1-delta)
    p <- p_dist * p_deg * p_comm
    p <- p / sum(p)
    
    # Add edges
    for (j in 1:m) {
      # Choose target node
      target <- sample(1:(i-1), 1, prob = p)
      
      # Choose distance for new edge
      d_new <- runif(1) * r
      
      

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      generate_graph <- function(N, beta, alpha, delta, p, m) {
        
        # Step 1: Generate nodes on a unit square and distribute them spatially using quadtree
        library(quadprog)
        
        quadtree <- function(xmin, ymin, xmax, ymax, max_points = 1) {
          list(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax, points = NULL, 
               TL = NULL, TR = NULL, BL = NULL, BR = NULL)
        }
        
        subdivide <- function(node) {
          mx <- (node$xmin + node$xmax) / 2
          my <- (node$ymin + node$ymax) / 2
          node$TL <- quadtree(node$xmin, node$ymin, mx, my, max_points = max_points)
          node$TR <- quadtree(mx, node$ymin, node$xmax, my, max_points = max_points)
          node$BL <- quadtree(node$xmin, my, mx, node$ymax, max_points = max_points)
          node$BR <- quadtree(mx, my, node$xmax, node$ymax, max_points = max_points)
          for (p in node$points) insert_point(node, p)
          node$points <- NULL
        }
        
        insert_point <- function(node, point) {
          if (is.null(node$TL)) {
            node$points <- rbind(node$points, point)
            if (nrow(node$points) > max_points) subdivide(node)
          } else {
            if (point$x < (node$xmin + node$xmax) / 2) {
              if (point$y < (node$ymin + node$ymax) / 2) insert_point(node$TL, point)
              else insert_point(node$BL, point)
            } else {
              if (point$y < (node$ymin + node$ymax) / 2) insert_point(node$TR, point)
              else insert_point(node$BR, point)
            }
          }
        }
        
        max_points <- 10
        root <- quadtree(0, 0, 1, 1, max_points = max_points)
        
        for (i in 1:N) {
          x <- runif(1)
          y <- runif(1)
          point <- c(x, y)
          insert_point(root, point)
          
          
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          generate_spatial_scalefree_expander <- function(N, beta, alpha, delta, p, m, r, max_dist) {
            
            library(data.table)
            library(Matrix)
            library(ggplot2)
            library(RANN)
            
            # Function to calculate probability of attachment to existing nodes
            attachment_prob <- function(dists, degrees, beta) {
              return ((dists ^ (-beta)) * degrees)
            }
            
            # Create initial node set
            nodes <- data.table(x = runif(N), y = runif(N), community = 1)
            
            # Compute distance matrix and minimum distance
            dist_mat <- RANN::nn2(cbind(nodes$x, nodes$y), k = 2, searchtype = "default")$nn.dist
            min_dist <- min(dist_mat[dist_mat > 0])
            
            # Create adjacency matrix based on distance matrix
            adj_mat <- ifelse(dist_mat <= r, 1, 0)
            
            # Compute degree vector
            degree <- rowSums(adj_mat)
            
            # Initialize graph object
            graph_mat <- Matrix(adj_mat, sparse = TRUE)
            
            # Loop over remaining nodes
            for (i in seq(N)[-1]) {
              
              # Compute attachment probability for each existing node
              prob_vec <- attachment_prob(dist_mat[i, ], degree, beta)
              
              # Normalize probabilities
              prob_vec <- prob_vec / sum(prob_vec)
              
              # Select m nodes to attach to
              attach_nodes <- sample.int(N - 1, size = m, replace = TRUE, prob = prob_vec)
              
              # Add edges to selected nodes
              graph_mat[i, attach_nodes] <- 1
              graph_mat[attach_nodes, i] <- 1
              
              # Update degree vector
              degree[i] <- m
              degree[attach_nodes] <- degree[attach_nodes] + 1
              
              # Assign node to community
              community_probs <- degree[attach_nodes] ^ alpha
              community_probs <- community_probs / sum(community_probs)
              nodes$community[i] <- sample.int(delta, size = 1, replace = TRUE, prob = community_probs)
              
              # Rewire edges with probability p
              if (p > 0) {
                
                # Identify edges to potentially rewire
                potential_edges <- which(graph_mat[i, ] == 1 & degree > 1)
                
                # Loop over potential edges
                for (j in potential_edges) {
                  
                  # Compute distance to other node
                  dist_ij <- sqrt((nodes$x[i] - nodes$x[j]) ^ 2 + (nodes$y[i] - nodes$y[j]) ^ 2)
                  
                  # Randomly rewire with probability p
                  if (runif(1) < p && dist_ij <= max_dist) {
                    
                    # Select new node to attach to
                    rewire_node <- sample.int(N, size = 1, replace = TRUE)
                    
                    # Check if edge already exists or if node is too far away
                    if (graph_mat[i, rewire_node] == 0 && dist_mat[i, rewire_node] <= r) {
                      
                      # Rewire edge
                      graph_mat[i, j] <- 0
                      graph_mat[j, i] <- 0
                      graph_mat[i, rewire_node] <- 1
                      graph_mat[rewire_node, i] <- 1
                      
                      # Update degree vector
                      degree
                      
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                      # Required packages
                      library(Matrix)
                      library(quadprog)
                      library(RANN)
                      library(ggplot2)
                      
                      # Function to generate spatially distributed points on a unit square
                      generate_points <- function(n) {
                        # Generate random x,y coordinates
                        x <- runif(n)
                        y <- runif(n)
                        # Return a matrix with x and y coordinates
                        return(cbind(x, y))
                      }
                      
                      # Function to calculate the minimum distance between pairs of nodes
                      min_distance <- function(points) {
                        # Calculate distance matrix using RANN package
                        dist_matrix <- RANN::nn2(points)$nn.dist
                        # Set diagonal elements to a large value so they won't be chosen as minimum
                        diag(dist_matrix) <- Inf
                        # Return the minimum distance
                        return(min(dist_matrix))
                      }
                      
                      # Function to create the initial adjacency matrix
                      create_adj_matrix <- function(points, r) {
                        # Calculate distance matrix
                        dist_matrix <- as.matrix(dist(points))
                        # Set adjacency matrix element to 1 if dist_{ij} <= 'r'
                        adj_matrix <- ifelse(dist_matrix <= r, 1, 0)
                        # Return adjacency matrix
                        return(adj_matrix)
                      }
                      
                      # Function to calculate the degree of each node in the graph
                      get_degree <- function(adj_matrix) {
                        # Calculate row sums of adjacency matrix
                        degree <- rowSums(adj_matrix)
                        # Return degree vector
                        return(degree)
                      }
                      
                      # Function to calculate the probability of attaching to each existing node
                      get_probabilities <- function(points, adj_matrix, degree, beta, alpha, delta) {
                        # Number of nodes in the graph
                        n <- nrow(points)
                        # Calculate spatial distance between new node and existing nodes
                        dist_spatial <- as.matrix(dist(points[(n+1),], points))
                        # Calculate community structure probability based on delta parameter
                        comm_prob <- 1 + (1-delta)*abs(1:n - (n+1))
                        comm_prob <- comm_prob/sum(comm_prob)
                        # Calculate degree probability based on beta and alpha parameters
                        degree_prob <- (degree+alpha)^beta
                        degree_prob <- degree_prob/sum(degree_prob)
                        # Calculate total probability of attachment
                        total_prob <- comm_prob*degree_prob*1/dist_spatial
                        # Set probabilities of nodes that are already connected to 0
                        total_prob[which(adj_matrix[(n+1),] == 1)] <- 0
                        # Return total probability vector
                        return(total_prob)
                      }
                      
                      # Function to add a new node to the graph
                      add_node <- function(points, adj_matrix, degree, beta, alpha, delta) {
                        # Number of nodes in the graph
                        n <- nrow(points)
                        # Calculate probability of attaching to each existing node
                        probabilities <- get_probabilities(points, adj_matrix, degree, beta, alpha, delta)
                        # Choose an existing node to attach to based on probability
                        chosen_node <- sample(1:n, size=1, prob=probabilities)
                        # Add new row and column to adjacency matrix
                        adj_matrix[(n+1),] <- 0
                        adj_matrix[, (n+1)] <- 0
                        # Connect new node to chosen node
                        adj_matrix[(n+1), chosen_node] <- 1
                        adj_matrix[chosen_node, (n+1)] <- 1
                        # Return updated adjacency matrix
                        return(adj_matrix)
                      }
                      
                      # Function to randomly rewire edges
                      

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                      # Function to generate a spatial scale-free expander graph
                      # with specified parameters
                      #
                      # Arguments:
                      # - N: number of nodes to generate
                      # - beta: power law exponent for preferential attachment
                      # - alpha: strength of degree effect
                      # - delta: probability of attachment to nodes in the same community
                      # - p: rewiring probability for small world effect
                      # - r: distance cutoff for spatial distance effect
                      #
                      # Returns:
                      # - a list with elements:
                      #   - 'graph': the adjacency matrix of the generated graph
                      #   - 'pos': the positions of the nodes in the unit square
                      #   - 'communities': a vector indicating the community of each node
                      #   - 'laplacian': the Laplacian matrix of the graph
                      #
                      generate_spatial_scalefree_expander <- function(N, beta, alpha, delta, p, r) {
                        # Generate initial random graph with two nodes and one edge
                        graph <- matrix(0, nrow = 2, ncol = 2)
                        graph[1,2] <- 1
                        
                        # Initialize node positions and communities
                        pos <- matrix(runif(N*2), nrow = N, ncol = 2)
                        qt <- create_quadtree(pos)
                        communities <- rep(1, N)
                        
                        # Grow graph by adding one node at a time
                        for (i in 3:N) {
                          # Compute distances to existing nodes and find nearest neighbor
                          dists <- dist(pos[1:(i-1),])
                          nn <- which.min(dists[i-1,])
                          nn_dist <- dists[i-1, nn]
                          
                          # Compute probability of attachment to each node
                          prob <- rep(0, i-1)
                          for (j in 1:(i-1)) {
                            if (j == nn) {
                              # Preferential attachment to nearest neighbor
                              prob[j] <- (alpha + degree(graph, j))^beta
                            } else if (communities[i] == communities[j]) {
                              # Preferential attachment to nodes in the same community
                              prob[j] <- (delta + degree(graph, j))^beta
                            } else {
                              # Uniform attachment to other nodes
                              prob[j] <- degree(graph, j)^beta
                            }
                          }
                          prob <- prob / sum(prob)
                          
                          # Sample attachment node based on probabilities
                          attach_node <- sample(1:(i-1), size = 1, prob = prob)
                          
                          # Add edge to attachment node
                          graph[i, attach_node] <- 1
                          graph[attach_node, i] <- 1
                          
                          # Rewire edge with probability p
                          if (runif(1) < p) {
                            # Find candidate rewiring nodes within distance r
                            candidates <- which(dists[i-1,] > nn_dist & dists[i-1,] <= r)
                            if (length(candidates) > 0) {
                              # Choose random candidate and rewire edge
                              rewire_node <- sample(candidates, size = 1)
                              graph[i, attach_node] <- 0
                              graph[attach_node, i] <- 0
                              graph[i, rewire_node] <- 1
                              graph[rewire_node, i] <- 1
                            }
                          }
                          
                          # Update node position and community
                          pos[i,] <- sample_quadtree(qt, pos[i-1,], r/2)
                          communities[i] <- communities[attach_node]
                        }
                        
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        # Function to grow a spatial scale-free expander graph
                        # with given properties
                        # Args:
                        #   - N: Number of nodes in the graph
                        #   - r: Cutoff distance for spatial adjacency matrix
                        #   - beta: Power law parameter for scale-free preferential attachment
                        #   - alpha: Parameter controlling strength of degree effect
                        #   - delta: Probability of attachment to nodes within same community
                        #   - p: Rewiring probability for small-world effect
                        # Returns:
                        #   - A list with following elements:
                        #     - adj: Adjacency matrix of the graph
                        #     - coords: Coordinates of the nodes
                        #     - laplacian: Laplacian matrix of the graph
                        
                        grow_spatial_sf_expander <- function(N, r, beta, alpha, delta, p) {
                          # Generate initial node coordinates
                          coords <- matrix(runif(N*2), ncol=2)
                          
                          # Calculate distance matrix
                          dist_matrix <- as.matrix(dist(coords))
                          
                          # Create spatial adjacency matrix
                          A_spatial <- as.matrix(dist_matrix <= r)
                          
                          # Initialize adjacency matrix
                          adj <- matrix(0, nrow=N, ncol=N)
                          
                          # Initialize degree vector
                          deg <- rep(0, N)
                          
                          # Add first node
                          adj[1,1] <- 1
                          
                          # Loop over remaining nodes
                          for (i in 2:N) {
                            # Calculate probabilities for preferential attachment
                            comm <- rep(0, N)
                            comm[i] <- 1
                            dist_prob <- 1/(dist_matrix[i,]^alpha + 1e-10)
                            comm_prob <- delta*comm + (1-delta)*rep(1,N)
                            degree_prob <- (deg + 1)^beta
                            
                            # Calculate total probability for each node
                            prob <- dist_prob*comm_prob*degree_prob
                            prob[i] <- 0
                            prob <- prob/sum(prob)
                            
                            # Choose node to attach to
                            attach_to <- sample(1:N, size=1, prob=prob)
                            
                            # Add edge to selected node
                            adj[i, attach_to] <- 1
                            adj[attach_to, i] <- 1
                            
                            # Update degree vector
                            deg[i] <- deg[i] + 1
                            deg[attach_to] <- deg[attach_to] + 1
                            
                            # Rewire edges with small-world probability
                            if (runif(1) < p) {
                              # Choose edge to rewire
                              edges <- which(adj[i,] == 1)
                              rewire_to <- sample(setdiff(1:N, c(i, edges)), size=1)
                              
                              # Rewire edge
                              adj[i, edges] <- 0
                              adj[edges, i] <- 0
                              adj[i, rewire_to] <- 1
                              adj[rewire_to, i] <- 1
                            }
                          }
                          
                          # Calculate Laplacian matrix
                          deg_matrix <- diag(deg)
                          laplacian <- deg_matrix - adj
                          
                          # Return results
                          result <- list(adj=adj, coords=coords, laplacian=laplacian)
                          return(result)
                        }
                        
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        # Function to generate a spatial scale-free expander graph
                        spatial_scalefree_expander_graph <- function(N, beta, alpha, delta, p, r, edge_weight = FALSE, node_attr = FALSE) {
                          # N: Number of nodes in the graph
                          # beta: Power law parameter for the degree distribution
                          # alpha: Strength of degree effect
                          # delta: Probability of attachment to other nodes for community structure
                          # p: Rewiring probability for small world effect
                          # r: Cutoff distance for adjacency matrix based on spatial distance
                          
                          # Generate N nodes with spatial coordinates
                          nodes <- data.frame(x = runif(N), y = runif(N))
                          
                          # Calculate the distance matrix between all pairs of nodes
                          dist_matrix <- as.matrix(dist(nodes))
                          
                          # Create adjacency matrix based on the distance matrix
                          adj_matrix <- ifelse(dist_matrix <= r, 1, 0)
                          
                          # Generate the initial graph with m0 nodes
                          m0 <- 2
                          edges <- matrix(0, nrow = N, ncol = N)
                          edges[1:m0, 1:m0] <- adj_matrix[1:m0, 1:m0]
                          
                          # Initialize the degrees and distances of the nodes
                          degrees <- rowSums(edges)
                          distances <- dist_matrix[1,]
                          
                          # Add new nodes to the graph one at a time
                          for (i in (m0 + 1):N) {
                            # Calculate the probability of attachment to existing nodes
                            prob_attach <- (distances^(-alpha) + delta*degrees[i-1]^(-beta)) / sum(distances^(-alpha) + delta*degrees[i-1]^(-beta))
                            
                            # Choose the nodes to attach to based on the probability of attachment
                            attach_nodes <- sample(1:(i-1), m, prob = prob_attach[1:(i-1)])
                            
                            # Add edges to the attach nodes
                            edges[i, attach_nodes] <- 1
                            edges[attach_nodes, i] <- 1
                            
                            # Update the degrees and distances of the nodes
                            degrees <- rowSums(edges)
                            distances[i,] <- dist_matrix[i,]
                            
                            # Rewire edges with probability p
                            for (j in 1:(i-1)) {
                              if (edges[i,j] == 1 & runif(1) < p) {
                                new_attach_node <- sample(setdiff(1:N, c(i,j)), 1)
                                edges[i,j] <- 0
                                edges[j,i] <- 0
                                edges[i,new_attach_node] <- 1
                                edges[new_attach_node,i] <- 1
                              }
                            }
                          }
                          
                          # Add edge weights if requested
                          if (edge_weight) {
                            edges <- edges*runif(N*N)
                          }
                          
                          # Add node attributes if requested
                          if (node_attr) {
                            nodes$attr <- rnorm(N)
                          }
                          
                          # Compute the Laplacian matrix
                          degree_matrix <- diag(degrees)
                          laplacian_matrix <- degree_matrix - edges
                          
                          # Return the graph as a list
                          return(list(nodes = nodes, edges = edges, laplacian = laplacian_matrix))
                        }
                        
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        generate_spatial_scalefree_expander <- function(N, beta, alpha, delta, p, r, distance_metric = "euclidean") {
                          # Generate random points within a unit square
                          coords <- matrix(runif(N*2), ncol = 2)
                          
                          # Calculate distance matrix
                          if (distance_metric == "euclidean") {
                            dist_matrix <- as.matrix(dist(coords, method = "euclidean"))
                          } else if (distance_metric == "manhattan") {
                            dist_matrix <- as.matrix(dist(coords, method = "manhattan"))
                          }
                          
                          # Generate adjacency matrix based on distance matrix and cutoff distance
                          adj_matrix <- ifelse(dist_matrix <= r, 1, 0)
                          
                          # Loop through each new node to add to the graph
                          for (i in 1:(N-1)) {
                            # Calculate the probability of attachment to each existing node
                            distances <- dist_matrix[i, 1:i]
                            probabilities <- distances^(-beta) * (degree(adj_matrix[1:i, 1:i]) + alpha)^delta
                            probabilities <- probabilities / sum(probabilities)
                            
                            # Choose the nodes to attach to based on the probabilities
                            attachment_nodes <- sample(1:i, size = 1, prob = probabilities)
                            
                            # Add edges between the new node and the chosen attachment nodes
                            adj_matrix[i+1, attachment_nodes] <- 1
                            adj_matrix[attachment_nodes, i+1] <- 1
                            
                            # Apply small-world rewiring with probability p
                            if (p > 0) {
                              for (j in 1:i) {
                                if (adj_matrix[i+1, j] == 1 && runif(1) < p) {
                                  # Rewire the edge to a random node within a certain distance
                                  eligible_nodes <- which(dist_matrix[j,] > r & adj_matrix[i+1,] == 0)
                                  if (length(eligible_nodes) > 0) {
                                    new_attachment_node <- sample(eligible_nodes, size = 1)
                                    adj_matrix[i+1, j] <- 0
                                    adj_matrix[j, i+1] <- 0
                                    adj_matrix[i+1, new_attachment_node] <- 1
                                    adj_matrix[new_attachment_node, i+1] <- 1
                                  }
                                }
                              }
                            }
                          }
                          
                          # Compute Laplacian matrix and check eigenvalue gap
                          deg_matrix <- diag(rowSums(adj_matrix))
                          lap_matrix <- deg_matrix - adj_matrix
                          eigenvals <- eigen(lap_matrix)$values
                          eigenvals <- eigenvals[order(eigenvals)]
                          eigen_gap <- eigenvals[2] - eigenvals[1]
                          print(paste0("Eigenvalue gap: ", eigen_gap))
                          
                          # Return adjacency matrix
                          return(adj_matrix)
                        }
                        
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        library(Matrix)  # for sparse matrix operations
                        
                        # function to generate a spatial scale-free expander graph
                        # inputs: 
                        # - N: number of nodes to add to the graph
                        # - r: cutoff distance for spatial adjacency matrix
                        # - beta: power law parameter for preferential attachment
                        # - alpha: parameter for degree effect strength
                        # - delta: community structure probability
                        # - p: rewiring probability for small world effect
                        # - weighted: whether to generate edge weights or not
                        # - node_attr: list of node attributes (optional)
                        # output:
                        # - igraph object representing the generated graph
                        generate_spatial_scalefree_graph <- function(N, r=0.1, beta=1, alpha=0, delta=0, p=0, weighted=FALSE, node_attr=NULL) {
                          # create initial random graph with 3 nodes
                          x <- runif(3)
                          y <- runif(3)
                          A <- matrix(0, 3, 3)
                          for (i in 1:2) {
                            for (j in (i+1):3) {
                              if (sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2) <= r) {
                                A[i,j] <- 1
                                A[j,i] <- 1
                              }
                            }
                          }
                          deg <- rowSums(A)
                          
                          # loop to add new nodes and edges
                          for (i in 4:N) {
                            # add new node with random position in unit square
                            x[i] <- runif(1)
                            y[i] <- runif(1)
                            
                            # calculate distances to existing nodes
                            dist <- sqrt((x[1:(i-1)] - x[i])^2 + (y[1:(i-1)] - y[i])^2)
                            
                            # create adjacency matrix based on distance and community structure
                            A_new <- matrix(0, i, i)
                            A_new[1:(i-1),1:(i-1)] <- A
                            A_new[i,1:(i-1)] <- (dist <= r) * (1 + alpha*deg[1:(i-1)]) * (1 + delta*(runif(i-1) < delta))
                            A_new[1:(i-1),i] <- (dist <= r) * (1 + alpha*deg[1:(i-1)]) * (1 + delta*(runif(i-1) < delta))
                            
                            # calculate degree distribution for preferential attachment
                            deg_new <- rowSums(A_new)
                            prob_new <- (deg_new^(beta+1)) / sum(deg_new^(beta+1))
                            
                            # add edges based on preferential attachment
                            for (j in 1:(i-1)) {
                              if (runif(1) < prob_new[j]) {
                                A_new[i,j] <- 1
                                A_new[j,i] <- 1
                              }
                            }
                            
                            # add small world rewiring
                            if (p > 0) {
                              for (j in 1:(i-1)) {
                                if (A_new[i,j] == 0 && runif(1) < p) {
                                  dist_ij <- sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2)
                                  if (dist_ij <= r) {
                                    A_new[i,j] <- 1
                                    A_new[j,i] <- 1
                                  } else {
                                    dist <-
                                      
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                      spatial_scale_free_graph <- function(N, r, delta, beta, alpha, p, weighted = FALSE) {
                                        # N: number of nodes
                                        # r: cutoff distance for spatial adjacency matrix
                                        # delta: probability of attaching to a node in the same community
                                        # beta: power law parameter for preferential attachment
                                        # alpha: parameter controlling the strength of the degree effect
                                        # p: rewiring probability for small world effect
                                        # weighted: if TRUE, edge weights are assigned as 1/distance
                                        
                                        # generate random spatial coordinates for nodes on a unit square
                                        coords <- matrix(runif(N*2), ncol = 2)
                                        
                                        # calculate distance matrix between all pairs of points
                                        dist_matrix <- as.matrix(dist(coords))
                                        
                                  
                                      # create adjacency matrix based on distance matrix
                                        adj_matrix <- ifelse(dist_matrix <= r, 1, 0)
                                        
                                        # initialize graph with a single node and no edges
                                        graph <- list()
                                        graph$adj_matrix <- matrix(0, nrow = N, ncol = N)
                                        graph$adj_matrix[1, 1] <- 1
                                        graph$degrees <- rep(0, N)
                                        
                                        # add nodes one at a time
                                        for (i in 2:N) {
                                          # calculate probability of attaching to each existing node
                                          probs <- rep(0, i-1)
                                          for (j in 1:(i-1)) {
                                            # spatial distance effect
                                            d_ij <- dist(coords[c(i, j), ])
                                            spatial_prob <- ifelse(d_ij <= r, 1/1+d_ij^alpha, 0)
                                            
                                            # community structure effect
                                            community_prob <- ifelse(runif(1) < delta, adj_matrix[i, j], 0)
                                            
                                            # preferential attachment effect
                                            degree_effect <- graph$degrees[j]^beta
                                            
                                            # total probability of attaching to node j
                                            probs[j] <- (spatial_prob + community_prob + degree_effect) / sum(spatial_prob + community_prob + degree_effect)
                                          
                                          
                                          # randomly select an existing node to attach to based on the probabilities
                                          attach_node <- sample(1:(i-1), size = 1, prob = probs)
                                          
                                          # add edge between new node and attach_node
                                          graph$adj_matrix[i, attach_node] <- 1}}
                                          graph$adj_matrix[attach_node, i] <- 1
                                          graph$degrees[i] <- graph$degrees[i] + 1
                                          graph$degrees[attach_node] <- graph$degrees[attach_node] + 1
                                          
                                          # small world effect: randomly rewire edge with probability p
                                          if (runif(1) < p) {
                                            possible_rewire_nodes <- setdiff(1:N, c(i, attach_node, which(graph$adj_matrix[i, ] == 1)))
                                            if (length(possible_rewire_nodes) > 0) {
                                              rewire_node <- sample(possible_rewire_nodes, size = 1)
                                              graph$adj_matrix[i, attach_node] <- 0
                                              graph$adj_matrix[attach_node, i] <- 0
                                              graph$adj_matrix[i, rewire_node] <- 1
                                              graph$adj_matrix[rewire_node, i] <- 1
                                              graph$degrees[i] <- graph$degrees[i] - 1
                                              graph$degrees[attach_node] <- graph$degrees[attach_node] - 1
                                              graph$degrees[rewire_node] <- graph

                                              
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                              library(Matrix)
                                              
                                              grow_spatial_scalefree_expander_graph <- function(N, beta, alpha, delta, r, p, distance_func = "euclidean", rewiring_distance = Inf, seed = NULL) {
                                                # Set seed if provided
                                                if (!is.null(seed)) {
                                                  set.seed(seed)
                                                }
                                                
                                                # Generate initial set of nodes
                                                coords <- matrix(runif(N*2), ncol = 2)
                                                dist_mat <- as.matrix(dist(coords, method = distance_func))
                                                
                                                # Create adjacency matrix based on distance matrix
                                                adj_mat <- ifelse(dist_mat <= r, 1, 0)
                                                
                                                # Initialize Laplacian matrix
                                                lap_mat <- diag(rowSums(adj_mat)) - adj_mat
                                                
                                                # Generate new nodes and attach edges
                                                for (i in seq(N)[-1]) {
                                                  # Calculate probability of attachment to existing nodes
                                                  dist_to_i <- dist_mat[i, ]
                                                  degree_of_i <- rowSums(adj_mat[i, ])
                                                  prob_of_attach <- (dist_to_i ^ (-beta)) + (alpha * degree_of_i) + (delta * (lap_mat[i,] / max(lap_mat[i,])))
                                                  prob_of_attach[i] <- 0
                                                  prob_of_attach <- prob_of_attach / sum(prob_of_attach)
                                                  
                                                  # Choose node to attach to
                                                  attach_to <- sample(seq(N)[-i], size = 1, prob = prob_of_attach[-i])
                                                  
                                                  # Add edge to adjacency matrix
                                                  adj_mat[i, attach_to] <- 1
                                                  adj_mat[attach_to, i] <- 1
                                                  
                                                  # Update Laplacian matrix
                                                  lap_mat[i, i] <- lap_mat[i, i] + 1
                                                  lap_mat[attach_to, attach_to] <- lap_mat[attach_to, attach_to] + 1
                                                  lap_mat[i, attach_to] <- -1
                                                  lap_mat[attach_to, i] <- -1
                                                  
                                                  # Rewire edges with probability p
                                                  if (p > 0) {
                                                    rewiring_probs <- adj_mat[i, ] * (dist_to_i <= rewiring_distance)
                                                    rewiring_probs[i] <- 0
                                                    if (sum(rewiring_probs) > 0) {
                                                      rewire_to <- sample(seq(N), size = 1, prob = rewiring_probs / sum(rewiring_probs))
                                                      adj_mat[i, rewire_to] <- 1
                                                      adj_mat[rewire_to, i] <- 1
                                                      
                                                      # Update Laplacian matrix
                                                      lap_mat[i, i] <- lap_mat[i, i] + 1
                                                      lap_mat[rewire_to, rewire_to] <- lap_mat[rewire_to, rewire_to] + 1
                                                      lap_mat[i, rewire_to] <- -1
                                                      lap_mat[rewire_to, i] <- -1
                                                      lap_mat[attach_to, i] <- lap_mat[attach_to, i] - 1
                                                      lap_mat[i, attach_to] <- lap_mat[i, attach_to] - 1
                                                    }
                                                  }
                                                  
                                                  # Update distance matrix
                                                  dist_mat <- min(dist_mat, as.matrix(dist(coords[i,], coords, method = distance_func)))
                                                }
                                                
                                                # Create weighted adjacency matrix
                                                edge_weights <- ifelse(adj_mat == 1, runif(N*(N-1)/2, min = 0.1, max = 1), 0)
                                                wadj_mat <- t(edge_weights)
                                                
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                              # Define function to grow the spatial scale-free expander graph
                                              grow_spatial_scale_free_graph <- function(N, r, beta, alpha, delta, p, max_attempts = 100) {
                                                # Generate nodes uniformly at random within a unit square
                                                coords <- matrix(runif(2 * N), ncol = 2)
                                                
                                                # Calculate the distance matrix between all pairs of points
                                                dist_matrix <- as.matrix(dist(coords))
                                                
                                                # Create an adjacency matrix based on the distance matrix and the cutoff distance 'r'
                                                adj_matrix <- ifelse(dist_matrix <= r, 1, 0)
                                                diag(adj_matrix) <- 0 # no self-loops
                                                
                                                # Initialize the graph with the first two nodes connected
                                                graph <- matrix(0, ncol = N, nrow = N)
                                                graph[1, 2] <- 1
                                                graph[2, 1] <- 1
                                                
                                                # Initialize the node degrees and community assignments
                                                degrees <- rep(2, N)
                                                communities <- rep(1, N)
                                                num_communities <- 1
                                                
                                                # Loop through the remaining nodes and connect them to the graph
                                                for (i in 3:N) {
                                                  attempts <- 0
                                                  while (attempts < max_attempts) {
                                                    # Calculate the probability of attachment to each existing node
                                                    # based on spatial distance, community structure, and degree
                                                    dist_probs <- exp(-dist_matrix[i, ]^beta)
                                                    comm_probs <- ifelse(communities == communities[i], delta, 1)
                                                    degree_probs <- degrees^alpha
                                                    attach_probs <- dist_probs * comm_probs * degree_probs
                                                    attach_probs[i] <- 0 # exclude self-attachment
                                                    
                                                    # Choose an existing node to attach to with probability proportional to attachment probability
                                                    if (sum(attach_probs) == 0) break # no suitable attachment node found
                                                    attach_probs <- attach_probs / sum(attach_probs)
                                                    attach_node <- sample(1:N, size = 1, prob = attach_probs)
                                                    
                                                    # Rewire edge with probability 'p' to a random node within a certain distance
                                                    if (runif(1) < p) {
                                                      # Find candidate rewiring nodes
                                                      cand_nodes <- which(dist_matrix[attach_node, ] > r & 
                                                                          communities != communities[i] & 
                                                                          degrees > 0)
                                                      if (length(cand_nodes) > 0) {
                                                        # Choose rewiring node at random
                                                        rewire_node <- sample(cand_nodes, size = 1)
                                                        # Rewire edge from i to attach_node to i to rewire_node
                                                        graph[i, attach_node] <- 0
                                                        graph[attach_node, i] <- 0
                                                        graph[i, rewire_node] <- 1
                                                        graph[rewire_node, i] <- 1
                                                        # Update degrees and communities
                                                        degrees[attach_node] <- degrees[attach_node] - 1
                                                        degrees[rewire_node] <- degrees[rewire_node] + 1
                                                        communities[i] <- communities[rewire_node]
                                                        break
                                                      }
                                                    }
                                                    
                                                    # Attach new node to chosen node
                                                    if (graph[i, attach_node] == 0) {
                                                      graph[i, attach_node] <- 1
                                                      graph[attach_node, i] <- 1
                                                      degrees[i] <- degrees[i] + 1
                                                      degrees[attach_node] <- degrees[attach_node] + 1
                                                      communities[i] <-
                                                        
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+sample code
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Define a function to generate a spatial graph on a d-dimensional torus
torus_graph <- function(num_nodes, dimensions=2, min_distance=0.1, max_distance=0.4) {
  
  # Define a quadtree node class to represent spatial subdivisions
  quadtree_node <- function(xmin, xmax, ymin, ymax, zmin=NULL, zmax=NULL, level=0) {
    node <- list(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, level=level)
    if (dimensions >= 3) {
      node$zmin <- zmin
      node$zmax <- zmax
    }
    node$children <- list()
    node$points <- list()
    node$centroid <- NULL
    return(node)
  }
  
  # Define a function to insert a point into a quadtree
  insert_point <- function(node, point) {
    if (length(node$children) == 0) {
      # If this node is a leaf, add the point to its list of points
      node$points[[length(node$points) + 1]] <- point
      return(node)
    } else {
      # Otherwise, insert the point into the appropriate child node
      child_index <- get_child_index(node, point)
      node$children[[child_index]] <- insert_point(node$children[[child_index]], point)
      return(node)
    }
  }
  
  # Define a function to retrieve the index of the child node that contains a given point
  get_child_index <- function(node, point) {
    xmid <- (node$xmin + node$xmax) / 2
    ymid <- (node$ymin + node$ymax) / 2
    if (dimensions >= 3) {
      zmid <- (node$zmin + node$zmax) / 2
    }
    if (point$x < xmid) {
      if (point$y < ymid) {
        if (dimensions >= 3 && point$z < zmid) {
          return(1)
        } else {
          return(1)
        }
      } else {
        if (dimensions >= 3 && point$z < zmid) {
          return(3)
        } else {
          return(2)
        }
      }
    } else {
      if (point$y < ymid) {
        if (dimensions >= 3 && point$z < zmid) {
          return(2)
        } else {
          return(3)
        }
      } else {
        if (dimensions >= 3 && point$z < zmid) {
          return(4)
        } else {
          return(4)
        }
      }
    }
  }
  
  # Define a function to calculate the centroid of a quadtree node
  calculate_centroid <- function(node) {
    num_points <- length(node$points)
    if (num_points == 0) {
      return(NULL)
    } else {
      x <- sapply(node$points, function(p) p$x)
      y <- sapply(node$points, function(p) p$y)
      if (dimensions >= 3) {
        z <- sapply(node$points, function(p) p$z)
        centroid <- list(x=mean(x), y=mean(y), z=mean(z))
      } else {
        centroid <- list(x=mean(x), y=mean(y))
      }
      return(centroid)
    }
  }
  
  # Generate an initial set of random points on the tor
  

# library(sp)
# library(rgeos)
# 
# quadtree <- quadtree(nodes)
# dist_matrix <- as.matrix(gDistance(quadtree, byid = TRUE))
# adj_matrix <- ifelse(dist_matrix <= r, 1, 0)