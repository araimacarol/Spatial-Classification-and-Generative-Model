library(future)
library(parallel)
library(future.apply)
library(doParallel)
library(progressr)
library(igraph)
library(dplyr)
library(randomForest)

###--Simulation of Pathogen----###
simPathE <- function(G, nTicks=100, beta=0.5, gamma=0.1, propInfected=0.1, initialState=NULL, nInfected=1, useProportion=F) {
  nVertices = gorder(G)
  infectionState = as.data.frame(matrix(0, ncol = nVertices, nrow = nTicks+1))
  
  # Set up the initial state: all vertices are either exposed (state = 0) or infected (state = 1); none can be recovered yet (2).
  if (is.null(initialState)) {
    if (useProportion==T) {
      infectionState[1,] <- rbinom(nVertices, 1, propInfected) # set initial state of each nodes if known (note, no recovered node at initial state)
    } else {
      infected <- rep(1, nInfected) # just create a vector of the right number of 1s
      exposed <- rep(0, (nVertices - nInfected))
      infectionState[1,] <- sample(c(infected, exposed), nVertices, replace=FALSE)
    }
  } else {
    if (length(initialState) != nVertices) {
      return ("Initial state and order of Graph (number of vertices) are incompatible.  Check the sizes of your input.")
    }
    infectionState <- initialState # initial existing state.
  }
  
  adjacencyList <- as_adj_list(G) # this is a list of which vertices are adjacent to each vertex i
  
  # Now do the simulation through time:
  for (t in 1:nTicks) {
    # FIRST phase: transmission: S -> I
    for (i in which(infectionState[t,] == 0)) { # for all susceptible nodes (denoted 0) in previous time step
      infectionState[t+1,i] <- 0 # node remains as susceptible, since not all contact leads to an infection.
      numDodgyNeighbours = length(which(infectionState[t,adjacencyList[[i]]] == 1))
      if ((runif(1) < 1.0 - (1.0 - beta)^numDodgyNeighbours)) {
        infectionState[t+1,i] <- 1;
      }
      # for (j in adjacencyList[[i]]) { # for all neighbours of i
      #   if (infectionState[t,j][1] == 1) { # vertex j is infectious
      #     if ((runif(1)*G[i,j]) <= beta) { # ... and passes it on!
      #       infectionState[t+1,i] <- 1;
      #       break # assign node as infected if above condition is met, and break out of loop: we don't need
      #       # to check any more adjacent vertices.
      #     }
      #   }
      # }
    }
    
    # SECOND phase: recovery: I -> R
    for (i in which(infectionState[t,] == 1)) { # for all infected nodes (denoted 1) in previous time step
      if (runif(1) <= gamma) { # compares a randomly generated uniform number to recovery rate
        infectionState[t+1,i] <- 2 # node is recovered
      } else {
        infectionState[t+1,i] <- 1 # node remains infected
      }
    }
    
    # THIRD phase: recovered stays recovered:
    for (i in which(infectionState[t,] == 2)) { # for all recovered nodes (denoted 2) in previous time step
      infectionState[t+1,i] <- 2 # node stays recovered
    }
  }
  rownames(infectionState) = 0:nTicks
  return(infectionState)
}

##########################################################################################################
# Count the number of times each value (state) occurs per row.
##########################################################################################################
countStates <- function(DF=NULL, states=NULL) {
  nr = dim(DF)[1]
  counts <- as.data.frame(matrix(NA,ncol=length(states), nrow=nr))
  colnames(counts) = as.character(states)
  for (i in 1:nr) {
    col = 1
    for (state in states) {
      counts[i,col] = length(which(DF[i,]==state))
      col = col+1
    }
  }
  return(counts)
}



#############--------------MAKING OF SPATIAL GRAPH-------------------####################


# Helper functions to convert between (row, column) pairs and the index in a list of Cells.
cellCoordsToIndex <- function(i, j, size) {
  return((i-1)*size+j)
}
indexToCellCoords <- function(idx, size) {
  j <- (idx-1) %% size + 1
  i <- (idx + size-1) %/% size
  return(c(i,j))
}

# - end of Helper functions

# 
# for (i in 1:4) {
#   for (j in 1:4) {
#     idx <- cellCoordsToIndex(i,j,4)
#     print(paste("i =", i, "; j = ", j , "index = " , idx ,
#                 "; reverse = ", indexToCellCoords(idx, 4)[1], indexToCellCoords(idx, 4)[2]))
#     print("")
#   }
# }


###Function to create Spatial network ##
fastSpatialNetwork <- function(n=100,r=0.1, makeConnected=FALSE, keepCellsSeparate=FALSE) {
  # Divide the grid into cells of diagonal length r (so side length r/sqrt(2))
  # All pairs of points in the same cell are within distance r of each other
  # Look around each cell at the points in the 8 adjacent cells
  # Test each for Euclidean distance.
  # Connect all that are close enough; note those that aren't for later addition
  
  # Set up the coordinates:
  v <- data.frame(matrix(nrow=n))
  v$x <- runif(n)
  v$y <- runif(n)
  v[,1] <- 1:n
  colnames(v) <- c("id","x","y")
  
  # put points into bins:
  cellSize <- r/sqrt(2)
  R2 <- r*r
  gridSize <- ceiling(sqrt(2)/r)
  # Create a grid of cells and assign vertices to them based on their coordinates:
  Cell <- vector(mode="list", gridSize*gridSize)
  for (i in 1:n) {
    #		rowNum <- floor(v$y[i]/cellSize) + 1
    #		colNum <- floor(v$x[i]/cellSize) + 1
    idx <- (floor(v$y[i]/cellSize))*gridSize + floor(v$x[i]/cellSize) + 1
    Cell[[idx]] <- cbind(Cell[[idx]], i)
  }
  E <- matrix(nrow=0, ncol=2, byrow=TRUE)
  
  # join all vertices that are in the same cell since they must be at most distance r apart
  for (cell in Cell) {
    if (length(cell) > 1) {
      E <- rbind(E, t(combn(cell, 2)))
    }
  }
  proximalPairs <- data.frame(matrix(5,ncol=3,nrow=1))
  colnames(proximalPairs) <- c("r","i","j")
  
  if (keepCellsSeparate == FALSE) {
    for (i in 2:gridSize) {
      # first column:
      # join any points in [ ]
      #                    [ ] that are close enough
      thisIdx <- cellCoordsToIndex(1,i, gridSize)
      aboveIdx <- cellCoordsToIndex(1,i-1, gridSize)
      for (a in Cell[[thisIdx]]) {
        for (b in Cell[[aboveIdx]]) {
          d <- (v[a,2]-v[b,2])^2 + (v[a,3]-v[b,3])^2
          if (d < R2) {
            E <- as.matrix(rbind(E,c(a,b)))
          } else {
            proximalPairs <- rbind(proximalPairs, c(d, a, b))
          }
        }
      }
      # first row:
      # join any points in [ ][ ] that are close enough.
      thisIdx <- cellCoordsToIndex(i,1, gridSize)
      leftIdx <- cellCoordsToIndex(i-1,1, gridSize)
      for (a in Cell[[thisIdx]]) {
        for (b in Cell[[leftIdx]]) {
          d <- (v[a,2]-v[b,2])^2 + (v[a,3]-v[b,3])^2
          if (d < R2) {
            E <- as.matrix(rbind(E,c(a,b)))
          } else {
            proximalPairs <- rbind(proximalPairs, c(d, a, b))
          }
        }
      }
      
      for (j in 2:gridSize) {
        # check all neighbours above and to the left
        #
        # aboveLeftIdx      aboveIdx
        #             XX       ||
        #     leftIdx    ==  thisIdx
        thisIdx <- cellCoordsToIndex(i,j,gridSize)
        aboveIdx <- cellCoordsToIndex(i,j-1,gridSize)
        
        # this + above:
        for (a in Cell[[thisIdx]]) {
          for (b in Cell[[aboveIdx]]) {
            d <- (v[a,2]-v[b,2])^2 + (v[a,3]-v[b,3])^2
            if (d < R2) {
              E <- as.matrix(rbind(E,c(a,b)))
            } else {
              proximalPairs <- rbind(proximalPairs, c(d, a, b))
            }
          }
        }
        
        # this + left:
        leftIdx <- cellCoordsToIndex(i-1,j,gridSize)
        for (a in Cell[[thisIdx]]) {
          for (b in Cell[[leftIdx]]) {
            d <- (v[a,2]-v[b,2])^2 + (v[a,3]-v[b,3])^2
            if (d < R2) {
              E <- as.matrix(rbind(E,c(a,b)))
            } else {
              proximalPairs <- rbind(proximalPairs, c(d, a, b))
            }
          }
        }
        
        # this + aboveLeft
        aboveLeftIdx <- cellCoordsToIndex(i-1,j-1,gridSize)
        for (a in Cell[[thisIdx]]) {
          for (b in Cell[[aboveLeftIdx]]) {
            d <- (v[a,2]-v[b,2])^2 + (v[a,3]-v[b,3])^2
            if (d < R2) {
              E <- as.matrix(rbind(E,c(a,b)))
            } else {
              proximalPairs <- rbind(proximalPairs, c(d, a, b))
            }
          }
        }
        
        # left and above:
        for (a in Cell[[leftIdx]]) {
          for (b in Cell[[aboveIdx]]) {
            d <- (v[a,2]-v[b,2])^2 + (v[a,3]-v[b,3])^2
            if (d < R2) {
              E <- as.matrix(rbind(E,c(a,b)))
            } else {
              proximalPairs <- rbind(proximalPairs, c(d, a, b))
            }
          }
        }
      }
    }
  }
  
  G = make_empty_graph(n, directed=FALSE)
  G <- add_edges(G, t(E), directed=FALSE)
  vertex_attr(G)$id <- 1:n
  G$layout <- as.matrix(v[,2:3])
  G$name <- "Spatial"
  G$type <- "Spatial"
  G$id <- "1"
  if (makeConnected == F) {
    return(G)
  }
  if (is.connected(G)) {
    return(G)
  }
  # sort the pairs of coordinates in ascending order of distance apart
  proximalPairs <- proximalPairs[-1,]
  rownames(proximalPairs) <- seq(1,dim(proximalPairs)[1])
  indx <- order(proximalPairs$r)
  
  comps <- components(G, mode="weak")
  
  for (nu in indx) { # in increasing order of edge length
    i <- proximalPairs$i[nu]
    j <- proximalPairs$j[nu]
    
    # Consider vertices i and j for joining:
    if (comps$membership[i] != comps$membership[j]) { # if two vertices are not in the same component
      G <- G + edge(i, j, directed=F) # then join them
      comps <- components(G, mode="weak")
    } else {
      # i and j are in the same component.  Moving on...
    }
    if (is.connected(G) == T) {
      break # we are done!
    }
  }
  return(G)
}

induced.graph=function(Graphs){
initgraph=NULL;G=NULL
for (i in 1:length(Graphs)){  
  initgraph=Graphs[[i]]
  components = igraph::clusters(initgraph , mode="weak")
  biggest_cluster_id = which.max(components$csize)
  # # ids
  vert_ids = V(initgraph)[components$membership== biggest_cluster_id]
  # # subgraph
  G[[i]]=igraph::induced_subgraph(initgraph, vert_ids)
}
  return(initgraph)
}

#induced.graph(x)
#induced.graph(x[1])


##########################################################################################################
# Create a list of graphs for simulation experiment
##########################################################################################################
makeSpatialGraphs<- function(node.size=25,Radius=0.8) {
  Graphs = list()#;G=list() # set up the list of Graphs first
  initgraph= list()
  i= 1
   print("Creating fastspatial networks")
   Graphs[[i]] = fastSpatialNetwork(n = node.size, r = Radius, makeConnected=T,keepCellsSeparate=FALSE)
   #Graphs[[i]]=induced.graph(G[i])
    # initgraph=Graphs[[i]]
    # components = igraph::clusters(initgraph , mode="weak")
    # biggest_cluster_id = which.max(components$csize)
    # # # ids
   # vert_ids = V(initgraph)[components$membership== biggest_cluster_id]
    # # subgraph
  #  Graphs[[i]]=igraph::induced_subgraph(initgraph, vert_ids)
    Graphs[[i]]$type = "Spatial"
    Graphs[[i]]$id = "1"
    Graphs[[i]]$name="Spatial"
    i <- i+1
  return(Graphs)
}


##########################################################################################################
# Input an igraph object from a file, treating it as a static graph
##########################################################################################################
getGraphFromFile <- function(file, simplify=TRUE, useBiggestComponent=TRUE, asUndirected=TRUE) {
  
  dat <- read.table(file) # just read static graph: ignoring third column
  G <- graph_from_data_frame(dat)
  if (asUndirected==TRUE) {
    G <- as.undirected(G, "collapse")
  }
  
  #  g_names <- gsub(".edges","",networks[i]) # one edge for each pair of connect vertices (not sure what this is for)
  
  if (useBiggestComponent==TRUE) {
    netclust <- clusters(G) #look for subgraphs
    gcc <- V(G)[netclust$membership == which.max(netclust$csize)]#select vertices from the largest sub-graph
    G <- induced.subgraph(G, gcc) #make it a igraph object.
  }
  if (simplify==TRUE) {
    G <- igraph::simplify(G, remove.multiple = TRUE, remove.loops = TRUE)
  }
  return(G)
}



calcGraphFeatures <- function(Graphs=NULL) {
  
  features <- c(
    "order",                    # number of vertices
    "edges",                     # number of edges
    "connected",                # True / False
    "max_component",            # maximum component size (=order iff the graph is connected)
    "minDegree",                # minimum degree of any vertex
    "maxDegree",                # maximum degree of any vertex
    "mean_degree",                # average degree of any vertex
    "minCut",                   # minimum cut weight of the graph (might take a while to compute)
    "FiedlerValue",             # second-highest eigenvalue of the Laplacian matrix
    "Normalized_FiedlerValue",   # second-highest eigenvalue of the Normaized Laplacian matrix
    "closeness",                # average inverse of distance between any pair of vertices
    "modularity",               # DEFINITION REQUIRED
    "diameter",                 # maximum distance between any two vertices (NAN if not connected)
    "betweenness",              # max_{v} proportion of shortest paths going through vertex v
    "transitivity",             # aka Clustering Coefficient, is proportion of connected triples that form triangles: e.g., (a--b--c--a) when (a--b--c) is present.
    "threshold",                 # 1/max(eigen value of A)
    "spectral_radius"         # max (eigen value of A)
    
  )
  
  df <- as.data.frame(matrix(ncol=length(features),nrow=length(Graphs)))
  colnames(df)=features
  
  # Stuff that is simple to apply and needs no interim components:
  
  df$order = base::as.numeric(lapply(Graphs, gorder))
  df$edges = base::as.numeric(lapply(Graphs, gsize))
  df$connected = base::as.numeric(lapply(Graphs, is.connected))
  df$minCut = base::as.numeric(lapply(Graphs, min_cut))
  df$diameter = base::as.numeric(lapply(Graphs, diameter))
  df$transitivity = base::as.numeric(lapply(Graphs, transitivity))
  
  # stuff that needs interim things:
  degrees = lapply(Graphs,igraph::degree )
  df$minDegree = base::as.numeric(lapply(degrees, min))
  df$maxDegree = base::as.numeric(lapply(degrees, max))
  df$mean_degree = base::as.numeric(lapply(degrees, mean))
  
  # stuff that lapply doesn't like so has to be done in a loop:
  communities <-lapply(Graphs, cluster_walktrap) #lapply(Graphs, cluster_leading_eigen)
  Adj <- lapply(Graphs, as_adjacency_matrix)
  L <- lapply(Graphs, laplacian_matrix)
  Norm_Lap<-lapply(Graphs, normalized_laplacian)
  Fiedler.value=NULL
  norm.fiedler.value=NULL
  
  for (i in 1:length(Graphs)) {
    if (is.null(Graphs[[i]]$type)) { Graphs[[i]]$type = "untyped" }
    df$modularity[i] <- modularity(communities[[i]])
    df$spectral_radius[i] <- eigen(Adj[[i]], symmetric=TRUE, only.values=TRUE)$values[1]
    
    Fiedler.value[[i]]=eigen(L[[i]], symmetric=TRUE, only.values=TRUE)$values
    
    df$FiedlerValue[i] <- Fiedler.value[[i]][length(Fiedler.value[[i]])-1]
    
    norm.fiedler.value[[i]]=eigen(Norm_Lap[[i]], symmetric=TRUE, only.values=TRUE)$values
    
    df$Normalized_FiedlerValue[i] <- norm.fiedler.value[[i]][length(norm.fiedler.value[[i]])-1]
    
    df$centrality_eigen[i] <- centr_eigen(Graphs[[i]])$centralization
    
    df$betweenness[i] <- centr_betw(Graphs[[i]])$centralization
    
    df$max_component[i] <- max(components(Graphs[[i]])$csize)
    
    df$threshold[i] <- 1/(df$spectral_radius[i])
    
    if (df$connected[i]==TRUE) {
      df$closeness[i] = mean(closeness(Graphs[[i]]))
    } else { # handle the case where G isn't connected
      df$closeness[i] = -1
    }
  }
  return (df)
}

calcGraphFeatures(list(x))

# calcGraphFeatures <- function(Graphs=NULL) {
#   
#   features <- c(
#     "order",                    # number of vertices
#     "size",                     # number of edges
#     "connected",                # True / False
#     "max_component",            # maximum component size (=order iff the graph is connected)
#     "minDegree",                # minimum degree of any vertex
#     "maxDegree",                # maximum degree of any vertex
#     "aveDegree",                # average degree of any vertex
#     "minCut",                   # minimum cut weight of the graph (might take a while to compute)
#     "FiedlerValue",             # second-highest eigenvalue of the Laplacian matrix
#     "Normalized_FiedlerValue",   # second-highest eigenvalue of the Normaized Laplacian matrix
#     "closeness",                # average inverse of distance between any pair of vertices
#     "modularity",               # DEFINITION REQUIRED
#     "diameter",                 # maximum distance between any two vertices (NAN if not connected)
#     #    "betweenness",              # max_{v} proportion of shortest paths going through vertex v
#     "transitivity",             # aka Clustering Coefficient, is proportion of connected triples that form triangles: e.g., (a--b--c--a) when (a--b--c) is present.
#     "threshold",                 # 1/max(eigen value of A)
#     "spectral_radius"         # max (eigen value of A)
#     
#   )
#   
#   df <- as.data.frame(matrix(ncol=length(features),nrow=length(Graphs)))
#   colnames(df)=features
#   
#   # Stuff that is simple to apply and needs no interim components:
#   
#   df$order = as.numeric(lapply(Graphs, gorder))
#   df$size = as.numeric(lapply(Graphs, gsize))
#   df$connected = as.numeric(lapply(Graphs, is.connected))
#   df$minCut = as.numeric(lapply(Graphs, min_cut))
#   df$diameter = as.numeric(lapply(Graphs, diameter))
#   df$transitivity = as.numeric(lapply(Graphs, transitivity))
#   
#   # stuff that needs interim things:
#   degrees = lapply(Graphs, degree)
#   df$minDegree = as.numeric(lapply(degrees, min))
#   df$maxDegree = as.numeric(lapply(degrees, max))
#   df$aveDegree = as.numeric(lapply(degrees, mean))
#   
#   # stuff that lapply doesn't like so has to be done in a loop:
#   communities <- lapply(Graphs, cluster_walktrap)#finds community mebership/structures (densely connected sub groups) via random walk
#   
#   
#   Adj <- lapply(Graphs, as_adjacency_matrix)
#   
#   L <- lapply(Graphs, laplacian_matrix)
#   
#   Norm_Lap<-lapply(Graphs, normalized_laplacian)
#   
#   Fiedler.Value=NULL
#   Normalized.FiedlerValue=NULL
#   
#   for (i in 1:length(Graphs)) {
#     if (is.null(Graphs[[i]]$type)) { Graphs[[i]]$type = "untyped" }
#     
#     df$modularity[i] <- modularity(communities[[i]])#calculates the modularity of the graph using the memberships
#     
#     f2.adj <- function(x, extra=NULL) { cat("."); as.vector(Adj[[i]] %*% x) }
#     
#     df$spectral_radius[i] <- arpack(f2.adj, sym=TRUE, options=list(n=vcount(Graphs[[i]]), nev=3, ncv=5,
#                                                                    which="LA",maxiter=3000000))$values[1]
#     #eigen(Adj[[i]], symmetric=TRUE, only.values=TRUE)$values[1]#RSpectra::eigs(Adj[[i]], k, opts = list(retvec = F))$values[[1]] 
#     
#     f2.lap <- function(x, extra=NULL) { cat("."); as.vector(L[[i]] %*% x) }
#     
#     Fiedler.Value[[i]] <- arpack(f2.lap, sym=TRUE, options=list(n=vcount(Graphs[[i]]), nev=3, ncv=5,
#                                                                 which="SA",maxiter=3000000))$values
#     
#     #eigen(L[[i]], symmetric=TRUE, only.values=TRUE)#RSpectra::eigs(L[[i]], k, opts = list(retvec = F))$values[[1]]#
#     
#     df$FiedlerValue[i] = Fiedler.Value[[i]][length(Fiedler.Value[[i]])-1]
#     
#     
#     f2.normlap <- function(x, extra=NULL) { cat("."); as.vector(Norm_Lap[[i]] %*% x) }
#     
#     Normalized.FiedlerValue[[i]]<- arpack(f2.normlap, sym=TRUE, options=list(n=vcount(Graphs[[i]]), nev=3, ncv=5,
#                                                                              which="SA",maxiter=3000000))$values
#     
#     #eigen(Norm_Lap[[i]], symmetric=TRUE, only.values=TRUE)#RSpectra::eigs(Norm_Lap[[i]], k, opts = list(retvec = FALSE))#
#     
#     df$Normalized_FiedlerValue[i]=Normalized.FiedlerValue[[i]][length(Normalized.FiedlerValue[[i]])-1]
#     
#     df$centrality_eigen[i] <- centr_eigen(Graphs[[i]])$centralization
#     
#     df$betweenness[i] <- centr_betw(Graphs[[i]])$centralization
#     
#     df$max_component[i] <- max(components(Graphs[[i]])$csize)
#     
#     df$threshold[i] <- 1/(df$spectral_radius[i])
#     
#     if (df$connected[i]==TRUE) {
#       df$closeness[i] = mean(closeness(Graphs[[i]]))
#     } else { # handle the case where G isn't connected
#       df$closeness[i] = -1
#     }
#   }
#   return (df)
# }



##----------- Graph Features----------#####
## Used to perform runs on multiple simulated graphs on any given network
RunSimOnGraphFeatures<-function(Graphs, nreps=nreps,output_file="GraphFeatures.csv", seed=-1) {
  set.seed(1)
  # ### Definition and initialization of parameters for graphfeatures
  graphProperties=list()
  # ### Definition and initialization of Graph Prefix
  graphid=list();graphreplicate=list(); graphname=list(); GraphPrefix=list(); analysis=list()
  for (g in 1:length(Graphs)){
    for (reps in 1:nreps) {
      ### Calculate the graph features for each simulated graph of all the synthetic networks
      print(paste("Calculating graph features on", Graphs[[g]]$name))
      #graphProperties[[reps]] <- calcGraphFeatures(Graphs[g])
      graphProperties[[reps]] <- calcGraphFeatures(Graphs[g])
    }
    graphname[[g]]=Graphs[[g]]$type
    graphid[[g]]=Graphs[[g]]$id
    graphreplicate[[g]]=c(1:nreps)
    GraphPrefix=cbind(graphname[[g]],graphid[[g]],graphreplicate[[g]])
    colnames(GraphPrefix)=c("GraphName","GraphID","GraphReplicate")
    analysis[[g]]=as.data.frame(cbind(GraphPrefix,graphProperties[[reps]]))
    row.names(analysis[[g]])=1:nreps
  }
  All_results=do.call(rbind,analysis)
  write.csv(All_results, file=output_file)
  return( All_results)
}

# Create a list of graphs for simulation experiment
##########################################################################################################
makeTheoGraphs <- function(nSamples=10, order=25,edges=45,dim.sw=2,nei.sw=2,
                           p.sw=0.2,power.sf=2,m.sf=10,r.sp=0.345,dim.lat=2,nei.lat=1) {
  
  graphTypes <- c("Erdos-Renyi", "Small World", "Scale Free", "Square Lattice", "Spatial")
  
  Graphs <- list(mode="any", nSamples * length(graphTypes)) # set up the list of Graphs first
  initgraph=list(mode="any", nSamples * length(graphTypes))
  
  print("Creating Erdos-Renyi graphs")
  i <- 1
  for (j in 1:nSamples) {
    Graphs[[i+j-1]] <-igraph::erdos.renyi.game(order,edges, type = "gnm" , directed = F , loops = F)
    while (is.connected(Graphs[[i+j-1]]) == F) {
      Graphs[[i+j-1]] <- igraph::erdos.renyi.game(order,edges , type = "gnm" , directed = F , loops = F)
    }
    initgraph=Graphs[[i+j-1]]
    components<- igraph::clusters(initgraph , mode="weak")
    biggest_cluster_id<- which.max(components$csize)
    # # ids
    vert_ids<- V(initgraph)[components$membership== biggest_cluster_id]
    # # subgraph
    Graphs[[i+j-1]]=igraph::induced_subgraph(initgraph, vert_ids)
    Graphs[[i+j-1]]$type <- "ER"
    Graphs[[i+j-1]]$id <- j
  }
  i <- i+nSamples
  
  print("Creating small-world graphs using the Watts-Strogatz method")
  for (j in 1:nSamples) {
    Graphs[[i+j-1]] <-igraph::sample_smallworld(dim=dim.sw, size=round(sqrt(order)),nei=nei.sw,p=p.sw) 
    initgraph=Graphs[[i+j-1]]
    components<- igraph::clusters(initgraph , mode="weak")
    biggest_cluster_id<- which.max(components$csize)
    # # ids
    vert_ids<- V(initgraph)[components$membership== biggest_cluster_id]
    # # subgraph
    Graphs[[i+j-1]]=igraph::induced_subgraph(initgraph, vert_ids)
    Graphs[[i+j-1]]$type <- "SW"
    Graphs[[i+j-1]]$id <- j
  }
  i <- i+nSamples
  
  print("Creating scale-free graphs using the Barabasi-Albert method")
  for (j in 1:nSamples) {
    Graphs[[i+j-1]] <- igraph::sample_pa(n=order, power=power.sf, m=m.sf, directed=FALSE, algorithm="psumtree")
    initgraph=Graphs[[i+j-1]]
    components<- igraph::clusters(initgraph , mode="weak")
    biggest_cluster_id<- which.max(components$csize)
    # # ids
    vert_ids<- V(initgraph)[components$membership== biggest_cluster_id]
    # # subgraph
    Graphs[[i+j-1]]=igraph::induced_subgraph(initgraph, vert_ids)
    Graphs[[i+j-1]]$type <- "SF"
    Graphs[[i+j-1]]$id <- j
  }
  i <- i+nSamples
  
  print("Creating fastspatial networks")
  for (j in 1:nSamples) {
    Graphs[[i+j-1]] <- fastSpatialNetwork(n = order, r = r.sp, makeConnected=T)
    initgraph=Graphs[[i+j-1]]
    components<- igraph::clusters(initgraph , mode="weak")
    biggest_cluster_id<- which.max(components$csize)
    # # ids
    vert_ids<- V(initgraph)[components$membership== biggest_cluster_id]
    # # subgraph
    Graphs[[i+j-1]]=igraph::induced_subgraph(initgraph, vert_ids)
    Graphs[[i+j-1]]$type <- "SP"
    Graphs[[i+j-1]]$id <- as.character(j)
  }
  i <- i+nSamples
  
  print("Creating a single square lattice graph (they're all the same!)")
  Graphs[[i]] <- igraph::make_lattice(length = ceiling(sqrt(order)),dim=dim.lat,nei = nei.lat)
  initgraph=Graphs[[i]]
  components<- igraph::clusters(initgraph , mode="weak")
  biggest_cluster_id<- which.max(components$csize)
  # # ids
  vert_ids<- V(initgraph)[components$membership== biggest_cluster_id]
  # # subgraph
  Graphs[[i]]=igraph::induced_subgraph(initgraph, vert_ids)
  Graphs[[i]]$type <- "Lat"
  Graphs[[i]]$id <- 1
  i <- i+1
  
  return(Graphs)
}


##----------- Multiple simulation of the given disease parameters(a single beta and gamma value)----------#####

### Note this function is used when we have only one value for the disease parameters beta and gamma values. Thus when length(beta)=1 and length(gamma)=1
### For multiple/range of beta and gamma values, we can use the 'ParalleEpicSimOnGraphs' function below

stderror <- function(x) sd(x)/sqrt(length(x))
## This experiment run/simulate a given epidemic process (eg SIR etc) on any given graph/network and calculate some 
## epidemic measures and any reported status at each timesteps
Run_EpicSim_And_Measures<-function(Graphs, nticks=10, beta=0.5,gamma=0.2,propInfected=0.1, initialState=NULL, nInfected=1, useProportion=F, nreps=2, seed=-1,epi_output_file="EpicMeasures.csv",report="i") {
  set.seed(1) 
  ### Definition and initialization of parameters for pathogen simulation
  populate <- list(data.frame(matrix(NA,nrow=(nticks+1),ncol=0)))# simulation results go here
  sickness=list()
  states=list()
  statesOfInterest <- 0:2
  status=list()
  
  #### Definition of final states
  
  Total_susc=list()
  Total_inf=list()
  Total_rec=list()
  Total_states=list()
  
  ### Definition and initialization of parameters for calculating the Epidemic Measures
  TimetoMaxOfReportedState=list()
  StdOfReportedState=list()
  MaxNumOfReportedState=list()
  MeanNumOfReportedState=list()
  QuantileOfReportedState_1st=list()
  QuantileOfReportedState_2nd=list()
  QuantileOfReportedState_3rd=list()
  QuantileOfReportedState_4th=list()
  MedianOfReportedState=list()
  StandarderrorOfReportedState=list()
  Epi_measures=list()
  
  analysis=list()
  status_prefix=list()
  #results=data.frame(matrix(nrow=nreps*length(Graphs),ncol=0))
  
  
  #### Definition and initialization of Graph Prefix
  graphid=list();graphreplicate=list(); graphname=list(); GraphPrefix=list(); analysis=list()
  epi_analysis=list()
  
  m=1
  for (g in Graphs) {
    for (idx in 1:nreps) {
      
      ### Simulate pathogen on a given graph
      print(paste("Simulating pathogen on", g$name))
      sickness<- simPathE(g,nTicks=nticks,beta = beta,gamma=gamma,propInfected=propInfected,initialState=initialState,nInfected=nInfected,useProportion=useProportion)
      GraphColnames <- c(paste(g$type,"s", idx , sep="-"), paste(g$type,"i",idx, sep="-"), paste(g$type,"r",idx, sep="-"))
      column_names <- c(colnames(populate),GraphColnames)
      states[[idx]]=countStates(sickness, states=statesOfInterest)
      colnames(states[[idx]]) <- column_names
      
      
      ## Total number in each compartment after the final timestep
      Total_susc[[idx]]=as.vector(sapply(states[[idx]], tail, 1))[1]
      Total_rec[[idx]]=as.vector(sapply(states[[idx]], tail, 1))[3]
      Total_inf[[idx]]=as.vector(sapply(states[[idx]], tail, 1))[2]+Total_rec[[idx]]
      
    } 
    
    ### Combine by column all the states of the simulated graphs
    populate[[m]]<- do.call(cbind,states)
    
    
    ### Select only the reported state columns needed from each simulated graph
    status[[m]]=populate[[m]][grep(report, names(populate[[m]]))] 
    
    status_prefix[[m]]=t(status[[m]])
    status_prefix.col.names=paste("t", 0:nticks, sep='')
    status_prefix.row.names=1:nreps
    colnames(status_prefix[[m]])=status_prefix.col.names
    row.names(status_prefix[[m]])=status_prefix.row.names
    
    
    #### reporting initial values in each state
    # Init_states[[m]]=cbind(Init_susc,Init_inf,Init_rec)
    # row.names(Init_states)=NULL
    # colnames(Init_states[[m]])=c("Init_susc","Init_inf","Init_rec")
    
    #### reporting total values in each state
    Total_states[[m]]=cbind(Total_susc,Total_inf,Total_rec)
    row.names(Total_states)=NULL
    colnames(Total_states[[m]])=c("Total_susc","Total_inf","Total_rec")
    
    #       m=m+1
    #     }
    #     return(Total_states)    
    # }    
    ### Some Epidemic Measures to calculate on the simulated graphs based on the selected reported status
    print(paste("Calculating epidemic measures on", g$name))
    
    TimetoMaxOfReportedState[[m]]=apply(status[[m]],2,which.max)
    StdOfReportedState[[m]]=apply(status[[m]],2,sd)
    MeanNumOfReportedState[[m]]=apply(status[[m]],2,mean)
    MedianOfReportedState[[m]]=apply(status[[m]],2,median)
    StandarderrorOfReportedState[[m]]=apply(status[[m]],2,stderror)
    MaxNumOfReportedState[[m]]=apply(status[[m]],2,max)
    QuantileOfReportedState_1st[[m]]=apply(status[[m]],2,quantile, probs=c(0.25), na.rm=TRUE)
    QuantileOfReportedState_2nd[[m]]=apply(status[[m]],2,quantile, probs=c(0.50), na.rm=TRUE)
    QuantileOfReportedState_3rd[[m]]=apply(status[[m]],2,quantile, probs=c(0.75), na.rm=TRUE)
    QuantileOfReportedState_4th[[m]]=apply(status[[m]],2,quantile, probs=c(1), na.rm=TRUE)
    Epi_measures=cbind(TimetoMaxOfReportedState[[m]],StdOfReportedState[[m]],MeanNumOfReportedState[[m]], MaxNumOfReportedState[[m]],MedianOfReportedState[[m]],StandarderrorOfReportedState[[m]],QuantileOfReportedState_1st[[m]],QuantileOfReportedState_2nd[[m]],QuantileOfReportedState_3rd[[m]],QuantileOfReportedState_4th[[m]])
    Epi_measures.col.names=c("TimetoMax","StdOf","MeanNumOf","MaxNumOf","MedianNumOf","StderrorOf","25th_quantile","50th_quantile","75th_quantile","100th_quantile")
    colnames(Epi_measures)=paste(Epi_measures.col.names,report,sep = "-")
    
    
    ## Graph Prefix for each of the simulated graphs
    graphname[[m]]=Graphs[[m]]$type
    graphid[[m]]=Graphs[[m]]$id
    graphreplicate[[m]]=c(1:nreps)
    GraphPrefix=cbind(graphname[[m]],graphid[[m]],graphreplicate[[m]])
    colnames(GraphPrefix)=c("GraphName","GraphID","GraphReplicate")
    
    ### Data frame of results for each simulated graph
    # graph_analysis[[m]]=as.data.frame(cbind(GraphPrefix,graphProperties[[idx]]))
    # row.names(graph_analysis[[m]])=1:nreps
    
    
    ### Data frame of results for each simulated disease
    epi_analysis[[m]]=as.data.frame(cbind(GraphPrefix,beta,gamma,Epi_measures,Total_states[[m]],report,status_prefix[[m]]))
    row.names(epi_analysis[[m]])=1:nreps
    
    m=m+1
  }
  
  
  ## epidemic measures 
  epi_results = do.call(rbind,epi_analysis)
  write.csv(as.matrix(epi_results), file=epi_output_file)
  
  return(epi_results)
}


##----------- Parallel multiple simulation of each combination of the disease parameters (multiple beta and gamma values)----------#####

### Note this function is used when we have a range of beta and gamma values. thus length(beta)>1 and length(gamma)>1
### For a single beta and gamma value, we can use the 'Run_EpicSim_And_Measures' function above


#cl=makeCluster(detectCores()-1)
#cl <- parallel::makeCluster(detectCores()-1)
#registerDoParallel(cl)
#plan(cluster, workers=cl)
plan(multisession, workers = 6)

#plan(multisession,workers=20) 


# # registerDoParallel(cl)
#plan(multicore) ## parallelize on local computer
results=list()
ListOfResults=list()
All_results=list()
#ListOfResults=list(data.frame(matrix(NA,nrow=length(BETA)*length(GAMMA),ncol=0)))
start <- Sys.time()

ParalleEpicSimOnGraphs<-function(Graphs, betaVals=betaVals,gammaVals=gammaVals, nreps=nreps,output_file="EpicSimOnGraphs.csv",report="s",nticks=10) {
  #  p = progressor(along=Graphs)
  pathogn=expand.grid(data.frame(betaVals,gammaVals))
  
  for (m in 1:nrow(pathogn)){## Counter to account for the cartesian product of all beta and gamma values
    results[[m]] = future_lapply(seq_along(Graphs), function(j) { 
      # Sys.sleep(6.0-j)
      #  p(sprintf("j=%g", j))
      data.frame(Run_EpicSim_And_Measures(Graphs[j], beta=pathogn$betaVals[m],gamma=pathogn$gammaVals[m],
                                          nticks=nticks,report=report,nreps = nreps))} 
      , future.seed = 0xBEEF, future.chunk.size=1)
    
  }
  
  ListOfResults=do.call(rbind,results)
  All_results=do.call(rbind,ListOfResults)#length of these data frame=length(beta)*length(gamma)*nreps*nSample
  All_results=as.matrix(All_results)
  # ## Produces a .CSV file of the final result
  write.csv(All_results, file=output_file)
  return(All_results)
}





####################---For Simulations---####################################################
########### Plot function for the observed networks with 50 nodes ######  
set.seed(82626)
nsim=10
plotfunc_obs<-function(Name="ER",ID=1,betaVal=.01,gammaVal=.2,nreps=nsim,Data="observed_networks_50_avrgdeg_4.csv",plot_title="",nticks=100,net_size=10){
  data_obs=read.csv(Data,header = T, sep = ",")
  df_obs=data_obs%>% filter(GraphName==Name,GraphID==ID, beta==betaVal,gamma==gammaVal)
  df_obs=df_obs%>% select(c(beta,gamma,t0:paste("t",nticks,sep = "")))
  
  p_obs=data.frame(0:(length(df_obs[grep("t", names(df_obs))])-1),t(df_obs[grep("t", names(df_obs))]))
  colnames(p_obs)= c("Timesteps",paste("Infecteds",1:nreps,sep = "_"))
  df_obs_long_format <- melt(p_obs, id="Timesteps")  # convert to long format
  
  
  plot_obs<- ggplot(data= df_obs_long_format,
                    aes(x=Timesteps, y=value/net_size, colour=variable)) +
    geom_line(show.legend = FALSE)+theme_classic()+labs(title="",tag=Name )+
    xlab("Timesteps (days)")+ylab("Prop-infected")+
    scale_y_continuous(limits = c(0,1))+theme(text = element_text(size = 26),
                                              axis.title = element_text(size = 26),
                                              axis.text.y = element_text(size = 26),
                                              axis.text.x = element_text(size = 26),
                                              legend.position = "none",
                                              plot.title.position = "plot",
                                              plot.tag.position = c(0.1, 0.98))
  #theme(text = element_text(size = 22))    
  #scale_y_continuous(limits = c(0,net_size))
  return(plot_obs)   
}

############----------AVERAGE PROPORTION OF INFECTED FUNCTION------###############
AvrgPropOf<-function(Name="ER",ID=1,betaVal=.01,gammaVal=.2,nreps=nsim,Data="observed_networks_50_avrgdeg_4.csv",plot_title="",nticks=100,net_size=10){
  df_obs=data.frame(read.csv(Data,header = T, sep = ","))
  df_obs=df_obs%>% filter(GraphName==Name,GraphID==ID, beta==betaVal,gamma==gammaVal)
  df_obs=df_obs%>%select(c(GraphName,t0:paste("t",nticks,sep = "")))
  #y=f%>%group_by(GraphName)%>%summarize(colMeans(f[sapply(f, is.numeric)]))
  df=data.frame(colMeans(df_obs[sapply(df_obs, is.numeric)]))
  colnames(df)=c("value")
  df$Timesteps=1:(nticks+1)
  df$GraphName=Name
  df=df%>%select(GraphName, Timesteps,value)

  return(df)   
}


###--Normalized Laplacina function--##
normalized_laplacian=function(Graphs){
  laplacian_matrix(Graphs,normalized = T)
}



###########----------- SUMMARY OF GRAPH FEATURES SYNTHETIC GRAPHS AND DOLPHIN/HYENA--------------

Graphfeatures=function(Name="ER",data="data.csv"){
  df_feat=data.frame(read.csv(data,header = T, sep = ","))
  df_feat=df_feat%>% filter(GraphName==Name)
  return(df_feat)  
}



##############------COMPARING REAL AND THEORETICAL NETWORKs###################
####################---GRAPH DRAWING-----------##########

spectral_drawing<-function(graph_type=G1[[10]]){
  laplac_graph_type<-laplacian_matrix(graph_type,sparse = FALSE)
  Spec_graph_type<-eigen(laplac_graph_type)
  #----Note-on--Laplacian--spectrum--
  #(1)The closer lambda_2 is to 0 the more disconnected the graph is
  #(2)Eigenvectors are orthogonal
  #(3)lambda_1(fiedler-value) is >0 for connected graph and has a multiplicity of exactly 1
  Fiedler_graph_type<-Spec_graph_type$values[length(Spec_graph_type$value)-1]
  FiedlerVec_graph_type<- Spec_graph_type$vectors[,length(Spec_graph_type$values)-1]
  
  #(12)Visualization with smallest eigenvectors
  #if  eigenvalue is greater than 0 graph is connected
  FiedlerVect_1<-FiedlerVec_graph_type
  FiedlerVect_2<- Spec_graph_type$vectors[,length(Spec_graph_type$values)-2]
  # round(FiedlerVect_ER1,digit=2)
  # round(FiedlerVect_ER2,digit=2)
  #Showing the first and second non-trivial eigen vectors
  # par(mfrow=c(1,2))
  # plot(FiedlerVect_ER1);plot(FiedlerVect_ER2)
  
  
  #Reduce arrow size and set arrow mode for clearer plot
  set.seed(2523)#to enable reproducibility of layout
  coords_graph_type<- as.matrix(cbind(round(FiedlerVect_1,digit=2),round(FiedlerVect_2,digit=2)))
  graph_type$layout=coords_graph_type
  #Vertex Options: Color
  #v_col<- rep("grey80", vcount(graph_type))
  V(graph_type)[degree(graph_type, mode="in")>=10]$color <- "red"  #Destinguishing High Degree Nodes as red
  V(graph_type)$size=degree(graph_type, mode = "in")/.15 #because we have wide range, I am dividing by 8 to keep the high in-degree nodes from overshadowing everything else.
  fig=plot(graph_type,vertex.label.font=0.1,
           vertex.label.font=0.1,edge.color="gray",
           vertex.label = NA,vertex.color="orange",vertex.size=5,
           layout=coords_graph_type)
  
  return(fig)          
}


# Helper functions to convert between (row, column) pairs and the index in a list of Cells.
cellCoordsToIndex <- function(i, j, size) {
  return((i-1)*size+j)
}
indexToCellCoords <- function(idx, size) {
  j <- (idx-1) %% size + 1
  i <- (idx + size-1) %/% size
  return(c(i,j))
}

# - end of Helper functions


# for (i in 1:4) {
#   for (j in 1:4) {
#     idx <- cellCoordsToIndex(i,j,4)
#     print(paste("i =", i, "; j = ", j , "index = " , idx ,
#                 "; reverse = ", indexToCellCoords(idx, 4)[1], indexToCellCoords(idx, 4)[2]))
#     print("")
#   }
# }


###Function to create Spatial network ##
Sp.hybrid.net<- function(n=100,r=0.1, makeConnected=TRUE, keepCellsSeparate=FALSE,prob=0.3) {
  # Divide the grid into cells of diagonal length r (so side length r/sqrt(2))
  # All pairs of points in the same cell are within distance r of each other
  # Look around each cell at the points in the 8 adjacent cells
  # Test each for Euclidean distance.
  # Connect all that are close enough; note those that aren't for later addition
  
  # Set up the coordinates:
  v <- data.frame(matrix(nrow=n))
  v$x <- runif(n)# generate x coordinate
  v$y <- runif(n)# generate y coordinate
  v[,1] <- 1:n # first column is the total number of nodes
  colnames(v) <- c("id","x","y")
  
  # put points into bins:
  cellSize <- r/sqrt(2)
  R2 <- r*r
  # create grid size to distribute the x,y points in with the radius. eg gridSize=5
  gridSize <- ceiling(sqrt(2)/r) 
  
  # # Create a vector of grid of cells eg. Cell = 5x5 if gridSize=5 and
  # # assign vertices/nodes to them based on their coordinates: this is done by a method
  Cell = vector(mode="list", gridSize*gridSize)
  ### 1.we divide each point in y and x by the cellsize and 
  ## 2. multiply the resulted y points by grid size and add the resulted x points by 1
  ## 3. Assign the ith node to the jth cell with this method
  ## eg if i=1 then idx=9 and we assign node 1 to cell 9 etc
  for (i in 1:n) {
    idx <- (floor(v$y[i]/cellSize))*gridSize + floor(v$x[i]/cellSize) + 1
    Cell[[idx]] <- cbind(Cell[[idx]], i)
  }
  E <- matrix(nrow=0, ncol=2, byrow=TRUE) # create an empty matrix
  
  # # join all vertices/nodes that are in the same cell since they must be at most distance r apart
  ## Thus if the number of nodes in each cell in >1, then create all possible pairs of nodes with each nodes
  ## Thus (n,2) combinations. since egdes are formed between nodes
  ## eg if cell 22 has node 16,22,26 and 43 the we have 16-22,16-26,16-43,22-26,22-43,26-43 
  for (cell in Cell) {
    if (length(cell) > 1) {
      E <- rbind(E, t(combn(cell, 2)))# form pairs of nodes (edges) with every node in a given cell
    }
  }
  ## Create data frame of three columns with radius,i and j as names
  proximalPairs <- data.frame(matrix(5,ncol=3,nrow=1))
  colnames(proximalPairs) <- c("r","i","j")
  
  if (keepCellsSeparate == FALSE) {
    for (i in 2:gridSize) {
      # first column:
      # join any points in [ ]
      #                    [ ] that are close enough
      thisIdx <- cellCoordsToIndex(1,i, gridSize)
      aboveIdx <- cellCoordsToIndex(1,i-1, gridSize)
      for (a in Cell[[thisIdx]]) {
        for (b in Cell[[aboveIdx]]) {
          d <- (v[a,2]-v[b,2])^2 + (v[a,3]-v[b,3])^2
          if (d < R2) {
            E <- as.matrix(rbind(E,c(a,b)))
          } else {
            proximalPairs <- rbind(proximalPairs, c(d, a, b))
          }
        }
      }
      # first row:
      # join any points in [ ][ ] that are close enough.
      thisIdx <- cellCoordsToIndex(i,1, gridSize)
      leftIdx <- cellCoordsToIndex(i-1,1, gridSize)
      for (a in Cell[[thisIdx]]) {
        for (b in Cell[[leftIdx]]) {
          d <- (v[a,2]-v[b,2])^2 + (v[a,3]-v[b,3])^2
          if (d < R2) {
            E <- as.matrix(rbind(E,c(a,b)))
          } else {
            proximalPairs <- rbind(proximalPairs, c(d, a, b))
          }
        }
      }
      
      for (j in 2:gridSize) {
        # check all neighbours above and to the left
        #
        # aboveLeftIdx      aboveIdx
        #             XX       ||
        #     leftIdx    ==  thisIdx
        thisIdx <- cellCoordsToIndex(i,j,gridSize)
        aboveIdx <- cellCoordsToIndex(i,j-1,gridSize)
        
        # this + above:
        for (a in Cell[[thisIdx]]) {
          for (b in Cell[[aboveIdx]]) {
            d <- (v[a,2]-v[b,2])^2 + (v[a,3]-v[b,3])^2
            if (d < R2) {
              E <- as.matrix(rbind(E,c(a,b)))
            } else {
              proximalPairs <- rbind(proximalPairs, c(d, a, b))
            }
          }
        }
        
        # this + left:
        leftIdx <- cellCoordsToIndex(i-1,j,gridSize)
        for (a in Cell[[thisIdx]]) {
          for (b in Cell[[leftIdx]]) {
            d <- (v[a,2]-v[b,2])^2 + (v[a,3]-v[b,3])^2
            if (d < R2) {
              E <- as.matrix(rbind(E,c(a,b)))
            } else {
              proximalPairs <- rbind(proximalPairs, c(d, a, b))
            }
          }
        }
        
        # this + aboveLeft
        aboveLeftIdx <- cellCoordsToIndex(i-1,j-1,gridSize)
        for (a in Cell[[thisIdx]]) {
          for (b in Cell[[aboveLeftIdx]]) {
            d <- (v[a,2]-v[b,2])^2 + (v[a,3]-v[b,3])^2
            if (d < R2) {
              E <- as.matrix(rbind(E,c(a,b)))
            } else {
              proximalPairs <- rbind(proximalPairs, c(d, a, b))
            }
          }
        }
        
        # left and above:
        for (a in Cell[[leftIdx]]) {
          for (b in Cell[[aboveIdx]]) {
            d <- (v[a,2]-v[b,2])^2 + (v[a,3]-v[b,3])^2
            if (d < R2) {
              E <- as.matrix(rbind(E,c(a,b)))
            } else {
              proximalPairs <- rbind(proximalPairs, c(d, a, b))
            }
          }
        }
      }
    }
  }
  ### make the spatial gragh
  G = igraph::make_empty_graph(n, directed=FALSE)
  G <- igraph::add_edges(G, t(E), directed=FALSE) # addind edges based on spatial construction
  #   ### rewiring edges based on some probabibility
  G %>% igraph::rewire(each_edge(p = prob, loops = FALSE,multiple = FALSE))
  vertex_attr(G)$id <- 1:n
  G$layout <- as.matrix(v[,2:3])
  G$name <- "Spatial"
  G$type <- "Spatial"
  G$id <- "1"
  if (makeConnected == T) {
    return(G)
  }
  if (is.connected(G)) {
    return(G)
  }
  # sort the pairs of coordinates in ascending order of distance apart
  proximalPairs <- proximalPairs[-1,]
  rownames(proximalPairs) <- seq(1,dim(proximalPairs)[1])
  indx <- order(proximalPairs$r)
  
  comps <- components(G, mode="weak")
  
  for (nu in indx) { # in increasing order of edge length
    i <- proximalPairs$i[nu]
    j <- proximalPairs$j[nu]
    
    # Consider vertices i and j for joining:
    if (comps$membership[i] != comps$membership[j]) { # if two vertices are not in the same component
      G <- G + edge(i, j, directed=F) # then join them
      comps <- components(G, mode="weak")
    } else {
      # i and j are in the same component.  Moving on...
    }
    if (is.connected(G) == T) {
      break # we are done!
    }
  }
  return(idx)
}

#Sp.hybrid.net(n=60,r=0.5,p=0.1, makeConnected=FALSE, keepCellsSeparate=FALSE)

# Rewire.graph=function(G,p=0.2){
#   ### rewiring edges based on some probabibility
#   G %>% rewire(each_edge(p = p, loops = FALSE,multiple = FALSE))
# }
#Sp.hybrid.net(n=10,r=0.1, makeConnected=FALSE, keepCellsSeparate=FALSE,prob=0.3)
#makeSpatialHybrid(node.size=15,Radius=0.8,prob=0.5)
##########################################################################################################
# Create a list of graphs for simulation experiment
##########################################################################################################
makeSpatialHybrid<- function(node.size=25,Radius=0.8,prob=0.5) {
  Graphs = NULL # set up the list of Graphs first
  G=NULL
  #initgraph= list()
  i= 1
  print("Creating Spatial-Hybrid networks")
  G[[i]]=Sp.hybrid.net(n=node.size,r=Radius,p=prob, makeConnected=TRUE, keepCellsSeparate=FALSE)
  Graphs[[i]]=induced.graph(G[i])
  Graphs[[i]]$type = "Spatial-hybrid"
  Graphs[[i]]$id = "1"
  Graphs[[i]]$name="Spatial-hydrid"
  i <- i+1
  return(G)
}

#U=makeSpatialHybrid(node.size=300,Radius=0.05,prob=0.2)
#U


###---Simulate--spatial--networks---##
simulate.spatial<-function(N=50,radius=0.4,nsim=100){
  spatial.graph=NULL
  for (i in 1:nsim){
    spatial.graph[[i]]=makeSpatialGraphs(node.size=N,Radius=radius)
  }
  return(spatial.graph)
}

###---Simulate spatial--hybrid--networks---##
simulate.spatialhybrid<-function(N=50,radius=0.4,p=0.3,nsim=100){
  spatial.graph=NULL
  for (i in 1:nsim){
    spatial.graph[[i]]=makeSpatialHybrid(node.size=N,Radius=radius,prob=p)
  }
  return(spatial.graph)
}


####----Data--Generation----###
# networks.hybrid=function(R=c(0.3,0.4,0.6,0.05),p=0.2,n=100){
#   net=NULL;data=NULL
#   k=1
#   for (i in 1:length(R)){
#     net[i]=makeSpatialHybrid(node.size =n,Radius = R[i],prob=p)
#     data=RunSimOnGraphFeatures(net,nreps = 1)
#     k=k+1
#   }
#   df=cbind(R,p,data)
#   return(df)
# }

# m15.100=networks.hybrid(R=rep(0.15,100),p=rep(0.8,100),n=100)
# write.csv(m15.100,"m15-100.csv")
# 
# m25.100=networks.hybrid(R=rep(0.2,100),p=rep(0.8,100),n=100)
# write.csv(m25.100,"m25-100.csv")
# 
# m31.100=networks.hybrid(R=rep(0.3,100),p=rep(0.2,100),n=100)
# write.csv(m31.100,"m31-100.csv")
# 
# m41.100=networks.hybrid(R=rep(0.4,100),p=rep(0.2,100),n=100)
# write.csv(m41.100,"m41-100.csv")
# 
# m51.100=networks.hybrid(R=rep(0.5,100),p=rep(0.2,100),n=100)
# write.csv(m51.100,"m51-100.csv")
# 
# m61.100=networks.hybrid(R=rep(0.6,100),p=rep(0.2,100),n=100)
# write.csv(m61.100,"m61-100.csv")
# 
# m71.100=networks.hybrid(R=rep(0.7,100),p=rep(0.2,100),n=100)
# write.csv(m71.100,"m71-100.csv")
# 
# m81.100=networks.hybrid(R=rep(0.8,100),p=rep(0.2,100),n=100)
# write.csv(m81.100,"m81-100.csv")
# 
# m91.100=networks.hybrid(R=rep(0.9,100),p=rep(0.2,100),n=100)
# write.csv(m91.100,"m91-100.csv")

# m10.5.100=networks.hybrid(R=rep(1,100),p=rep(0.8,100),n=100)
# write.csv(m10.5.100,"m10-5-100.csv")


# source("SPATIAL-PIPELINE.R")
# 
# set.seed(123)
# r1 <- round(runif(100, 0.01, 0.5),digit=3)
# r2 <- round(runif(100, 0.5, 1),digits = 3)
# r3 <- round(runif(300, 0.07, 0.3),digits=3)
# r4 <- round(runif(300, 0.4, 0.8),digits = 3)
# r5 <- round(runif(500, 0.02, 0.05),digits=3)
# r6 <- round(runif(500, 0.4, 0.9),digits=3)
# r7 <- round(runif(750, 0.05, 0.3),digits = 3)
# r8 <- round(runif(750, 0.02, 0.4),digits = 3)
# r9 <- round(runif(1000, 0.01, 0.03),digits = 3)
# r10 <-round(runif(1000, 0.4, 0.8),digits = 3)
# 
# 
# R=r3
# node.num=300
# 
# 
# net=NULL
# for (i in 1:length(R)){
#   # net=fastSpatialNetwork(n=node.num,r=R[[i]],makeConnected=TRUE, keepCellsSeparate=FALSE)
#   net[i]=makeSpatialGraphs(node.size =node.num,Radius=R[i])
# }
# 
# del_vertices=function(G){
#   Isolated = which(degree(G)==0)
#   Net=igraph::delete.vertices(G, Isolated)
#   return(Net)  
# } 
# 
# net.1=lapply(net,del_vertices)
# 
# Data=RunSimOnGraphFeatures(net.1,nreps = 1)
# 
# DF.DATA=cbind(R,Data)
# 
# 
# df100.1=write.csv(DF.DATA,"r.100.1.csv")
# df100.2=write.csv(DF.DATA,"r.100.2.csv")
# df300.1=write.csv(DF.DATA,"r.300.1.csv")
# df300.2=write.csv(DF.DATA,"r.300.2.csv")
# df500.1=write.csv(DF.DATA,"r.500.1.csv")  
# df500.2=write.csv(DF.DATA,"r.500.2.csv")
# df750.1=write.csv(DF.DATA,"r.750.1.csv")
# df750.2=write.csv(DF.DATA,"r.750.2.csv")
# df1000.1=write.csv(DF.DATA,"r.1000.1.csv")
# df1000.2=write.csv(DF.DATA,"r.1000.2.csv")
# 
# 
# head(df2)
# 
# data=df2[sample(1:nrow(df2)),]
#
# node.num=100
# # Drop variables
# df.data<-data %>%dplyr::filter(order==node.num)%>% 
#   dplyr::select(-c(GraphID,order, GraphReplicate, order,max_component,minDegree,
#                    threshold,maxDegree,GraphName,connected,size))
# 
# dim(df.data) 
# head(df.data)


# #######--Feature--Processing---######
# feature_data<- data_df %>%
#   dplyr::select(-c(GraphID,order, GraphReplicate, order,max_component,minDegree,
#                    threshold,maxDegree,X,GraphName,connected,maxDegree,minDegree,
#                    max_component,order,size)) 
# 
# 
# ####--Correlation---coefficient---####
# features.data=cor(feature_data, method = "pearson", use = "complete.obs")
# round(features.data, 2)
# 
# ####--Correlation---matrix---with--significance--(p)--level####
# #rcorr(features.data, type = c("pearson","spearman"))
# feature.matrix=rcorr(as.matrix(features.data))
# 
# # Extract the correlation coefficients
# feature.matrix$r
# # Extract p-values
# feature.matrix$P
# 
# # ++++++++++++++++++++++++++++
# # flattenCorrelationMatrix
# # ++++++++++++++++++++++++++++
# # cormat : matrix of the correlation coefficients
# # pmat : matrix of the correlation p-values
# flattenCorrMatrix <- function(cormat, pmat) {
#   ut <- upper.tri(cormat)
#   data.frame(
#     row = rownames(cormat)[row(cormat)[ut]],
#     column = rownames(cormat)[col(cormat)[ut]],
#     cor  =(cormat)[ut],
#     p = pmat[ut]
#   )
# }
# 
# 
# feature.matrix.flatten<-feature.matrix
# flattenCorrMatrix(feature.matrix.flatten$r, feature.matrix.flatten$P)
# 
# # Visualize correlation matrix
# # There are different ways for visualizing a correlation matrix in R software :
# # 
# # symnum() function
# # corrplot() function to plot a correlogram
# # scatter plots
# # heatmap
# # symnum(x, cutpoints = c(0.3, 0.6, 0.8, 0.9, 0.95),
# #        symbols = c(" ", ".", ",", "+", "*", "B"),
# #        abbr.colnames = TRUE)
# 
# symnum(features.data, abbr.colnames = FALSE)
# 
# ####-------Correlogram------
# corrplot(features.data, type = "upper", order = "hclust", 
#          tl.col = "black", tl.srt = 45)
# 
# # Insignificant correlation are crossed
# corrplot(feature.matrix$r, type="upper", order="hclust", 
#          p.mat = feature.matrix$P, sig.level = 0.01, insig = "cross")
# # Insignificant correlations are leaved blank
# corrplot(feature.matrix$r, type="upper", order="hclust", 
#          p.mat = feature.matrix$P, sig.level = 0.01, insig = "blank")
# 
# #####------chart of a correlation matrix.
# chart.Correlation(feature_data, histogram=TRUE, pch=19)
# 
# ####------a chart of a correlation matrix--##
# col<- colorRampPalette(c("blue", "white", "red"))(20)
# heatmap(x = features.data, col = col, symm = TRUE)
# 
# 
# ###+++++++++++=+++Feature---Processing--with---Boruta+++++++++++++++
# # feat_data1<- data_df %>%
# #   dplyr::select(-c(GraphID,order, GraphReplicate, order,max_component,minDegree,
# #             threshold,maxDegree,X,GraphName,connected,maxDegree,minDegree,
# #             max_component,order,size,diameter,minCut)) 
# # 
# # boruta=Boruta(r~., data=feat_data1, doTrace=2,maxRuns=500)
# # plot(boruta,las=2, cex.axis=0.7)
# # plotImpHistory(boruta)
# # attStats(boruta)
# # 
# # 
# # Ten.feat=TentativeRoughFix(boruta)
# # plot(Ten.feat)
# ```
# 
# **PCA**
#   ```{r pressure, echo=FALSE, include=FALSE}
# 
# df.data=read.csv("r.100.2.csv",header = T, sep=",")
# data<- df.data[sample(1:nrow(df.data)), ]
# #data$R=r1
# head(data)
# 
# node.num=100
# # Drop variables
# df.data<-data %>%dplyr::filter(order==node.num)%>% 
#   dplyr::select(-c(GraphID,order, GraphReplicate, order,max_component,minDegree,
#                    threshold,maxDegree,GraphName,connected,X))
# 
# 
# head(df.data)
# dim(df.data) 
# 
# ### find correlations
# corel.data=df.data.corr()
# 
# PCA.1=princomp(df.data[,c(-1)],cor = T)
# sum_pca=summary(PCA, loadings=T)
# 
# df.pca=write.csv(sum_pca$loadings,"PCA.csv")
# 
# PCA.2=prcomp(df.data[,-1], scale=T)
# 
# PCA.2
# 
# summary(PCA.2)
# 
# plot(PCA.2, type="l")
# 
# biplot(PCA.2,scale=0)
# 
# str(PCA.2)
# 
# PCA.2$x
# 
# df=cbind(df.data,PCA.2$x[,1:2])
# 
# head(df)
# 
# ggplot(df, aes(PC1,PC2, col=R,fill=R))+
#   stat_ellipse(geom="polygon",col="black",alpha=0.5)+
#   geom_point(shape=21,col="black")
# 
# 
# cor(df.data[,-1],df[,14:15])
