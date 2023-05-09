if(!require(igraph)) {
  install.packages("igraph")
  library(igraph)
}
if(!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}


if(!require(tidyverse)) {
  install.packages("tidyverse")
  library(tidyverse)
}


if(!require(EpiModel)) {
  install.packages("EpiModel")
  library(EpiModel)
}


if(!require(dplyr)) {
  install.packages("dplyr")
  library(dplyr)
}


if(!require(future.apply)) {
  install.packages("future.apply")
  library(future.apply)
}


if(!require(doParallel)) {
  install.packages("doParallel")
  library(doParallel)
}



if(!require(progressr)) {
  install.packages("progressr")
  library(progressr)
}


if(!require(plotly)) {
  install.packages("plotly")
  library(plotly)
}


if(!require(reshape2)) {
  install.packages("reshape2")
  library(reshape2)
}


if(!require(ggpubr)) {
  install.packages("ggpubr")
  library(ggpubr)
}


library(dplyr)
library(reshape2)
library(plotly)


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




for (i in 1:4) {
  for (j in 1:4) {
    idx <- cellCoordsToIndex(i,j,4)
    print(paste("i =", i, "; j = ", j , "index = " , idx ,
                "; reverse = ", indexToCellCoords(idx, 4)[1], indexToCellCoords(idx, 4)[2]))
    print("")
  }
}

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
  G$name <- "Spatial graph"
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



##########################################################################################################
# Create a list of graphs for simulation experiment
##########################################################################################################
makeGraphs1 <- function(nSamples=10, order=25) {
  
  graphTypes <- c("Erdos-Renyi", "Small World", "Scale Free", "Square Lattice", "Spatial")
  
  Graphs <- list(mode="any", nSamples * length(graphTypes)) # set up the list of Graphs first
  initgraph= list(mode="any", nSamples * length(graphTypes))
  
  print("Creating Erdos-Renyi graphs")
  i <- 1
  for (j in 1:nSamples) {
    Graphs[[i+j-1]] <- sample_gnm(n=order, m=(2*order), directed=FALSE, loops=FALSE)
    initgraph=Graphs[[i+j-1]]
    components<- igraph::clusters(initgraph, mode="weak")
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
    Graphs[[i+j-1]] <- simplify(sample_smallworld(dim=2, size=round(sqrt(order)), nei=round(log(order)/4), p=0.2))
    initgraph=Graphs[[i+j-1]]
    components<- igraph::clusters(initgraph, mode="weak")
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
    Graphs[[i+j-1]] <- sample_pa(n=order, power=1, m=2, directed=FALSE, algorithm="psumtree")
    initgraph=Graphs[[i+j-1]]
    components<- igraph::clusters(initgraph, mode="weak")
    biggest_cluster_id<- which.max(components$csize)
    # # ids
    vert_ids<- V(initgraph)[components$membership== biggest_cluster_id]
    # # subgraph
    Graphs[[i+j-1]]=igraph::induced_subgraph(initgraph, vert_ids)
    Graphs[[i+j-1]]$type <- "SF"
    Graphs[[i+j-1]]$id <- j
  }
  i <- i+nSamples
  # 
  # print("Creating spatial networks using naive algorithm")
  # for (j in 1:nSamples) {
  #   Graphs[[i+j-1]] <- makeSpatialNetwork(n = order, r = 0.1, makeConnected=T)
  #   Graphs[[i+j-1]]$type <- "Sp"
  #   Graphs[[i+j-1]]$id <- as.character(j)
  # }
  # i <- i+nSamples
  
  print("Creating a single square lattice graph (they're all the same!)")
  Graphs[[i]] <- make_lattice(length = floor(sqrt(order)),dim=2)
  Graphs[[i]]$type <- "Lat"
  Graphs[[i]]$id <- 1
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

##########################################################################################################
# Calculate a range of graph values from diameter, degree distribution, through to 
# modularilty via estimating community structure and spectral analysis.
##########################################################################################################
# library(RSpectra)
# library(Rcpp)

normalized_laplacian=function(Graphs){
    laplacian_matrix(Graphs,normalized = T)
}

calcGraphFeatures1 <- function(Graphs=NULL) {
  
  features <- c(
    "order",                    # number of vertices
    "size",                     # number of edges
    "connected",                # True / False
    "max_component",            # maximum component size (=order iff the graph is connected)
    "minDegree",                # minimum degree of any vertex
    "maxDegree",                # maximum degree of any vertex
    "aveDegree",                # average degree of any vertex
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
  
  df$order = as.numeric(lapply(Graphs, gorder))
  df$size = as.numeric(lapply(Graphs, gsize))
  df$connected = as.numeric(lapply(Graphs, is.connected))
  df$minCut = as.numeric(lapply(Graphs, min_cut))
  df$diameter = as.numeric(lapply(Graphs, diameter))
  df$transitivity = as.numeric(lapply(Graphs, transitivity))
  
  # stuff that needs interim things:
  degrees = lapply(Graphs, degree)
  df$minDegree = as.numeric(lapply(degrees, min))
  df$maxDegree = as.numeric(lapply(degrees, max))
  df$aveDegree = as.numeric(lapply(degrees, mean))
  
  # stuff that lapply doesn't like so has to be done in a loop:
  communities <- lapply(Graphs, cluster_leading_eigen)
  Adj <- lapply(Graphs, as_adjacency_matrix)
  L <- lapply(Graphs, laplacian_matrix)
  Norm_Lap<-lapply(Graphs, normalized_laplacian)
  for (i in 1:length(Graphs)) {
    if (is.null(Graphs[[i]]$type)) { Graphs[[i]]$type = "untyped" }
    df$modularity[i] <- modularity(communities[[i]])
    df$spectral_radius[i] <- eigen(Adj[[i]], symmetric=TRUE, only.values=TRUE)$values[1]
    df$FiedlerValue[i] <- eigen(L[[i]], symmetric=TRUE, only.values=TRUE)$values[1]
    df$Normalized_FiedlerValue[i] <- eigen(Norm_Lap[[i]], symmetric=TRUE, only.values=TRUE)$values[1]
    
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
      graphProperties[[reps]] <- calcGraphFeatures1(Graphs[g])
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


##Test
#d=makeGraphs1(nSamples = 2,order=36)
RunSimOnGraphFeatures(h,nreps = 2)


##----------- Multiple simulation of the given disease parameters(a single beta and gamma value)----------#####

### Note this function is used when we have only one value for the disease parameters beta and gamma values. Thus when length(beta)=1 and length(gamma)=1
### For multiple/range of beta and gamma values, we can use the 'ParalleEpicSimOnGraphs' function below

#stderror <- function(x) sd(x)/sqrt(length(x))
## This experiment run/simulate a given epidemic process (eg SIR etc) on any given graph/network and calculate some 
## epidemic measures and any reported status at each timesteps

EpicMeasures3<-function(Graphs, nticks=10, beta=0.5,gamma=0.1,propInfected=0.1, initialState=NULL, nInfected=1, useProportion=F, nreps=2, seed=-1) {
  
  set.seed(1) 
  
  ### Definition and initialization of parameters for pathogen simulation
  graphProperties=list()
  populate <- list(data.frame(matrix(NA,nrow=(nticks+1),ncol=0)))# simulation results go here
  sickness=list()
  states=list()
  statesOfInterest <- 0:2
  Infecteds_only=list()
  
  ### Definition and initialization of parameters for calculating the Epidemic Measures
  TimetoMaxInfectn=list()
  StdOfInfectds=list()
  MaxNumOfInfectds=list()
  MeanNumOfInfectds=list()
  Epi_measures=list()
  
  ### Definition and initialization of Graph Prefix
  graphid=list()
  graphreplicate=list()
  graphname=list()
  GraphPrefix=list()
  
  ### Variable that contain data frame of analysis on each of the simulated graphs
  analysis=list()
  #results=data.frame(matrix(nrow=nreps*length(Graphs),ncol=0))
  
  for (g in 1:length(Graphs)) {
    for (idx in 1:nreps) {
      ### Calculate the graph features for each simulated graph of all the synthetic networks
      graphProperties[[idx]] <- calcGraphFeatures1(Graphs[g])
      
      ### Simulate pathogen on a given graph
      print(paste("Simulating pathogen on", Graphs[[g]]$name))
      sickness<- simPathE(Graphs[[g]],nTicks=nticks,beta = beta,gamma=gamma,propInfected=propInfected,initialState=initialState,nInfected=nInfected,useProportion=useProportion)
      GraphColnames <- c(paste(Graphs[[g]]$type,"S", idx , sep="-"), paste(Graphs[[g]]$type,"I",idx, sep="-"), paste(Graphs[[g]]$type,"R",idx, sep="-"))
      column_names <- c(colnames(populate),GraphColnames)
      states[[idx]]=countStates(sickness, states=statesOfInterest)
      colnames(states[[idx]]) <- column_names
      
      
    } 
    
    ### Combine by column all the states of the simulated graphs
    populate[[g]]<- do.call(cbind,states)
    
    ### Select only the infected columns from each simulated graph
    Infecteds_only[[g]]=populate[[g]][grep("I", names(populate[[g]]))] 
    
    print(paste("Calculating epidemic measures on", Graphs[[g]]$name))
    
    ### Some Epidemic Measures to calculate on the simulated graphs
    TimetoMaxInfectn[[g]]=apply(Infecteds_only[[g]],2,which.max)
    StdOfInfectds[[g]]=apply(Infecteds_only[[g]],2,sd)
    MeanNumOfInfectds[[g]]=apply(Infecteds_only[[g]],2,mean)
    MaxNumOfInfectds[[g]]=apply(Infecteds_only[[g]],2,max)
    Epi_measures=cbind(TimetoMaxInfectn[[g]],StdOfInfectds[[g]],MeanNumOfInfectds[[g]], MaxNumOfInfectds[[g]])
    colnames(Epi_measures)=c("TimetoMaxInfectn","StdOfInfectds","MeanNumOfInfectds","MaxNumOfInfectds")
    
    
    ### Graph Prefix for each of the simulated graphs
    graphname[[g]]=Graphs[[g]]$type
    graphid[[g]]=Graphs[[g]]$id
    graphreplicate[[g]]=c(1:nreps)
    GraphPrefix=cbind(graphname[[g]],graphid[[g]],graphreplicate[[g]])
    colnames(GraphPrefix)=c("GraphName","GraphID","GraphReplicate")
    
    ### Data frame of results for each simulated graph
    analysis[[g]]=as.data.frame(cbind(GraphPrefix,beta,gamma,Epi_measures,graphProperties[[idx]]))
    row.names(analysis[[g]])=1:nreps 
    
    
  }
  results = do.call(rbind,analysis)
  
  #write.csv(results, file=simulationOutputFile)
  
  return(results)
}



##################################################################################
### Epidemic Simulation on any BETA AND GAMMA unique combinations on Graphs
##################################################################################
#betaVals=c(.01,.1,.5)
#gammaVals=c(.1,.2,.5)

## This experiment run/simulate a given epidemic process (eg SIR etc) on any given graph/network and calculate some 
## epidemic measures and graph features for each of the simulates graphs with any range of values of 
##transmission rate (beta) and recovery rate (gamma)

results=list()
ListOfResults=list(data.frame(matrix(NA,nrow=length(BETA)*length(GAMMA),ncol=0)))

RunEpicSimOnGraphs<-function(Graphs, betaVals, gammaVals,nreps=nreps,output_file="EpicSimOnGraphs.csv") {
  
  for(g in 1:length(Graphs)){
    k=1 ## Counter to account for the cartesian product of all beta and gamma values
    for (Beta in betaVals){
      for (Gamma in gammaVals){
        results[[k]]= cbind(data.frame(EpicMeasures3(Graphs[g], beta=Beta,gamma=Gamma,nreps = nreps)))
        k = k+1
      }
    }
    ListOfResults[[g]]=do.call(rbind,results)
    All_results=do.call(rbind,ListOfResults)#length of these data frame=length(beta/gamma)*nreps*nSample
    
    ## Produces a .CSV file of the final result
    write.csv(All_results, file=output_file)
  }
  return( All_results)
  
}



############------Make---synthetic--Graphs to compare to dolphin----##########
makeSimGraphs_dolph<- function(nSamples=10, order=25) {
  
  graphTypes <- c("Erdos-Renyi", "Small World", "Scale Free", "Square Lattice", "Spatial")
  
  Graphs <- list(mode="any", nSamples * length(graphTypes)) # set up the list of Graphs first
  initgraph=list(mode="any", nSamples * length(graphTypes))
  
  print("Creating Erdos-Renyi graphs")
  i <- 1
  for (j in 1:nSamples) {
    Graphs[[i+j-1]] <-erdos.renyi.game(vcount(ggk1),ecount(ggk1) , type = "gnm" , directed = F , loops = F)
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
    Graphs[[i+j-1]] <-sample_smallworld(dim=2, size=round(sqrt(vcount(ggk1))),nei=1,p=0.2) 
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
    Graphs[[i+j-1]] <- sample_pa(n=vcount(ggk1), power=2, m=2, directed=FALSE, algorithm="psumtree")
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
    Graphs[[i+j-1]] <- fastSpatialNetwork(n = vcount(ggk1), r = 0.1795, makeConnected=T)
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
  Graphs[[i]] <- make_lattice(length = ceiling(sqrt(vcount(ggk1))),dim=2,nei = 1)
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


# A1=makeSimGraphs_dolph(order = vcount(ggk1), nSamples = 1)
# RunSimOnGraphFeatures(A1,nreps = 1)

############------Make---synthetic--Graphs to compare to Hyena----##########
makeSimGraphs_hyena<- function(nSamples=10, order=25) {
  
  graphTypes <- c("Erdos-Renyi", "Small World", "Scale Free", "Square Lattice", "Spatial")
  
  Graphs <- list(mode="any", nSamples * length(graphTypes)) # set up the list of Graphs first
  initgraph=list(mode="any", nSamples * length(graphTypes))
  
  print("Creating Erdos-Renyi graphs")
  i <- 1
  for (j in 1:nSamples) {
    Graphs[[i+j-1]] <-erdos.renyi.game(vcount(gg1),ecount(gg1) , type = "gnm" , directed = F , loops = F)
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
    Graphs[[i+j-1]] <-sample_smallworld(dim=2, size=round(sqrt(vcount(gg1))),nei=4,p=0.1) 
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
    Graphs[[i+j-1]] <- sample_pa(n=vcount(gg1), power=2, m=22, directed=FALSE, algorithm="psumtree")
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
    Graphs[[i+j-1]] <- fastSpatialNetwork(n = vcount(gg1), r = 0.8, makeConnected=T)
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
  Graphs[[i]] <- make_lattice(length = ceiling(sqrt(vcount(gg1))),dim=2,nei = 6)
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



# A2=makeSimGraphs_hyena(order = vcount(gg1), nSamples = 1)
# RunSimOnGraphFeatures(A2,nreps = 1)

############------Dolphin----------######
kg <- read.table("mammalia-dolphin-social.edges")
kg1=graph_from_data_frame(as.matrix(kg),directed=FALSE)
ggk1=simplify(kg1)
ggk1$type="Dolphin"
ggk1$id="1"
G=list(ggk1)

ecount(ggk1)
vcount(ggk1)

#######-----------------------Hyena----------------##############
dd <- read.table("mammalia-hyena-networkc.edges")
gg=graph_from_data_frame(as.matrix(dd),directed=FALSE)
gg1=simplify(gg)
gg1$type="Hyena"
gg1$id="1"


ecount(gg1)
vcount(gg1)









#################----Spectra based distances----------############################

set.seed(120)


library(RSpectra)
# k=50
# b1=eigs(Adj_matrix_of_Graphs[[1]],k)$values
# 
# b2=eigs(Adj_matrix_of_Graphs[[2]],k)$values
# 
# b3=sqrt(sum(b1-b2)^2)
normalized_laplacian=function(Graphs){
  laplacian_matrix(Graphs,normalized = T)
}

spectral_distance=function(Graphs=Graphs_for_spectral_distance,name="dolphin",k=50){
  Adj_matrix_of_Graphs <- lapply(Graphs, as_adjacency_matrix)
  Lap_matrix_of_Graphs <- lapply(Graphs, laplacian_matrix)
  Norm_Lap_matrix_of_Graphs<-lapply(Graphs, normalized_laplacian)

 adj_graph=Adj_matrix_of_Graphs
 k=50;m=1;a3=list()
  for (i in adj_graph){
    for (j in adj_graph){
      a1=eigs(i,k)$values
      a2=eigs(j,k)$values
      a3[[m]]=sqrt(sum(a1-a2)^2)
      m=m+1
    }
    adj_spec_dist=do.call(rbind,a3)
    adj_spec_dist=unique( adj_spec_dist)
    adj_spec_dist= adj_spec_dist[-1]
 }
 
 
 lap_graph=Lap_matrix_of_Graphs
 m=1;lap3=list()
 for (i in lap_graph){
   for (j in lap_graph){
     lap1=eigs(i,k)$values
     lap2=eigs(j,k)$values
     lap3[[m]]=sqrt(sum(lap1-lap2)^2)
     m=m+1
   }
   lap_spec_dist=do.call(rbind,lap3)
   lap_spec_dist=unique(lap_spec_dist)
   lap_spec_dist=lap_spec_dist[-1]
 }
 
 
 
 norm_lap_graph=Norm_Lap_matrix_of_Graphs
 m=1;norm_lap3=list()
 for (i in norm_lap_graph){
   for (j in norm_lap_graph){
     norm_lap1=eigs(i,k)$values
     norm_lap2=eigs(j,k)$values
     norm_lap3[[m]]=sqrt(sum(norm_lap1-norm_lap2)^2)
     m=m+1
   }
   norm_lap_spec_dist=do.call(rbind,norm_lap3)
   norm_lap_spec_dist=unique( norm_lap_spec_dist)
   norm_lap_spec_dist=norm_lap_spec_dist[-1]
   
 }
 
 all_spectral_distance=as.data.frame(cbind(adj_spec_dist,lap_spec_dist,norm_lap_spec_dist))
 all_spectral_distance$Compared_networks=c("ER-SW","ER-SF","ER-SP","ER-LAT",
                                         paste("ER",sep = "-",name),"SW-SF","SW-SP","SW-LAT",
                                         paste("SW",sep = "-",name),"SF-SP","SF-LAT",paste("SF",sep = "-",name),
                                         "SP-LAT", paste("SP",sep = "-",name),paste("LAT",sep = "-",name))
 all_spectral_distance <- all_spectral_distance[, c("Compared_networks","adj_spec_dist","lap_spec_dist",
                                                    "norm_lap_spec_dist")]
 

 return(all_spectral_distance)
}

###---DOlPHIN----####
set.seed(120)
dolph11=makeSimGraphs_dolph(nSamples = 50,order = vcount(ggk1))
dolph11$"real"=ggk1
Graphs_for_spectral_distance_dolphin=list(dolph11[[2]],dolph11[[51]],dolph11[[102]],dolph11[[185]],
                                  dolph11[[201]],dolph11[[202]])


Dolphinspectral=spectral_distance(Graphs=Graphs_for_spectral_distance_dolphin,name="DOLPHIN",k=50)
data_spect=write.csv(Dolphinspectral,'C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi/DolphinSpectralFeat.csv')



###---HYENA----####
set.seed(120)
h11=makeSimGraphs_hyena(nSamples = 50,order = vcount(gg1))
h11$"real"=gg1
Graphs_for_spectral_distance_hyena=list(h11[[2]],h11[[51]],h11[[102]],h11[[196]],
                                        h11[[201]],h11[[202]])

Hyenaspectral=spectral_distance(Graphs=Graphs_for_spectral_distance_hyena,name="HYENA",k=30)

data_spect=write.csv(all_spectral_distance,'C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi/HyenaSpectralFeat.csv')





Graphs_for_spectral_distance_hyena=list(h11[[2]],h11[[51]],h11[[102]],h11[[196]],
                                        h11[[201]],h11[[202]])

Adj_matrix_of_Graphs_hyena <- lapply(Graphs_for_spectral_distance_hyena, as_adjacency_matrix)
Lap_matrix_of_Graphs_hyena <- lapply(Graphs_for_spectral_distance_hyena, laplacian_matrix)
Norm_Lap_matrix_of_Graphs_hyena<-lapply(Graphs_for_spectral_distance_hyena, normalized_laplacian)

name="HYENA"
adj_graph_hyena=Adj_matrix_of_Graphs_hyena
k=30;m=1;a3=list()
for (i in adj_graph_hyena){
  for (j in adj_graph_hyena){
    a1=eigs(i,k)$values
    a2=eigs(j,k)$values
    a3[[m]]=sqrt(sum(a1-a2)^2)
    m=m+1
  }
  adj_spec_dist=do.call(rbind,a3)
  adj_spec_dist=unique( adj_spec_dist)
  adj_spec_dist= adj_spec_dist[-1]
}

lap_graph_hyena=Lap_matrix_of_Graphs_hyena
m=1;lap3=list()
for (i in lap_graph_hyena){
  for (j in lap_graph_hyena){
    lap1=eigs(i,k)$values
    lap2=eigs(j,k)$values
    lap3[[m]]=sqrt(sum(lap1-lap2)^2)
    m=m+1
  }
  lap_spec_dist=do.call(rbind,lap3)
  lap_spec_dist=unique(lap_spec_dist)
  lap_spec_dist=lap_spec_dist[-1]
}


########################----DOLPHIN AND ITS SYNTHETIC NETWORKS------------#################
Graphs_for_spectral_distance_dolphin=list(dolph11[[2]],dolph11[[51]],dolph11[[102]],dolph11[[185]],
                                          dolph11[[201]],dolph11[[202]])

Adj_matrix_of_Graphs_dolphin <- lapply(Graphs_for_spectral_distance_dolphin, as_adjacency_matrix)
Lap_matrix_of_Graphs_dolphin <- lapply(Graphs_for_spectral_distance_dolphin, laplacian_matrix)
Norm_Lap_matrix_of_Graphs_dolphin<-lapply(Graphs_for_spectral_distance_dolphin, normalized_laplacian)

adj_graph_dolphin=Adj_matrix_of_Graphs_dolphin
k=50
#####----Ajacency eigenvalues boxplot ----######
ER=eigs(adj_graph_dolphin[[1]],k)$values
SW=eigs(adj_graph_dolphin[[2]],k)$values
SF=eigs(adj_graph_dolphin[[3]],k)$values
SP=eigs(adj_graph_dolphin[[4]],k)$values
LAT=eigs(adj_graph_dolphin[[5]],k)$values
DOLPHIN=eigs(adj_graph_dolphin[[6]],k)$values
df_adj_dolphin=data.frame(ER,SW,SF,SP,LAT,DOLPHIN)

#melt(df, id.vars=c("ER","SW","SF","SP","LAT","dolphin"))
plot1_dolph=gather(df_adj_dolphin, GraphName, Values,ER,SW,SF,SP,LAT,DOLPHIN, factor_key=TRUE)%>%
  ggplot(aes(x = GraphName, y = Values, group_by=GraphName, color = GraphName)) +
  geom_boxplot(show.legend = T) +ylab("Values")+labs(title = "")+
  #scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.text.x=element_text(size = 15),
        legend.position = "none",plot.title.position = "plot",
        plot.title = element_text(hjust = 0.6,vjust = -0.6))

plot1_dolph

D1=plot1_dolph+annotate("text", x = "ER", y =7,
                  label = expression(
                    ~ lambda ["k,er"]^"A"))

D2=D1+annotate("text", x = "SW", y =6,
               label = expression(
                 ~ lambda ["k,sw"]^"A"))

D3=D2+annotate("text", x = "SF", y =11,
               label = expression(
                 ~ lambda ["k,sf"]^"A"))                  

D4=D3+annotate("text", x = "SP", y =8,
               label = expression(
                 ~ lambda ["k,sp"]^"A"))  

D5=D4+annotate("text", x = "LAT", y =5.2,
               label = expression(
                 ~ lambda ["k,lat"]^"A"))
D6=D5+annotate("text", x = "DOLPHIN", y =8,
               label = expression(
                 ~ lambda ["k,dolphin"]^"A"))
D6    

ggsave("Adjacency_eigenval_dolphin.png", width = 8, height = 8)




#####----laplacian eigenvalues boxplot ----######
lap_graph_dolphin=Lap_matrix_of_Graphs_dolphin
ER=eigs(lap_graph_dolphin[[1]],k)$values
SW=eigs(lap_graph_dolphin[[2]],k)$values
SF=eigs(lap_graph_dolphin[[3]],k)$values
SP=eigs(lap_graph_dolphin[[4]],k)$values
LAT=eigs(lap_graph_dolphin[[5]],k)$values
DOLPHIN=eigs(lap_graph_dolphin[[6]],k)$values
df_lap_dolphin=data.frame(ER,SW,SF,SP,LAT,DOLPHIN)


plot2_dolph=gather(df_lap_dolphin, GraphName, Values,ER,SW,SF,SP,LAT,DOLPHIN, factor_key=TRUE)%>%
  ggplot(aes(x = GraphName, y = Values, group_by=GraphName, color = GraphName)) +
  geom_boxplot(show.legend = T) +ylab("Values")+labs(title = "")+
  # scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.text.x=element_text(size = 15),
        legend.position = "none",plot.title.position = "plot",
        plot.title = element_text(hjust = 0.6,vjust = -0.6))

plot2_dolph

D2=plot2_dolph+annotate("text", x = "ER", y =15,
                        label = expression(
                          ~ lambda ["k,er"]^"L"))

D3=D2+annotate("text", x = "SW", y =13,
               label = expression(
                 ~ lambda ["k,sw"]^"L"))

D4=D3+annotate("text", x = "SF", y =54,
               label = expression(
                 ~ lambda ["k,sf"]^"L"))                  

D5=D4+annotate("text", x = "SP", y =14,
               label = expression(
                 ~ lambda ["k,sp"]^"L"))  

D6=D5+annotate("text", x = "LAT", y =10,
               label = expression(
                 ~ lambda ["k,lat"]^"L"))
D7=D6+annotate("text", x = "DOLPHIN", y =16.5,
               label = expression(
                 ~ lambda ["k,dolphin"]^"L"))
D7    

ggsave("laplacian_eigenval_dolphin.png", width = 8, height = 8)




######Normalized laplacian boxplot#########
k=50
norm_lap_graph_dolphin=Norm_Lap_matrix_of_Graphs_dolphin
ER=eigs(norm_lap_graph_dolphin[[1]],k)$values
SW=eigs(norm_lap_graph_dolphin[[2]],k)$values
SF=eigs(norm_lap_graph_dolphin[[3]],k)$values
SP=eigs(norm_lap_graph_dolphin[[4]],k)$values
LAT=eigs(norm_lap_graph_dolphin[[5]],k)$values
DOLPHIN=eigs(norm_lap_graph_dolphin[[6]],k)$values
df_lap_dolphin=data.frame(ER,SW,SF,SP,LAT,DOLPHIN)


plot3_dolph=gather(df_lap_dolphin, GraphName, Values,ER,SW,SF,SP,LAT,DOLPHIN, factor_key=TRUE)%>%
  ggplot(aes(x = GraphName, y = Values, group_by=GraphName, color = GraphName)) +
  geom_boxplot(show.legend = T) +ylab("Values")+labs(title = "")+
  # scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.text.x=element_text(size = 15),
        legend.position = "none",plot.title.position = "plot",
        plot.title = element_text(hjust = 0.6,vjust = -0.6))

plot3_dolph

D3=plot3_dolph+annotate("text", x = "ER", y =1.9,
                        label = expression(
                          ~ lambda ["k,er"]^"|L|"))

D4=D3+annotate("text", x = "SW", y =2,
               label = expression(
                 ~ lambda ["k,sw"]^"|L|"))

D5=D4+annotate("text", x = "SF", y =2,
               label = expression(
                 ~ lambda ["k,sf"]^"|L|"))                  

D6=D5+annotate("text", x = "SP", y =1.8,
               label = expression(
                 ~ lambda ["k,sp"]^"|L|"))  

D7=D6+annotate("text", x = "LAT", y =2.1,
               label = expression(
                 ~ lambda ["k,lat"]^"|L|"))
D8=D7+annotate("text", x = "DOLPHIN", y =1.9,
               label = expression(
                 ~ lambda ["k,dolphin"]^"|L|"))
D8    

ggsave("normalized_laplacian_eigenval_dolphin.png", width = 8, height = 8)








############## HYENA NETWORK ###########################
norm_lap_graph_hyena=Norm_Lap_matrix_of_Graphs_hyena
m=1;norm_lap3=list()
for (i in norm_lap_graph_hyena){
  for (j in norm_lap_graph_hyena){
    norm_lap1=eigs(i,k)$values
    norm_lap2=eigs(j,k)$values
    norm_lap3[[m]]=sqrt(sum(norm_lap1-norm_lap2)^2)
    m=m+1
  }
  norm_lap_spec_dist=do.call(rbind,norm_lap3)
  norm_lap_spec_dist=unique( norm_lap_spec_dist)
  norm_lap_spec_dist=norm_lap_spec_dist[-1]
  
}

all_spectral_distance=as.data.frame(cbind(adj_spec_dist,lap_spec_dist,norm_lap_spec_dist))
all_spectral_distance$Compared_networks=c("ER-SW","ER-SF","ER-SP","ER-LAT",
                                          paste("ER",sep = "-",name),"SW-SF","SW-SP","SW-LAT",
                                          paste("SW",sep = "-",name),"SF-SP","SF-LAT",paste("SF",sep = "-",name),
                                          "SP-LAT", paste("SP",sep = "-",name),paste("LAT",sep = "-",name))
all_spectral_distance_hyena <- all_spectral_distance[, c("Compared_networks","adj_spec_dist","lap_spec_dist",
                                                   "norm_lap_spec_dist")]

data_spect_hyena=write.csv(all_spectral_distance_hyena ,'C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi/HyenaSpectralFeat.csv')







########################----HYENA AND ITS SYNTHETIC NETWORKS------------#################
#####----Ajacency eigenvalues boxplot ----######
ER=eigs(adj_graph_hyena[[1]],k=30)$values
SW=eigs(adj_graph_hyena[[2]],k=30)$values
SF=eigs(adj_graph_hyena[[3]],k=30)$values
SP=eigs(adj_graph_hyena[[4]],k=30)$values
LAT=eigs(adj_graph_hyena[[5]],k=30)$values
HYENA=eigs(adj_graph_hyena[[6]],k=30)$values
df=data.frame(ER,SW,SF,SP,LAT,HYENA)

#melt(df, id.vars=c("ER","SW","SF","SP","LAT","HYENA"))
plot1=gather(df, GraphName, Values,ER,SW,SF,SP,LAT,HYENA, factor_key=TRUE)%>%
ggplot(aes(x = GraphName, y = Values, group_by=GraphName, color = GraphName)) +
  geom_boxplot(show.legend = T) +ylab("Values")+labs(title = "")+
  #scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.text.x=element_text(size = 15),
        legend.position = "none",plot.title.position = "plot",
        plot.title = element_text(hjust = 0.6,vjust = -0.6))
plot1

p1=plot1+annotate("text", x = "ER", y =31.5,
                  label = expression(
                    ~ lambda ["k,er"]^"A"))

p2=p1+annotate("text", x = "SW", y =32,
               label = expression(
                 ~ lambda ["k,sw"]^"A"))

p3=p2+annotate("text", x = "SF", y =32,
               label = expression(
                 ~ lambda ["k,sf"]^"A"))                  

p4=p3+annotate("text", x = "SP", y =32,
               label = expression(
                 ~ lambda ["k,sp"]^"A"))  

p5=p4+annotate("text", x = "LAT", y =34,
               label = expression(
                 ~ lambda ["k,lat"]^"A"))
p6=p5+annotate("text", x = "HYENA", y =32,
               label = expression(
                 ~ lambda ["k,hyena"]^"A"))
p6    

ggsave("Adjacency_eigenval.png", width = 8, height = 8)



#####----laplacian eigenvalues boxplot ----######
ER=eigs(lap_graph_hyena[[1]],k=30)$values
SW=eigs(lap_graph_hyena[[2]],k=30)$values
SF=eigs(lap_graph_hyena[[3]],k=30)$values
SP=eigs(lap_graph_hyena[[4]],k=30)$values
LAT=eigs(lap_graph_hyena[[5]],k=30)$values
HYENA=eigs(lap_graph_hyena[[6]],k=30)$values
df_lap_hyena=data.frame(ER,SW,SF,SP,LAT,HYENA)


plot2=gather(df_lap_hyena, GraphName, Values,ER,SW,SF,SP,LAT,HYENA, factor_key=TRUE)%>%
  ggplot(aes(x = GraphName, y = Values, group_by=GraphName, color = GraphName)) +
  geom_boxplot(show.legend = T) +ylab("Values")+labs(title = "")+
  # scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.text.x=element_text(size = 15),
        legend.position = "none",plot.title.position = "plot",
        plot.title = element_text(hjust = 0.6,vjust = -0.6))
plot2


p2=plot2+annotate("text", x = "ER", y =34.3,
                  label = expression(
                    ~ lambda ["k,er"]^"L"))

p3=p2+annotate("text", x = "SW", y =35.3,
               label = expression(
                 ~ lambda ["k,sw"]^"L"))

p4=p3+annotate("text", x = "SF", y =35.3,
               label = expression(
                 ~ lambda ["k,sf"]^"L"))                  

p5=p4+annotate("text", x = "SP", y =35.3,
               label = expression(
                 ~ lambda ["k,sp"]^"L"))  

p6=p5+annotate("text", x = "LAT", y =36.5,
               label = expression(
                 ~ lambda ["k,lat"]^"L"))
p7=p6+annotate("text", x = "HYENA", y =35.5,
               label = expression(
                 ~ lambda ["k,hyena"]^"L"))
p7    

ggsave("laplacian_eigenval.png", width = 8, height = 8)


######Normalized normalizedlaplacian boxplot#########
ER=eigs(norm_lap_graph_hyena[[1]],k=30)$values
SW=eigs(norm_lap_graph_hyena[[2]],k=30)$values
SF=eigs(norm_lap_graph_hyena[[3]],k=30)$values
SP=eigs(norm_lap_graph_hyena[[4]],k=30)$values
LAT=eigs(norm_lap_graph_hyena[[5]],k=30)$values
HYENA=eigs(norm_lap_graph_hyena[[6]],k=30)$values
df_lap_hyena=data.frame(ER,SW,SF,SP,LAT,HYENA)


plot3=gather(df_lap_hyena, GraphName, Values,ER,SW,SF,SP,LAT,HYENA, factor_key=TRUE)%>%
  ggplot(aes(x = GraphName, y = Values, group_by=GraphName, color = GraphName)) +
  geom_boxplot(show.legend = T) +ylab("Values")+labs(title = "")+
  # scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.text.x=element_text(size = 15),
        legend.position = "none",plot.title.position = "plot",
        plot.title = element_text(hjust = 0.6,vjust = -0.6))
plot3


p4=plot3+annotate("text", x = "ER", y =1.18,
               label = expression(
                 ~ lambda ["k,er"]^"|L|"))

p5=p4+annotate("text", x = "SW", y =1.185,
                 label = expression(
                   ~ lambda ["k,sw"]^"|L|"))

p6=p5+annotate("text", x = "SF", y =1.19,
               label = expression(
                 ~ lambda ["k,sf"]^"|L|"))                  

p7=p6+annotate("text", x = "SP", y =1.18,
               label = expression(
                 ~ lambda ["k,sp"]^"|L|"))  

p8=p7+annotate("text", x = "LAT", y =1.21,
               label = expression(
                 ~ lambda ["k,lat"]^"|L|"))
p9=p8+annotate("text", x = "HYENA", y =1.18,
               label = expression(
                 ~ lambda ["k,hyena"]^"|L|"))
p9    

ggsave("normalized_laplacian_eigenval.png", width = 8, height = 8)



######################------------GRAPH FEATURES OF DOLPHIN AND HYENA-----------#######################
#################----Feature based distances----------############################
set.seed(120)
dolph1=makeSimGraphs_dolph(nSamples = 50,order = vcount(ggk1))
dolph1$"real"=ggk1

h=makeSimGraphs_hyena(nSamples = 50,order = vcount(gg1))
h$"real"=gg1


nsim=10 
nreps=10
nticks=100 

nsim=100

# G=makeSimGraphs_dolph(nSamples = 1,order=5)
# betaVals=c(0.01) 
# gammaVals=c(0.033)

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=0.143
# gammaVals=c(0.033,0.143,0.25,0.5)

Hyena_vs_synthetic_networks=RunEpicSimOnGraphs(h, betaVals, gammaVals, nreps=nreps,output_file="hyena_vs_synthetic_network_new4.csv")

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=0.143
#gammaVals=c(0.033,0.143,0.25,0.5)

Dolphin_vs_synthetic_networks=RunEpicSimOnGraphs(dolph1, betaVals, gammaVals, nreps=nreps,output_file="Dolphine_vs_synthetic_network_new4.csv")


##########-------Quantiles of features-----------#################
Graphfeatures=function(Name="ER"){
  df_feat=data.frame(read.csv(data,header = T, sep = ","))
  df_feat=df_feat%>% filter(GraphName==Name)
  return(df_feat)  
}


df_dolph=read.csv("Dolphine_vs_synthetic_network_new4.csv")

####---Beta=0.01, Gamma=0.143----#########
asd=df_dolph%>%filter(beta==0.01, gamma==0.143)
View(asd)

QuantilesAnalysis=function(data="Dolphine_vs_synthetic_network_new4.csv",betaval=0.01,gammaval=0.143,realdata="Dolphin"){
  df_dolph=data.frame(read.csv(data,header = T, sep = ","))
  df_dolph%>%filter(beta==betaval, gamma==gammaval)
  
  df_ER=df_dolph%>%filter(GraphName=="ER")
  ERDATA=data.frame(min(df_ER$minDegree),max(df_ER$maxDegree),
               quantile(df_ER$FiedlerValue,probs = c(0.25, 0.75)),
               quantile(df_ER$Normalized_FiedlerValue,probs = c(0.25, 0.75)),
               quantile(df_ER$spectral_radius,probs = c(0.25, 0.75)),
               quantile(df_ER$modularity,probs = c(0.25, 0.75)),
               quantile(df_ER$betweenness,probs = c(0.25, 0.75)),
               quantile(df_ER$transitivity,probs = c(0.25, 0.75)),
               quantile(df_ER$diameter,probs = c(0.25, 0.75)))    
  
  colnames(ERDATA)=c("minDegree","maxDegree","FiedlerQuantiles","NormFiedlerQuantiles","SpectralRadiusQuantiles",
                     "modularityQuantiles","BetweenessQuantiles","TransitivityQuantiles","Diameter")
  ERDATA$GraphName="ER"                   
  
  
  df_SW=df_dolph%>%filter(GraphName=="SW")
  SWDATA=data.frame(min(df_SW$minDegree),max(df_SW$maxDegree),
                          quantile(df_SW$FiedlerValue,probs = c(0.25, 0.75)),
                    quantile(df_SW$Normalized_FiedlerValue,probs = c(0.25, 0.75)),
                          quantile(df_SW$spectral_radius,probs = c(0.25, 0.75)),
                          quantile(df_SW$modularity,probs = c(0.25, 0.75)),
                          quantile(df_SW$betweenness,probs = c(0.25, 0.75)),
                          quantile(df_ER$transitivity,probs = c(0.25, 0.75)),
                          quantile(df_ER$diameter,probs = c(0.25, 0.75)))
  
  colnames(SWDATA)=c("minDegree","maxDegree","FiedlerQuantiles","NormFiedlerQuantiles","SpectralRadiusQuantiles",
                     "modularityQuantiles","BetweenessQuantiles","TransitivityQuantiles","Diameter")
  SWDATA$GraphName="SW" 
  
  df_SF=df_dolph%>%filter(GraphName=="SF")
  SFDATA=data.frame(min(df_SF$minDegree),max(df_SF$maxDegree),
                    quantile(df_SF$FiedlerValue,probs = c(0.25, 0.75)),
                    quantile(df_SF$Normalized_FiedlerValue,probs = c(0.25, 0.75)),
                    quantile(df_SF$spectral_radius,probs = c(0.25, 0.75)),
                    quantile(df_SF$modularity,probs = c(0.25, 0.75)),
                    quantile(df_SF$betweenness,probs = c(0.25, 0.75)),
                    quantile(df_ER$transitivity,probs = c(0.25, 0.75)),
                    quantile(df_ER$diameter,probs = c(0.25, 0.75)))
  
  colnames(SFDATA)=c("minDegree","maxDegree","FiedlerQuantiles","NormFiedlerQuantiles","SpectralRadiusQuantiles",
                     "modularityQuantiles","BetweenessQuantiles","TransitivityQuantiles","Diameter")
  SFDATA$GraphName="SF" 
  
  
  df_SP=df_dolph%>%filter(GraphName=="SP")
  SPDATA=data.frame(min(df_SP$minDegree),max(df_SP$maxDegree),
                    quantile(df_SP$FiedlerValue,probs = c(0.25, 0.75)),
                    quantile(df_SP$Normalized_FiedlerValue,probs = c(0.25, 0.75)),
                    quantile(df_SP$spectral_radius,probs = c(0.25, 0.75)),
                    quantile(df_SP$modularity,probs = c(0.25, 0.75)),
                    quantile(df_SP$betweenness,probs = c(0.25, 0.75)),
                    quantile(df_ER$transitivity,probs = c(0.25, 0.75)),
                    quantile(df_ER$diameter,probs = c(0.25, 0.75)))
  
  colnames(SPDATA)=c("minDegree","maxDegree","FiedlerQuantiles","NormFiedlerQuantiles","SpectralRadiusQuantiles",
                     "modularityQuantiles","BetweenessQuantiles","TransitivityQuantiles","Diameter")
  SPDATA$GraphName="SP" 
  
  df_Lat=df_dolph%>%filter(GraphName=="Lat")
  LatDATA=data.frame(min(df_Lat$minDegree),max(df_Lat$maxDegree),
                     quantile(df_Lat$FiedlerValue,probs = c(0.25, 0.75)),
                     quantile(df_Lat$Normalized_FiedlerValue,probs = c(0.25, 0.75)),
                     quantile(df_Lat$spectral_radius,probs = c(0.25, 0.75)),
                     quantile(df_Lat$modularity,probs = c(0.25, 0.75)),
                     quantile(df_Lat$betweenness,probs = c(0.25, 0.75)),
                     quantile(df_ER$transitivity,probs = c(0.25, 0.75)),
                     quantile(df_ER$diameter,probs = c(0.25, 0.75)))
  
  colnames(LatDATA)=c("minDegree","maxDegree","FiedlerQuantiles","NormFiedlerQuantiles","SpectralRadiusQuantiles",
                      "modularityQuantiles","BetweenessQuantiles","TransitivityQuantiles","Diameter")
  LatDATA$GraphName="Lat" 
  
  
  df_realdata=df_dolph%>%filter(GraphName==realdata)
  DATA=data.frame(min(df_realdata$minDegree),max(df_realdata$maxDegree),
                         quantile(df_realdata$FiedlerValue,probs = c(0.25, 0.75)),
                  quantile(df_realdata$Normalized_FiedlerValue,probs = c(0.25, 0.75)),
                         quantile(df_realdata$spectral_radius,probs = c(0.25, 0.75)),
                         quantile(df_realdata$modularity,probs = c(0.25, 0.75)),
                         quantile(df_realdata$betweenness,probs = c(0.25, 0.75)),
                         quantile(df_ER$transitivity,probs = c(0.25, 0.75)),
                         quantile(df_ER$diameter,probs = c(0.25, 0.75)))
  
  colnames(DATA)=c("minDegree","maxDegree","FiedlerQuantiles","NormFiedlerQuantiles","SpectralRadiusQuantiles",
                          "modularityQuantiles","BetweenessQuantiles","TransitivityQuantiles","Diameter")
  DATA$GraphName=realdata 
  
  finaldata=rbind(ERDATA,SWDATA,SFDATA,SPDATA,LatDATA,DATA)
  return(finaldata)
}

###################--------DOLPHIN---------##############################
#### beta=0.01, gamma=0.143
z1=QuantilesAnalysis(data="Dolphine_vs_synthetic_network_new4.csv",betaval=0.01,gammaval=0.143,realdata="Dolphin")
z1
#### beta=0.1, gamma=0.143
z2=QuantilesAnalysis(data="Dolphine_vs_synthetic_network_new4.csv",betaval=0.1,gammaval=0.143,realdata="Dolphin")
z2
#### beta=0.33, gamma=0.143
z3=QuantilesAnalysis(data="Dolphine_vs_synthetic_network_new4.csv",betaval=0.33,gammaval=0.143,realdata="Dolphin")
z3
#### beta=0.5, gamma=0.143
z4=QuantilesAnalysis(data="Dolphine_vs_synthetic_network_new4.csv",betaval=0.5,gammaval=0.143,realdata="Dolphin")
z4



###################--------HYENA---------##############################
#### beta=0.01, gamma=0.143
y1=QuantilesAnalysis(data="hyena_vs_synthetic_network_new4.csv",betaval=0.01,gammaval=0.143,realdata="Hyena")
y1
#### beta=0.1, gamma=0.143
y2=QuantilesAnalysis(data="hyena_vs_synthetic_network_new4.csv",betaval=0.1,gammaval=0.143)
y2
#### beta=0.33, gamma=0.143
y3=QuantilesAnalysis(data="hyena_vs_synthetic_network_new4.csv",betaval=0.33,gammaval=0.143)
y3
#### beta=0.5, gamma=0.143
y4=QuantilesAnalysis(data="hyena_vs_synthetic_network_new4.csv",betaval=0.5,gammaval=0.143)
y4




# x=Graphfeatures(Name="ER",data)
# min(x$minDegree)
# max(x$maxDegree)
# quantile(x$FiedlerValue,probs = c(0.25, 0.75))
# quantile(x$modularity,probs = c(0.25, 0.75))
# quantile(x$betweenness,probs = c(0.25, 0.75))
# quantile(x$transitivity,probs = c(0.25, 0.75))
# 
# ##SW
# x=Graphfeatures(Name="SW",data)
# min(x$minDegree)
# max(x$maxDegree)
# quantile(x$FiedlerValue,probs = c(0.25, 0.75))
# quantile(x$modularity,probs = c(0.25, 0.75))
# quantile(x$betweenness,probs = c(0.25, 0.75))
# quantile(x$transitivity,probs = c(0.25, 0.75))
# 
# ##SF
# x=Graphfeatures(Name="SF",data)
# min(x$minDegree)
# max(x$maxDegree)
# quantile(x$FiedlerValue,probs = c(0.25, 0.75))
# quantile(x$modularity,probs = c(0.25, 0.75))
# quantile(x$betweenness,probs = c(0.25, 0.75))
# quantile(x$transitivity,probs = c(0.25, 0.75))
# 
# 
# ##SP
# x=Graphfeatures(Name="SP",data)
# min(x$minDegree)
# max(x$maxDegree)
# quantile(x$FiedlerValue,probs = c(0.25, 0.75))
# quantile(x$modularity,probs = c(0.25, 0.75))
# quantile(x$betweenness,probs = c(0.25, 0.75))
# quantile(x$transitivity,probs = c(0.25, 0.75))
# 
# 
# ##Lat
# x=Graphfeatures(Name="Lat",data)
# min(x$minDegree)
# max(x$maxDegree)
# quantile(x$FiedlerValue,probs = c(0.25, 0.75))
# quantile(x$modularity,probs = c(0.25, 0.75))
# quantile(x$betweenness,probs = c(0.25, 0.75))
# quantile(x$transitivity,probs = c(0.25, 0.75))
# 
# ##Hyena
# x=Graphfeatures(Name="Hyena",data)
# min(x$minDegree)
# max(x$maxDegree)
# quantile(x$FiedlerValue,probs = c(0.25, 0.75))
# quantile(x$modularity,probs = c(0.25, 0.75))
# quantile(x$betweenness,probs = c(0.25, 0.75))
# quantile(x$transitivity,probs = c(0.25, 0.75))


# simPathE <- function(G, nTicks=100, beta=0.5, gamma=0.1, propInfected=0.1, initialState=NULL, nInfected=1, useProportion=F) {
#   nVertices = gorder(G)
#   infectionState = as.data.frame(matrix(0, ncol = nVertices, nrow = nTicks+1))
#   
#   # Set up the initial state: all vertices are either exposed (state = 0) or infected (state = 1); none can be recovered yet (2).
#   if (is.null(initialState)) {
#     if (useProportion==T) {
#       infectionState[1,] <- rbinom(nVertices, 1, propInfected) # set initial state of each nodes if known (note, no recovered node at initial state)
#     } else {
#       infected <- rep(1, nInfected) # just create a vector of the right number of 1s
#       exposed <- rep(0, (nVertices - nInfected))
#       infectionState[1,] <- sample(c(infected, exposed), nVertices, replace=FALSE)
#     }
#   } else {
#     if (length(initialState) != nVertices) {
#       return ("Initial state and order of Graph (number of vertices) are incompatible.  Check the sizes of your input.")
#     }
#     infectionState <- initialState # initial existing state.
#   }
#   
#   adjacencyList <- as_adj_list(G) # this is a list of which vertices are adjacent to each vertex i
#   
#   # Now do the simulation through time:
#   for (t in 1:nTicks) {
#     # FIRST phase: transmission: S -> I
#     for (i in which(infectionState[t,] == 0)) { # for all susceptible nodes (denoted 0) in previous time step
#       infectionState[t+1,i] <- 0 # node remains as susceptible, since not all contact leads to an infection.
#       for (j in adjacencyList[[i]]) { # for all neighbours of i
#         if (infectionState[t,j][1] == 1) { # vertex j is infectious
#           if ((runif(1)*G[i,j]) <= beta) { # ... and passes it on!
#             infectionState[t+1,i] <- 1;
#             break # assign node as infected if above condition is met, and break out of loop: we don't need
#             # to check any more adjacent vertices.
#           }
#         }
#       }
#     }
#     
#     # SECOND phase: recovery: I -> R
#     for (i in which(infectionState[t,] == 1)) { # for all infected nodes (denoted 1) in previous time step
#       if (runif(1) <= gamma) { # compares a randomly generated uniform number to recovery rate
#         infectionState[t+1,i] <- 2 # node is recovered
#       } else {
#         infectionState[t+1,i] <- 1 # node remains infected 
#       }
#     }
#     
#     # THIRD phase: recovered stays recovered:
#     for (i in which(infectionState[t,] == 2)) { # for all recovered nodes (denoted 2) in previous time step
#       infectionState[t+1,i] <- 2 # node stays recovered
#     }
#   }
#   rownames(infectionState) = 0:nTicks 
#   return(infectionState)
# }

# 
# ############----------AVERAGE PROPORTION OF INFECTED FUNCTION------###############
# AvrgPropOf<-function(Name="ER",ID=1,betaVal=.01,gammaVal=.2,nreps=nsim,Data="observed_networks_50_avrgdeg_4.csv",plot_title="",nticks=100,net_size=10){
#   df_obs=data.frame(read.csv(Data,header = T, sep = ","))
#   df_obs=df_obs%>% filter(GraphName==Name,GraphID==ID, beta==betaVal,gamma==gammaVal)
#   df_obs=df_obs%>%select(c(GraphName,t0:paste("t",nticks,sep = "")))
#   #y=f%>%group_by(GraphName)%>%summarize(colMeans(f[sapply(f, is.numeric)]))
#   df=data.frame(colMeans(df_obs[sapply(df_obs, is.numeric)]))
#   colnames(df)=c("value")
#   df$Timesteps=1:(nticks+1)
#   df$GraphName=Name
#   df=df%>%select(GraphName, Timesteps,value)
#   
#   return(df)   
# }
# 
# plotfunc_obs<-function(Name="ER",ID=1,betaVal=.01,gammaVal=.2,Xlabel="Timesteps (days)",Ylabel="Prop-Infected",nreps=nsim,Data="observed_networks_50_avrgdeg_4.csv",plot_title="",nticks=100,net_size=10){
#   data_obs=read.csv(Data,header = T, sep = ",")
#   df_obs=data_obs%>% filter(GraphName==Name,GraphID==ID, beta==betaVal,gamma==gammaVal)
#   df_obs=df_obs%>% select(c(beta,gamma,t0:paste("t",nticks,sep = "")))
#   
#   p_obs=data.frame(0:(length(df_obs[grep("t", names(df_obs))])-1),t(df_obs[grep("t", names(df_obs))]))
#   colnames(p_obs)= c("Timesteps",paste("Infecteds",1:nreps,sep = "_"))
#   df_obs_long_format <- melt(p_obs, id="Timesteps")  # convert to long format
#   
#   
#   plot_obs<- ggplot(data= df_obs_long_format,
#                     aes(x=Timesteps, y=value/net_size, colour=variable)) +
#     geom_line(show.legend = FALSE)+theme_classic()+labs(title="",tag=Name )+
#     xlab(Xlabel)+ylab(Ylabel)+
#     scale_y_continuous(limits = c(0,1))+theme(text = element_text(size = 34,face = "italic"),
#                                               axis.title = element_text(size = 34),
#                                               axis.text.y = element_text(size = 34),
#                                               axis.text.x = element_text(size = 34),
#                                               legend.position = "none",
#                                               plot.title.position = "plot",
#                                               plot.tag.position = c(0.1, 0.98))
#   #theme(text = element_text(size = 22))    
#   #scale_y_continuous(limits = c(0,net_size))
#   return(plot_obs)   
# }
# 
# 
# 
# Complete_Graph<-function(beta=0.1,gamma=0.4,avrdg=4,sim=10,ntime=100,net_size=10,Xlabel="Timesteps (days)",Ylabel="Prop-Infected"){
#   param <- param.icm(inf.prob = beta, act.rate = avrdg, rec.rate = gamma,
#                      a.rate = 0, ds.rate = 0, di.rate = 0,
#                      dr.rate = 0)
#   init <- init.icm(s.num = net_size, i.num = 1, r.num = 0)
#   control <- control.icm(type = "SIR", nsteps = ntime, nsims = sim)
#   sim <- icm(param, init, control)
#   #plot(sim, legend=FALSE)
#   #plot(sim, y = "i.num", sim.lines = TRUE, mean.smooth = FALSE, qnts.smooth = FALSE,main = "CG-Infected-dynamics-beta:0.5,gamma=0.033")
#   f1=as.data.frame(sim, out = "vals")
#   f1$GraphName="CG"
#   
#   CGPlot=ggplot(data= f1,
#                 aes(x=time, y=i.num/num, colour=as.factor(sim))) +
#     geom_line(show.legend = FALSE)+theme_classic()+labs(title="",tag = "CG")+
#     xlab(Xlabel)+ylab(Ylabel)+
#     scale_y_continuous(limits = c(0,1))+
#     theme(text = element_text(size = 34,face = "italic"),
#           axis.title = element_text(size = 34),
#           axis.text.y = element_text(size = 34),
#           axis.text.x = element_text(size = 34),
#           legend.position = "none",
#           plot.title.position = "plot",
#           plot.tag.position = c(0.1, 0.98))
#   
#   
#   return(CGPlot)
#   
# }
