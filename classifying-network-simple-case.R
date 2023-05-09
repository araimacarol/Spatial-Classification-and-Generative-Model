setwd("C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi-New")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Building the Classification model
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#This project aim to build a simple classifcation model to classify empirical networks
#by the network models and the spatial network and investigate which categorical class
# the empirical networks belongss to and with which probability from the feature vectors


#++++++++++++++++++++++++++++++++
#+ Loading packages
#+++++++++++++++++++++++++++++++

suppressPackageStartupMessages({
  library(Matrix)
  library(stats)
  library(igraph)
  library(tidymodels)
  library(tidyverse)
  library(ggplot2)
  library(janitor)
  library(vip)
  ##--Libraries for EDA and Data Wrangling
  library(visdat)
  library(gt)
  library(plotly)
  library(skimr)
  library(GGally)
  library("corrplot")
  library(psych)
  library(tidyverse)
  library(lubridate)
  library(recipeselectors)  
  
  ##--Libraries for Feature Engineering
  library(Boruta)
  
  
  
  ##--Libraries for Models
  #library(mrIML)
  library(tidymodels)
  tidymodels_prefer()
  library(embed)
  library(textrecipes)
  library(vip)
  library("randomForest")
  library("gbm")
  library("tidyverse")
  library("MLmetrics")
  library("readxl")
  library("varImp")
  library(ranger)
  
  ##--Libraries (Others)
  library(igraph)
  library(purrr)
  library(beans)
  library(ggforce)
  library(patchwork)
  library(bestNormalize)
  library(Metrics)
  library(future)
  library("future.apply")
  library(usemodels)
  library(here)
  library(parsnip)
})



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Data generation for the machine learning model
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++
#+ Helper functions
#+++++++++++++++++++++++++++++++
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

###--Normalized Laplacina function--##
normalized_laplacian=function(Graphs){
  laplacian_matrix(Graphs,normalized = T)
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Calculate graph features
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
    "closeness_centr",                # average inverse of distance between any pair of vertices
    "modularity",               # DEFINITION REQUIRED
    "diameter",                 # maximum distance between any two vertices (NAN if not connected)
    "betw_centr",              # max_{v} proportion of shortest paths going through vertex v
    "transitivity",             # aka Clustering Coefficient, is proportion of connected triples that form triangles: e.g., (a--b--c--a) when (a--b--c) is present.
    "threshold",                 # 1/max(eigen value of A)
    "spectral_radius"         # max (eigen value of A)
    
  )
  
  df <- as.data.frame(matrix(ncol=length(features),nrow=length(Graphs)))
  colnames(df)=features
  
  # Stuff that is simple to apply and needs no interim components:
  
  df$order = base::as.numeric(lapply(Graphs, gorder))
  df$edges = base::as.numeric(lapply(Graphs, gsize))
  df$connected = base::as.numeric(lapply(Graphs, is_connected))
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
    
    df$eigen_centr[i] <- centr_eigen(Graphs[[i]])$centralization
    df$deg_centr[i] <- centr_degree(Graphs[[i]])$centralization
    
    df$betw_centr[i] <- centr_betw(Graphs[[i]])$centralization
    
    df$max_component[i] <- max(components(Graphs[[i]])$csize)
    df$mean_eccentr[i]<-mean(eccentricity(Graphs[[i]]))
    df$radius[i]<-radius(Graphs[[i]])
    df$mean_path_length[i]<-average.path.length(Graphs[[i]])
    #df$trace[i]<-sum(diag(Adj[[i]]))
    df$graph_energy[i]<-sum(abs(eigen(Adj[[i]], symmetric=TRUE, only.values=TRUE)$values))
    df$min_triangle[i]= min(count_triangles(Graphs[[i]]))
    df$mean_triangle[i]= mean(count_triangles(Graphs[[i]]))
    df$sd_triangle[i]= sd(count_triangles(Graphs[[i]]))
    df$max_triangle[i]= max(count_triangles(Graphs[[i]]))
    df$num_triangle[i]= sum(count_triangles(Graphs[[i]]))
    df$deg_assort_coef[i]=assortativity_degree(Graphs[[i]])
    
    df$threshold[i] <- 1/(df$spectral_radius[i])
    
    if (df$connected[i]==TRUE) {
      df$closeness_centr[i] = mean(closeness(Graphs[[i]]))
    } else { # handle the case where G isn't connected
      df$closeness_centr[i] = -1
    }
  }
  return (df)
}

#calcGraphFeatures(x)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Multiple Calculation of graph features
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Used to perform runs on multiple simulated graphs on any given network
RunSimOnGraphFeatures<-function(Graphs, nreps=nreps,output_file=NULL, seed=-1) {
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
  # write.csv(All_results, file=output_file)
  return( All_results)
}




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Generative Graph models
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++
#+ #----MAKING OF SPATIAL GRAPH---
#+++++++++++++++++++++++++++++++++++++++++++++++

# Helper functions for spatial graphs to convert between (row, column) pairs and the index in a list of Cells.
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


####---Induce graph function---######
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


#+++++++++++++++++++++++++++++++++++++++++++
# Function to make spatial network
#+++++++++++++++++++++++++++++++++++++++++++++
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

#induced.graph(x)
#induced.graph(x[1])

###########################################################
# Create spatial graphs for simulation experiment
###########################################################
makeSpatialGraphs<- function(n=25,r=0.8) {
  Graphs = list()#;G=list() # set up the list of Graphs first
  initgraph= list()
  i= 1
  print("Creating fastspatial networks")
  Graphs[[i]] = fastSpatialNetwork(n = n, r = r, makeConnected=T,keepCellsSeparate=FALSE)
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

sim.sp.net=function(r=c(0.3,0.4,0.6,0.05),n=100){
  net=NULL;data=NULL
  net=as.list(mapply(FUN =makeSpatialGraphs,n,r))
  data= future.apply::future_lapply(list(net),RunSimOnGraphFeatures,nreps=1,future.seed = 0xBEEF)
  
  df=cbind(data.frame(r),data)  
  return(df)
}

#sim.sp.net(r=rep(0.1,2),n=20)

plan(multisession, workers = 32)
set.seed(0xBEEF)

# sim.sp.net=function(r=c(0.3,0.4,0.6,0.05),n=100){
#   net=NULL;data=NULL
#   k=1
#   for (i in 1:length(r)){
#     net[i]= makeSpatialGraphs(n,r[i])
#     data=RunSimOnGraphFeatures(net,nreps = 1)
#     k=k+1
#   }
#   df=cbind(r,data)
#   return(df)
# }
#sim.sp.net(r=rep(0.1,2),n=50)
### 50 nodes 
SP1.50=sim.sp.net(r=rep(0.1,250),n=50)
write.csv(SP1.50,"SP1.50.csv")
SP2.50=sim.sp.net(r=rep(0.2,250),n=50)
write.csv(SP2.50,"SP2.50.csv")
SP3.50=sim.sp.net(r=rep(0.3,250),n=50)
write.csv(SP3.50,"SP3.50.csv")
SP4.50=sim.sp.net(r=rep(0.4,250),n=50)
write.csv(SP4.50,"SP4.50.csv")
SP5.50=sim.sp.net(r=rep(0.5,250),n=50)
write.csv(SP5.50,"SP5.50.csv")
SP6.50=sim.sp.net(r=rep(0.6,250),n=50)
write.csv(SP6.50,"SP6.50.csv")
SP7.50=sim.sp.net(r=rep(0.7,250),n=50)
write.csv(SP7.50,"SP7.50.csv")
SP8.50=sim.sp.net(r=rep(0.8,250),n=50)
write.csv(SP8.50,"SP8.50.csv")
SP9.50=sim.sp.net(r=rep(0.9,250),n=50)
write.csv(SP9.50,"SP9.50.csv")
######## Saving all SP graphs 50 nodes data to csv
SP1=read.csv("SP1.50.csv");SP2=read.csv("SP2.50.csv");SP3=read.csv("SP3.50.csv")
SP4=read.csv("SP4.50.csv");SP5=read.csv("SP5.50.csv");SP6=read.csv("SP6.50.csv")
SP7=read.csv("SP7.50.csv");SP8=read.csv("SP8.50.csv");SP9=read.csv("SP9.50.csv")

df.SP.50=rbind(SP1,SP2,SP3,SP4,SP5,SP6,SP7,SP8,SP9)
write.csv(df.SP.50,"df.SP.50.csv")

### 100 nodes 
SP1.100=sim.sp.net(r=rep(0.1,250),n=100)
write.csv(SP1.100,"SP1.100.csv")
SP2.100=sim.sp.net(r=rep(0.2,250),n=100)
write.csv(SP2.100,"SP2.100.csv")
SP3.100=sim.sp.net(r=rep(0.3,250),n=100)
write.csv(SP3.100,"SP3.100.csv")
SP4.100=sim.sp.net(r=rep(0.4,250),n=100)
write.csv(SP4.100,"SP4.100.csv")
SP5.100=sim.sp.net(r=rep(0.5,250),n=100)
write.csv(SP5.100,"SP5.100.csv")
SP6.100=sim.sp.net(r=rep(0.6,250),n=100)
write.csv(SP6.100,"SP6.100.csv")
SP7.100=sim.sp.net(r=rep(0.7,250),n=100)
write.csv(SP7.100,"SP7.100.csv")
SP8.100=sim.sp.net(r=rep(0.8,250),n=100)
write.csv(SP8.100,"SP8.100.csv")
SP9.100=sim.sp.net(r=rep(0.9,250),n=100)
write.csv(SP9.100,"SP9.100.csv")
######## Saving all SP graphs 100 nodes data to csv
SP1=read.csv("SP1.100.csv");SP2=read.csv("SP2.100.csv");SP3=read.csv("SP3.100.csv")
SP4=read.csv("SP4.100.csv");SP5=read.csv("SP5.100.csv");SP6=read.csv("SP6.100.csv")
SP7=read.csv("SP7.100.csv");SP8=read.csv("SP8.100.csv");SP9=read.csv("SP9.100.csv")

df.SP.100=rbind(SP1,SP2,SP3,SP4,SP5,SP6,SP7,SP8,SP9)
write.csv(df.SP.100,"df.SP.100.csv")

### 150 nodes 
SP1.150=sim.sp.net(r=rep(0.1,250),n=150)
write.csv(SP1.150,"SP1.150.csv")
SP2.150=sim.sp.net(r=rep(0.2,250),n=150)
write.csv(SP2.150,"SP2.150.csv")
SP3.150=sim.sp.net(r=rep(0.3,250),n=150)
write.csv(SP3.150,"SP3.150.csv")
SP4.150=sim.sp.net(r=rep(0.4,250),n=150)
write.csv(SP4.150,"SP4.150.csv")
SP5.150=sim.sp.net(r=rep(0.5,250),n=150)
write.csv(SP5.150,"SP5.150.csv")
SP6.150=sim.sp.net(r=rep(0.6,250),n=150)
write.csv(SP6.150,"SP6.150.csv")
SP7.150=sim.sp.net(r=rep(0.7,250),n=150)
write.csv(SP7.150,"SP7.150.csv")
SP8.150=sim.sp.net(r=rep(0.8,250),n=150)
write.csv(SP8.150,"SP8.150.csv")
SP9.150=sim.sp.net(r=rep(0.9,250),n=150)
write.csv(SP9.150,"SP9.150.csv")
######## Saving all SP graphs 150 nodes data to csv
SP1=read.csv("SP1.150.csv");SP2=read.csv("SP2.150.csv");SP3=read.csv("SP3.150.csv")
SP4=read.csv("SP4.150.csv");SP5=read.csv("SP5.150.csv");SP6=read.csv("SP6.150.csv")
SP7=read.csv("SP7.150.csv");SP8=read.csv("SP8.150.csv");SP9=read.csv("SP9.150.csv")

df.SP.150=rbind(SP1,SP2,SP3,SP4,SP5,SP6,SP7,SP8,SP9)
write.csv(df.SP.150,"df.SP.150.csv")

### 200 nodes 
SP1.200=sim.sp.net(r=rep(0.1,250),n=200)
write.csv(SP1.200,"SP1.200.csv")
SP2.200=sim.sp.net(r=rep(0.2,250),n=200)
write.csv(SP2.200,"SP2.200.csv")
SP3.200=sim.sp.net(r=rep(0.3,250),n=200)
write.csv(SP3.200,"SP3.200.csv")
SP4.200=sim.sp.net(r=rep(0.4,250),n=200)
write.csv(SP4.200,"SP4.200.csv")
SP5.200=sim.sp.net(r=rep(0.5,250),n=200)
write.csv(SP5.200,"SP5.200.csv")
SP6.200=sim.sp.net(r=rep(0.6,250),n=200)
write.csv(SP6.200,"SP6.200.csv")
SP7.200=sim.sp.net(r=rep(0.7,250),n=200)
write.csv(SP7.200,"SP7.200.csv")
SP8.200=sim.sp.net(r=rep(0.8,250),n=200)
write.csv(SP8.200,"SP8.200.csv")
SP9.200=sim.sp.net(r=rep(0.9,250),n=200)
write.csv(SP9.200,"SP9.200.csv")
######## Saving all SP graphs 200 nodes data to csv
SP1=read.csv("SP1.200.csv");SP2=read.csv("SP2.200.csv");SP3=read.csv("SP3.200.csv")
SP4=read.csv("SP4.200.csv");SP5=read.csv("SP5.200.csv");SP6=read.csv("SP6.200.csv")
SP7=read.csv("SP7.200.csv");SP8=read.csv("SP8.200.csv");SP9=read.csv("SP9.200.csv")

df.SP.200=rbind(SP1,SP2,SP3,SP4,SP5,SP6,SP7,SP8,SP9)
write.csv(df.SP.200,"df.SP.200.csv")

### 250 nodes 
SP1.250=sim.sp.net(r=rep(0.1,250),n=250)
write.csv(SP1.250,"SP1.250.csv")
SP2.250=sim.sp.net(r=rep(0.2,250),n=250)
write.csv(SP2.250,"SP2.250.csv")
SP3.250=sim.sp.net(r=rep(0.3,250),n=250)
write.csv(SP3.250,"SP3.250.csv")
SP4.250=sim.sp.net(r=rep(0.4,250),n=250)
write.csv(SP4.250,"SP4.250.csv")
SP5.250=sim.sp.net(r=rep(0.5,250),n=250)
write.csv(SP5.250,"SP5.250.csv")
SP6.250=sim.sp.net(r=rep(0.6,250),n=250)
write.csv(SP6.250,"SP6.250.csv")
SP7.250=sim.sp.net(r=rep(0.7,250),n=250)
write.csv(SP7.250,"SP7.250.csv")
SP8.250=sim.sp.net(r=rep(0.8,250),n=250)
write.csv(SP8.250,"SP8.250.csv")
SP9.250=sim.sp.net(r=rep(0.9,250),n=250)
write.csv(SP9.250,"SP9.250.csv")
######## Saving all SP graphs 250 nodes data to csv
SP1=read.csv("SP1.250.csv");SP2=read.csv("SP2.250.csv");SP3=read.csv("SP3.250.csv")
SP4=read.csv("SP4.250.csv");SP5=read.csv("SP5.250.csv");SP6=read.csv("SP6.250.csv")
SP7=read.csv("SP7.250.csv");SP8=read.csv("SP8.250.csv");SP9=read.csv("SP9.250.csv")

df.SP.250=rbind(SP1,SP2,SP3,SP4,SP5,SP6,SP7,SP8,SP9)
write.csv(df.SP.250,"df.SP.250.csv")

### 300 nodes 
SP1.300=sim.sp.net(r=rep(0.1,250),n=300)
write.csv(SP1.300,"SP1.300.csv")
SP2.300=sim.sp.net(r=rep(0.2,250),n=300)
write.csv(SP2.300,"SP2.300.csv")
SP3.300=sim.sp.net(r=rep(0.3,250),n=300)
write.csv(SP3.300,"SP3.300.csv")
SP4.300=sim.sp.net(r=rep(0.4,250),n=300)
write.csv(SP4.300,"SP4.300.csv")
SP5.300=sim.sp.net(r=rep(0.5,250),n=300)
write.csv(SP5.300,"SP5.300.csv")
SP6.300=sim.sp.net(r=rep(0.6,250),n=300)
write.csv(SP6.300,"SP6.300.csv")
SP7.300=sim.sp.net(r=rep(0.7,250),n=300)
write.csv(SP7.300,"SP7.300.csv")
SP8.300=sim.sp.net(r=rep(0.8,250),n=300)
write.csv(SP8.300,"SP8.300.csv")
SP9.300=sim.sp.net(r=rep(0.9,250),n=300)
write.csv(SP9.300,"SP9.300.csv")
######## Saving all SP graphs 300 nodes data to csv
SP1=read.csv("SP1.300.csv");SP2=read.csv("SP2.300.csv");SP3=read.csv("SP3.300.csv")
SP4=read.csv("SP4.300.csv");SP5=read.csv("SP5.300.csv");SP6=read.csv("SP6.300.csv")
SP7=read.csv("SP7.300.csv");SP8=read.csv("SP8.300.csv");SP9=read.csv("SP9.300.csv")

df.SP.300=rbind(SP1,SP2,SP3,SP4,SP5,SP6,SP7,SP8,SP9)
write.csv(df.SP.300,"df.SP.300.csv")

### 350 nodes 
SP1.350=sim.sp.net(r=rep(0.1,250),n=350)
write.csv(SP1.350,"SP1.350.csv")
SP2.350=sim.sp.net(r=rep(0.2,250),n=350)
write.csv(SP2.350,"SP2.350.csv")
SP3.350=sim.sp.net(r=rep(0.3,250),n=350)
write.csv(SP3.350,"SP3.350.csv")
SP4.350=sim.sp.net(r=rep(0.4,250),n=350)
write.csv(SP4.350,"SP4.350.csv")
SP5.350=sim.sp.net(r=rep(0.5,250),n=350)
write.csv(SP5.350,"SP5.350.csv")
SP6.350=sim.sp.net(r=rep(0.6,250),n=350)
write.csv(SP6.350,"SP6.350.csv")
SP7.350=sim.sp.net(r=rep(0.7,250),n=350)
write.csv(SP7.350,"SP7.350.csv")
SP8.350=sim.sp.net(r=rep(0.8,250),n=350)
write.csv(SP8.350,"SP8.350.csv")
SP9.350=sim.sp.net(r=rep(0.9,250),n=350)
write.csv(SP9.350,"SP9.350.csv")
######## Saving all SP graphs 350 nodes data to csv
SP1=read.csv("SP1.350.csv");SP2=read.csv("SP2.350.csv");SP3=read.csv("SP3.350.csv")
SP4=read.csv("SP4.350.csv");SP5=read.csv("SP5.350.csv");SP6=read.csv("SP6.350.csv")
SP7=read.csv("SP7.350.csv");SP8=read.csv("SP8.350.csv");SP9=read.csv("SP9.350.csv")

df.SP.350=rbind(SP1,SP2,SP3,SP4,SP5,SP6,SP7,SP8,SP9)
write.csv(df.SP.350,"df.SP.350.csv")

### 400 nodes 
SP1.400=sim.sp.net(r=rep(0.1,250),n=400)
write.csv(SP1.400,"SP1.400.csv")
SP2.400=sim.sp.net(r=rep(0.2,250),n=400)
write.csv(SP2.400,"SP2.400.csv")
SP3.400=sim.sp.net(r=rep(0.3,250),n=400)
write.csv(SP3.400,"SP3.400.csv")
SP4.400=sim.sp.net(r=rep(0.4,250),n=400)
write.csv(SP4.400,"SP4.400.csv")
SP5.400=sim.sp.net(r=rep(0.5,250),n=400)
write.csv(SP5.400,"SP5.400.csv")
SP6.400=sim.sp.net(r=rep(0.6,250),n=400)
write.csv(SP6.400,"SP6.400.csv")
SP7.400=sim.sp.net(r=rep(0.7,250),n=400)
write.csv(SP7.400,"SP7.400.csv")
SP8.400=sim.sp.net(r=rep(0.8,250),n=400)
write.csv(SP8.400,"SP8.400.csv")
SP9.400=sim.sp.net(r=rep(0.9,250),n=400)
write.csv(SP9.400,"SP9.400.csv")
######## Saving all SP graphs 400 nodes data to csv
SP1=read.csv("SP1.400.csv");SP2=read.csv("SP2.400.csv");SP3=read.csv("SP3.400.csv")
SP4=read.csv("SP4.400.csv");SP5=read.csv("SP5.400.csv");SP6=read.csv("SP6.400.csv")
SP7=read.csv("SP7.400.csv");SP8=read.csv("SP8.400.csv");SP9=read.csv("SP9.400.csv")

df.SP.400=rbind(SP1,SP2,SP3,SP4,SP5,SP6,SP7,SP8,SP9)
write.csv(df.SP.400,"df.SP.400.csv")

### 500 nodes 
SP1.500=sim.sp.net(r=rep(0.1,250),n=500)
write.csv(SP1.500,"SP1.500.csv")
SP2.500=sim.sp.net(r=rep(0.2,250),n=500)
write.csv(SP2.500,"SP2.500.csv")
SP3.500=sim.sp.net(r=rep(0.3,250),n=500)
write.csv(SP3.500,"SP3.500.csv")
SP4.500=sim.sp.net(r=rep(0.4,250),n=500)
write.csv(SP4.500,"SP4.500.csv")
SP5.500=sim.sp.net(r=rep(0.5,250),n=500)
write.csv(SP5.500,"SP5.500.csv")
SP6.500=sim.sp.net(r=rep(0.6,250),n=500)
write.csv(SP6.500,"SP6.500.csv")
SP7.500=sim.sp.net(r=rep(0.7,250),n=500)
write.csv(SP7.500,"SP7.500.csv")
SP8.500=sim.sp.net(r=rep(0.8,250),n=500)
write.csv(SP8.500,"SP8.500.csv")
SP9.500=sim.sp.net(r=rep(0.9,250),n=500)
write.csv(SP9.500,"SP9.500.csv")
######## Saving all SP graphs 500 nodes data to csv
SP1=read.csv("SP1.500.csv");SP2=read.csv("SP2.500.csv");SP3=read.csv("SP3.500.csv")
SP4=read.csv("SP4.500.csv");SP5=read.csv("SP5.500.csv");SP6=read.csv("SP6.500.csv")
SP7=read.csv("SP7.500.csv");SP8=read.csv("SP8.500.csv");SP9=read.csv("SP9.500.csv")

df.SP.500=rbind(SP1,SP2,SP3,SP4,SP5,SP6,SP7,SP8,SP9)
write.csv(df.SP.500,"df.SP.500.csv")

### 750 nodes 
SP1.750=sim.sp.net(r=rep(0.1,250),n=750)
write.csv(SP1.750,"SP1.750.csv")
SP2.750=sim.sp.net(r=rep(0.2,250),n=750)
write.csv(SP2.750,"SP2.750.csv")
SP3.750=sim.sp.net(r=rep(0.3,250),n=750)
write.csv(SP3.750,"SP3.750.csv")
SP4.750=sim.sp.net(r=rep(0.4,250),n=750)
write.csv(SP4.750,"SP4.750.csv")
SP5.750=sim.sp.net(r=rep(0.5,250),n=750)
write.csv(SP5.750,"SP5.750.csv")
SP6.750=sim.sp.net(r=rep(0.6,250),n=750)
write.csv(SP6.750,"SP6.750.csv")
SP7.750=sim.sp.net(r=rep(0.7,250),n=750)
write.csv(SP7.750,"SP7.750.csv")
SP8.750=sim.sp.net(r=rep(0.8,250),n=750)
write.csv(SP8.750,"SP8.750.csv")
SP9.750=sim.sp.net(r=rep(0.9,250),n=750)
write.csv(SP9.750,"SP9.750.csv")
######## Saving all SP graphs 750 nodes data to csv
SP1=read.csv("SP1.750.csv");SP2=read.csv("SP2.750.csv");SP3=read.csv("SP3.750.csv")
SP4=read.csv("SP4.750.csv");SP5=read.csv("SP5.750.csv");SP6=read.csv("SP6.750.csv")
SP7=read.csv("SP7.750.csv");SP8=read.csv("SP8.750.csv");SP9=read.csv("SP9.750.csv")

df.SP.750=rbind(SP1,SP2,SP3,SP4,SP5,SP6,SP7,SP8,SP9)
write.csv(df.SP.750,"df.SP.750.csv")

### 1000 nodes 
SP1.1000=sim.sp.net(r=rep(0.1,250),n=1000)
write.csv(SP1.1000,"SP1.1000.csv")
SP2.1000=sim.sp.net(r=rep(0.2,250),n=1000)
write.csv(SP2.1000,"SP2.1000.csv")
SP3.1000=sim.sp.net(r=rep(0.3,250),n=1000)
write.csv(SP3.1000,"SP3.1000.csv")
SP4.1000=sim.sp.net(r=rep(0.4,250),n=1000)
write.csv(SP4.1000,"SP4.1000.csv")
SP5.1000=sim.sp.net(r=rep(0.5,250),n=1000)
write.csv(SP5.1000,"SP5.1000.csv")
SP6.1000=sim.sp.net(r=rep(0.6,250),n=1000)
write.csv(SP6.1000,"SP6.1000.csv")
SP7.1000=sim.sp.net(r=rep(0.7,250),n=1000)
write.csv(SP7.1000,"SP7.1000.csv")
SP8.1000=sim.sp.net(r=rep(0.8,250),n=1000)
write.csv(SP8.1000,"SP8.1000.csv")
SP9.1000=sim.sp.net(r=rep(0.9,250),n=1000)
write.csv(SP9.1000,"SP9.1000.csv")
######## Saving all SP graphs 1000 nodes data to csv
SP1=read.csv("SP1.1000.csv");SP2=read.csv("SP2.1000.csv");SP3=read.csv("SP3.1000.csv")
SP4=read.csv("SP4.1000.csv");SP5=read.csv("SP5.1000.csv");SP6=read.csv("SP6.1000.csv")
SP7=read.csv("SP7.1000.csv");SP8=read.csv("SP8.1000.csv");SP9=read.csv("SP9.1000.csv")

df.SP.1000=rbind(SP1,SP2,SP3,SP4,SP5,SP6,SP7,SP8,SP9)

write.csv(df.SP.1000,"df.SP.1000.csv")


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Random graph
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+# (1) Function to generate a single Random graph

erdos<-function(n,p){
  erdos_graph=sample_gnp(n, p, directed = FALSE, loops = FALSE)
}

##---Test case----
#f1=erdos(10,0.6)
#plot(f1)

# #Function for Replicating the random network N times
# set.seed(1)
# erdos_graph=list()
# Erdos<-function(n,p,N){
#   for (i in 1:N){
#     erdos_graph[[i]]=sample_gnp(n,p, directed = FALSE, loops =FALSE)
#     components(erdos_graph[[i]], mode = c( "strong"))
#   }
#   return(erdos_graph)
# }

MKERGraphs<- function(n=25,p=0.8) {
  Graphs = list()#;G=list() # set up the list of Graphs first
  initgraph= list()
  i= 1
  print("Creating ER Graphs")
  Graphs[[i]] = erdos(n,p)
  
  Graphs[[i]]$type = "ER"
  Graphs[[i]]$id = "1"
  Graphs[[i]]$name="ER"
  i <- i+1
  return(Graphs)
}


sim.erdos.net=function(p=c(0.3,0.4,0.6,0.05),n=100){
  net=NULL;data=NULL
  k=1
  for (i in 1:length(p)){
    net[i]= MKERGraphs(n,p[i])
    data=RunSimOnGraphFeatures(net,nreps = 1)
    k=k+1
  }
  df=cbind(p,data)
  return(df)
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Scale-Free
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Function to generate a single scale-free graph
scale_free<-function(n,pw,m){
  scale_free_graph=sample_pa(n,pw,m, directed= F)
}

#Function for Replicating the scale-free network N times

# set.seed(1)
# scale_free_net=list()
# Scale_Free<-function(n,pw,m,N){
#   for (i in 1:N){
#     scale_free_net[[i]]=sample_pa(n,pw,m, directed= F)
#   }
#   return(scale_free_net)
# }
MKSFGraphs<- function(n=50, power=2, m=4) {
  Graphs = list()#;G=list() # set up the list of Graphs first
  initgraph= list()
  i= 1
  print("Creating SF Graphs")
  Graphs[[i]] = sample_pa(n, power, m, directed=FALSE, algorithm="psumtree")
  
  Graphs[[i]]$type = "SF"
  Graphs[[i]]$id = "1"
  Graphs[[i]]$name="SF"
  i <- i+1
  return(Graphs)
}


sim.sf.net=function(power=powerVals,m=mVals,n=500,nreps=1){
  net=NULL;data=NULL
  dt=expand.grid(powerVals,mVals)
  colnames(dt)=c("powerVals","mVals")
  for (i in 1:nrow(dt)){
    net[i]= MKSFGraphs(n,power=dt$powerVals[i],m=dt$mVals[i])
    data=RunSimOnGraphFeatures(net,nreps)
  }
  df=cbind(power,m,data)  
  return(df)
}
powerVals=1:10
mVals=1:100

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Small world
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Function to generate a single small world graph
sm_world<-function(d,s,n,p){
  sm_graph=sample_smallworld(d, s, n, p, loops = F, multiple = F)
}

#Function for Replicating the small world network N times
# set.seed(1)
# sm_graph=list()
# SM_World<-function(d,s,n,p,N){
#   for (i in 1:N){
#     sm_graph[[i]]=sample_smallworld(d, s, n, p, loops = F, multiple = F)
#   }
#   return(sm_graph)
# }

MKSWGraphs<- function(n=50,nei=6,p=0.2) {
  Graphs = list()#;G=list() # set up the list of Graphs first
  initgraph= list()
  i= 1
  print("Creating SW Graphs")
  Graphs[[i]] = sample_smallworld(dim=1, size=n, nei=nei, p=p)
  
  Graphs[[i]]$type = "SW"
  Graphs[[i]]$id = "1"
  Graphs[[i]]$name="SW"
  i <- i+1
  return(Graphs)
}


sim.sw.net=function(p=pVals,nei=neiVals,n=50,nreps=1){
  net=NULL;data=NULL
  dat=expand.grid(pVals,neiVals)
  colnames(dat)=c("prewireVals","neiVals")
  for (i in 1:nrow(dat)){
    net[i]= MKSWGraphs(n,p=dat$prewireVals[i],nei=dat$neiVals[i])
    data=RunSimOnGraphFeatures(net,nreps)
  }
  df=cbind(p,nei,data)  
  return(df)
}

pVals=c(0.001,0.005,0.01,0.03,0.05,0.07,0.09,0.1,0.2,0.3)
neiVals=1:10

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Spatial graph
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Function for a sngle  Spatial network
#Spatial graph of x,y coordinate join in space
# spatial <- function(n,r) {
#   v <- NULL
#   for (i in 1:n) {
#     x <- runif(1)
#     y <- runif(1)
#     v <- rbind(v,c(i,x,y))
#     coords<<-rbind(coords, c(x,y))
#   }
#   
#   E <- NULL
#   R2 <- r*r;
#   for (i in seq(1,n-1)) {
#     for (j in seq(i+1,n)) {
#       if (((v[i,2]-v[j,2])^2 + (v[i,3]-v[j,3])^2) < R2) {
#         E <- as.data.frame(rbind(E,c(i,j)))
#       }
#     }
#   }
#   return(E)
# }
# 
# coords=NULL
# 
# E = spatial(10,0.3)
# 
# G=graph_from_data_frame(E, directed=FALSE)
# 
# V(G)$color<-"grey"
# G$layout=coords
# 
# ##--Test case---
# plot(G)

## (2) Replicating the Spatial network N times
# set.seed(1)
# spatial_net=list()
# SPATIAL<-function(n,r,N){
#   for (i in 1:N){
#     spatial_net[[i]]=spatial(n,r)
#   }
#   return(spatial_net)
# }
# h1=SPATIAL(20,.4,5)
#h2=graph_from_data_frame(h1[[5]], directed=FALSE)
#plot(h2)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Lattice graph
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Function to generate a single lattice graph
#make_lattice(length = ceiling(sqrt(order)),dim=dim.lat,nei = nei.lat)
lattice_net<-function(n,nei){
  lattice_graph=make_lattice(length=floor(sqrt(n)), dim = 2,nei = nei,directed= F, mutual= T, circular = F)
}

#Function for Replicating the lattice network N times

 # lattice_graph=list()
 # Lattice_net<-function(d,n,nei,N){
 #   for (i in 1:N){
 #     lattice_graph[[i]]=lattice_net(d,n,nei)
 #   }
 #   return(lattice_graph)
 # }

##---Test case----
#f3=Lattice(2,5)
#plot(f3)

MKLATGraphs<- function(n=50,nei=3) {
   Graphs = list()#;G=list() # set up the list of Graphs first
   initgraph= list()
   i= 1
   print("Creating LAT Graphs they are all the same")
   Graphs[[i]] = lattice_net(n,nei)
   Graphs[[i]]$type = "LAT"
   Graphs[[i]]$id = "1"
   Graphs[[i]]$name="LAT"
   i <- i+1
   return(Graphs)
 }


sim.lat.net=function(n=30,nei=c(1,2,4)){
  net=NULL;data=NULL
  for (i in 1:length(nei)){
    net[i]= MKLATGraphs(n,nei=nei[i])
    data=RunSimOnGraphFeatures(net,nreps=1)
  }
  df=cbind(nei,data)  
  return(df)
}

nei=1:20
### 50 nodes 
lat.50=sim.lat.net(nei,n=50)

write.csv(lat.50,"lat.50.csv")

### 100 nodes 
lat.100=sim.lat.net(nei,n=100)

write.csv(lat.100,"lat.100.csv")

### 150 nodes 
lat.150=sim.lat.net(nei,n=150)
write.csv(lat.150,"lat.150.csv")

### 200 nodes 
lat.200=sim.lat.net(nei,n=200)
write.csv(lat.200,"lat.200.csv")

### 250 nodes 
lat.250=sim.lat.net(nei,n=250)
write.csv(lat.250,"lat.250.csv")

### 300 nodes 
lat.300=sim.lat.net(nei,n=300)
write.csv(lat.300,"lat.300.csv")

### 350 nodes 
lat.350=sim.lat.net(nei,n=350)
write.csv(lat.350,"lat.350.csv")

### 400 nodes 
lat.400=sim.lat.net(nei,n=400)
write.csv(lat.400,"lat.400.csv")

### 500 nodes 
lat.500=sim.lat.net(nei,n=500)
write.csv(lat.500,"lat.500.csv")

### 750 nodes 
lat.750=sim.lat.net(nei,n=750)
write.csv(lat.750,"lat.750.csv")

### 1000 nodes 
lat.1000=sim.lat.net(nei,n=1000)
write.csv(lat.1000,"lat.1000.csv")

### 1500 nodes 
lat.1500=sim.lat.net(nei,n=1500)
write.csv(lat.1500,"lat.1500.csv")

### 2000 nodes 
lat.2000=sim.lat.net(nei,n=2000)
write.csv(lat.2000,"lat.2000.csv")

######## Saving all lat gradhs to csv
lat.50=read.csv("lat.50.csv");lat.100=read.csv("lat.100.csv");
lat.150=read.csv("lat.150.csv")
lat.200=read.csv("lat.100.csv");
lat.250=read.csv("lat.250.csv");lat.300=read.csv("lat.300.csv")
lat.350=read.csv("lat.350.csv");lat.400=read.csv("lat.400.csv")
lat.500=read.csv("lat.500.csv");lat.750=read.csv("lat.750.csv")
lat.1000=read.csv("lat.1000.csv");lat.1500=read.csv("lat.1500.csv")
lat.2000=read.csv("lat.2000.csv")

df.lat=rbind(lat.50,lat.150,lat.200,lat.250,lat.300,lat.350,
            lat.400, lat.500,lat.750,lat.1000,lat.1500,
            lat.2000)

write.csv(df.lat,"df.lat.csv")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Stochastic block model
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Two groups with not only few connection between groups
#pm <- cbind(c(.05, .002), c(.002, .05))
sbm_net <- function(pm,n){
  sbm.graph=sample_sbm(n, pref.matrix = pm, block.sizes = c(0.4*n, 0.6*n))
}

#plot(g)

MKSBMGraphs<- function(pm,n) {
  Graphs = list()#;G=list() # set up the list of Graphs first
  initgraph= list()
  i= 1
  print("Creating Stochastic block matrix")
  Graphs[[i]] = sbm_net(pm,n)
  Graphs[[i]]$type = "sbm"
  Graphs[[i]]$id = "1"
  Graphs[[i]]$name="sbm"
  i <- i+1
  return(Graphs)
}


sim.sbm.net=function(pm,n=50){
    net= lapply(pm, MKSBMGraphs,n)
    data=lapply(net,RunSimOnGraphFeatures,nreps=1)
    df <-  as.data.frame(do.call(rbind, data))
  return(df)
}

numofComm=2
pm=list();mat=list()
for (i in 1:10000){
  mat[[i]]<- matrix(runif(numofComm^2), numofComm, numofComm)
  pm[[i]]=mat[[i]]*t(mat[[i]])
}


### 50 nodes 
sbm.50=sim.sbm.net(pm,n=50)

write.csv(sbm.50,"sbm.50.csv")

### 100 nodes 
sbm.100=sim.sbm.net(pm,n=100)
write.csv(sbm.100,"sbm.100.csv")

### 150 nodes 
sbm.150=sim.sbm.net(pm,n=150)
write.csv(sbm.150,"sbm.150.csv")

### 200 nodes 
sbm.200=sim.sbm.net(pm,n=200)
write.csv(sbm.200,"sbm.200.csv")

### 250 nodes 
sbm.250=sim.sbm.net(pm,n=250)
write.csv(sbm.250,"sbm.250.csv")

### 300 nodes 
sbm.300=sim.sbm.net(pm,n=300)
write.csv(sbm.300,"sbm.300.csv")

### 350 nodes 
sbm.350=sim.sbm.net(pm,n=350)
write.csv(sbm.350,"sbm.350.csv")

### 400 nodes 
sbm.400=sim.sbm.net(pm,n=400)
write.csv(sbm.400,"sbm.400.csv")

### 500 nodes 
sbm.500=sim.sbm.net(pm,n=500)
write.csv(sbm.500,"sbm.500.csv")

### 750 nodes 
sbm.750=sim.sbm.net(pm,n=750)
write.csv(sbm.750,"sbm.750.csv")

### 1000 nodes 
sbm.1000=sim.sbm.net(pm,n=1000)
write.csv(sbm.1000,"sbm.1000.csv")

### 1500 nodes 
sbm.1500=sim.sbm.net(pm,n=1500)
write.csv(sbm.1500,"sbm.1500.csv")

### 2000 nodes 
sbm.2000=sim.sbm.net(pm,n=2000)
write.csv(sbm.2000,"sbm.2000.csv")

######## Saving all sbm gradhs to csv
sbm.50=read.csv("sbm.50.csv");sbm.100=read.csv("sbm.100.csv");
sbm.150=read.csv("sbm.150.csv")
sbm.200=read.csv("sbm.100.csv");
sbm.250=read.csv("sbm.250.csv");sbm.300=read.csv("sbm.300.csv")
sbm.350=read.csv("sbm.350.csv");sbm.400=read.csv("sbm.400.csv")
sbm.500=read.csv("sbm.500.csv");sbm.750=read.csv("sbm.750.csv")
sbm.1000=read.csv("sbm.1000.csv");sbm.1500=read.csv("sbm.1500.csv")
sbm.2000=read.csv("sbm.2000.csv")

df.sbm=rbind(sbm.50,sbm.100,sbm.150,sbm.200,sbm.250,sbm.300,sbm.350,
             sbm.400, sbm.500,sbm.750,sbm.1000,sbm.1500,
             sbm.2000)

write.csv(df.sbm,"df.sbm.csv")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Forest fire
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# g <- sample_forestfire(100, fw.prob = 0.37, bw.factor = 0.32 / 0.37,directed = F)
# dd1 <- degree_distribution(g, mode = "in")
# dd2 <- degree_distribution(g, mode = "out")
# plot(seq(along.with = dd1) - 1, dd1, log = "xy")
# points(seq(along.with = dd2) - 1, dd2, col = 2, pch = 2)
# 
# plot(g)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Data processing and cleaning
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
df1=read.csv("df.lat.csv")%>%
  filter(order<=1000)


df1=df1%>%dplyr::select(c(GraphName,order,edges,
                          mean_eccentr,mean_path_length,graph_energy,
                          modularity,diameter,betw_centr,transitivity,
                          spectral_radius,eigen_centr,deg_centr,
                          mean_degree,minCut,FiedlerValue,Normalized_FiedlerValue,
                          closeness_centr))%>%
  mutate_if(is.character,factor)

# df2=rbind(read.csv("sbm.50.csv",sep = ",", header = T),
#           read.csv("sbm.100.csv",sep = ",", header = T),
#           read.csv("sbm.150.csv",sep = ",", header = T),
#           read.csv("sbm.150.csv",sep = ",", header = T))

df2=read.csv("df.sbm.csv")%>%
  filter(order<=1000)


df2=df2%>%dplyr::select(c(GraphName,order,edges,
                          mean_eccentr,mean_path_length,graph_energy,
                          modularity,diameter,betw_centr,transitivity,
                          spectral_radius,eigen_centr,deg_centr,
                          mean_degree,minCut,FiedlerValue,Normalized_FiedlerValue,
                          closeness_centr))%>%
  mutate_if(is.character,factor)

df3=rbind(read.csv("df.er.50.csv",sep = ",", header = T),
          read.csv("df.er.100.csv",sep = ",", header = T),
          read.csv("df.er.150.csv",sep = ",", header = T),
          read.csv("df.er.200.csv",sep = ",", header = T),
          read.csv("df.er.250.csv",sep = ",", header = T),
          read.csv("df.er.350.csv",sep = ",", header = T),
          read.csv("df.er.400.csv",sep = ",", header = T),
          read.csv("df.er.500.csv",sep = ",", header = T),
          read.csv("df.er.750.csv",sep = ",", header = T))


df3=df3%>%dplyr::select(c(GraphName,order,edges,
                          mean_eccentr,mean_path_length,graph_energy,
                          modularity,diameter,betw_centr,transitivity,
                          spectral_radius,eigen_centr,deg_centr,
                          mean_degree,minCut,FiedlerValue,Normalized_FiedlerValue,
                          closeness_centr))%>%
  mutate_if(is.character,factor)


df4=rbind(read.csv("df.SP.50.csv",sep = ",", header = T),
          read.csv("df.SP.100.csv",sep = ",", header = T),
          read.csv("df.SP.150.csv",sep = ",", header = T),
          read.csv("df.SP.200.csv",sep = ",", header = T),
          read.csv("df.SP.250.csv",sep = ",", header = T),
          read.csv("df.SP.300.csv",sep = ",", header = T),
          read.csv("df.SP.350.csv",sep = ",", header = T))

df4=df4%>%dplyr::select(c(GraphName,order,edges,
                          mean_eccentr,mean_path_length,graph_energy,
                          modularity,diameter,betw_centr,transitivity,
                          spectral_radius,eigen_centr,deg_centr,
                          mean_degree,minCut,FiedlerValue,Normalized_FiedlerValue,
                          closeness_centr))%>%
  mutate_if(is.character,factor)

df5=rbind(read.csv("SF1.50.csv",sep = ",", header = T),
          read.csv("SF1.100.csv",sep = ",", header = T),
          read.csv("SF1.150.csv",sep = ",", header = T),
          read.csv("SF1.200.csv",sep = ",", header = T))

df5=df5%>%rename("power"="prewireVals","mvals"="neiVals")

df5=df5%>%
  filter(power %in% c("2") & GraphReplicate <= 250)

df5=df5%>%dplyr::select(c(GraphName,order,edges,
                          mean_eccentr,mean_path_length,graph_energy,
                          modularity,diameter,betw_centr,transitivity,
                          spectral_radius,eigen_centr,deg_centr,
                          mean_degree,minCut,FiedlerValue,Normalized_FiedlerValue,
                          closeness_centr))%>%
  mutate_if(is.character,factor)


df6=rbind(read.csv("SW.50.csv",sep = ",", header = T),
          read.csv("SW.100.csv",sep = ",", header = T),
          read.csv("SW.150.csv",sep = ",", header = T),
          read.csv("SW.200.csv",sep = ",", header = T),
          read.csv("SW.250.csv",sep = ",", header = T),
          read.csv("SW.300.csv",sep = ",", header = T))


df6=df6%>%
  filter(prewireVals %in% c("0.2") & GraphReplicate <= 250)
#df6[df6$GraphReplicate <= 250,]
df6=df6%>%dplyr::select(c(GraphName,order,edges,
                          mean_eccentr,mean_path_length,graph_energy,
                          modularity,diameter,betw_centr,transitivity,
                          spectral_radius,eigen_centr,deg_centr,
                          mean_degree,minCut,FiedlerValue,Normalized_FiedlerValue,
                          closeness_centr))%>%
  mutate_if(is.character,factor)




Data=rbind(df1,df2,df3,df4,df5,df6)

write.csv(Data,"Data.CSV")
##--Shuffle--data
Data = Data[sample(1:nrow(Data)), ]##shuffle row indices and randomly re order 
Data%>%na.omit()
  #replace(is.na(.), 0)#change NA's in deg_assort_coef to 0

# Change data frame to tibble (similar to data-frame)
df=as_tibble(Data)
df%>%count(GraphName,sort=T)

# Show first n row of the data
df%>%
  slice_head(n=10)
#df%>%View()

###--------Load the janitor package to clean the data---------###
df<-df%>%
#mutate(GraphName=factor(GraphName,levels=c("ER","sbm","LAT")))%>%
  clean_names()

df%>%
  glimpse()

##--Visualizing all data structure
vis_dat(df)

df=df%>%na.omit()

##--Visualizing all data structure
vis_dat(df)

##--Count the outcome variable
df %>%
  dplyr::count(graph_name)

#or
skim(df)

##---The function ggscatmat from the package GGally creates a matrix with scatterplots, densities and correlations for numeric columns
#ggscatmat(df, alpha=0.2)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                Find correlations between variables(predictors)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##----Correlation plot
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# corell.plot=pairs.panels(df[,-1],
#                          gap = 0,
#                          bg = c("red", "yellow", "blue")[df$graph_name],
#                          smooth = TRUE,
#                          stars = TRUE, 
#                          ci=TRUE,
#                          pch=24)
# corell.plot

###----Convert data from wide to long format----for----sub plots####
#we use facet-wrap() function for subplot of data and pivot_longer() to change data structure
theme_set(theme_light())
df_long<-df%>%
  pivot_longer(!graph_name, names_to="graph_features",
               values_to = "values")

###----Box plot for each predictor variable----####
df_long%>%
  ggplot(mapping = aes(x=graph_name,y=values,fill=graph_features))+
  geom_boxplot()+
  facet_wrap(~graph_features,scales = "free",ncol = 6)+
  scale_colour_viridis_d(option = "plasma",end = .9)+
  theme(legend.position = "none")


#########################################################################################
# Feature Engineering
#########################################################################################
df=df%>%na.omit()

##Preparing data
df.rec <- recipe(graph_name ~ ., data = df) %>%
#  recipes::step_nzv(all_nominal()) %>%
  recipes::step_zv(all_predictors()) %>%
  recipes::step_normalize(all_numeric_predictors())%>%
  step_smote(graph_name) 
  #step_dummy(all_predictors())
  
  
 
df.prep=df.rec%>%prep()

df.juice=df.prep%>%juice()%>%clean_names()



##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#----Data Partitioning----###
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
set.seed(384)
##sample and split
df_split=rsample::initial_split(df.juice, strata = "graph_name",prop = 0.7)

df.train <- rsample::training(df_split)
df.test <- rsample::testing(df_split)

#train.smote = DMwR::SMOTE(graph_name ~ ., df.train, perc.over = 100, perc.under=200)

###----Print number of training and testing set----###
cat("training set:", nrow(df.train), "\n",
    "testing set :", nrow(df.test), "\n")


#####################################################
# Boruta feature scaling/selection
#####################################################

## first method via tidymodels
df.rec.boruta <- 
  recipe(graph_name ~ ., data = df)%>%
  #update_role(order,new_role = "nodes")%>%
  #step_impute_mean(all_numeric(),-all_outcomes())%>%
  #step_nzv(all_numeric(),-all_outcomes())%>%
  #step_corr(all_numeric(),-all_outcomes(), threshold = 0.7)%>%
  step_select_boruta(all_predictors(), outcome = "graph_name")


preproc_boruta_data=df.rec.boruta  %>% 
  prep() %>% 
  juice()

preproc_boruta_data%>%
  clean_names()


## second method
set.seed(9532)
boruta.train <- Boruta(graph_name~., data = df.train, doTrace = 2)
print(boruta.train)

#Plot boruta object showing var importance
plot(boruta.train, xlab = "", xaxt = "n")

lz<-lapply(1:ncol(boruta.train$ImpHistory),function(i)
  boruta.train$ImpHistory[is.finite(boruta.train$ImpHistory[,i]),i])
names(lz) <- colnames(boruta.train$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(boruta.train$ImpHistory), cex.axis = 0.7)

#Dealing with the tentative attributes
final.boruta <- TentativeRoughFix(boruta.train)
print(final.boruta)


# list of confirmed selected attributes
selected.attr=getSelectedAttributes(final.boruta, withTentative = F)
selected.attr

# final imortance score of each predictor variable with boruta
boruta.df <- attStats(final.boruta)
#boruta.df$selected_features <- row.names(boruta.df)   

# Selecting top N predictors values by group
boruta.df.top <- boruta.df %>%
  select(-c(normHits,decision)) %>% 
  arrange(desc(meanImp),desc(medianImp),desc(minImp),desc(maxImp))%>%
  slice(1:6)

boruta.df.top

#########################################################################################
#+ The models
#########################################################################################

## Create recipe with selected features for the model
boruta.final.data <- 
  recipe(graph_name ~ ., data = df)%>%# all feature selected by boruta
  prep() %>% 
  juice() %>% 
  clean_names()

boruta.final.data

################################################################################
#+ Creating final data for partitioning
################################################################################
set.seed(384)
##sample and split
df_split=rsample::initial_split(df, strata = "graph_name",prop = 0.7)

df.train <- rsample::training(df_split)
df.test <- rsample::testing(df_split)

###----Print number of training and testing set----###
cat("training set:", nrow(df.train), "\n",
    "testing set :", nrow(df.test), "\n")

##----Creating Cross Validation set
df.cv.splits <- vfold_cv(df.train, v = 10)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#----Build a Neural net classification model----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The model has one penalty tuning parameter which is the amount of regularization
## Multinormial regression via neural net
model_spec=multinom_reg(
  penalty = 1,
  engine = "nnet",
  mode = "classification"
)

##--Receipe
model_rec=recipe(graph_name~., data = df.train)%>%
  step_normalize(all_numeric_predictors())

##--Workflow
model_wkflw=workflow(preprocessor = model_rec,spec = model_spec)

##--FitWorkflow
model_fit=model_wkflw%>%
  fit(data = df.train)

##--PrintWkflw object
model_fit

###----Evaluate Model----####
#Prediction on test set
pred_result=augment(model_fit,df.test)

predVals=pred_result%>%
  slice_head(n=10)

#View(predVals)

##Confusion Matrix
pred_result%>%
  conf_mat(truth = graph_name,estimate =.pred_class)
#%>%summary()

#-Visualize conf matrix
update_geom_defaults(geom = "tile", new=list(color="black",alpha=0.5))

pred_result%>%
  conf_mat(truth = graph_name,estimate =.pred_class)%>%
  autoplot(type="heatmap")

###----Statistical summary of confusion matrix----###
Accuracy=yardstick::accuracy(data=pred_result,truth = graph_name,estimate =.pred_class)
Precision=yardstick::precision(data=pred_result,truth = graph_name,estimate =.pred_class)
Recall=yardstick::recall(data=pred_result,truth = graph_name,estimate =.pred_class)

Accuracy%>%
  bind_rows(Precision)%>%
  bind_rows(Recall)

##--ROC CURVE
nnet.roc=pred_result%>%
  roc_curve(graph_name,c(.pred_LAT,.pred_sbm,.pred_ER,.pred_Spatial,.pred_SF,.pred_SW))%>%
  ggplot(aes(x=1-specificity,y=sensitivity, color=.level))+
  geom_abline(lty=2,color="gray80",size=.9)+
  geom_path(show.legend = T,alpha=0.6,size=1.2)+
  coord_equal()

nnet.roc

##ROC accuracy score
pred_result%>%
  roc_auc(graph_name,c(.pred_LAT,.pred_sbm,.pred_ER,.pred_Spatial,.pred_SF,.pred_SW))

###----Save Trained----Workflow----###
library(here)
#save trained classifier
saveRDS(model_fit,"nnet_fitted_class_model.rds")

#To make predictions on new data set, we call the save model eg
load_fitted_model_nnet=readRDS("nnet_fitted_class_model.rds")


#########################################################################################
#+ ALE ANALYSIS FOR NNET
#########################################################################################

#ALE for top six feaures selected by boruta
fitted.model=load_fitted_model_nnet
DAT=boruta.final.data%>%
  select(-c(graph_name))
yhat <- function(X.model, newdata) as.numeric(predict(X.model, newdata))
## Calculate and plot the ALE main and second-order interaction effects of x1, x2, x3
par(mfrow = c(2,3))
ALE.1=DAT%>%
  ALEPlot(fitted.model, pred.fun=yhat, J=1, K=50, NA.plot = TRUE)



ALE.2=ALEPlot(DAT[,2:4], fitted.model, pred.fun=yhat, J=2, K=50, NA.plot = TRUE)
ALE.3=ALEPlot(DAT[,2:4], fitted.model, pred.fun=yhat, J=3, K=50, NA.plot = TRUE)
ALE.12=ALEPlot(DAT[,2:4], fitted.model, pred.fun=yhat, J=c(1,2), K=20, NA.plot = TRUE)
ALE.13=ALEPlot(DAT[,2:4], fitted.model, pred.fun=yhat, J=c(1,3), K=20, NA.plot = TRUE)
ALE.23=ALEPlot(DAT[,2:4], fitted.model, pred.fun=yhat, J=c(2,3), K=20, NA.plot = TRUE)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#######====PREDICTION ON NEW DATA====#####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##----Empirical Data----##

emp.data=read.csv("GraphFeaturesOnAnimalNetworks.csv")

emp.data=
  emp.data%>%
  rename("graph_name"="GraphNames")%>%
dplyr::select(c(graph_name,order,edges,
                mean_eccentr,mean_path_length,graph_energy,
                modularity,diameter,betw_centr,transitivity,
                spectral_radius,eigen_centr,deg_centr,
                mean_degree,minCut,FiedlerValue,Normalized_FiedlerValue,
                closeness_centr))%>%
  clean_names()%>%
  as_tibble()%>%
  mutate_if(is.character,factor)%>%
  na.omit()
  #replace(is.na(.), 0)#change NA's in deg_assort_coef to 0



#new_preds=new_preds%>%
# # slice_head(n=10)

#View(new_preds)

##----Data----frame----for newly predicted data----##

new_preds=load_fitted_model_nnet%>%
  augment(new_data=emp.data)

nnet.df.frame=data.frame(new_preds$graph_name,new_preds$.pred_class)
colnames(nnet.df.frame)=c("networks","predicted_classes")
nnet.df.frame

###----Stack bar plot----###
nnet.plot <- ggplot(nnet.df.frame, aes(x = predicted_classes,fill = networks))+
  geom_bar()+
  theme_classic()+
theme(text = element_text(size = 18),
  # legend.key.height= unit(2, 'cm'),
  #     legend.key.width= unit(4, 'cm'),
      legend.position="none")

nnet.plot


#ggsave("stack.png", width =32, height = 24)

# geom_bar(position = "fill") + ylab("proportion") +
#   stat_count(geom = "text", 
#              aes(label = stat(count)),
#              position=position_fill(vjust=0.5), colour="white")
# p


# ###----Statistical summary of confusion matrix----###
# Accuracy=accuracy(data=pred_result,truth = graph_name,estimate =.pred_class)
# Precision=precision(data=pred_result,truth = graph_name,estimate =.pred_class)
# Recall=recall(data=pred_result,truth = graph_name,estimate =.pred_class)
# 
# Accuracy%>%
#   bind_rows(Precision)%>%
#   bind_rows(Recall)
# 
# ##--ROC CURVE
# pred_result%>%
#   roc_curve(graph_name,c(.pred_LAT,.pred_sbm,.pred_ER,.pred_Spatial))%>%
#   ggplot(aes(x=1-specificity,y=sensitivity, color=.level))+
#   geom_abline(lty=2,color="gray80",size=.9)+
#   geom_path(show.legend = T,alpha=0.6,size=1.2)+
#   coord_equal()
# 
# ##ROC
# pred_result%>%
#   roc_auc(graph_name,c(.pred_LAT,.pred_sbm,.pred_ER,.pred_Spatial))


####++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++############
#----Random----Forest
####++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++############
cores <- parallel::detectCores()
cores

#Setting up model with tunning parameters
rf.model.spec <- 
  rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>% 
  set_engine("ranger", num.threads = cores) %>% 
  set_mode("classification")

#Creating recipe object
rf.rec <- 
  recipe(graph_name ~ ., data = train) 

rf.prep <- prep(rf.rec)
rf.juiced <- juice(rf.prep)

#Creating workflow object
rf.workflow <- 
  workflow() %>% 
  add_model(rf.model.spec) %>% 
  add_recipe(rf.rec)

# show what will be tuned
rf.model.spec
extract_parameter_set_dials(rf.model.spec)

#replacing any Na's in train and test set of deg_assort_coef column with 0
 train=train%>%
   na.omit()
 #  replace(is.na(.), 0)
 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Train RF hyperparameters
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 test=test%>%
   na.omit()
  # replace(is.na(.), 0)

# Create validation set
validatn.set <- validation_split(train, 
                            strata = graph_name, 
                            prop = 0.80)

rf.fold=vfold_cv(train)
rf.resample=bootstraps(train)


#library(doParallel)
#doParallel::registerDoParallel()

#all_cores <- parallel::detectCores(logical = FALSE)
 
# cl <- makePSOCKcluster(all_cores)
# registerDoParallel(cl)
#Result

set.seed(345)
rf.result <- 
  rf.workflow%>%
  tune_grid(rf.fold,
            grid = 25,
            control = control_grid(save_pred = TRUE),
            metrics = metric_set(roc_auc)
    )


library(here)
#save trained RF classifier
saveRDS(rf.result,"rf.result.rds")
#Show best metric
rf.result %>% 
  show_best(metric = "roc_auc")

#Shows best estimated tuned parameters
autoplot(rf.result)

#Another way to show tuned parameters

rf.result %>%
   collect_metrics() %>%
   filter(.metric == "roc_auc") %>%
   select(mean, min_n, mtry) %>%
   pivot_longer(min_n:mtry,
                values_to = "value",
                names_to = "parameter"
   ) %>%
   ggplot(aes(value, mean, color = parameter)) +
   geom_point(show.legend = FALSE) +
   facet_wrap(~parameter, scales = "free_x") +
   labs(x = NULL, y = "AUC")

##----Tune again with range of best tuned parameter values 
# rf.grid <- grid_regular(
#   mtry(range = c(2, 8)),
#   min_n(range = c(2, 13)),
#   levels = 5
# )
# 
# rf.grid



#We tune one more time, but this time in a more targeted way with this rf_grid.

# set.seed(456)
# regular.result <- tune_grid(
#   rf.workflow,
#   resamples = rf.fold,
#   grid = rf.grid
# )
# 
# regular.result


#Viewing best results 
rf.result %>%
  collect_metrics() %>%
  filter(.metric == "roc_auc") %>%
  mutate(min_n = factor(min_n)) %>%
  ggplot(aes(mtry, mean, color = min_n)) +
  geom_line(alpha = 0.5, size = 1.5) +
  geom_point() +
  labs(y = "AUC")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Best Random Forest Model
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


####----Select best tuned parameters for final model-----####
best_auc <- select_best(rf.result, "roc_auc")
best_auc
rf.final.model <- finalize_model(
  rf.model.spec,
  best_auc
)

rf.final.model

# rf.best.model <- 
#   rf.result %>% 
#   select_best(metric = "roc_auc")
# rf.best.params

####----Variable importance of final model----####
set.seed(6224)
library(vip)

rf.final.model%>%
  set_engine("ranger", importance = "permutation") %>%
  fit(graph_name ~ .,
      data = juice(rf.prep)
  ) %>%
  vip(geom = "point")


#The last workflow
rf.final.wkflw <- workflow() %>%
  add_recipe(rf.rec) %>%
  add_model(rf.final.model)

#The last fitted model
rf.final.fit <- 
  rf.final.wkflw %>%
  last_fit(df_split)

###----Save----Trained----RF----last----model----###
#save trained RF classifier
saveRDS(rf.final.fit ,"RfFinalMModel.rds")

#collect results
rf.final.fit %>%
  collect_metrics()

#Collect prediction
rf.final.fit %>% 
  collect_predictions()

rf.final.fit %>% 
  collect_predictions() %>% 
  roc_curve(graph_name, c(.pred_LAT,.pred_sbm,.pred_ER,.pred_Spatial)) %>% 
  autoplot()


#Filter results for our best RF model with the best parameters
NNET.roc=pred_result%>%
  roc_curve(graph_name,c(.pred_LAT,.pred_sbm,.pred_ER,.pred_Spatial))%>%
  mutate(model = "Neural-net")

rf.auc <- 
  rf.result %>% 
  collect_predictions(parameters = rf.final.model) %>% 
  roc_curve(graph_name, c(.pred_LAT,.pred_sbm,.pred_ER,.pred_Spatial)) %>% 
  mutate(model = "Random Forest")

bind_rows(rf.auc, NNET.roc) %>% 
  ggplot(aes(x = 1 - specificity, y = sensitivity, col = model)) + 
  geom_path(linewidth = 1.5, alpha = 0.8) +
  geom_abline(lty = 3) + 
  coord_equal() + 
  scale_color_viridis_d(option = "plasma", end = .6)

##----Lats RF model with best fitted parameters
# We’ll start by building our parsnip model object again from scratch.
# We take our best hyperparameter values from our random forest model. 
# When we set the engine, we add a new argument: importance = "impurity".
# This will provide variable importance scores for this last model, which 
# gives some insight into which predictors drive model performance.

#----The last model
# last.rf.model <- 
#   rand_forest(mtry = 7, min_n = 7, trees = 1000) %>% 
#   set_engine("ranger", num.threads = cores, importance = "impurity") %>% 
#   set_mode("classification")

#----The last workflow
# last.rf.workflow <- 
#   rf.workflow %>% 
#   update_model(last.rf.model)

#----The last fitted model
# set.seed(345)
# last.rf.fit <- 
#   load_fitted_model_rf %>% 
#   last_fit(df_split)
# 
# last.rf.fit



#To make predictions on new data set, we call the saved RF model eg
load_fitted_rf_model=readRDS("RfFinalMModel.rds")
#extract workflow fro fitted model
final.fr.model=extract_workflow(load_fitted_rf_model)

#----Show performance of RF model
# rf.final.fit%>% 
#   collect_metrics()

#----variable importance
# rf.final.fit %>% 
#   extract_fit_parsnip() %>% 
#   vip(num_features = 22)

#----Final Roc curve for the predicted class
# rf.final.fit %>% 
#   collect_predictions() %>% 
#   roc_curve(graph_name, c(.pred_LAT,.pred_sbm,.pred_ER,.pred_Spatial)) %>% 
#   autoplot()

# rf_new_preds=final.fr.model%>%
#   augment.last_fit(emp.data)

###----Prediction on empirical data
#The last fitted model
#set.seed(345)
df.emp=emp.data[sample(1:nrow(emp.data)), ]%>%
  na.omit()
  #replace(is.na(.), 0)#%>%
 # initial_split

###----Evaluate Model----####
#Predictions on empirical data 

rf.emp.pred.result=predict(final.fr.model,df.emp)

#emp.pred.result=pred_result%>%
#  slice_head(n=10)
rf.df.frame.emp=data.frame(df.emp$graph_name,rf.emp.pred.result)
colnames(rf.df.frame.emp)=c("networks","predicted_classes")
rf.df.frame.emp

###----Stack bar plot----###
rf.bar.plot <- ggplot(rf.df.frame.emp, aes(x = predicted_classes,fill = networks))+
  geom_bar()+
  theme_classic()+
  theme(text = element_text(size = 8),
        #       legend.key.height= unit(2, 'cm'),
        #       legend.key.width= unit(4, 'cm'),
        legend.position="none")

rf.bar.plot

##Confusion Matrix
emp.pred.result%>%
  conf_mat(truth = graph_name,estimate =.pred_class)


##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##+XGBOOST MODEL
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
xgb.fold=vfold_cv(train)

xgb.model.spec <- boost_tree(
  trees = 1000,
  tree_depth = tune(), min_n = tune(),
  loss_reduction = tune(),                     ## first three: model complexity
  sample_size = tune(), mtry = tune(),         ## randomness
  learn_rate = tune()                          ## step size
) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

xgb.model.spec



xgb.grid <- grid_latin_hypercube(
  tree_depth(),
  min_n(),
  loss_reduction(),
  sample_size = sample_prop(),
  finalize(mtry(), xgb.fold),
  learn_rate(),
  size = 30
)


xgb.wkflow <- workflow() %>%
  add_formula(graph_name ~ .) %>%
  add_model(xgb.model.spec)

xgb.wkflow

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Train XGBOOST hyperparameters
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
doParallel::registerDoParallel()

set.seed(234)
xgb.res <- tune_grid(
  xgb.wkflow,
  resamples = xgb.fold,
  grid = xgb.grid,
  control = control_grid(save_pred = TRUE)
)

xgb.res


collect_metrics(xgb.res)
###----Explore results----###
xgb.res %>%
  collect_metrics() %>%
  filter(.metric == "roc_auc") %>%
  select(mean, mtry:sample_size) %>%
  pivot_longer(mtry:sample_size,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "AUC")



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Best XGBOOST  Model
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#What are the best performing sets of parameters?
show_best(xgb.res, "roc_auc")


best_xgb_auc <- select_best(xgb.res, "roc_auc")
best_xgb_auc

#let’s finalize our tuneable workflow with these parameter values.

final.xgb.wkfl <- finalize_workflow(
  xgb.wkflow,
  best_xgb_auc
)

final.xgb.wkfl

###----Variable importance for xgboost
library(vip)

final.xgb.wkfl %>%
  fit(data = train) %>%
  pull_workflow_fit() %>%
  vip(geom = "point")


xgb.final.fit <- last_fit(final.xgb.wkfl, df_split)
###----Save----Trained----XGB----last----model----###
#save trained RF classifier
saveRDS(xgb.final.fit  ,"xgb.final.fit.rds")

xgb.final.fit %>%
collect_metrics()



xgb.final.fit%>%
  collect_predictions() %>%
  roc_curve(graph_name, c(.pred_LAT,.pred_sbm,.pred_ER,.pred_Spatial)) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(size = 1.5, color = "midnightblue") +
  geom_abline(
    lty = 2, alpha = 0.5,
    color = "gray50",
    size = 1.2
  )


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#######====PREDICTION ON NEW DATA====#####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##----Empirical Data----##

emp.data=read.csv("GraphFeaturesOnAnimalNetworks.csv")

emp.data=
  emp.data%>%
  rename("graph_name"="GraphNames")%>%
  dplyr::select(c(graph_name,order,edges,mean_eccentr,
                  radius,mean_path_length,graph_energy,min_triangle,mean_triangle,
                  sd_triangle,modularity,diameter,betw_centr,transitivity,threshold,
                  spectral_radius,eigen_centr,deg_centr,
                  mean_degree,minCut,FiedlerValue,Normalized_FiedlerValue,
                  closeness_centr),max_triangle,num_triangle,deg_assort_coef)%>%
  clean_names()%>%
  as_tibble()%>%
  mutate_if(is.character,factor)%>%
  replace(is.na(.), 0)#change NA's in deg_assort_coef to 0

#set.seed(345)
df.emp=emp.data[sample(1:nrow(emp.data)), ]%>%
  na.omit()


#new_preds=load_fitted_model%>%
 # augment(new_data=emp.data)

#To make predictions on new data set, we call the saved RF model eg
load_fitted_xgb_model=readRDS("xgb.final.fit.rds")
#extract workflow fro fitted model
final.xgb.model=extract_workflow(load_fitted_xgb_model)

xgb.emp.pred.result=predict(final.xgb.model,df.emp)


xgb.df.frame.emp=data.frame(df.emp$graph_name,xgb.emp.pred.result)
colnames(xgb.df.frame.emp)=c("networks","predicted_classes")
xgb.df.frame.emp



###----Stack bar plot----###
p <- ggplot(xgb.df.frame.emp, aes(x = predicted_classes,fill = networks))+
  geom_bar()+
  theme_classic()+
   theme(text = element_text(size = 8),
  #       legend.key.height= unit(2, 'cm'),
  #       legend.key.width= unit(4, 'cm'),
        legend.position ="none")

p


#ggsave("stack.png", width =32, height = 24)
