library(Matrix)
library(stats)
#library(spatstat)
library(igraph)
library(tidymodels)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Data generation for the machine learning model
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###---Simulate--spatial--networks---##
# simulate.spatial<-function(N=50,radius=0.4,nsim=100){
#   spatial.graph=NULL
#   for (i in 1:nsim){
#     spatial.graph[[i]]=makeSpatialGraphs(node.size=N,Radius=radius)
#   }
#   return(spatial.graph)
# }

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
##----------- Graph Features----------#####
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
### 50 nodes 
SF1.50=sim.sf.net(power=powerVals,m=mVals,nreps = 100,n=50)
write.csv(SF1.50,"SF1.50.csv")

### 100 nodes 
SF1.100=sim.sf.net(power=powerVals,m=mVals,nreps = 100,n=100)
write.csv(SF1.100,"SF1.100.csv")

### 150 nodes 
SF1.150=sim.sf.net(power=powerVals,m=mVals,nreps = 100,n=150)
write.csv(SF1.150,"SF1.150.csv")

### 200 nodes 
SF1.200=sim.sf.net(power=powerVals,m=mVals,nreps = 100,n=200)
write.csv(SF1.200,"SF1.200.csv")

### 250 nodes 
SF1.250=sim.sf.net(power=powerVals,m=mVals,nreps = 100,n=250)
write.csv(SF1.250,"SF1.250.csv")

### 300 nodes 
SF1.300=sim.sf.net(power=powerVals,m=mVals,nreps = 100,n=300)
write.csv(SF1.300,"SF1.300.csv")

### 350 nodes 
SF1.350=sim.sf.net(power=powerVals,m=mVals,nreps = 100,n=350)
write.csv(SF1.350,"SF1.350.csv")

### 400 nodes 
SF1.400=sim.sf.net(power=powerVals,m=mVals,nreps = 100,n=400)
write.csv(SF1.400,"SF1.400.csv")

### 500 nodes 
SF1.500=sim.sf.net(power=powerVals,m=mVals,nreps = 100,n=500)
write.csv(SF1.500,"SF1.500.csv")

### 750 nodes 
SF1.750=sim.sf.net(power=powerVals,m=mVals,nreps = 100,n=750)
write.csv(SF1.750,"SF1.750.csv")

### 1000 nodes 
SF1.1000=sim.sf.net(power=powerVals,m=mVals,nreps = 100,n=1000)
write.csv(SF1.1000,"SF1.1000.csv")

### 1500 nodes 
SF1.1500=sim.sf.net(power=powerVals,m=mVals,nreps = 100,n=1500)
write.csv(SF1.1500,"SF1.1500.csv")

### 2000 nodes 
SF1.2000=sim.sf.net(power=powerVals,m=mVals,nreps = 100,n=2000)
write.csv(SF1.2000,"SF1.2000.csv")

######## Saving all SF graphs to csv
SF.50=read.csv("SF1.50.csv");SF.100=read.csv("SF1.100.csv");
SF.150=read.csv("SF1.150.csv")
SF.200=read.csv("SF1.100.csv");
SF.250=read.csv("SF1.250.csv");SF.300=read.csv("SF1.300.csv")
SF.350=read.csv("SF1.350.csv");SF.400=read.csv("SF1.400.csv")
SF.500=read.csv("SF1.500.csv");SF.750=read.csv("SF1.750.csv")
SF.1000=read.csv("SF1.1000.csv");SF.1500=read.csv("SF1.1500.csv")
SF.2000=read.csv("SF1.2000.csv")
df.sf=rbind(SF.50,SF.150,SF.200,SF.250,SF.300,SF.350,
            SF.400, SF.500,SF.750,SF.1000,SF.1500,
            SF.2000)

write.csv(df.sf,"df.sf.csv")

n <- 50
m <- 520

# Work out how theye divide into each other
rem <- m %% n
div <- m %/% n

set.seed(123)
if(rem != 0) {
  g <- sample_smallworld(1, n, div+1, p = 0.001)
  # Randomly delete the unwanted edges. Should be quite homegenous
  g <- delete_edges(g, sample(1:gsize(g), size = gsize(g) - m))
} else {
  g <- sample_smallworld(1, n, div, p = 0.001)
}



lat=make_lattice(dimvector=1000, nei=2, circular=T)
plot(lat)

g <- sample_smallworld(dim=1, size=172, nei=6, p=0.05)
g.new=rewire(make_lattice(dimvector=500, nei=3, circular=T), with=each_edge(p=0.9, loops=F))
plot(g.new)
