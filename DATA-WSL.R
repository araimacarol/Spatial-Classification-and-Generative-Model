setwd("C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi-New")
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


 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                    Erdos Renyi  data

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MKERGraphs<- function(n=25,p=0.8) {
  Graphs = list()#;G=list() # set up the list of Graphs first
  initgraph= list()
  i= 1
  print("Creating ER Graphs")
  Graphs[[i]] = erdos.renyi.game(n,p,type = c("gnp"),directed = F,loops = F)
  
  Graphs[[i]]$type = "ER"
  Graphs[[i]]$id = i
  Graphs[[i]]$name="ER"
  i <- i+1
  return(Graphs)
}

sim.erdos.net=function(p=c(0.3,0.4,0.6,0.05),n=100){
net=NULL;data=NULL
net=as.list(mapply(FUN =MKERGraphs,n,p))
data= future.apply::future_lapply(list(net),RunSimOnGraphFeatures,nreps=1,future.seed = 0xBEEF)

df=cbind(data.frame(p),data)  
return(df)
}
# sim.erdos.net=function(p=c(0.3,0.4,0.6,0.05),n=100){
#   net=NULL;data=NULL
#   k=1
#   for (i in 1:length(p)){
#     net[i]= MKERGraphs(n,p[i])
#     data=RunSimOnGraphFeatures(net,nreps = 1)
#     k=k+1
#   }
#   df=cbind(p,data)
#   return(df)
# }
plan(multisession, workers = 32)
set.seed(0xBEEF)
### 50 nodes 
ER1.50=sim.erdos.net(p=rep(0.1,250),n=50)
write.csv(ER1.50,"ER1.50.csv")
ER2.50=sim.erdos.net(p=rep(0.2,250),n=50)
write.csv(ER2.50,"ER2.50.csv")
ER3.50=sim.erdos.net(p=rep(0.3,250),n=50)
write.csv(ER3.50,"ER3.50.csv")
ER4.50=sim.erdos.net(p=rep(0.4,250),n=50)
write.csv(ER4.50,"ER4.50.csv")
ER5.50=sim.erdos.net(p=rep(0.5,250),n=50)
write.csv(ER5.50,"ER5.50.csv")
ER6.50=sim.erdos.net(p=rep(0.6,250),n=50)
write.csv(ER6.50,"ER6.50.csv")
ER7.50=sim.erdos.net(p=rep(0.7,250),n=50)
write.csv(ER7.50,"ER7.50.csv")
ER8.50=sim.erdos.net(p=rep(0.8,250),n=50)
write.csv(ER8.50,"ER8.50.csv")
ER9.50=sim.erdos.net(p=rep(0.9,250),n=50)
write.csv(ER9.50,"ER9.50.csv")
######## Saving all ER graphs 50 nodes data to csv
ER1=read.csv("ER1.50.csv");ER2=read.csv("ER2.50.csv");ER3=read.csv("ER3.50.csv")
ER4=read.csv("ER4.50.csv");ER5=read.csv("ER5.50.csv");ER6=read.csv("ER6.50.csv")
ER7=read.csv("ER7.50.csv");ER8=read.csv("ER8.50.csv");ER9=read.csv("ER9.50.csv")

df.er.50=rbind(ER1,ER2,ER3,ER4,ER5,ER6,ER7,ER8,ER9)
write.csv(df.er.50,"df.er.50.csv")

### 100 nodes 
ER1.100=sim.erdos.net(p=rep(0.1,250),n=100)
write.csv(ER1.100,"ER1.100.csv")
ER2.100=sim.erdos.net(p=rep(0.2,250),n=100)
write.csv(ER2.100,"ER2.100.csv")
ER3.100=sim.erdos.net(p=rep(0.3,250),n=100)
write.csv(ER3.100,"ER3.100.csv")
ER4.100=sim.erdos.net(p=rep(0.4,250),n=100)
write.csv(ER4.100,"ER4.100.csv")
ER5.100=sim.erdos.net(p=rep(0.5,250),n=100)
write.csv(ER5.100,"ER5.100.csv")
ER6.100=sim.erdos.net(p=rep(0.6,250),n=100)
write.csv(ER6.100,"ER6.100.csv")
ER7.100=sim.erdos.net(p=rep(0.7,250),n=100)
write.csv(ER7.100,"ER7.100.csv")
ER8.100=sim.erdos.net(p=rep(0.8,250),n=100)
write.csv(ER8.100,"ER8.100.csv")
ER9.100=sim.erdos.net(p=rep(0.9,250),n=100)
write.csv(ER9.100,"ER9.100.csv")
######## Saving all ER graphs 100 nodes data to csv
ER1=read.csv("ER1.100.csv");ER2=read.csv("ER2.100.csv");ER3=read.csv("ER3.100.csv")
ER4=read.csv("ER4.100.csv");ER5=read.csv("ER5.100.csv");ER6=read.csv("ER6.100.csv")
ER7=read.csv("ER7.100.csv");ER8=read.csv("ER8.100.csv");ER9=read.csv("ER9.100.csv")

df.er.100=rbind(ER1,ER2,ER3,ER4,ER5,ER6,ER7,ER8,ER9)
write.csv(df.er.100,"df.er.100.csv")

### 150 nodes 
ER1.150=sim.erdos.net(p=rep(0.1,250),n=150)
write.csv(ER1.150,"ER1.150.csv")
ER2.150=sim.erdos.net(p=rep(0.2,250),n=150)
write.csv(ER2.150,"ER2.150.csv")
ER3.150=sim.erdos.net(p=rep(0.3,250),n=150)
write.csv(ER3.150,"ER3.150.csv")
ER4.150=sim.erdos.net(p=rep(0.4,250),n=150)
write.csv(ER4.150,"ER4.150.csv")
ER5.150=sim.erdos.net(p=rep(0.5,250),n=150)
write.csv(ER5.150,"ER5.150.csv")
ER6.150=sim.erdos.net(p=rep(0.6,250),n=150)
write.csv(ER6.150,"ER6.150.csv")
ER7.150=sim.erdos.net(p=rep(0.7,250),n=150)
write.csv(ER7.150,"ER7.150.csv")
ER8.150=sim.erdos.net(p=rep(0.8,250),n=150)
write.csv(ER8.150,"ER8.150.csv")
ER9.150=sim.erdos.net(p=rep(0.9,250),n=150)
write.csv(ER9.150,"ER9.150.csv")
######## Saving all ER graphs 150 nodes data to csv
ER1=read.csv("ER1.150.csv");ER2=read.csv("ER2.150.csv");ER3=read.csv("ER3.150.csv")
ER4=read.csv("ER4.150.csv");ER5=read.csv("ER5.150.csv");ER6=read.csv("ER6.150.csv")
ER7=read.csv("ER7.150.csv");ER8=read.csv("ER8.150.csv");ER9=read.csv("ER9.150.csv")

df.er.150=rbind(ER1,ER2,ER3,ER4,ER5,ER6,ER7,ER8,ER9)
write.csv(df.er.150,"df.er.150.csv")

### 200 nodes 
ER1.200=sim.erdos.net(p=rep(0.1,250),n=200)
write.csv(ER1.200,"ER1.200.csv")
ER2.200=sim.erdos.net(p=rep(0.2,250),n=200)
write.csv(ER2.200,"ER2.200.csv")
ER3.200=sim.erdos.net(p=rep(0.3,250),n=200)
write.csv(ER3.200,"ER3.200.csv")
ER4.200=sim.erdos.net(p=rep(0.4,250),n=200)
write.csv(ER4.200,"ER4.200.csv")
ER5.200=sim.erdos.net(p=rep(0.5,250),n=200)
write.csv(ER5.200,"ER5.200.csv")
ER6.200=sim.erdos.net(p=rep(0.6,250),n=200)
write.csv(ER6.200,"ER6.200.csv")
ER7.200=sim.erdos.net(p=rep(0.7,250),n=200)
write.csv(ER7.200,"ER7.200.csv")
ER8.200=sim.erdos.net(p=rep(0.8,250),n=200)
write.csv(ER8.200,"ER8.200.csv")
ER9.200=sim.erdos.net(p=rep(0.9,250),n=200)
write.csv(ER9.200,"ER9.200.csv")
######## Saving all ER graphs 200 nodes data to csv
ER1=read.csv("ER1.200.csv");ER2=read.csv("ER2.200.csv");ER3=read.csv("ER3.200.csv")
ER4=read.csv("ER4.200.csv");ER5=read.csv("ER5.200.csv");ER6=read.csv("ER6.200.csv")
ER7=read.csv("ER7.200.csv");ER8=read.csv("ER8.200.csv");ER9=read.csv("ER9.200.csv")

df.er.200=rbind(ER1,ER2,ER3,ER4,ER5,ER6,ER7,ER8,ER9)
write.csv(df.er.200,"df.er.200.csv")

### 250 nodes 
ER1.250=sim.erdos.net(p=rep(0.1,250),n=250)
write.csv(ER1.250,"ER1.250.csv")
ER2.250=sim.erdos.net(p=rep(0.2,250),n=250)
write.csv(ER2.250,"ER2.250.csv")
ER3.250=sim.erdos.net(p=rep(0.3,250),n=250)
write.csv(ER3.250,"ER3.250.csv")
ER4.250=sim.erdos.net(p=rep(0.4,250),n=250)
write.csv(ER4.250,"ER4.250.csv")
ER5.250=sim.erdos.net(p=rep(0.5,250),n=250)
write.csv(ER5.250,"ER5.250.csv")
ER6.250=sim.erdos.net(p=rep(0.6,250),n=250)
write.csv(ER6.250,"ER6.250.csv")
ER7.250=sim.erdos.net(p=rep(0.7,250),n=250)
write.csv(ER7.250,"ER7.250.csv")
ER8.250=sim.erdos.net(p=rep(0.8,250),n=250)
write.csv(ER8.250,"ER8.250.csv")
ER9.250=sim.erdos.net(p=rep(0.9,250),n=250)
write.csv(ER9.250,"ER9.250.csv")
######## Saving all ER graphs 250 nodes data to csv
ER1=read.csv("ER1.250.csv");ER2=read.csv("ER2.250.csv");ER3=read.csv("ER3.250.csv")
ER4=read.csv("ER4.250.csv");ER5=read.csv("ER5.250.csv");ER6=read.csv("ER6.250.csv")
ER7=read.csv("ER7.250.csv");ER8=read.csv("ER8.250.csv");ER9=read.csv("ER9.250.csv")

df.er.250=rbind(ER1,ER2,ER3,ER4,ER5,ER6,ER7,ER8,ER9)
write.csv(df.er.250,"df.er.250.csv")

### 300 nodes 
ER1.300=sim.erdos.net(p=rep(0.1,250),n=300)
write.csv(ER1.300,"ER1.300.csv")
ER2.300=sim.erdos.net(p=rep(0.2,250),n=300)
write.csv(ER2.300,"ER2.300.csv")
ER3.300=sim.erdos.net(p=rep(0.3,250),n=300)
write.csv(ER3.300,"ER3.300.csv")
ER4.300=sim.erdos.net(p=rep(0.4,250),n=300)
write.csv(ER4.300,"ER4.300.csv")
ER5.300=sim.erdos.net(p=rep(0.5,250),n=300)
write.csv(ER5.300,"ER5.300.csv")
ER6.300=sim.erdos.net(p=rep(0.6,250),n=300)
write.csv(ER6.300,"ER6.300.csv")
ER7.300=sim.erdos.net(p=rep(0.7,250),n=300)
write.csv(ER7.300,"ER7.300.csv")
ER8.300=sim.erdos.net(p=rep(0.8,250),n=300)
write.csv(ER8.300,"ER8.300.csv")
ER9.300=sim.erdos.net(p=rep(0.9,250),n=300)
write.csv(ER9.300,"ER9.300.csv")
######## Saving all ER graphs 300 nodes data to csv
ER1=read.csv("ER1.300.csv");ER2=read.csv("ER2.300.csv");ER3=read.csv("ER3.300.csv")
ER4=read.csv("ER4.300.csv");ER5=read.csv("ER5.300.csv");ER6=read.csv("ER6.300.csv")
ER7=read.csv("ER7.300.csv");ER8=read.csv("ER8.300.csv");ER9=read.csv("ER9.300.csv")

df.er.300=rbind(ER1,ER2,ER3,ER4,ER5,ER6,ER7,ER8,ER9)
write.csv(df.er.300,"df.er.300.csv")

### 350 nodes 
ER1.350=sim.erdos.net(p=rep(0.1,250),n=350)
write.csv(ER1.350,"ER1.350.csv")
ER2.350=sim.erdos.net(p=rep(0.2,250),n=350)
write.csv(ER2.350,"ER2.350.csv")
ER3.350=sim.erdos.net(p=rep(0.3,250),n=350)
write.csv(ER3.350,"ER3.350.csv")
ER4.350=sim.erdos.net(p=rep(0.4,250),n=350)
write.csv(ER4.350,"ER4.350.csv")
ER5.350=sim.erdos.net(p=rep(0.5,250),n=350)
write.csv(ER5.350,"ER5.350.csv")
ER6.350=sim.erdos.net(p=rep(0.6,250),n=350)
write.csv(ER6.350,"ER6.350.csv")
ER7.350=sim.erdos.net(p=rep(0.7,250),n=350)
write.csv(ER7.350,"ER7.350.csv")
ER8.350=sim.erdos.net(p=rep(0.8,250),n=350)
write.csv(ER8.350,"ER8.350.csv")
ER9.350=sim.erdos.net(p=rep(0.9,250),n=350)
write.csv(ER9.350,"ER9.350.csv")
######## Saving all ER graphs 350 nodes data to csv
ER1=read.csv("ER1.350.csv");ER2=read.csv("ER2.350.csv");ER3=read.csv("ER3.350.csv")
ER4=read.csv("ER4.350.csv");ER5=read.csv("ER5.350.csv");ER6=read.csv("ER6.350.csv")
ER7=read.csv("ER7.350.csv");ER8=read.csv("ER8.350.csv");ER9=read.csv("ER9.350.csv")

df.er.350=rbind(ER1,ER2,ER3,ER4,ER5,ER6,ER7,ER8,ER9)
write.csv(df.er.350,"df.er.350.csv")

### 400 nodes 
ER1.400=sim.erdos.net(p=rep(0.1,250),n=400)
write.csv(ER1.400,"ER1.400.csv")
ER2.400=sim.erdos.net(p=rep(0.2,250),n=400)
write.csv(ER2.400,"ER2.400.csv")
ER3.400=sim.erdos.net(p=rep(0.3,250),n=400)
write.csv(ER3.400,"ER3.400.csv")
ER4.400=sim.erdos.net(p=rep(0.4,250),n=400)
write.csv(ER4.400,"ER4.400.csv")
ER5.400=sim.erdos.net(p=rep(0.5,250),n=400)
write.csv(ER5.400,"ER5.400.csv")
ER6.400=sim.erdos.net(p=rep(0.6,250),n=400)
write.csv(ER6.400,"ER6.400.csv")
ER7.400=sim.erdos.net(p=rep(0.7,250),n=400)
write.csv(ER7.400,"ER7.400.csv")
ER8.400=sim.erdos.net(p=rep(0.8,250),n=400)
write.csv(ER8.400,"ER8.400.csv")
ER9.400=sim.erdos.net(p=rep(0.9,250),n=400)
write.csv(ER9.400,"ER9.400.csv")
######## Saving all ER graphs 400 nodes data to csv
ER1=read.csv("ER1.400.csv");ER2=read.csv("ER2.400.csv");ER3=read.csv("ER3.400.csv")
ER4=read.csv("ER4.400.csv");ER5=read.csv("ER5.400.csv");ER6=read.csv("ER6.400.csv")
ER7=read.csv("ER7.400.csv");ER8=read.csv("ER8.400.csv");ER9=read.csv("ER9.400.csv")

df.er.400=rbind(ER1,ER2,ER3,ER4,ER5,ER6,ER7,ER8,ER9)
write.csv(df.er.400,"df.er.400.csv")

### 500 nodes 
ER1.500=sim.erdos.net(p=rep(0.1,250),n=500)
write.csv(ER1.500,"ER1.500.csv")
ER2.500=sim.erdos.net(p=rep(0.2,250),n=500)
write.csv(ER2.500,"ER2.500.csv")
ER3.500=sim.erdos.net(p=rep(0.3,250),n=500)
write.csv(ER3.500,"ER3.500.csv")
ER4.500=sim.erdos.net(p=rep(0.4,250),n=500)
write.csv(ER4.500,"ER4.500.csv")
ER5.500=sim.erdos.net(p=rep(0.5,250),n=500)
write.csv(ER5.500,"ER5.500.csv")
ER6.500=sim.erdos.net(p=rep(0.6,250),n=500)
write.csv(ER6.500,"ER6.500.csv")
ER7.500=sim.erdos.net(p=rep(0.7,250),n=500)
write.csv(ER7.500,"ER7.500.csv")
ER8.500=sim.erdos.net(p=rep(0.8,250),n=500)
write.csv(ER8.500,"ER8.500.csv")
ER9.500=sim.erdos.net(p=rep(0.9,250),n=500)
write.csv(ER9.500,"ER9.500.csv")
######## Saving all ER graphs 500 nodes data to csv
ER1=read.csv("ER1.500.csv");ER2=read.csv("ER2.500.csv");ER3=read.csv("ER3.500.csv")
ER4=read.csv("ER4.500.csv");ER5=read.csv("ER5.500.csv");ER6=read.csv("ER6.500.csv")
ER7=read.csv("ER7.500.csv");ER8=read.csv("ER8.500.csv");ER9=read.csv("ER9.500.csv")

df.er.500=rbind(ER1,ER2,ER3,ER4,ER5,ER6,ER7,ER8,ER9)
write.csv(df.er.500,"df.er.500.csv")

### 750 nodes 
ER1.750=sim.erdos.net(p=rep(0.1,250),n=750)
write.csv(ER1.750,"ER1.750.csv")
ER2.750=sim.erdos.net(p=rep(0.2,250),n=750)
write.csv(ER2.750,"ER2.750.csv")
ER3.750=sim.erdos.net(p=rep(0.3,250),n=750)
write.csv(ER3.750,"ER3.750.csv")
ER4.750=sim.erdos.net(p=rep(0.4,250),n=750)
write.csv(ER4.750,"ER4.750.csv")
ER5.750=sim.erdos.net(p=rep(0.5,250),n=750)
write.csv(ER5.750,"ER5.750.csv")
ER6.750=sim.erdos.net(p=rep(0.6,250),n=750)
write.csv(ER6.750,"ER6.750.csv")
ER7.750=sim.erdos.net(p=rep(0.7,250),n=750)
write.csv(ER7.750,"ER7.750.csv")
ER8.750=sim.erdos.net(p=rep(0.8,250),n=750)
write.csv(ER8.750,"ER8.750.csv")
ER9.750=sim.erdos.net(p=rep(0.9,250),n=750)
write.csv(ER9.750,"ER9.750.csv")
######## Saving all ER graphs 750 nodes data to csv
ER1=read.csv("ER1.750.csv");ER2=read.csv("ER2.750.csv");ER3=read.csv("ER3.750.csv")
ER4=read.csv("ER4.750.csv");ER5=read.csv("ER5.750.csv");ER6=read.csv("ER6.750.csv")
ER7=read.csv("ER7.750.csv");ER8=read.csv("ER8.750.csv");ER9=read.csv("ER9.750.csv")

df.er.750=rbind(ER1,ER2,ER3,ER4,ER5,ER6,ER7,ER8,ER9)
write.csv(df.er.750,"df.er.750.csv")

### 1000 nodes 
ER1.1000=sim.erdos.net(p=rep(0.1,250),n=1000)
write.csv(ER1.1000,"ER1.1000.csv")
ER2.1000=sim.erdos.net(p=rep(0.2,250),n=1000)
write.csv(ER2.1000,"ER2.1000.csv")
ER3.1000=sim.erdos.net(p=rep(0.3,250),n=1000)
write.csv(ER3.1000,"ER3.1000.csv")
ER4.1000=sim.erdos.net(p=rep(0.4,250),n=1000)
write.csv(ER4.1000,"ER4.1000.csv")
ER5.1000=sim.erdos.net(p=rep(0.5,250),n=1000)
write.csv(ER5.1000,"ER5.1000.csv")
ER6.1000=sim.erdos.net(p=rep(0.6,250),n=1000)
write.csv(ER6.1000,"ER6.1000.csv")
ER7.1000=sim.erdos.net(p=rep(0.7,250),n=1000)
write.csv(ER7.1000,"ER7.1000.csv")
ER8.1000=sim.erdos.net(p=rep(0.8,250),n=1000)
write.csv(ER8.1000,"ER8.1000.csv")
ER9.1000=sim.erdos.net(p=rep(0.9,250),n=1000)
write.csv(ER9.1000,"ER9.1000.csv")
######## Saving all ER graphs 1000 nodes data to csv
ER1=read.csv("ER1.1000.csv");ER2=read.csv("ER2.1000.csv");ER3=read.csv("ER3.1000.csv")
ER4=read.csv("ER4.1000.csv");ER5=read.csv("ER5.1000.csv");ER6=read.csv("ER6.1000.csv")
ER7=read.csv("ER7.1000.csv");ER8=read.csv("ER8.1000.csv");ER9=read.csv("ER9.1000.csv")

df.er.1000=rbind(ER1,ER2,ER3,ER4,ER5,ER6,ER7,ER8,ER9)
write.csv(df.er.1000,"df.er.1000.csv")

### 1500 nodes 
ER1.1500=sim.erdos.net(p=rep(0.1,250),n=1500)
write.csv(ER1.1500,"ER1.1500.csv")
ER2.1500=sim.erdos.net(p=rep(0.2,250),n=1500)
write.csv(ER2.1500,"ER2.1500.csv")
ER3.1500=sim.erdos.net(p=rep(0.3,250),n=1500)
write.csv(ER3.1500,"ER3.1500.csv")
ER4.1500=sim.erdos.net(p=rep(0.4,250),n=1500)
write.csv(ER4.1500,"ER4.1500.csv")
ER5.1500=sim.erdos.net(p=rep(0.5,250),n=1500)
write.csv(ER5.1500,"ER5.1500.csv")
ER6.1500=sim.erdos.net(p=rep(0.6,250),n=1500)
write.csv(ER6.1500,"ER6.1500.csv")
ER7.1500=sim.erdos.net(p=rep(0.7,250),n=1500)
write.csv(ER7.1500,"ER7.1500.csv")
ER8.1500=sim.erdos.net(p=rep(0.8,250),n=1500)
write.csv(ER8.1500,"ER8.1500.csv")
ER9.1500=sim.erdos.net(p=rep(0.9,250),n=1500)
write.csv(ER9.1500,"ER9.1500.csv")
######## Saving all ER graphs 1500 nodes data to csv
ER1=read.csv("ER1.1500.csv");ER2=read.csv("ER2.1500.csv");ER3=read.csv("ER3.1500.csv")
ER4=read.csv("ER4.1500.csv");ER5=read.csv("ER5.1500.csv");ER6=read.csv("ER6.1500.csv")
ER7=read.csv("ER7.1500.csv");ER8=read.csv("ER8.1500.csv");ER9=read.csv("ER9.1500.csv")

df.er.1500=rbind(ER1,ER2,ER3,ER4,ER5,ER6,ER7,ER8,ER9)
write.csv(df.er.1500,"df.er.1500.csv")

### 2000 nodes 
ER1.2000=sim.erdos.net(p=rep(0.1,250),n=2000)
write.csv(ER1.2000,"ER1.2000.csv")
ER2.2000=sim.erdos.net(p=rep(0.2,250),n=2000)
write.csv(ER2.2000,"ER2.2000.csv")
ER3.2000=sim.erdos.net(p=rep(0.3,250),n=2000)
write.csv(ER3.2000,"ER3.2000.csv")
ER4.2000=sim.erdos.net(p=rep(0.4,250),n=2000)
write.csv(ER4.2000,"ER4.2000.csv")
ER5.2000=sim.erdos.net(p=rep(0.5,250),n=2000)
write.csv(ER5.2000,"ER5.2000.csv")
ER6.2000=sim.erdos.net(p=rep(0.6,250),n=2000)
write.csv(ER6.2000,"ER6.2000.csv")
ER7.2000=sim.erdos.net(p=rep(0.7,250),n=2000)
write.csv(ER7.2000,"ER7.2000.csv")
ER8.2000=sim.erdos.net(p=rep(0.8,250),n=2000)
write.csv(ER8.2000,"ER8.2000.csv")
ER9.2000=sim.erdos.net(p=rep(0.9,250),n=2000)
write.csv(ER9.2000,"ER9.2000.csv")
######## Saving all ER graphs 2000 nodes data to csv
ER1=read.csv("ER1.2000.csv");ER2=read.csv("ER2.2000.csv");ER3=read.csv("ER3.2000.csv")
ER4=read.csv("ER4.2000.csv");ER5=read.csv("ER5.2000.csv");ER6=read.csv("ER6.2000.csv")
ER7=read.csv("ER7.2000.csv");ER8=read.csv("ER8.2000.csv");ER9=read.csv("ER9.2000.csv")

df.er.2000=rbind(ER1,ER2,ER3,ER4,ER5,ER6,ER7,ER8,ER9)
write.csv(df.er.2000,"df.er.2000.csv")


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                    Scale free  data

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MKSFGraphs<- function(n=50, power=2, m=4) {
  Graphs = list()#;G=list() # set up the list of Graphs first
  initgraph= list()
  i= 1
  print("Creating SF Graphs")
  Graphs[[i]] = sample_pa(n=n, power=power, m=m, directed=FALSE, algorithm="psumtree")
  
  Graphs[[i]]$type = "SF"
  Graphs[[i]]$id = i
  Graphs[[i]]$name="SF"
  i <- i+1
  return(Graphs)
}


# ParalleEpicSimOnGraphs<-function(Graphs, betaVals=betaVals,gammaVals=gammaVals, nreps=nreps,output_file="EpicSimOnGraphs.csv",report="s",nticks=10) {
#   #  p = progressor(along=Graphs)
#   pathogn=expand.grid(data.frame(betaVals,gammaVals))
#   
#   for (m in 1:nrow(pathogn)){## Counter to account for the cartesian product of all beta and gamma values
#     results[[m]] = future_lapply(seq_along(Graphs), function(j) { 
#       # Sys.sleep(6.0-j)
#       #  p(sprintf("j=%g", j))
#       data.frame(Run_EpicSim_And_Measures(Graphs[j], beta=pathogn$betaVals[m],gamma=pathogn$gammaVals[m],
#                                           nticks=nticks,report=report,nreps = nreps))} 
#       , future.seed = 0xBEEF, future.chunk.size=1)
#     
#   }}


sim.sf.net=function(power=powerVals,m=mVals,n=50,nreps=1){
  net=NULL;data=NULL
  dt=expand.grid(n,powerVals,mVals)
  colnames(dt)=c("nodes","prewireVals","neiVals")
  #dt.split=setNames(split(dt, seq(nrow(dt))), rownames(dtname))
  net=as.list(mapply(FUN =MKSFGraphs,n=dt$nodes,power=dt$prewireVals,m=dt$neiVals))
  #k=1
 # start <- Sys.time()
  #net[i]= MKSFGraphs(n,power=dt$powerVals[i],m=dt$mVals[i])
  data= future.apply::future_lapply(list(net),RunSimOnGraphFeatures,nreps,future.seed = 0xBEEF)
  #k=k+1
    
  df=cbind(dt,data)  
  
  return(df)
}

plan(multisession, workers = 6)
set.seed(0xBEEF)
#powerVals=1:14
#mVals=5:15
#x=sim.sf.net(power=powerVals,m=mVals,nreps = 50,n=50)
powerVals=1:5
mVals=1:50
### 50 nodes 
SF1.50=sim.sf.net(power=powerVals,m=mVals,nreps = 1200,n=50)
write.csv(SF1.50,"SF1.50.csv")

### 100 nodes 
SF1.100=sim.sf.net(power=powerVals,m=mVals,nreps = 1200,n=100)
write.csv(SF1.100,"SF1.100.csv")

### 150 nodes 
SF1.150=sim.sf.net(power=powerVals,m=mVals,nreps = 1200,n=150)
write.csv(SF1.150,"SF1.150.csv")

### 200 nodes 
SF1.200=sim.sf.net(power=powerVals,m=mVals,nreps = 1200,n=200)
write.csv(SF1.200,"SF1.200.csv")

### 250 nodes 
SF1.250=sim.sf.net(power=powerVals,m=mVals,nreps = 1200,n=250)
write.csv(SF1.250,"SF1.250.csv")

### 300 nodes 
SF1.300=sim.sf.net(power=powerVals,m=mVals,nreps = 1200,n=300)
write.csv(SF1.300,"SF1.300.csv")

### 350 nodes 
SF1.350=sim.sf.net(power=powerVals,m=mVals,nreps = 1200,n=350)
write.csv(SF1.350,"SF1.350.csv")

### 400 nodes 
SF1.400=sim.sf.net(power=powerVals,m=mVals,nreps = 1200,n=400)
write.csv(SF1.400,"SF1.400.csv")

### 500 nodes 
SF1.500=sim.sf.net(power=powerVals,m=mVals,nreps = 1200,n=500)
write.csv(SF1.500,"SF1.500.csv")

### 750 nodes 
SF1.750=sim.sf.net(power=powerVals,m=mVals,nreps = 1200,n=750)
write.csv(SF1.750,"SF1.750.csv")

### 1000 nodes 
SF1.1000=sim.sf.net(power=powerVals,m=mVals,nreps = 1200,n=1000)
write.csv(SF1.1000,"SF1.1000.csv")

### 1500 nodes 
SF1.1500=sim.sf.net(power=powerVals,m=mVals,nreps = 1200,n=1500)
write.csv(SF1.1500,"SF1.1500.csv")

### 2000 nodes 
SF1.2000=sim.sf.net(power=powerVals,m=mVals,nreps = 1200,n=2000)
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


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                    Small world  data

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
  dt=expand.grid(n,pVals,neiVals)
  colnames(dt)=c("nodes","prewireVals","neiVals")
  net=as.list(mapply(FUN =MKSWGraphs,n=dt$nodes,p=dt$prewireVals,nei=dt$neiVals))
  data= future.apply::future_lapply(list(net),RunSimOnGraphFeatures,nreps,future.seed = 0xBEEF)
 
  df=cbind(dt,data)  
  return(df)
}

plan(multisession, workers = 32)
set.seed(0xBEEF)

pVals=c(0.01,0.03,0.05,0.07,0.09,0.1,0.2,0.3,0.4)
neiVals=1:20

### 50 nodes 
SW.50=sim.sw.net(p=pVals,nei=neiVals,n=50,nreps=1200)
write.csv(SW.50,"SW.50.csv")

### 100 nodes 
SW.100=sim.sw.net(p=pVals,nei=neiVals,n=100,nreps=1200)
write.csv(SW.100,"SW.100.csv")

### 150 nodes 
SW.150=sim.sw.net(p=pVals,nei=neiVals,n=150,nreps=1200)
write.csv(SW.150,"SW.150.csv")

### 200 nodes 
SW.200=sim.sw.net(p=pVals,nei=neiVals,n=200,nreps=1200)
write.csv(SW.200,"SW.200.csv")

### 250 nodes 
SW.250=sim.sw.net(p=pVals,nei=neiVals,n=250,nreps=1200)
write.csv(SW.250,"SW.250.csv")

### 300 nodes 
SW.300=sim.sw.net(p=pVals,nei=neiVals,n=300,nreps=1200)
write.csv(SW.300,"SW.300.csv")

### 350 nodes 
SW.350=sim.sw.net(p=pVals,nei=neiVals,n=350,nreps=1200)
write.csv(SW.350,"SW.350.csv")

### 400 nodes 
SW.400=sim.sw.net(p=pVals,nei=neiVals,n=400,nreps=1200)
write.csv(SW.400,"SW.400.csv")

### 500 nodes 
SW.500=sim.sw.net(p=pVals,nei=neiVals,n=500,nreps=1200)
write.csv(SW.500,"SW.500.csv")

### 750 nodes 
SW.750=sim.sw.net(p=pVals,nei=neiVals,n=750,nreps=1200)
write.csv(SW.750,"SW.750.csv")

### 1000 nodes 
SW.1000=sim.sw.net(p=pVals,nei=neiVals,n=1000,nreps=1200)
write.csv(SW.1000,"SW.1000.csv")

### 1500 nodes 
SW.1500=sim.sw.net(p=pVals,nei=neiVals,n=1500,nreps=1200)
write.csv(SW.1500,"SW.1500.csv")

### 2000 nodes 
SW.2000=sim.sw.net(p=pVals,nei=neiVals,n=2000,nreps=1200)
write.csv(SW.2000,"SW.2000.csv")

######## Saving all SW graphs to csv
SW.50=read.csv("SW.50.csv");SW.100=read.csv("SW.100.csv");
SW.150=read.csv("SW.150.csv")
SW.200=read.csv("SW.100.csv");
SW.250=read.csv("SW.250.csv");SW.300=read.csv("SW.300.csv")
SW.350=read.csv("SW.350.csv");SW.400=read.csv("SW.400.csv")
SW.500=read.csv("SW.500.csv");SW.750=read.csv("SW.750.csv")
SW.1000=read.csv("SW.1000.csv");SW.1500=read.csv("SW.1500.csv")
SW.2000=read.csv("SW.2000.csv")

df.sw=rbind(SW.50,SW.150,SW.200,SW.250,SW.300,SW.350,
            SW.400, SW.500,SW.750,SW.1000,SW.1500,
            SW.2000)

write.csv(df.sw,"df.sw.csv")


# library(Matrix)
# library(stats)
# #library(spatstat)
# library(igraph)
# library(tidymodels)
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# #+ Data generation for the machine learning model
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ###---Simulate--spatial--networks---##
# # simulate.spatial<-function(N=50,radius=0.4,nsim=100){
# #   spatial.graph=NULL
# #   for (i in 1:nsim){
# #     spatial.graph[[i]]=makeSpatialGraphs(node.size=N,Radius=radius)
# #   }
# #   return(spatial.graph)
# # }
# 
# ##########################################################################################################
# # Input an igraph object from a file, treating it as a static graph
# ##########################################################################################################
# getGraphFromFile <- function(file, simplify=TRUE, useBiggestComponent=TRUE, asUndirected=TRUE) {
#   
#   dat <- read.table(file) # just read static graph: ignoring third column
#   G <- graph_from_data_frame(dat)
#   if (asUndirected==TRUE) {
#     G <- as.undirected(G, "collapse")
#   }
#   
#   #  g_names <- gsub(".edges","",networks[i]) # one edge for each pair of connect vertices (not sure what this is for)
#   
#   if (useBiggestComponent==TRUE) {
#     netclust <- clusters(G) #look for subgraphs
#     gcc <- V(G)[netclust$membership == which.max(netclust$csize)]#select vertices from the largest sub-graph
#     G <- induced.subgraph(G, gcc) #make it a igraph object.
#   }
#   if (simplify==TRUE) {
#     G <- igraph::simplify(G, remove.multiple = TRUE, remove.loops = TRUE)
#   }
#   return(G)
# }
# 
# ###--Normalized Laplacina function--##
# normalized_laplacian=function(Graphs){
#   laplacian_matrix(Graphs,normalized = T)
# }
# 
# 
# 
# calcGraphFeatures <- function(Graphs=NULL) {
#   
#   features <- c(
#     "order",                    # number of vertices
#     "edges",                     # number of edges
#     "connected",                # True / False
#     "max_component",            # maximum component size (=order iff the graph is connected)
#     "minDegree",                # minimum degree of any vertex
#     "maxDegree",                # maximum degree of any vertex
#     "mean_degree",                # average degree of any vertex
#     "minCut",                   # minimum cut weight of the graph (might take a while to compute)
#     "FiedlerValue",             # second-highest eigenvalue of the Laplacian matrix
#     "Normalized_FiedlerValue",   # second-highest eigenvalue of the Normaized Laplacian matrix
#     "closeness_centr",                # average inverse of distance between any pair of vertices
#     "modularity",               # DEFINITION REQUIRED
#     "diameter",                 # maximum distance between any two vertices (NAN if not connected)
#     "betw_centr",              # max_{v} proportion of shortest paths going through vertex v
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
#   df$order = base::as.numeric(lapply(Graphs, gorder))
#   df$edges = base::as.numeric(lapply(Graphs, gsize))
#   df$connected = base::as.numeric(lapply(Graphs, is_connected))
#   df$minCut = base::as.numeric(lapply(Graphs, min_cut))
#   df$diameter = base::as.numeric(lapply(Graphs, diameter))
#   df$transitivity = base::as.numeric(lapply(Graphs, transitivity))
#   
#   # stuff that needs interim things:
#   degrees = lapply(Graphs,igraph::degree )
#   df$minDegree = base::as.numeric(lapply(degrees, min))
#   df$maxDegree = base::as.numeric(lapply(degrees, max))
#   df$mean_degree = base::as.numeric(lapply(degrees, mean))
#   
#   # stuff that lapply doesn't like so has to be done in a loop:
#   communities <-lapply(Graphs, cluster_walktrap) #lapply(Graphs, cluster_leading_eigen)
#   Adj <- lapply(Graphs, as_adjacency_matrix)
#   L <- lapply(Graphs, laplacian_matrix)
#   Norm_Lap<-lapply(Graphs, normalized_laplacian)
#   Fiedler.value=NULL
#   norm.fiedler.value=NULL
#   
#   for (i in 1:length(Graphs)) {
#     if (is.null(Graphs[[i]]$type)) { Graphs[[i]]$type = "untyped" }
#     df$modularity[i] <- modularity(communities[[i]])
#     df$spectral_radius[i] <- eigen(Adj[[i]], symmetric=TRUE, only.values=TRUE)$values[1]
#     
#     Fiedler.value[[i]]=eigen(L[[i]], symmetric=TRUE, only.values=TRUE)$values
#     
#     df$FiedlerValue[i] <- Fiedler.value[[i]][length(Fiedler.value[[i]])-1]
#     
#     norm.fiedler.value[[i]]=eigen(Norm_Lap[[i]], symmetric=TRUE, only.values=TRUE)$values
#     
#     df$Normalized_FiedlerValue[i] <- norm.fiedler.value[[i]][length(norm.fiedler.value[[i]])-1]
#     
#     df$eigen_centr[i] <- centr_eigen(Graphs[[i]])$centralization
#     df$deg_centr[i] <- centr_degree(Graphs[[i]])$centralization
#     
#     df$betw_centr[i] <- centr_betw(Graphs[[i]])$centralization
#     
#     df$max_component[i] <- max(components(Graphs[[i]])$csize)
#     df$mean_eccentr[i]<-mean(eccentricity(Graphs[[i]]))
#     df$radius[i]<-radius(Graphs[[i]])
#     df$mean_path_length[i]<-average.path.length(Graphs[[i]])
#     #df$trace[i]<-sum(diag(Adj[[i]]))
#     df$graph_energy[i]<-sum(abs(eigen(Adj[[i]], symmetric=TRUE, only.values=TRUE)$values))
#     df$min_triangle[i]= min(count_triangles(Graphs[[i]]))
#     df$mean_triangle[i]= mean(count_triangles(Graphs[[i]]))
#     df$sd_triangle[i]= sd(count_triangles(Graphs[[i]]))
#     df$max_triangle[i]= max(count_triangles(Graphs[[i]]))
#     df$num_triangle[i]= sum(count_triangles(Graphs[[i]]))
#     df$deg_assort_coef[i]=assortativity_degree(Graphs[[i]])
#     
#     df$threshold[i] <- 1/(df$spectral_radius[i])
#     
#     if (df$connected[i]==TRUE) {
#       df$closeness_centr[i] = mean(closeness(Graphs[[i]]))
#     } else { # handle the case where G isn't connected
#       df$closeness_centr[i] = -1
#     }
#   }
#   return (df)
# }
# 
# #calcGraphFeatures(x)
# ##----------- Graph Features----------#####
# ## Used to perform runs on multiple simulated graphs on any given network
# RunSimOnGraphFeatures<-function(Graphs, nreps=nreps,output_file=NULL, seed=-1) {
#   set.seed(1)
#   # ### Definition and initialization of parameters for graphfeatures
#   graphProperties=list()
#   # ### Definition and initialization of Graph Prefix
#   graphid=list();graphreplicate=list(); graphname=list(); GraphPrefix=list(); analysis=list()
#   for (g in 1:length(Graphs)){
#     for (reps in 1:nreps) {
#       ### Calculate the graph features for each simulated graph of all the synthetic networks
#       print(paste("Calculating graph features on", Graphs[[g]]$name))
#       #graphProperties[[reps]] <- calcGraphFeatures(Graphs[g])
#       graphProperties[[reps]] <- calcGraphFeatures(Graphs[g])
#     }
#     graphname[[g]]=Graphs[[g]]$type
#     graphid[[g]]=Graphs[[g]]$id
#     graphreplicate[[g]]=c(1:nreps)
#     GraphPrefix=cbind(graphname[[g]],graphid[[g]],graphreplicate[[g]])
#     colnames(GraphPrefix)=c("GraphName","GraphID","GraphReplicate")
#     analysis[[g]]=as.data.frame(cbind(GraphPrefix,graphProperties[[reps]]))
#     row.names(analysis[[g]])=1:nreps
#   }
#   All_results=do.call(rbind,analysis)
#   # write.csv(All_results, file=output_file)
#   return( All_results)
# }
# 
# 
# MKSFGraphs<- function(n=50, power=2, m=4) {
#   Graphs = list()#;G=list() # set up the list of Graphs first
#   initgraph= list()
#   i= 1
#   print("Creating SF Graphs")
#   Graphs[[i]] = sample_pa(n, power, m, directed=FALSE, algorithm="psumtree")
#   
#   Graphs[[i]]$type = "SF"
#   Graphs[[i]]$id = "1"
#   Graphs[[i]]$name="SF"
#   i <- i+1
#   return(Graphs)
# }
# 
# 
# sim.sf.net=function(power=powerVals,m=mVals,n=500,nreps=1){
#   net=NULL;data=NULL
#   dt=expand.grid(powerVals,mVals)
#   colnames(dt)=c("powerVals","mVals")
#   for (i in 1:nrow(dt)){
#     net[i]= MKSFGraphs(n,power=dt$powerVals[i],m=dt$mVals[i])
#     data=RunSimOnGraphFeatures(net,nreps)
#   }
#   df=cbind(power,m,data)  
#   return(df)
# }
# powerVals=1:10
# mVals=1:100
# ### 50 nodes 
# SF1.50=sim.sf.net(power=powerVals,m=mVals,nreps = 100,n=50)
# write.csv(SF1.50,"SF1.50.csv")
# 
# ### 100 nodes 
# SF1.100=sim.sf.net(power=powerVals,m=mVals,nreps = 100,n=100)
# write.csv(SF1.100,"SF1.100.csv")
# 
# ### 150 nodes 
# SF1.150=sim.sf.net(power=powerVals,m=mVals,nreps = 100,n=150)
# write.csv(SF1.150,"SF1.150.csv")
# 
# ### 200 nodes 
# SF1.200=sim.sf.net(power=powerVals,m=mVals,nreps = 100,n=200)
# write.csv(SF1.200,"SF1.200.csv")
# 
# ### 250 nodes 
# SF1.250=sim.sf.net(power=powerVals,m=mVals,nreps = 100,n=250)
# write.csv(SF1.250,"SF1.250.csv")
# 
# ### 300 nodes 
# SF1.300=sim.sf.net(power=powerVals,m=mVals,nreps = 100,n=300)
# write.csv(SF1.300,"SF1.300.csv")
# 
# ### 350 nodes 
# SF1.350=sim.sf.net(power=powerVals,m=mVals,nreps = 100,n=350)
# write.csv(SF1.350,"SF1.350.csv")
# 
# ### 400 nodes 
# SF1.400=sim.sf.net(power=powerVals,m=mVals,nreps = 100,n=400)
# write.csv(SF1.400,"SF1.400.csv")
# 
# ### 500 nodes 
# SF1.500=sim.sf.net(power=powerVals,m=mVals,nreps = 100,n=500)
# write.csv(SF1.500,"SF1.500.csv")
# 
# ### 750 nodes 
# SF1.750=sim.sf.net(power=powerVals,m=mVals,nreps = 100,n=750)
# write.csv(SF1.750,"SF1.750.csv")
# 
# ### 1000 nodes 
# SF1.1000=sim.sf.net(power=powerVals,m=mVals,nreps = 100,n=1000)
# write.csv(SF1.1000,"SF1.1000.csv")
# 
# ### 1500 nodes 
# SF1.1500=sim.sf.net(power=powerVals,m=mVals,nreps = 100,n=1500)
# write.csv(SF1.1500,"SF1.1500.csv")
# 
# ### 2000 nodes 
# SF1.2000=sim.sf.net(power=powerVals,m=mVals,nreps = 100,n=2000)
# write.csv(SF1.2000,"SF1.2000.csv")
# 
# ######## Saving all SF graphs to csv
# SF.50=read.csv("SF1.50.csv");SF.100=read.csv("SF1.100.csv");
# SF.150=read.csv("SF1.150.csv")
# SF.200=read.csv("SF1.100.csv");
# SF.250=read.csv("SF1.250.csv");SF.300=read.csv("SF1.300.csv")
# SF.350=read.csv("SF1.350.csv");SF.400=read.csv("SF1.400.csv")
# SF.500=read.csv("SF1.500.csv");SF.750=read.csv("SF1.750.csv")
# SF.1000=read.csv("SF1.1000.csv");SF.1500=read.csv("SF1.1500.csv")
# SF.2000=read.csv("SF1.2000.csv")
# df.sf=rbind(SF.50,SF.150,SF.200,SF.250,SF.300,SF.350,
#             SF.400, SF.500,SF.750,SF.1000,SF.1500,
#             SF.2000)
# 
# write.csv(df.sf,"df.sf.csv")
# 
# n <- 50
# m <- 520
# 
# # Work out how theye divide into each other
# rem <- m %% n
# div <- m %/% n
# 
# set.seed(123)
# if(rem != 0) {
#   g <- sample_smallworld(1, n, div+1, p = 0.001)
#   # Randomly delete the unwanted edges. Should be quite homegenous
#   g <- delete_edges(g, sample(1:gsize(g), size = gsize(g) - m))
# } else {
#   g <- sample_smallworld(1, n, div, p = 0.001)
# }
# 
# 
# 
# lat=make_lattice(dimvector=1000, nei=2, circular=T)
# plot(lat)
# 
# g <- sample_smallworld(dim=1, size=172, nei=6, p=0.05)
# g.new=rewire(make_lattice(dimvector=500, nei=3, circular=T), with=each_edge(p=0.9, loops=F))
# plot(g.new)
