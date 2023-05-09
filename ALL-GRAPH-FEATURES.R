#+++++++++++++++++++++++++++++++++++++++++++++++++++++
# Network properties for multiple graphs
#++++++++++++++++++++++++++++++++++++++++++++++++++

library(igraph)

# del_vertices=function(G){
#     Isolated = which(igraph::degree(G)==0)
#     Net=igraph::delete.vertices(G, Isolated)
#     
#     return(Net)
#   }

largest_comp=function(df){
  G=graph_from_data_frame(as.matrix(df),directed=FALSE)
  df.graph=igraph::simplify(G,remove.multiple = T,remove.loops = T)
  Isolated = which(igraph::degree(df.graph)==0)
  Net=igraph::delete.vertices(df.graph, Isolated)
  components = igraph::clusters(Net, mode="weak")
  biggest_cluster_id = which.max(components$csize)
  vert_ids = V(df.graph)[components$membership== biggest_cluster_id]
  graph=igraph::induced_subgraph(df.graph, vert_ids)
}
#create graph object
#FolderName = c('C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi/Networks')

Network.Summary <- function(FolderName){
  
  Net.lists = list.files(FolderName)
  Net.Feat.Summary = data.frame()
 
     g = lapply(Net.lists,read.table)
     GraphNames = cbind(lapply(Net.lists,function(x) gsub(".edges","",x)))
     
     G = lapply(g, largest_comp)
     
     g.features=calcGraphFeatures(G)
     all_graphs=cbind(GraphNames, g.features)
  return(all_graphs)
}

FolderName = c('networks')
z=Network.Summary(FolderName )
z
#calcGraphFeatures(Graphs=m)


x=makeSpatialGraphs(496,0.06)#m=4,f=0.017,s=9, e=984
y=calcGraphFeatures(x)
y


r1=runif(1000,0.74,0.78);n1=35
r2=runif(1000,0.44,0.48);n2=17
r3=runif(1000,0.145,0.15);n3=117
r4=runif(1000,0.24,0.25);n4=20
r5=runif(1000,0.17,0.18);n5=62
r6=runif(1000,0.055,0.06);n6=496

R=r6
node.num=n6
 
 
 net=NULL
 for (i in 1:length(R)){
   # net=fastSpatialNetwork(n=node.num,r=R[[i]],makeConnected=TRUE, keepCellsSeparate=FALSE)
   net[i]=makeSpatialGraphs(node.size =node.num,Radius=R[i])
 }

 
 Data=RunSimOnGraphFeatures(net,nreps = 1)
 
 DF.DATA=cbind(R,Data)
 x=write.csv(DF.DATA,"node496-r1000sample.csv")
 