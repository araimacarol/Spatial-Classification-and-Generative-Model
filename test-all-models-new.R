#x.1=all.nets(model=lasso.pred.bestLambda,data="mammalia-hyena-networkb.edges",graphname="hyena")
library(intergraph);library(ergm)
data="mammalia-dolphin-social.edges"
graphname="dolphin"

df <- read.table(data)  
df_new=graph_from_data_frame(as.matrix(df),directed=FALSE)
df.graph=igraph::simplify(df_new,remove.multiple = T,remove.loops = T)
components = igraph::clusters(df.graph , mode="weak")
biggest_cluster_id = which.max(components$csize)
vert_ids = V(df.graph)[components$membership== biggest_cluster_id]
graph=igraph::induced_subgraph(df.graph, vert_ids)

graph$name=graphname
graph$type=graphname
graph$id="1"
G=list(graph)

###---Graph--Features--###
data=RunSimOnGraphFeatures(G,nreps = 1)
data

newdata=data %>%dplyr::select(-c(GraphID,order, GraphReplicate,max_component,minDegree,threshold,maxDegree,GraphName,connected))

#newdata=data %>%dplyr::select(c(FiedlerValue,Normalized_FiedlerValue,spectral_radius,transitivity,modularity))

#  newdata=data %>%dplyr::select(c(edges,minCut,FiedlerValue,closeness, 
# modularity,diameter,betweenness,centrality_eigen))

predicted_radius.lasso=predict(lasso.mod, s = bestlam, newx = as.matrix(newdata))
predicted_radius.ridge=predict(ridge.mod , s = bestlam, newx = as.matrix(newdata))
predicted_radius.rfmodel2=predict(rfmodel2,newdata)
predicted_radius.glm=predict(lm.model.1,newdata)




###--Predictiv--models
# spatial.net.lasso=makeSpatialGraphs(node.size=vcount(df.graph),Radius=predicted_radius.lasso)
# spatial.net.ridge=makeSpatialGraphs(node.size=vcount(df.graph),Radius=predicted_radius.ridge)
# spatial.net.rfmodel2=makeSpatialGraphs(node.size=vcount(df.graph),Radius=predicted_radius.rfmodel2)
# spatial.net.glm=makeSpatialGraphs(node.size=vcount(df.graph),Radius=predicted_radius.glm)

###---Theoretical--models
z=makeTheoGraphs(nSamples=1, order=62,edges=159,dim.sw=2,nei.sw=1,
                 p.sw=0.2,power.sf=4,m.sf=2,r.sp=0.19,dim.lat=2,nei.lat=1)
z

th.models=z

##--ERGMs---
set.seed(569)

ergm_net <- asNetwork(graph)
summary(ergm_net ~ edges)
ergm.model <- ergm(ergm_net ~ edges) #1st fitted ERGM
summary(ergm.model)
gf1=gof(ergm.model)
ergm.sim.model <- simulate(ergm.model,nsim=100)

ergm.sim.net=lapply(ergm.sim.model,asIgraph)
edgeCount_ergm=c(lapply(ergm.sim.net, ecount))
edgeCount_ergm

ergm.graph=ergm.sim.net[69]

#fastSpatialNetwork(n=vcount(graph),r=predicted_radius,makeConnected=TRUE, keepCellsSeparate=FALSE)

###--Graph--Features--for--empirical--and--spatial--net
# spatial.net$name="Spatial"
# spatial.net$type="Spatial"
# spatial.net$id="1"
set.seed(7349)
net=c(G,ergm.graph,spatial.net.lasso,spatial.net.ridge,spatial.net.rfmodel2,spatial.net.glm,th.models)
feat=RunSimOnGraphFeatures(net,nreps = 1)
feat$GraphName=c(graphname,"ERGM","SPNet.Lasso","SPNet.Ridge","SPNet.RF","SPNet.GLM","ER",
                 "SW","SF","SP","LAT")
feat

write.csv(feat,"allnets3.csv")
#new.list <- list(spatial.net[[1]],graph,feat)
#x=c(G,spatial.net)
#x$summary=feat
#return(x)


hist(igraph::degree(graph),main = "hyena")
hist(igraph::degree(spatial.net.lasso[[1]]),main="splasso")
hist(igraph::degree(spatial.net.ridge[[1]]), main="ridge")
hist(igraph::degree(spatial.net.rfmodel2[[1]]),main="rf")
hist(igraph::degree(spatial.net.glm[[1]]),main="glm")
hist(igraph::degree(ergm.graph[[1]]),main="ergm")

par(mar=c(2,2,2,2))
par(mfrow=c(1,3))
g1=th.models[4]#spatial.net.ridge#spatial.net.lasso
g2=spatial.net.lasso
g3=graph

plot(g1[[1]]);plot(g2[[1]]);plot(g3)
