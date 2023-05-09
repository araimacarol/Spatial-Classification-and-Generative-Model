setwd("C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi-New")
source("SPATIAL-PIPELINE-NEW.R")
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                        libraries
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(rsample)      # data splitting 
library(gbm)          # basic implementation
library(xgboost)      # a faster implementation of gbm
library(caret)        # an aggregator package for performing many machine learning models
#library(h2o)          # a java-based platform
library(pdp)          # model visualization
library(ggplot2)      # model visualization
library(lime)         # model visualization
library(randomForest)
library(Matrix)
library(magrittr)
library(dplyr)
library(glmnet)
library(caret)
library(leaps)
library("mgcv")
library(elasticnet)
library("gam")
library("brnn")
library(intergraph)
library("nnet")
library("clusterSim")
library("randomForest")
library("gbm")
library("plyr")
library(ggbiplot)
library("rpart")
library("mlbench")
library("tidyverse")
library("MLmetrics")
library("readxl")
library("varImp")
library(caTools)
library(mrIML)
library(tidymodels)
library(future)
library("future.apply")
library("Hmisc")
library(corrplot)
library("PerformanceAnalytics")
library(mlbench)
library(igraph)
library(purrr)
library(RSpectra)
library(tidymodels)
tidymodels_prefer()
library(beans)
library(learntidymodels)
library(ggforce)
library(patchwork)
library(bestNormalize)
library(Metrics)
library(MASS)
library(ergm)
library(psych)
library(betareg)
library(lmtest)
library(emmeans)
library(car)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# library(igraph)
# library(randomForest)

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
#function(x)read.table(x,header=TRUE,sep="\t")
library(data.table)
Network.Summary <- function(FolderName){
  
  Net.lists = list.files(FolderName, full.names=TRUE)
  Net.Feat.Summary = data.frame()
  
  g = lapply(Net.lists,fread)
  GraphName = cbind(lapply(Net.lists,function(x) gsub(".edges","",x)))
  
  G = lapply(g, largest_comp)
  
  g.features=calcGraphFeatures(G)
  all_graphs=cbind(GraphName, g.features)
  results=list(G,all_graphs)
  return(results)
}

FolderName = c('Animal-Social-Networks')
empdata=Network.Summary(FolderName)
empdata[[1]]
empdata[[2]]
# FolderName = c('C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi-New/networks')

#####++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++########
#    Making 100 instances of each of the synthetic graphs 
#####++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++########

# Create a list of graphs for simulation experiment
##########################################################################################################
makeSynthGraphs <- function(nSamples=10, order=25,edges=45,dim.sw=2,nei.sw=2,
                           p.sw=0.2,power.sf1=2,m.sf1=10,power.sf2=2,m.sf2=10,r.sp=0.345) {
  
  graphTypes <- c("Erdos-Renyi", "Small World", "Scale Free", "Square Lattice", "Spatial")
  
  Graphs <- list(mode="any", nSamples * length(graphTypes)) # set up the list of Graphs first
  initgraph=list(mode="any", nSamples * length(graphTypes))
  
  print("Creating Erdos-Renyi graphs")
  i <- 1
  for (j in 1:nSamples) {
    Graphs[[i+j-1]] <-igraph::erdos.renyi.game(order,edges, type = "gnm" , directed = F , loops = F)
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
    #initgraph=simplify(Graphs[[i+j-1]],remove.multiple = T,remove.loops = T)
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
  
  print("Creating scale-free linear attachment graphs using the Barabasi-Albert method")
  for (j in 1:nSamples) {
    Graphs[[i+j-1]] <- igraph::sample_pa(n=order, power=power.sf1, m=m.sf1, directed=FALSE, algorithm="psumtree")
    initgraph=Graphs[[i+j-1]]
    components<- igraph::clusters(initgraph , mode="weak")
    biggest_cluster_id<- which.max(components$csize)
    # # ids
    vert_ids<- V(initgraph)[components$membership== biggest_cluster_id]
    # # subgraph
    Graphs[[i+j-1]]=igraph::induced_subgraph(initgraph, vert_ids)
    Graphs[[i+j-1]]$type <- "SF-linear"
    Graphs[[i+j-1]]$name <- "SF-linear"
    Graphs[[i+j-1]]$id <- j
  }
  i <- i+nSamples
  
  print("Creating scale-free non-linear attachment graphs using the Barabasi-Albert method")
  for (j in 1:nSamples) {
    Graphs[[i+j-1]] <- igraph::sample_pa(n=order, power=power.sf2, m=m.sf2, directed=FALSE, algorithm="psumtree")
    initgraph=Graphs[[i+j-1]]
    components<- igraph::clusters(initgraph , mode="weak")
    biggest_cluster_id<- which.max(components$csize)
    # # ids
    vert_ids<- V(initgraph)[components$membership== biggest_cluster_id]
    # # subgraph
    Graphs[[i+j-1]]=igraph::induced_subgraph(initgraph, vert_ids)
    Graphs[[i+j-1]]$type <- "SF-non-linear"
    Graphs[[i+j-1]]$name <- "SF-non-linear"
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
  
  return(Graphs)
}

#####---Wildbird--synthetic--equivalent--Networks---
set.seed(2748)
AvesWbirdSynth.Net1=makeSynthGraphs(nSamples=100, order=131,edges=1444,dim.sw=2,nei.sw=3,
                            p.sw=0.2,power.sf1=1,m.sf1=12,power.sf2=2,m.sf2=12,r.sp=0.27) 
  
AvesWbirdSynth.Net2=makeSynthGraphs(nSamples=100, order=82,edges=1170,dim.sw=2,nei.sw=3.5,
                                    p.sw=0.2,power.sf1=1,m.sf1=16,power.sf2=2,m.sf2=16,
                                    r.sp=0.42)


AvesWbirdSynth.Net3=makeSynthGraphs(nSamples=100, order=126,edges=1615,dim.sw=2,nei.sw=3.8,
                                    p.sw=0.2,power.sf1=1,m.sf1=14,power.sf2=2,m.sf2=14,
                                    r.sp=0.30)


AvesWbirdSynth.Net4=makeSynthGraphs(nSamples=100, order=93,edges=1689,dim.sw=2,nei.sw=4,
                                    p.sw=0.2,power.sf1=1,m.sf1=20,power.sf2=2,m.sf2=20,
                                    r.sp=0.43)

AvesWbirdSynth.Net5=makeSynthGraphs(nSamples=100, order=145,edges=2512,dim.sw=2,nei.sw=4,
                                    p.sw=0.2,power.sf1=1,m.sf1=19,power.sf2=2,m.sf2=19,
                                    r.sp=0.33)


AvesWbirdSynth.Net6=makeSynthGraphs(nSamples=100, order=136,edges=2798,dim.sw=2,nei.sw=4,
                                    p.sw=0.2,power.sf1=1,m.sf1=23,power.sf2=2,m.sf2=23,
                                    r.sp=0.38)

AvesWbirdSynth.Net7=makeSynthGraphs(nSamples=100, order=202,edges=4574,dim.sw=2,nei.sw=4,
                                    p.sw=0.2,power.sf1=1,m.sf1=24,power.sf2=2,m.sf2=24,
                                    r.sp=0.32)

Theoretical.nets.wildbird=c(AvesWbirdSynth.Net1,AvesWbirdSynth.Net2,AvesWbirdSynth.Net3,
                       AvesWbirdSynth.Net4,AvesWbirdSynth.Net5,AvesWbirdSynth.Net6,
                       AvesWbirdSynth.Net7)



#####---Hyena--synthetic--equivalent--Networks---
Hyena.Net1=makeSynthGraphs(nSamples=100, order=35,edges=521,dim.sw=2,nei.sw=4,
                                    p.sw=0.2,power.sf1=1,m.sf1=22,power.sf2=2,m.sf2=22,
                           r.sp=0.79) 

Hyena.Net2=makeSynthGraphs(nSamples=100, order=36,edges=585,dim.sw=2,nei.sw=5,
                           p.sw=0.2,power.sf1=1,m.sf1=26,power.sf2=2,m.sf2=26,
                           r.sp=0.86) 

Hyena.Net3=makeSynthGraphs(nSamples=100, order=35,edges=509,dim.sw=2,nei.sw=4,
                           p.sw=0.2,power.sf1=1,m.sf1=21,power.sf2=2,m.sf2=21,r.sp=0.8)

Theoretical.nets.Hyena=c(Hyena.Net1,Hyena.Net2,Hyena.Net3)

###--Ant--synthetic--equivalent--Networks
Ant.Net1=makeSynthGraphs(nSamples=100, order=113,edges=4550,dim.sw=2,nei.sw=6,
                        p.sw=0.2,power.sf1=1,m.sf1=53,power.sf2=2,m.sf2=53,r.sp=0.68)

Ant.Net2=makeSynthGraphs(nSamples=100, order=113,edges=4573,dim.sw=2,nei.sw=6,
                         p.sw=0.2,power.sf1=1,m.sf1=53,power.sf2=2,m.sf2=53,r.sp=0.682)

###--Bats--synthetic--equivalent--Networks
Bat.Net1=makeSynthGraphs(nSamples=100, order=43,edges=546,dim.sw=2,nei.sw=3,
                         p.sw=0.2,power.sf1=1,m.sf1=16,power.sf2=2,m.sf2=16,r.sp=0.6)

Bat.Net2=makeSynthGraphs(nSamples=100, order=21,edges=72,dim.sw=2,nei.sw=1.5,
                         p.sw=0.2,power.sf1=1,m.sf1=4,power.sf2=2,m.sf2=4,r.sp=0.4)

###--Tortoise--synthetic--equivalent--Networks
Tortoise.Net=makeSynthGraphs(nSamples=100, order=496,edges=984,dim.sw=2,nei.sw=1,
                         p.sw=0.2,power.sf1=1,m.sf1=1,power.sf2=2,m.sf2=1,r.sp=0.06)

###--Sparrow--Social--synthetic--equivalent--Networks
Sparrow.Net1=makeSynthGraphs(nSamples=100, order=31,edges=211,dim.sw=2,nei.sw=2,
                         p.sw=0.2,power.sf1=1,m.sf1=8,power.sf2=2,m.sf2=8,r.sp=0.47)

Sparrow.Net2=makeSynthGraphs(nSamples=100, order=40,edges=305,dim.sw=2,nei.sw=3,
                             p.sw=0.2,power.sf1=1,m.sf1=9,power.sf2=2,m.sf2=9,r.sp=0.45)

###--Sparrowlyon--synthetic--equivalent--Networks
Sparrowlyon.Net1=makeSynthGraphs(nSamples=100, order=46,edges=348,dim.sw=2,nei.sw=2.5,
                             p.sw=0.2,power.sf1=1,m.sf1=8,power.sf2=2,m.sf2=8,r.sp=0.43)

Sparrowlyon.Net2=makeSynthGraphs(nSamples=100, order=27,edges=163,dim.sw=2,nei.sw=2,
                                 p.sw=0.2,power.sf1=1,m.sf1=7,power.sf2=2,m.sf2=7,r.sp=0.48)

###--Aves--Weaver--social--synthetic--equivalent--Networks
Aves.WeaverSocial.Net=makeSynthGraphs(nSamples=100, order=117,edges=304,dim.sw=2,nei.sw=1,
                                 p.sw=0.2,power.sf1=1,m.sf1=3,power.sf2=2,m.sf2=3,r.sp=0.13)

Aves.SongbirdSocial.Net=makeSynthGraphs(nSamples=100, order=108,edges=1026,dim.sw=2,nei.sw=3,
                              p.sw=0.2,power.sf1=1,m.sf1=10,power.sf2=2,m.sf2=10,r.sp=0.282)

##--Cattle--synthetic--equivalent--Networks
CattleDom.Net=makeSynthGraphs(nSamples=100, order=28,edges=205,dim.sw=2,nei.sw=3,
                              p.sw=0.2,power.sf1=1,m.sf1=9,power.sf2=1,m.sf2=9,r.sp=0.6)

###--Dolphin--synthetic--equivalent--Networks
DolphinSoc.Net=makeSynthGraphs(nSamples=100, order=62,edges=159,dim.sw=2,nei.sw=1,
                              p.sw=0.2,power.sf1=1,m.sf1=3,power.sf2=1,m.sf2=3,r.sp=0.17)


#######----Making--Arbitrary--Graphs--of--Different--Sizes--for-each--Model
net50_1=makeSynthGraphs(nSamples=100, order=50,edges=256,dim.sw=2,nei.sw=2,
                p.sw=0.2,power.sf1=1,m.sf1=5,power.sf2=2,m.sf2=5,r.sp=0.3)

net50_2=makeSynthGraphs(nSamples=100, order=50,edges=589,dim.sw=2,nei.sw=3,
                        p.sw=0.2,power.sf1=1,m.sf1=14,power.sf2=2,m.sf2=14,r.sp=0.5)

net50_3=makeSynthGraphs(nSamples=100, order=50,edges=1107,dim.sw=2,nei.sw=5,
                        p.sw=0.2,power.sf1=1,m.sf1=28,power.sf2=2,m.sf2=28,r.sp=0.88)

net100_1=makeSynthGraphs(nSamples=100, order=100,edges=1214,dim.sw=2,nei.sw=3,
                        p.sw=0.2,power.sf1=1,m.sf1=13,power.sf2=2,m.sf2=13,r.sp=0.33)

net100_2=makeSynthGraphs(nSamples=100, order=100,edges=2674,dim.sw=2,nei.sw=5,
                         p.sw=0.2,power.sf1=1,m.sf1=32,power.sf2=2,m.sf2=32,r.sp=0.56)

net100_3=makeSynthGraphs(nSamples=100, order=100,edges=4218,dim.sw=2,nei.sw=7,
                         p.sw=0.2,power.sf1=1,m.sf1=62,power.sf2=2,m.sf2=62,r.sp=0.79)

######--Comibine--All--Theoretical--Networks
AllTheoreticalNetwork=c(Theoretical.nets.wildbird,Theoretical.nets.Hyena,Ant.Net1,
                        Ant.Net2,Bat.Net1,Bat.Net2,Tortoise.Net,
                        Sparrow.Net1,Sparrow.Net2,Sparrowlyon.Net1,Sparrowlyon.Net2,
                        Aves.WeaverSocial.Net,Aves.SongbirdSocial.Net,CattleDom.Net,
                        DolphinSoc.Net,net50_1,net50_2,net50_3,
                        net100_1,net100_2,net100_3)

# n=35:534
# e=100:5000
# 
# n=5:10
# e=7:10
# pos=expand.grid(n,e)
# list.er=NULL
# generate_er=function(N,x){
# x=pos[,1]
# y=pos[,2]
# k=k+1
# list.er=NULL
# for (i in 1: length(x)){
#   for (j in 1: length(y)){
#     list.er[[k]]=erdos.renyi.game(x[i],y[j], type = "gnm" , directed = F , loops = F)
# }
#   }
#   return(results)
# }
# 
# generate_er(5,x)
#####++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++########
#    Calculating Graph Features for the synthetic graphs 
#####++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++########
SynthGraphFeatures=RunSimOnGraphFeatures(AllTheoreticalNetwork, nreps=1,output_file="TrainingGraphFeatures.csv", seed=927)
head(SynthGraphFeatures)

#####--Making--GraphNames--as--factor
SynthGraphFeatures$GraphName <- as.factor(SynthGraphFeatures$GraphName)
#####++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++########
#    Splitting Data into Training and Testing
#####++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++####
#####+
df=read.csv("TrainingGraphFeatures.csv",sep = ",", header = T)

Data= df%>%dplyr::select(c(GraphName,order,edges,mean_degree,minCut,
                FiedlerValue,Normalized_FiedlerValue,closeness,modularity,diameter,
                betweenness,transitivity,spectral_radius,centrality_eigen))
##--Shuffle--data
df<- Data[sample(1:nrow(Data)), ]##shuffle row indices and randomly re order 
head(df)
dim(df)

Train_and_Test <- function(Data, size = 0.8, train = TRUE) {
  n_row = nrow(Data)
  df_split=size*n_row
  train_sample = 1:df_split
  if (train == TRUE) {
    return (Data[train_sample, ])
  } else {
    return (Data[-train_sample, ])
  }
}

train=Train_and_Test(df,0.8,train = T)
test=Train_and_Test(df,0.8,train = F)

dim(train)
dim(test)

head(train)
head(test)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Random forest Graph Classification Model 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#labels <- as.factor(SynthGraphFeatures$GraphName)
##--Tuning--the--number--of--mtrys
set.seed((6283))
tn=tuneRF(train[,-1],train$GraphNames,stepFactor=0.5,plot = TRUE,ntreeTry = 500,
          trace = TRUE,improve = 0.5) #works only if mtry
#> the number of variables(features). This is because mtry is the number of randomly sampled variable as candidate at each split. Base case use mtry=2  when this happens
train.control<- trainControl(method = "cv", number = 10,savePredictions = "all")
tn=as.data.frame(tn)
tn.min=tn$mtry[tn$OOBError== min(tn$OOBError)] 
set.seed(6748)
RfClassificationmodel <- randomForest::randomForest(GraphName~ ., data = train,
                                        trainControl=train.control,
                                        mtry=tn.min,ntree=500,proximity=TRUE,
                                        importance=TRUE)
print(RfClassificationmodel)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Predict labels for empirical networks
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
RealNet=empdata[[2]]
RealNet= RealNet%>%dplyr::select(c(GraphName,order,edges,mean_degree,minCut,
                                           FiedlerValue,Normalized_FiedlerValue,closeness,modularity,diameter,
                                           betweenness,transitivity,spectral_radius,centrality_eigen))

predicted_label <- predict(RfClassificationmodel, newdata = RealNet[,])
predicted_label

confusionMatrix(predicted_label, train$GraphName)

#####Avesaves-wildbird-net1 predicted as Erdos Renyi
er.aves1=AvesWbirdSynth.Net1 [[1]]
sw.aves1=AvesWbirdSynth.Net1[[110]]
sf.linear.aves1=AvesWbirdSynth.Net1[[210]]
sf.nonlinear.aves1=AvesWbirdSynth.Net1[[310]]
sp.aves1=AvesWbirdSynth.Net1[[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[9]],vertex.label=NA,vertex.size=2)
plot(er.aves1,vertex.label=NA,vertex.size=2)
plot(sw.aves1,vertex.label=NA,vertex.size=2)
plot(sf.linear.aves1,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.aves1,vertex.label=NA,vertex.size=2)
plot(sp.aves1,vertex.label=NA,vertex.size=2)


#####Aves wildbird net2 predicted as Spatial
er.aves2=AvesWbirdSynth.Net2 [[1]]
sw.aves2=AvesWbirdSynth.Net2[[110]]
sf.linear.aves2=AvesWbirdSynth.Net2[[210]]
sf.nonlinear.aves2=AvesWbirdSynth.Net2[[310]]
sp.aves2=AvesWbirdSynth.Net2[[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[10]],vertex.label=NA,vertex.size=2)
plot(er.aves2,vertex.label=NA,vertex.size=2)
plot(sw.aves2,vertex.label=NA,vertex.size=2)
plot(sf.linear.aves2,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.aves2,vertex.label=NA,vertex.size=2)
plot(sp.aves2,vertex.label=NA,vertex.size=2)

#####Tortoise predicted as scale-free with linear attachment
er.tortoise=Tortoise.Net[[1]]
sw.tortoise=Tortoise.Net[[110]]
sf.linear.tortoise=Tortoise.Net[[210]]
sf.nonlinear.tortoise=Tortoise.Net[[310]]
sp.tortoise=Tortoise.Net[[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[29]],vertex.label=NA,vertex.size=2)
plot(er.tortoise,vertex.label=NA,vertex.size=2)
plot(sw.tortoise,vertex.label=NA,vertex.size=2)
plot(sf.linear.tortoise,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.tortoise,vertex.label=NA,vertex.size=2)
plot(sp.tortoise,vertex.label=NA,vertex.size=2)

##### Bat roosting predicted as spatial
er.bats1=Bat.Net1 [[1]]
sw.bats1=Bat.Net1[[110]]
sf.linear.bats1=Bat.Net1[[210]]
sf.nonlinear.bats1=Bat.Net1[[310]]
sp.bats1=Bat.Net1[[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[19]],vertex.label=NA,vertex.size=2)
plot(er.bats1,vertex.label=NA,vertex.size=2)
plot(sw.bats1,vertex.label=NA,vertex.size=2)
plot(sf.linear.bats1,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.bats1,vertex.label=NA,vertex.size=2)
plot(sp.bats1,vertex.label=NA,vertex.size=2)

##### Vampire Bat foodsharing predicted as small world
er.bats2=Bat.Net2 [[1]]
sw.bats2=Bat.Net2[[110]]
sf.linear.bats2=Bat.Net2[[210]]
sf.nonlinear.bats2=Bat.Net2[[310]]
sp.bats2=Bat.Net2[[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[28]],vertex.label=NA,vertex.size=2)
plot(er.bats2,vertex.label=NA,vertex.size=2)
plot(sw.bats2,vertex.label=NA,vertex.size=2)
plot(sf.linear.bats2,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.bats2,vertex.label=NA,vertex.size=2)
plot(sp.bats2,vertex.label=NA,vertex.size=2)


#####Ant Colony net1 (day 1) predicted as spatial
er.ants1=Ant.Net1 [[1]]
sw.ants1=Ant.Net1[[110]]
sf.linear.ants1=Ant.Net1[[210]]
sf.nonlinear.ants1=Ant.Net1[[310]]
sp.ants1=Ant.Net1[[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[16]],vertex.label=NA,vertex.size=2)
plot(er.ants1,vertex.label=NA,vertex.size=2)
plot(sw.ants1,vertex.label=NA,vertex.size=2)
plot(sf.linear.ants1,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.ants1,vertex.label=NA,vertex.size=2)
plot(sp.ants1,vertex.label=NA,vertex.size=2)


#####Ant Colony Net2 (day 2) predicted as spatial
er.ants2=Ant.Net2 [[1]]
sw.ants2=Ant.Net2[[110]]
sf.linear.ants2=Ant.Net2[[210]]
sf.nonlinear.ants2=Ant.Net2[[310]]
sp.ants2=Ant.Net2[[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[17]],vertex.label=NA,vertex.size=2)
plot(er.ants2,vertex.label=NA,vertex.size=2)
plot(sw.ants2,vertex.label=NA,vertex.size=2)
plot(sf.linear.ants2,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.ants2,vertex.label=NA,vertex.size=2)
plot(sp.ants2,vertex.label=NA,vertex.size=2)

#####Hyena net1 predicted as Spatial
er.hyena1=Hyena.Net1 [[1]]
sw.hyena1=Hyena.Net1 [[110]]
sf.linear.hyena1=Hyena.Net1 [[210]]
sf.nonlinear.hyena1=Hyena.Net1 [[310]]
sp.hyena1=Hyena.Net1 [[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[22]],vertex.label=NA,vertex.size=2)
plot(er.hyena1,vertex.label=NA,vertex.size=2)
plot(sw.hyena1,vertex.label=NA,vertex.size=2)
plot(sf.linear.hyena1,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.hyena1,vertex.label=NA,vertex.size=2)
plot(sp.hyena1,vertex.label=NA,vertex.size=2)


#####Hyena net2 predicted as Spatial
er.hyena2=Hyena.Net2 [[1]]
sw.hyena2=Hyena.Net2 [[110]]
sf.linear.hyena2=Hyena.Net2 [[210]]
sf.nonlinear.hyena2=Hyena.Net2 [[310]]
sp.hyena2=Hyena.Net2 [[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[23]],vertex.label=NA,vertex.size=2)
plot(er.hyena2,vertex.label=NA,vertex.size=2)
plot(sw.hyena2,vertex.label=NA,vertex.size=2)
plot(sf.linear.hyena2,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.hyena2,vertex.label=NA,vertex.size=2)
plot(sp.hyena2,vertex.label=NA,vertex.size=2)


#####Hyena net3 predicted as Spatial
er.hyena3=Hyena.Net3 [[1]]
sw.hyena3=Hyena.Net3 [[110]]
sf.linear.hyena3=Hyena.Net3 [[210]]
sf.nonlinear.hyena3=Hyena.Net3 [[310]]
sp.hyena3=Hyena.Net3 [[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[24]],vertex.label=NA,vertex.size=2)
plot(er.hyena3,vertex.label=NA,vertex.size=2)
plot(sw.hyena3,vertex.label=NA,vertex.size=2)
plot(sf.linear.hyena3,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.hyena3,vertex.label=NA,vertex.size=2)
plot(sp.hyena3,vertex.label=NA,vertex.size=2)


#####Cattle dominance predicted as Scale-free with non linear attachment
er.cattle=CattleDom.Net [[1]]
sw.cattle=CattleDom.Net [[110]]
sf.linear.cattle=CattleDom.Net [[210]]
sf.nonlinear.cattle=CattleDom.Net [[310]]
sp.cattle=CattleDom.Net [[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[20]],vertex.label=NA,vertex.size=2)
plot(er.cattle,vertex.label=NA,vertex.size=2)
plot(sw.cattle,vertex.label=NA,vertex.size=2)
plot(sf.linear.cattle,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.cattle,vertex.label=NA,vertex.size=2)
plot(sp.cattle,vertex.label=NA,vertex.size=2)

#####Dolphine social predicted as Erdos Renyi
er.dolphine=DolphinSoc.Net [[1]]
sw.dolphine=DolphinSoc.Net [[110]]
sf.linear.dolphine=DolphinSoc.Net [[210]]
sf.nonlinear.dolphine=DolphinSoc.Net [[310]]
sp.dolphine=DolphinSoc.Net [[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[21]],vertex.label=NA,vertex.size=2)
plot(er.dolphine,vertex.label=NA,vertex.size=2)
plot(sw.dolphine,vertex.label=NA,vertex.size=2)
plot(sf.linear.dolphine,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.dolphine,vertex.label=NA,vertex.size=2)
plot(sp.dolphine,vertex.label=NA,vertex.size=2)

#####sparrowlyon net1 predicted spatial
er.sparrowlyon1=Sparrowlyon.Net1 [[1]]
sw.sparrowlyon1=Sparrowlyon.Net1 [[110]]
sf.linear.sparrowlyon1=Sparrowlyon.Net1 [[210]]
sf.nonlinear.sparrowlyon1=Sparrowlyon.Net1 [[310]]
sp.sparrowlyon1=Sparrowlyon.Net1 [[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[5]],vertex.label=NA,vertex.size=2)
plot(er.sparrowlyon1,vertex.label=NA,vertex.size=2)
plot(sw.sparrowlyon1,vertex.label=NA,vertex.size=2)
plot(sf.linear.sparrowlyon1,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.sparrowlyon1,vertex.label=NA,vertex.size=2)
plot(sp.sparrowlyon1,vertex.label=NA,vertex.size=2)

#####sparrowlyon net2 predicted spatial
er.sparrowlyon2=Sparrowlyon.Net2 [[1]]
sw.sparrowlyon2=Sparrowlyon.Net2 [[110]]
sf.linear.sparrowlyon2=Sparrowlyon.Net2 [[210]]
sf.nonlinear.sparrowlyon2=Sparrowlyon.Net2 [[310]]
sp.sparrowlyon2=Sparrowlyon.Net2 [[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[6]],vertex.label=NA,vertex.size=2)
plot(er.sparrowlyon2,vertex.label=NA,vertex.size=2)
plot(sw.sparrowlyon2,vertex.label=NA,vertex.size=2)
plot(sf.linear.sparrowlyon2,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.sparrowlyon2,vertex.label=NA,vertex.size=2)
plot(sp.sparrowlyon2,vertex.label=NA,vertex.size=2)
# spweaver=fastSpatialNetwork(n =vcount(empdata[[1]][[1]]), r = 0.14, makeConnected=T)
# spweaver
# sfweaver=sample_pa(n=vcount(empdata[[1]][[2]]),power=1,m=3,directed=FALSE)
# sfweaver
# par(mfrow=c(1,3))
# plot(empdata[[1]][[2]]);plot(spweaver);plot(sfweaver)#plot(Hyena.Net2[[410]])#plot(AvesWbirdSynth.Net5[[310]])

# ####--Compare--Graph--Features--for--empirical--and--the--synthetic--class
# emp.g=empdata[[1]][[7]]
# emp.g$name="wildbirdnet7"
# emp.g$type="wildbirdnet7"
# emp.g$id="1"
# 
# synth.g=AvesWbirdSynth.Net7[[401]]
# t1=list(emp.g,synth.g)
# RunSimOnGraphFeatures(t1,nreps = 1)



varImpPlot(RfClassificationmodel)


# #####--ERDOS--RENYI--#######
# er.simgraph=function(nodes=10,N=10,edges=200){
#   er.graphs=NULL
#   for (i in 1:N){
#     er.graphs[[i]]=erdos.renyi.game(nodes,edges, type = "gnm" , directed = F , loops = F)
#   }
#   return(er.graphs)
# }
# 
# 
# AvesWildbirdNet1.ER=er.simgraph(nodes=131,N=100,edges=1444)
# AvesWildbirdNet2.ER=er.simgraph(nodes=82,N=100,edges=1170)
# AvesWildbirdNet3.ER=er.simgraph(nodes=126,N=100,edges=1615)
# AvesWildbirdNet4.ER=er.simgraph(nodes=93,N=100,edges=1689)
# AvesWildbirdNet5.ER=er.simgraph(nodes=145,N=100,edges=2512)
# AvesWildbirdNet6.ER=er.simgraph(nodes=136,N=100,edges=2798)
# AvesWildbirdNet7.ER=er.simgraph(nodes=202,N=100,edges=4574)
# 
# ER.ALL=c(AvesWildbirdNet1.ER,AvesWildbirdNet2.ER,AvesWildbirdNet3.ER,
#              AvesWildbirdNet4.ER,AvesWildbirdNet5.ER,AvesWildbirdNet6.ER,
#              AvesWildbirdNet7.ER)
# length(ER.ALL)
# 
# #####--SMALL--WORLD--#######
# makeTheoGraphs() 



# n <- 100 # number of nodes
# m <- 4   # number of edges to add for each new node in Erdos-Rényi model
# p <- 0.1 # probability of adding a new edge between two nodes in Erdos-Rényi model
# 
# 
# # Generate theoretical networks with known labels
# 
# set.seed(123)
# theoretical_networks <- list(
#   small_world = sample_smallworld(dim = 1, size = n, nei = 2, p = 0.05),
#   scale_free = barabasi.game(n, power=2.1, directed=FALSE, algorithm="psumtree"), 
#   #sample_degseq(degree.sequence.powerlaw(100, exponent = 2.1), method = "vl")
#   spatial = sample_pa(n, power = 1, m = 2, directed = FALSE),
#   erdos_renyi = erdos.renyi.game(n, m, type = "gnp", p = p)
# )
# labels <- names(theoretical_networks)
# 
# # Generate empirical network to classify
# set.seed(456)
# empirical_network <-barabasi.game(n, power=1.2, directed=FALSE, algorithm="psumtree")
# #sample_degseq(degree.sequence.powerlaw(100, exponent = 2.1), method = "vl")
# 
# # Extract features from theoretical and empirical networks
# theoretical_features <- sapply(theoretical_networks, function(g) {
#   c(
#     mean(degree(g)),
#     diameter(g),
#     mean_distance(g),
#     transitivity(g)
#   )
# })
# colnames(theoretical_features) <- c("mean_degree", "diameter", "mean_distance", "transitivity")
# 
# emp_data=matrix(0,nrow=1,ncol = 4)
# empirical_features <- c(
#   mean(degree(empirical_network)),
#   diameter(empirical_network),
#   mean_distance(empirical_network),
#   transitivity(empirical_network)
# )
# empirical_features=data.frame(matrix(empirical_features,nrow=1,ncol = 4))
# colnames(empirical_features) <- colnames(theoretical_features)
# # Train random forest classifier to predict labels
# data <- rbind(
#   data.frame(label = labels, theoretical_features),
#   data.frame(label = "empirical", empirical_features)
# )
# 
# 
# data1=data[-nrow(data),]
# data1$label=as.factor(data1$label)
# 
# rf <- randomForest(label~., data = data1, importance = TRUE)
# 
# # Predict label for empirical network
# predicted_label <- predict(rf, newdata = data[nrow(data),])[1]
# predicted_label
