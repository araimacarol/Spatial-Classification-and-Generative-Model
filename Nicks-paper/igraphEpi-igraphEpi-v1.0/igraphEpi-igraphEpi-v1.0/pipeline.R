#' Vignette to run the EpicGraph package

#Load packages - keep this to a minimum
library(igraph)
library(EpiModel)
library(intergraph)
library(DiagrammeR)
library(tidyverse)

# load all function codes. This will disappear when we formally make this a function
# <<<<<<< HEAD
# #source("getGraphFeatures.R")
# =======
# source("getGraphFeatures.R")
# >>>>>>> 2efec6358e09a69a9b953184b13cbbb65f805e38


#Nicks experiment

#-----------------------------------------------------------------------------------
#Step 1: Summarise network features
#-----------------------------------------------------------------------------------
#folder_name = c('Networks') #.edges file need to be in the working directory too at this stage.

#NetSum <- getGraphFeatures(network_name = f,  multinetwork = TRUE)

#NetSum


#Raima's code

## sourcing idealised networks
source( "random_network.R")

source( "small-world_network.R")


source( "scale-free_network.R")


source( "spatial_network.R")

source( "lattice_network.R") #tbis isnt working for me




### Epidemic simulation functions
source( "epi_sim.R")

source( "epi-summary.R")

source( "network_summary.R")

source("Experiment-3.R")

source("Epic.R")

source("graph_features.R")



# one simulated random graph


# F1=erdos(100,0.6)
# s=list()
# s[[1]]=F1
# # n simulated random network 
# G_2=Erdos(10,.2,10)
# G=G_2
# 
# # multiple epidemic simulation on all simulated networks of a specific network type
# netsim<-episim(G,beta = 0.5,gamma=0.2,propInfected = 0.1,numInfected = 2,useProportion = T) 
# netsim
# 
# # multiple  epidemic simulation on each simulated networks of a specific network type
# #All_netsim=nets_episim(G)
# #All_netsim
# 
# # single epidemic summary on all network
# netsummary<-single_epi_summary(netsim) # epidemic summary on all simulated networks of a specific network type
# netsummary
# 
# #multiple epidemic summary on the networks
# All_netsummary=nets_epi_summary(netsim)
# All_netsummary
# 
# 
# 
# ###-------Experiment-1-------######
# n=10
# random_graph=erdos(n,.04)
# lattice=lattice_net(ceiling(sqrt(n)), floor(sqrt(n)))
# scalefree_graph=scale_free(n,1,2)
# spatial_graph=graph_from_data_frame(spatial(n,.12),directed=FALSE)
# 
# small_world=sm_world(1,n,2,.3)##needs attention



#gen_models <- list(random_graph, small_world, lattice, scalefree_graph, spatial_graph) 
#names(gen_models) <- c("random_graph","small_world","lattice","scale_free","spatial")


#data_frame_of_network_features=Network_Features(gen_models)
#data_frame_of_network_features
# =======
# F1=erdos(100,0.6)
# s=list()
# s[[1]]=F1
# # n simulated random network 
# G_2=Erdos(10,.2,10)
# G=G_2
# 
# # multiple epidemic simulation on all simulated networks of a specific network type
# netsim<-episim(G,beta = 0.5,gamma=0.2,propInfected = 0.1,numInfected = 2,useProportion = T) 
# netsim
# 
# # multiple  epidemic simulation on each simulated networks of a specific network type
# #All_netsim=nets_episim(G)
# #All_netsim
# 
# # single epidemic summary on all network
# netsummary<-single_epi_summary(netsim) # epidemic summary on all simulated networks of a specific network type
# netsummary
# 
# #multiple epidemic summary on the networks
# All_netsummary=nets_epi_summary(netsim)
# All_netsummary
# 
# 
# 
# ###-------Experiment-3-------######
# n=10
# random_graph=erdos(n,.04)
# 
# lat=lattice(ceiling(sqrt(n)), floor(sqrt(n))) #neede lattice? instead of lattice_net
# 
# scalefree_graph=scale_free(n,1) #Error in scale_free(n, 1, 2) : unused argument (2) #removed 2 and it worked?
# 
# spatial_graph=graph_from_data_frame(spatial(n,.12),directed=FALSE)
# 
# small_world=sm_world(1,n,2,.3)##needs attention #NFJ: What's wrong? Worked OK I think
# 

# 
# gen_models <- list(random_graph, small_world, lat, scalefree_graph, spatial_graph) 
# names(gen_models) <- c("random_graph","small_world","lattice","scale_free","spatial")


# data_frame_of_network_features=Network_Features(gen_models) # NFJ: this error: Error in vcount(g[[i]]) : Not a graph object 
# data_frame_of_network_features
# >>>>>>> 85c38d1944351d7baeb0339d358092860b4a6f91
# 



###-----Experiment 2 on random network-------
set.seed(25)
RandomGraphs=Erdos(20,.07,2)
number_of_replicate=10
GLOBAL_NETWORK_SUMMARY=Epic(Graphs = RandomGraphs,nrep = number_of_replicate, beta = 0.7,gamma = .4)
GLOBAL_NETWORK_SUMMARY

##-----------simulations for different network sizes--------

###---Small network size---####
# n=16
# r1=erdos(n,.2)
# s1=sm_world(2, sqrt(n), 5, 0.05)
# l1=Lattice(2,sqrt(n))
# sc1=scale_free(n,1)
# 
# epidemic_simuations=lapply(G,episim)
# 
# for (i in (1:length(epidemic_simuations))){
#   v[[i]]=apply(epidemic_simuations[[i]]==1, 1, sum)
#   plot(v[[i]], col=sample(rainbow(10)), type='l', 
#        main="",sub="",
#        xlab="Time",ylab="infected ")
#   par(new=TRUE)
# }
# 
# 
# 
# ###---medium network size---####
# r2=erdos(100,.2)
# s2=sm_world(1, 100, 5, 0.05)
# l2=Lattice(2,100)
# sc2=scale_free(100,1)
# 
# #episim(sc2,numInfected = 20)
# epidemic_simuations_2=lapply(sc2,episim)
# 
# for (i in (1:length(epidemic_simuations_2))){
#   v[[i]]=apply(epidemic_simuations_2[[i]]==1, 1, sum)
#   plot(v[[i]], col=sample(rainbow(10)), type='l', 
#        main="",sub="",
#        xlab="Time",ylab=" ",ylim=c(0,ncol(epidemic_simuations_2[[i]])))
#   par(new=TRUE)
# }
# 
# ###---large network size---####
# r3=erdos(1000,.2)
# s3=sm_world(1, 1000, 5, 0.05)
# l3=Lattice(2,1000)
# sc3=scale_free(1000,1)
# 
# epidemic_simuations_3=lapply(G,episim(G,numInfected = 200));
# 
# for (i in (1:length(epidemic_simuations_3))){
#   v[[i]]=apply(epidemic_simuations_3[[i]]==1, 1, sum)
#   plot(v[[i]], col=sample(rainbow(10)), type='l', 
#        main="",sub="",
#        xlab="Time",ylab=" ",ylim=c(0,ncol(epidemic_simuations_3[[i]])))
#   par(new=TRUE)
# }
# 
# 
# 


