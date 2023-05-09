library(igraph)
library(EpiModel)
library(intergraph)
library(DiagrammeR)
library(tidyverse)
library("dplyr")
library("purrr")



source( "epi_sim.R")

source( "epi-summary.R")

source( "network_summary.R")


Network_Features <- function(models,numb_simult){
  

  network_summary <- NULL
  
  g=list()
  network_size=list()
  g_names=list()
  avg_deg=list()
  
  for (i in 1:length(models))
  {
    g[[i]] <- models[[i]]
    
    g_names[[i]] <- names(models[i])
  
  
    network_size[[i]]<-vcount(g[[i]])
    
    avg_deg[[i]] <- mean(degree(g[[i]], mode="all"))


    network_summary[[i]] <- c(g_names[[i]], network_size[[i]], avg_deg[[i]]) 
    
  }
  
  
  network_summary <- do.call(rbind, network_summary)
  network_summary  <- as.data.frame(network_summary)
  
  colnames(network_summary) <- c('Network_type', 'Network_size','Average_degree')

  
  
  return(network_summary)
}





###---- Random network example----------
##--- Experiment with rando network

random_net_features <- function(rand_network,numb_simult){# input random network and number of simultns
  
  
  rand_network_summary <- NULL
  
  g=list()
  network_size=list()
  sim_number=list()
  avg_deg=list()
  num_infected=list()
  sim_number=list()
  
  for (i in 1:length(rand_network))
  {
    g[[i]] <- rand_network[[i]]#assigning each random network
    

    for (j in 1:numb_simult){
      sim_number[[j]] <- j
      network_size[[j]]<-vcount(g[[i]])
      avg_deg[[j]] <- mean(degree(g[[i]], mode="all"))
      num_infected[[j]]= single_epi_summary(episim(g))[[i]][][colnames(single_epi_summary(episim(g))[[i]][])=="I"] %>% select_if(is.numeric) %>% map_dbl(sum)
    }
    
    rand_network_summary[[i]] <- c( sim_number, network_size, avg_deg, num_infected) #summary of the network
    
  } 
  
   rand_network_summary <- do.call(rbind,  rand_network_summary)
   rand_network_summary  <- as.data.frame( rand_network_summary)
   
 #  colnames( rand_network_summary) <- c('simulation_number', 'Network_size','Average_degree','number_of_infected')
   
  return(rand_network_summary)
}







# x=All_netsummary[[10]][[4]]
# h[[1]][][colnames(h[[1]][])=="I"] %>%
#   select_if(is.numeric) %>%
#   map_dbl(sum)



# x=All_netsummary[[10]][[4]]
# x[colnames(x)=="I"] %>%
#   select_if(is.numeric) %>%
#   map_dbl(sum)






# g
# names(g) <- c("random_graph","small_world","lattice","scale_free","spatial")
# numb_simult=5
# g1= NULL
# for (j in 1:numb_simult){
#   if(names(g)[[i]]=="random_graph"){
#     g1[[j]]=erdos(4,.5)}
#   
#   if(names(g)[[i]]=="small_world"){
#     g1[[j]]=sm_world(1, 4, 2, 0.05)}
# }
