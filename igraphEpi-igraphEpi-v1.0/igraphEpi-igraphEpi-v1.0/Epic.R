library(igraph)
library(EpiModel)
library(intergraph)
library(DiagrammeR)
library(tidyverse)

## [1] simulating a single SIR process on each lists of igraph objects




library(igraph)
#Simple SIR disease simulation
# Disease status: 0=healthy/susceptible, 1=infected, 2=removed/dead/recovered

Epic<-function(Graphs=NULL, nticks=10, beta=0.1, gamma=0.2, propInfected=0.1, initState=NULL, numInfected=1, useProportion=F,nrep=NULL){
      
      Graph_features=GraphFeatures(Graphs,nrep)#Graph-features
      Epic_Sim=episim(Graphs=Graphs, nticks=nticks, beta=beta, gamma=gamma, propInfected=propInfected, initState=initState, numInfected=numInfected, useProportion=useProportion,nrep=nrep)#epidemic network simulation
      Epic_summary<-nets_epi_summary(Epic_Sim)#epidemic network summary
      
      Max_epid_size=vector(mode="list", length(Epic_summary))#peak/maximum size of infection/epidemic size
      time_Max_epid_size=vector(mode="list", length(Epic_summary))
      ##Calculate all other stuffs on inf_per_time   
  Global_network_summary=NULL 
  nGraphs=length(Epic_summary)
  
      for (idx in 1:nGraphs) { 
        for (sim_num in 1:nrep){
          Max_epid_size[[idx]][[sim_num]]<-c(max(Epic_summary[[idx]][[sim_num]]$I))
          time_Max_epid_size[[idx]][[sim_num]]<-which.max(Epic_summary[[idx]][[sim_num]]$I)-1#timestep start from 0
        
          Global_network_summary[[idx]]=as.data.frame(cbind(time_Max_epid_size[[idx]],Max_epid_size[[idx]]))
        
        
          
        colnames(Global_network_summary[[idx]])=c('Time_to_epidemic_size','Epidemic_size')
        
        
        
        } 
        Global_network_summary[[idx]]$"numb_of_rep"= 1:nrep
      }
      
      Global_network_summary=do.call(rbind,Global_network_summary)
  #    column_order=c('numb_graph','numb_of_rep','Time_to_epidemic_size','Epidemic_size')
  #    Global_network_summary=Global_network_summary[, column_order]
      Global_network_summary=cbind(Global_network_summary,Graph_features)
      graphNumbers=rep(1:nGraphs,each=nrep)
      Global_network_summary=cbind(Global_network_summary,graphNumbers)
      
return(Global_network_summary)
}
 









   # Max_epid_size=do.call(rbind,Max_epid_size)
    #  time_Max_epid_size=do.call(rbind,time_Max_epid_size)
    #    colnames(Max_epid_size[[idx]][[rep]])=c('Epidemic_size')
    #    colnames(time_Max_epid_size[[idx]][[rep]])=c('Time_to_epidemic_size')
    
# threshold_value=2
    
    
  #   Global_network_summary=cbind(rep_number,Max_epid_size,time_Max_epid_size,Graph_features)
  #   #graph xtics/rowise
  #   #append row to result
  # 
  #   
  # } 

  #return entire summary per simulation on a network/result
 
#}

# 
# Epic<- function(Graphs, nticks=10, beta=0.1, gamma=0.2, propInfected=0.1, initState=NULL, numInfected=1, useProportion=F,nrep=NULL){
#   
#   E=vector(mode = "list", length(Graphs))
#   n=vector(mode = "list", length(Graphs))
#   inf_per_time=vector(mode = "list", length(Graphs))
#   infected=vector(mode = "list", length(Graphs))
#   clear=vector(mode = "list", length(Graphs))
#   stat=vector(mode = "list", length(Graphs))
#   nextStat=vector(mode = "list", length(Graphs))
#   
#   for (idx in 1:length(Graphs)){
#     for(rep in 1:nrep){  
#       E[[idx]][[rep]]<- as.matrix(Graphs[[idx]][])
#       n[[idx]][[rep]]<- ncol(E[[idx]][[rep]]) #number of columns (nodes) in the list of network object
#       inf_per_time[[idx]][[rep]]<-matrix(data=0,nticks+1,n[[idx]][[rep]])#list of matrix related to the igraph objects with progression of infection/transition of state of each node for each time step
#       
#       
#       # Set initial state  
#       # initState is list of initial state of each vertex in G.  E.g., from (1,0,0...0) or (0,1,0,1,0,....)
#       
#       
#       
#       if (is.null(initState[[idx]])) {
#         if (useProportion==T) {
#           stat[[idx ]][[rep]] <- rbinom(n[[idx]][[rep]],1,propInfected)#set initial state of each nodes if non (note, no recovered node at initial state)
#         } else {
#           infected[[idx]][[rep]]<- rep(1,numInfected)
#           clear[[idx]][[rep]]<-rep(0,(n[[idx]][[rep]]-numInfected))
#           stat[[idx]][[rep]]<-sample(c(infected[[idx]][[rep]], clear[[idx]][[rep]]), n[[idx]][[rep]], replace=FALSE)
#         }
#         
#       } else {
#         stat[[idx]][[rep]]<- initState[[idx]][[rep]] #initial existing state.
#       }
#       #
#       
#       
#       current_stat=stat
#       
#       nextStat[[idx ]][[rep]]<-rep(0,length(stat[[idx]][[rep]]))
#       time=2 # since we have nticks+1
#       
#       #note: there is no recovered node at initial state only susceptible nodes and some infected nodes
#       #initial state vector must not contain 2, but only 0 and 1 for the first row
#       #set next state vector initially as initial state.
#       for (tick in 1:nticks) {
#         inf_per_time[[idx]][[rep]][1,]=current_stat[[idx]][[rep]]#time 0=initial state
#         # first transition phase: S -> I
#         for (i in which(stat[[idx]][[rep]]==0)) {#for all susceptible nodes (denoted 0) in initial state
#           for (j in (1:n[[idx]][[rep]])) { # for all other neighbor nodes
#             if ((E[[idx]][[rep]][i,j]>=1)&(stat[[idx]][[rep]][j]==1)&(runif(1)<=beta)) #each infected node infects its nearest neighbor/node it shares a link/edge with, with prob beta
#             {
#               nextStat[[idx]][[rep]][i] <- 1;break#assign node as infected if above condition is met else and break out of loop
#             }
#             else{
#               nextStat[[idx]][[rep]][i]<-0 #node remains as susceptible, since not all contact leads to an infection.
#             }
#           }
#         }
#         
#         # second transition phase: I -> R
#         for (i in which(stat[[idx]][[rep]]==1)) {#for all infected nodes in initial state
#           if (runif(1)<=gamma) { #compares a randomly generated uniform number to recovery rate
#             nextStat[[idx]][[rep]][i] <- 2 #assigns node as recovered if condition above is met else
#           }
#           else{
#             nextStat[[idx]][[rep]][i]<-1 #nodes remain infected 
#           }
#         }
#         
#         #for SIR, recovered node do not participate in further disease propagation (permanent immunity)    
#         
#         stat[[idx]][[rep]]<- nextStat[[idx]][[rep]]#assigns states of nodes for each time step as existing state 
#         
#         inf_per_time[[idx]][[rep]][time,] <- stat[[idx]][[rep]]
#         rownames(inf_per_time[[idx]][[rep]])=0:nticks#re-order row names from time step 0=initial state
#         colnames(inf_per_time[[idx]][[rep]])=NULL#make column names as null
#         
#         
#         time=time+1
#         
#       }
#     }  
    