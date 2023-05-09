## [1] simulating a single SIR process on each lists of igraph objects


library(igraph)

#Simple SIR disease simulation
# Disease status: 0=healthy/susceptible, 1=infected, 2=removed/dead/recovered

episim <- function(Graphs=NULL, nticks=10, beta=0.1, gamma=0.2, propInfected=0.1, initState=NULL, numInfected=1, useProportion=F,nrep=NULL){
  
 # graph_of_one=list() 
 # 
 # if (length(Graphs)==1){
 #   graph_of_one[[1]]=Graphs
 #    Graphs= graph_of_one
 # }
 
  # else if (length(Graphs)!=1){
  #    for (i in 1:length(Graphs)){
  #         Graphs[[i]]=Graphs[[i]]
  #  }
  # }
  # 
  
  # E=list()
  # n=list()
  # inf_per_time=list()
  # infected=list()
  # clear=list()
  # stat=list()
  # nextStat=list()
  # 

E=vector(mode = "list", length(Graphs))
n=vector(mode = "list", length(Graphs))
inf_per_time=vector(mode = "list", length(Graphs))
infected=vector(mode = "list", length(Graphs))
clear=vector(mode = "list", length(Graphs))
stat=vector(mode = "list", length(Graphs))
nextStat=vector(mode = "list", length(Graphs))

for (idx in 1:length(Graphs)){
  for(rep in 1:nrep){  
    E[[idx]][[rep]]<- as.matrix(Graphs[[idx]][])
    n[[idx]][[rep]]<- ncol(E[[idx]][[rep]]) #number of columns (nodes) in the list of network object
    inf_per_time[[idx]][[rep]]<-matrix(data=0,nticks+1,n[[idx]][[rep]])#list of matrix related to the igraph objects with progression of infection/transition of state of each node for each time step

 
# N=3
# y=NULL
# for(i in 1:N){
#   y[[i]]=as.matrix(G[[2]][])
# }
#  result=list()
 # graph_idx=1
  
#  for (G in Graphs){
    ## G is an igraph object
#    E <- get.adjacency(G,type = c("both")) #adjacency matrix of graph (network object)
 #   n <- ncol(E) #number of columns (nodes) in network object
#    inf_per_time<-matrix(data=0,nticks+1,n)#progression of infection/transition of state of each node for each time step
    
    # Set initial state  
    # initState is list of initial state of each vertex in G.  E.g., from (1,0,0...0) or (0,1,0,1,0,....)
  
  
  
  

      if (is.null(initState[[idx]])) {
        if (useProportion==T) {
          stat[[idx ]][[rep]] <- rbinom(n[[idx]][[rep]],1,propInfected)#set initial state of each nodes if non (note, no recovered node at initial state)
        } else {
          infected[[idx]][[rep]]<- rep(1,numInfected)
          clear[[idx]][[rep]]<-rep(0,(n[[idx]][[rep]]-numInfected))
          stat[[idx]][[rep]]<-sample(c(infected[[idx]][[rep]], clear[[idx]][[rep]]), n[[idx]][[rep]], replace=FALSE)
        }
  
       } else {
         stat[[idx]][[rep]]<- initState[[idx]][[rep]] #initial existing state.
       }
    #
  
    
       current_stat=stat
      
       nextStat[[idx ]][[rep]]<-rep(0,length(stat[[idx]][[rep]]))
       time=2 # since we have nticks+1
      
    #note: there is no recovered node at initial state only susceptible nodes and some infected nodes
    #initial state vector must not contain 2, but only 0 and 1 for the first row
    #set next state vector initially as initial state.
      for (tick in 1:nticks) {
        inf_per_time[[idx]][[rep]][1,]=current_stat[[idx]][[rep]]#time 0=initial state
        # first transition phase: S -> I
        for (i in which(stat[[idx]][[rep]]==0)) {#for all susceptible nodes (denoted 0) in initial state
          for (j in (1:n[[idx]][[rep]])) { # for all other neighbor nodes
            if ((E[[idx]][[rep]][i,j]>=1)&(stat[[idx]][[rep]][j]==1)&(runif(1)<=beta)) #each infected node infects its nearest neighbor/node it shares a link/edge with, with prob beta
            {
              nextStat[[idx]][[rep]][i] <- 1;break#assign node as infected if above condition is met else and break out of loop
            }
            else{
              nextStat[[idx]][[rep]][i]<-0 #node remains as susceptible, since not all contact leads to an infection.
            }
          }
        }
        
        # second transition phase: I -> R
        for (i in which(stat[[idx]][[rep]]==1)) {#for all infected nodes in initial state
          if (runif(1)<=gamma) { #compares a randomly generated uniform number to recovery rate
            nextStat[[idx]][[rep]][i] <- 2 #assigns node as recovered if condition above is met else
          }
          else{
            nextStat[[idx]][[rep]][i]<-1 #nodes remain infected 
          }
        }
      
      #for SIR, recovered node do not participate in further disease propagation (permanent immunity)    
      
        stat[[idx]][[rep]]<- nextStat[[idx]][[rep]]#assigns states of nodes for each time step as existing state 
        
        inf_per_time[[idx]][[rep]][time,] <- stat[[idx]][[rep]]
        rownames(inf_per_time[[idx]][[rep]])=0:nticks#re-order row names from time step 0=initial state
        colnames(inf_per_time[[idx]][[rep]])=NULL#make column names as null
        
        
          time=time+1
    }  
     }
    
   } 
  return(inf_per_time)
}


## [2] simulating multiple SIR process on each lists of igraph objects

# episim_all_graphs=NULL
# N=10
# nets_episim<-function(G){
#   for (i in 1:N){
#     episim_all_graphs[[i]]=episim(G)
#   }
#   return(episim_all_graphs)
#  }


# G=f1
# for (i in 1:length(G)){
#   if (is.igraph(G[[i]])){
#     G[[i]]=as.list(G[[i]])
#   }
#   return(G)
# }
# 
# E=list()
# n=list()
# for (k in (1:length(G))){
# E[[k]] <- as.matrix(get.adjacency(G[[k]]))
# n[[k]] <- ncol(E[[k]]) #number of columns (nodes) in network object
# }

# M=function(useprop=NULL){
#   if(useprop==T){
#     print("yeah")
#   }else{
#     print("naa")
#   }
# }