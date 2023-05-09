
#-------------------------------------------
# Calculating network properties and proportion infected from multiple graph simulations
#-------------------------------------------
#code by Nick Fountain-Jones May/June 2021

#library(igraph)
#library(EpiModel)
#library(intergraph)
#library(DiagrammeR)
#library(tidyverse)
#library(doParallel)
#library(foreach)
#source('simPathE.R')
#source('countStates.R ') 

#this function calculates prevalence ofan SIR pathogen network as well as the corresponding network characteristics
Network_sum <- function(folder_name){
  
  networks <- list.files(paste(folder_name))
  
  #Global_summary <- NULL
  
  #make sure all network object in the folder can load.
  
  netcheck <- function(networks ){ 
                        filename=paste(folder_name, networks, sep="/") #
                        
                        dat <- read.table(filename) 
                        options(show.error.messages = TRUE)
                        print(filename)
                        gtry <- try(igraph::graph_from_data_frame(dat) %>%  as.undirected( "collapse"))
                        if(is( gtry,"try-error")) {print('graph error')}
                        else { print('graph ok')}
                        }

 #g <- netcheck(networks[540])
  checkGraphTable <- networks %>%  map( netcheck)
  
  #make multicore
  cl <- makePSOCKcluster(5)
  registerDoParallel(cl)
  
  # removing mammalia-voles-kcs-trapping-01 as there is a strange eror with this one
  
  #will thow an eror if there is anything wrong with any graph in the folder
  }  
  #for (i in 1:length(networks))##
  Global_summary  <- foreach( i = 1:length(networks), combine=rbind,  .packages=c('tidyverse', 'igraph','intergraph', 'DiagrammeR')) %dopar% 
 
  {
    source('simPathE.R')
    source('countStates.R ') 
    
    print(networks[i])
    filename=paste(folder_name, networks[i], sep="/") #
    
    dat <- read.table(filename) #haven't dealt with the temporal nature of this yet. Third column appears to be a date...
    g <- igraph::graph_from_data_frame(dat) %>%  as.undirected( "collapse")
   
  
    #g <- read_graph(networks[4], format = c("edgelist")) %>% as.undirected( "collapse") #seem only to work on smaller graphs?
    #one edge for each pair of connect vertices
    g_names <- gsub(".edges","",networks[i]) ##
    
    netclust <- clusters(g) #look for subgraphs
    gcc <- V(g)[netclust$membership == which.max(netclust$csize)]#select vertices from the largest sub-graph
    #we could do this for other component graphs?
    gccnet <- induced.subgraph(g, gcc) #make it a igraph object.
    
    #need to remove loops etc
    gcc_s <- igraph::simplify(gccnet, remove.multiple = TRUE, remove.loops = TRUE)
    #plot(gcc_s,main=networks[1])##
    
    #-------------------------------------------
    # simulate SIR model over the network 100 times (for 100 time steps)
    #-------------------------------------------
    
    #beta =0.05 tests
    
   nsim = 100 #10 to start with
    
    inf_data_beta0.05 <- list()
    
       for(j in(1:nsim)){
      
          sm05 <- data.frame(simPathE(gcc_s, nTicks=100, beta=0.05, gamma=0.04, propInfected=0.1, initialState=NULL,nInfected=1, useProportion=T)) #r0 of 1.2 roughly first wave SARScoV2
    
           countEachState <- countStates(DF=sm05, states=c(0:2))  
    
           TwoPcInfected <-  sum(sm05[1,])+ceiling((vcount(gcc_s)/100)*2) #Sah et al use 2% threshold for pathogen invasion
    
          invasion_time <- which(countEachState[,2] >= TwoPcInfected )[1] #this selects the first time the number of infects
    # 2% of the population 
    #Need to exclude the first infection we simulate (hence the +1)
    
          prop_infected <- max(countEachState[,2])/vcount(gcc_s)
    
          inf_data_beta0.05[[j]] <-  c(invasion_time=invasion_time, prop_infected=prop_infected)
    
       }
    
    #collect data
    complete_data_beta0.05 <- as.data.frame(do.call(rbind,  inf_data_beta0.05))
    complete_data_noNA_beta0.05 <-   complete_data_beta0.05[complete.cases(complete_data_beta0.05 ),] #NA mean it never took off
    
    no_invasion_beta0.05 <- sum(is.na(complete_data_beta0.05$invasion_time)) #number of times there was no outbreak (>2% infected is interesting)
    
    avg_prop_infected_beta0.05 <- mean(complete_data_noNA_beta0.05$prop_infected)
    avg_invasion_time_beta0.05 <- mean(complete_data_noNA_beta0.05$invasion_time)
    
#------------------------------------------------------------------------------------------------------------
    #beta =0.1 tests (more infectious)
  
    
    inf_data_beta0.1 <- list()

       for(j in(1:nsim)){

        sm01 <- simPathE(gcc_s, nTicks=100, beta=0.1, gamma=0.4, propInfected=0.05, initialState=NULL,nInfected=1,  useProportion=T)
        countEachState <- countStates(DF=sm01, states=c(0:2))

         TwoPcInfected <-  sum(sm01[1,])+ceiling((vcount(gcc_s)/100)*2) #Sah et al use 2% threshold for pathogen invasion

        invasion_time <- which(countEachState[,2] >= TwoPcInfected )[1] #this selects the first time the number of infects
       # 2% of the population
        #Need to exclude the first infection we simulate (hence the +1)

        prop_infected <- max(countEachState[,2])/vcount(gcc_s)

        inf_data_beta0.1[[j]] <-  c(invasion_time=invasion_time, prop_infected=prop_infected)

    }

    #collect data
    complete_data_beta0.1 <- as.data.frame(do.call(rbind,  inf_data_beta0.1))
    complete_data_noNA_beta0.1 <-   complete_data_beta0.1[complete.cases(complete_data_beta0.1 ),] #NA mean it never took off
    
    no_invasion_beta0.1 <- sum(is.na(complete_data_beta0.1$invasion_time)) #number of times there was no outbreak (>2% infected is interesting)
    
    avg_prop_infected_beta0.1 <- mean(complete_data_noNA_beta0.1$prop_infected)
    avg_invasion_time_beta0.1 <- mean(complete_data_noNA_beta0.1$invasion_time)
      
  
    
    #------------------------------------------
    
    #network characteristics
    
    #convert to a laplacian matrix
    M <- laplacian_matrix(gcc_s,sparse = FALSE)
    
    #Adjacency matrix
    Adj <-as.matrix(gcc_s[])
    #Adjacency spectrum
    Adjacency_spectrum<-eigen(Adj)
    
    #calc eigenvectors
    spec <- eigen(M)
    #take the second smallest eigenvalue
    FiedlerValue<- spec$values[length(spec$value)-1]
    
    #Largest value of the adjacency
    ##This value controls the propagation of infection in SIS/SIR model
    ### The larger Adj_val the faster the growth of infection at the beginning time
    Adj_val<-Adjacency_spectrum$values[1]
    Eigen_central=Adjacency_spectrum$vectors[,ncol(Adjacency_spectrum$vectors)[1]]
    most_infected_node<-match(max(Eigen_central),Eigen_central)# highest node degree
    
    #Note:Probability of infection of nodes depends on V_1(first eigen vector of adjacency matrix)
    ## so the prob that a node is infected is proportional to its eigen vector centrality (V_1)
    ## since centrality depends on number of nodes its connected to and how well node are connected
    
    
    #Epidemic threshold
    ##Infection dies or survives if beta/gamma is < or > Rnot respectively
    Rnot<-1/Adj_val # spectral of adjacency for undirected network controls the epidemic threshold
    

    #Threshold comparism of 
    
    #Newman's relative Q (from Sah et al 2017) modularity
    
    #detect commnuities and calc relative Newman's Q
    wc <- cluster_walktrap(gcc_s) # This function tries to find densely connected subgraphs, 
    lc <- cluster_louvain(gcc_s) #this is the method used inSah et al 2017 - allows comparison? 
    #also called communities in a graph via random walks.
    
    mod <- modularity(gcc_s, membership(lc))
    #max_mod <- 
    #Qrel <- mod/max_mod #Sah et al normalize modularity this way, but Im not sure why maximum modulairity is the best way to do this
    #plot(lc,gccnet)  check out the plot if needed
    
    #other network characteristics
    deg <- mean(degree(gcc_s, mode="all", normalized=TRUE)) #normalized currently
    
    cent <-  eigen_centrality(gcc_s, directed=F, weights=NA)
    cent <- cent$value
    
    trans <- as.data.frame(transitivity(gcc_s, type="global"))
    
    di <- diameter(gcc_s, directed=F, weights=NA)
    network_size<-vcount(gcc_s)
    
    #Global_summary[1] <- bind_cols(g_names, propI2,propI5, FiedlerValue, mod, deg, cent, trans, di) %>%
    #set_names('g_names')
    
    data.frame(Network=g_names,  avg_invasion_time_beta0.1=avg_invasion_time_beta0.1,
                           avg_invasion_time_beta0.05=avg_invasion_time_beta0.05,  avg_prop_infected_beta0.1=avg_prop_infected_beta0.1, 
                            avg_prop_infected_beta0.05=avg_prop_infected_beta0.05,   no_invasion_beta0.1=no_invasion_beta0.1,   
                             no_invasion_beta0.05=no_invasion_beta0.05, Fiedler=  FiedlerValue, Adj_val=Adj_val, Rnot=Rnot,network_size=network_size, 
                             most_infected_node=most_infected_node, Modularity=mod, Mean_degree=deg, Centrality=cent,
                             Transivity=trans, Diameter=di)
  }
    
      #write.table(file="globalData.csv",append=TRUE,sep=",",col.names=TRUE,row.names=TRUE)
  #Global_summary[[i]] 
  
  #Global_summary_list 
    #write.table(Global_summary[,,i],file="globalData.i.csv",append=TRUE,sep=",",col.names=TRUE,row.names=TRUE) 
  
Global_summary <- do.call(rbind, Global_summary)


 # colnames(Global_summary) <- c('Network','Avg_infected_R02','Avg_infected_R05', 'Fiedler','Adj_val','Rnot','network_size','most_infected_node', 'Modularity', 'Mean_degree', 'Centrality', 'Transivity', 'Diameter')
# Global_summary  <- as.data.frame(Global_summary)
  
  return(Global_summary)
  
  parallel::stopCluster(cl)
  

}

#########################################################################################
#PIPELINE
#########################################################################################

#create graph object
folder_name = c('Networks') #.edges file need to be in the working directory too at this stage.

test <- Network_sum(folder_name)

str(test)
saveRDS(test,'graphdf' )


dfNumeric <- test %>% select(-Network)

write.csv(x=test, file='NetworkDataSim.csv')
DataExplorer::create_report(dfNumeric )