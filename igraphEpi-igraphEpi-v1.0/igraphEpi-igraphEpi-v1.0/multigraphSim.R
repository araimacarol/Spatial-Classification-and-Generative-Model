#'multigraphSim:  Wrapper to run episim on a collection of networks
#.
#'@param folder_name A \code{list} folder within the current working directory where the network objects are located

#'@details This function extracts a series of network properties from a list of graph objects that are suitable to be added 
#'to, for example, predictive models. The network properties include the Fiedler value','Adj_val','Rnot',network size,
#''most_infected_node', modularity, Mean degree, Centrality, Transivity, Diameter). output is a tidy dataframe.
#'#Raima need to add some details about the ones you included. 
#'
#'If the user wants to extract features from a single tree (multinetwork = FALSE) the network name is a single network object.
#'If the user wants to extract from multiple networks objects (multinetwork = TRUE) the network name is folder stroing the network objects

#' 
#' @example 
#' folder_name = c('Networks') #.edges file need to be in the working directory too at this stage.

#'  NetSimSum <- multigraphSim(folder_name)

#'  NetSimSum


multigraphSim <- function(network_name, nticks=10, beta=0.1, gamma=0.2, initState=NULL, propInfected=0.5){
  
    networks <- list.files(paste(folder_name))
    
    #Global_summary <- NULL
    
    net_num <- length(networks)
    
      netEpiStats <- lapply(seq(1, net_num ), function(i){
      
    #for (i in 1:length(networks)){
      
      dat <- read.table(networks[1]) #haven't dealt with the temporal nature of this yet. Third column appears to be a date...
      g <- graph_from_data_frame(dat)%>% as.undirected( "collapse")
      
      #g <- read_graph(networks[4], format = c("edgelist")) %>% as.undirected( "collapse") #seem only to work on smaller graphs?
      #one edge for each pair of connect vertices
      g_names <- gsub(".edges","",networks[1])
      
      netclust <- clusters(g) #look for subgraphs
      gcc <- V(g)[netclust$membership == which.max(netclust$csize)]#select vertices from the largest sub-graph
      #we could do this for other component graphs?
      gccnet <- induced.subgraph(g, gcc) #make it a igraph object.
      
      #need to remove loops etc
      gcc_s <- igraph::simplify(gccnet, remove.multiple = TRUE, remove.loops = TRUE)
      
       sm2 <- episim(gcc_s, nticks=nticks, beta=beta, gamma=gamma, initState=initState, propInfected=propInfected)
      
      # time_bins(sm2, middle = TRUE) #sets the time bins to calculate stats
      # time_bins(sm5, middle = TRUE)
       
       med2 <- as.data.frame(median(sm2))
       #med5 <- median(sm5)
       #med5$NI[is.na(med5$NI)] <- 0
       #make sure NAs are removed
       #plot(sm2)
       #plot(sm5)
       propI2 <- as.data.frame((mean(med2$NI))/length(V(gcc_s)))
       
       Global_summary<- c(g_names, propI2)
       #propI5 <- as.data.frame((mean(med5$NI))/length(V(gcc_s))) 
      })
      
}
      