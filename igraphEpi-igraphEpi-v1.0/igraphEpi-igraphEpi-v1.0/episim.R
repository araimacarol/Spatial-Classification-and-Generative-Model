##

getGraphFeatures <- function() {
  ## output: data frame
  
}

# converting downloaded networks into igraph objects
# use something liek this to convert:
#dat <- read.table(networks[i]) #haven't dealt with the temporal nature of this yet. Third column appears to be a date...
#g <- graph_from_data_frame(dat)%>% as.undirected( "collapse")

G <- graph_from_data_frame(dat)%>% as.undirected( "collapse")

# status: 0=healthy/susceptible, 1=infected, 2=recovered, 3=removed/dead/immune
episim <- function(G, nticks=100, beta=0.1, gamma=0.2, sigma=0, initState=NULL, propInfected=0.1) {
  ## G is igraph object
  ## extract edges and vertices of G
  # initState is list of initial state of each vertex in G.  E.g., from (1,0,0...0) or (0,1,0,1,0,....)
  E <- G[,] # from igraph -- wow.
<<<<<<< HEAD
  #  E <- as_edgelist(G)
  #print(E)
=======
#  E <- as_edgelist(G)
#  print(E)
>>>>>>> 2efec6358e09a69a9b953184b13cbbb65f805e38
  n <- ncol(E)
#  cat("n = ", n)
#  cat(sprintf("n=%d\n",n))
  if (is.null(initState)) {
    stat <- rbinom(n,1,propInfected)
  } else {
    stat <- initState
  }
  print(stat)
  nextStat <- rep(0,length(stat))
  for (tick in 1:nticks) {
    # first phase: S -> I
    for (i in which(stat==0)) {
      for (j in 1:ncol(E)) { # ncol(E) is the number of nodes (do we need to rexplude index i (with [-i] ?)
        if ((E[i,j]>=1) & (runif(1) < beta)) {
          nextStat[j] <- 1; break # don't need to check any more neighbours!
        }
      }
    }
    
    # can add in E state later: S->E->I->R
    
    # second phase: I -> R
    for (i in which(stat==1)) {
      if (runif(1) < gamma) {
        nextStat[j] <- 2
      }
    }
    
    # third phase: R -> S
    for (i in which(stat==2)) {
      if (runif(1) < sigma) {
        nextStat[j] <- 0
      }
    }
    
    stat <- nextStat
<<<<<<< HEAD
 #   print(length(which(stat==1)))
=======
#    print(length(which(stat==1)))
>>>>>>> 2efec6358e09a69a9b953184b13cbbb65f805e38
    print(stat)
  }
}

G <- erdos.renyi.game(20, .10,direct=F)    #random graph generated

episim(G)
 
