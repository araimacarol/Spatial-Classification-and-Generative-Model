###############################################################################################################################
# simPathE ("sympathy") SIMulates PATHogens in an Epidemic (on a graph)
###############################################################################################################################
simPathE <- function(G, nTicks=100, beta=0.5, gamma=0.2, propInfected=0.1, initialState=NULL, nInfected=1, useProportion=F) {
  nVertices = gorder(G)
  infectionState = as.data.frame(matrix(0, ncol = nVertices, nrow = nTicks+1))
  
  # Set up the initial state: all vertices are either exposed (state = 0) or infected (state = 1); none can be recovered yet (2).
  if (is.null(initialState)) {
    if (useProportion==T) {
      infectionState[1,] <- rbinom(nVertices, 1, propInfected) # set initial state of each nodes if known (note, no recovered node at initial state)
    } else {
      infected <- rep(1, nInfected) # just create a vector of the right number of 1s
      exposed <- rep(0, (nVertices - nInfected))
      infectionState[1,] <- sample(c(infected, exposed), nVertices, replace=FALSE)
    }
  } else {
    if (length(initialState) != nVertices) {
      return ("Initial state and order of Graph (number of vertices) are incompatible.  Check the sizes of your input.")
    }
    infectionState <- initialState # initial existing state.
  }
  
  adjacencyList <- as_adj_list(G) # this is a list of which vertices are adjacent to each vertex i
  
  # Now do the simulation through time:
  for (t in 1:nTicks) {
    # FIRST phase: transmission: S -> I
    for (i in which(infectionState[t,] == 0)) { # for all susceptible nodes (denoted 0) in previous time step
      infectionState[t+1,i] <- 0 # node remains as susceptible, since not all contact leads to an infection.
      for (j in adjacencyList[[i]]) { # for all neighbours of i
        if (infectionState[t,j][1] == 1) { # vertex j is infectious
          if ((runif(1)*G[i,j]) <= beta) { # ... and passes it on!
            infectionState[t+1,i] <- 1;
            break # assign node as infected if above condition is met, and break out of loop: we don't need
            # to check any more adjacent vertices.
          }
        }
      }
    }
    
    # SECOND phase: recovery: I -> R
    for (i in which(infectionState[t,] == 1)) { # for all infected nodes (denoted 1) in previous time step
      if (runif(1) <= gamma) { # compares a randomly generated uniform number to recovery rate
        infectionState[t+1,i] <- 2 # node is recovered
      } else {
        infectionState[t+1,i] <- 1 # node remains infected
      }
    }
    
    # THIRD phase: recovered stays recovered:
    for (i in which(infectionState[t,] == 2)) { # for all recovered nodes (denoted 2) in previous time step
      infectionState[t+1,i] <- 2 # node stays recovered
    }
  }
  rownames(infectionState) = 0:nTicks
  return(infectionState)
}