### Making epidemic summary on one network

# make state counts
# data frame for counts of state
# set number of rows and columns
# column names to SEIR
# epi_summary<-function(f,nticks,nrep){ 
#   #f=matrix of simulated sir process
#   #  if (is.data.frame(f)){
#   #    f=as.matrix(f)
#   #  }
#   s=vector(mode = "list",length(f))
#   i=vector(mode = "list",length(f))
#   r=vector(mode = "list",length(f))
#   
#   numrows=nticks+1
#   state_counts=vector(mode = "list",length(f))
#   
#  # colnames(state_counts)=c("S","I","R")
#   # inf_per_time[[idx]][[rep]][time,] <- stat[[idx]][[rep]]
#   # rownames(inf_per_time[[idx]][[rep]])=0:nticks#re-order row names from time step 0=initial state
#   # colnames(inf_per_time[[idx]][[rep]])=NULL
#   for (idx in 1:length(f))  
#     for(rep in 1:nrep){  
#       for (tick in 1:numrows){
#         s[[idx]][[rep]] = length(which(f[[idx]][[rep]][tick,]==0))
#         i[[idx]][[rep]] = length(which(f[[idx]][[rep]][tick,]==1))
#         r[[idx]][[rep]] = length(which(f[[idx]][[rep]][tick,]==2))
#         
#         state_counts[[rep]]<- cbind(s[[idx]][[rep]],i[[idx]][[rep]],r[[idx]][[rep]])
#     
#     
#   }
#   
#   
#     }  
#   state_counts <- do.call(rbind, state_counts)
#   state_counts<- as.data.frame(state_counts)
# #  rownames(state_counts[[idx]][[rep]])=0:(numrows-1)
#   state_counts
# }



epi_summary<-function(f){
  #  if (is.data.frame(f)){
  #    f=as.matrix(f)
  #  }
  numrows=dim(f)[1]
  state_counts=matrix(0,nrow=numrows,ncol=3)
  colnames(state_counts)=c("S","I","R")


  for (tick in 1:numrows){
    s = length(which(f[tick,]==0))
    i = length(which(f[tick,]==1))
    r = length(which(f[tick,]==2))

    state_counts[tick,] <- c(s,i,r)


  }

  #state_counts <- do.call(rbind, state_counts)
  state_counts <- as.data.frame(state_counts)
  rownames(state_counts)=0:(numrows-1)
  state_counts
}
