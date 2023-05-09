#[1] simulating a single SIR process on each simulated graph and extracting the summary

# single_epi_summary<-function(f){
#     summary_all_graphs=lapply(f,epi_summary)
# 
# return(summary_all_graphs)
# }
#[2] simulating multiple SIR processes on each simulated graph and extracting the summary

summary_all_graphs=NULL

nets_epi_summary<-function(f){
  for(idx in 1:length(f)){
    for (rep in 1:length(f)){
      summary_all_graphs[[rep]]=lapply(f[[idx]],epi_summary)
    }
  }
    return(summary_all_graphs)
}



#[3] Time to infection of 20% of individual
