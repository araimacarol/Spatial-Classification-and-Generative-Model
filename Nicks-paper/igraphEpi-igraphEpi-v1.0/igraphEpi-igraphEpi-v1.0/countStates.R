##########################################################################################################
# Count the number of times each value (state) occurs per row.
##########################################################################################################
countStates <- function(DF=NULL, states=NULL) {
  nr = dim(DF)[1]
  counts <- as.data.frame(matrix(0,ncol=length(states), nrow=nr))
  colnames(counts) = as.character(states)
  for (i in 1:nr) {
    col = 1
    for (state in states) {
      counts[i,col] = length(which(DF[i,]==state))
      col = col+1
    }
  }
  return(counts)
}