library(igraph)

############# Network model ########################
#runs an outbreak on a network where you specify the beta, gamma and sigma

g <- erdos.renyi.game(100, .10,direct=F)    #random graph generated
#very simple model; every possible edge is created with the same constant probability
#100 vertices
#probability of .10 of drawing an edge between two arbitrary vertices
ER<-get.adjacency(g) #make an adjacency network; convert graph to an adjacency matrix or an edge list
# #or you can import network:
# netwb <- as.matrix(read.table("raccoon_trial_matrix.txt",header=F)) 

netwb <- ER

#running disease simulation
time<-1
gp <- ncol(netwb)    # number of columns in network    
stat<-sample(c(numeric(gp-1),1));   #choose first node to get infected. This shuffles the order of one 1 and the rest zeros (to make bigger network, need to add more 0's)
#sample function can be used to return a random permutation of a vector
#numeric creates a real vector of the specified length. The elements of the vector are all 0.
#c combines values into a vector or list
m <- matrix(0,53,gp)  # make empty matrix for number of time steps long by number of individuals wide (this will be the progression of infection per time step per indiv.) 52 plus 1 more for time=0

beta<-0.2
gamma<-1
sigma<-1/(7.5)

write(c(0,stat),file="rabiesdata.txt",ncolumns=(gp+1),append=TRUE) #writes initial conditions to file

for (i in 1:52) {  
  print(i)
  #0 means susceptible; 1 means infectious; 3 means exposed; 2 means removed/dead
  
  # for each (weekly) time step
  statc <- stat  #assign stat to copy of stat so can change 'statc' while still working off of 'stat'
  for (i in (which(stat==0)))     #for every susceptible individual in original stat...
  {  for (j in (1:gp)[-i])          #for each j in 1 to number groups, with exclusion of itself (because that is always 0)
  {
    if ((netwb[i,j]>=1)&(stat[j]==1)&(rbinom(1,1,beta)==1))  #is there an edge to an infected individual, and does a transmission event take place? (By looking down column of that indiv in netw for a 1)  
    {
      statc[i] <- 1; break   #if so, assign 3 to the copy and get out of loop
    }
  }
  }
  
  # for (i in (which(stat==3)))     #for every exposed individual in original stat...
  # {
  #   if (rbinom(1,1,sigma)==1)  #is there an exposed individual, does it become infectious?
  #   {
  #     statc[i] <- 1   #if so, assign 1 to the copy and get out of loop
  #   }
  # }
  
  for (i in (which(stat==1)))     #for every infectious individual in original stat
  {
    if (rbinom(1,1,gamma)==1)
    {
      statc[i] <- 2
    }
  }
  stat <- statc
  m[time,] <- stat
  write(c(time,m[time,]),file="rabiesdata.txt",ncolumns=(gp+1),append=TRUE)
  #rbinom(n _number of observations,size _number of trials, prob of success on each trial)
  time<-(time+1)
}









#PLOTTING:     #(apply counts time so remove time column before plotting)
plot(apply(m==1, 1, sum), col="red", type='l', main="Simulation of disease spread on raccoon contact network",sub="",xlab="Time (weeks)",ylab="Population numbers",ylim=c(0,ncol(m)),xlim=c(0,50))   # ‘apply’ calculates the sum across the rows (middle value=1 vs 2 for columns) for all values of matrix m that create a TRUE when value =1.      INFECTIOUS
lines(apply(m==0, 1,sum), col="green")   #green is sum of #0  SUSCEPTIBLE
lines(apply(m==2, 1,sum), col="blue")     #blue is sum of #2  REMOVED
lines(apply(m==3, 1,sum), col="yellow")           #EXPOSED
legend("right", inset=.05,c("Susceptible","Exposed","Infectious","Removed"),col=c("green","yellow","red","blue"),lty=c(1,1,1,1))


epidemic <- data.frame(m)
epidemic2 <- epidemic[c(1:52),]
epidemic2$any.infx <- NA

### duration of epidemic ###
# 1 = infectious, so want the time when number of 1's = 0
for(i in 1:nrow(epidemic2)){
  print(i)
  temp <- epidemic[i,]
  epidemic2$any.infx[i] <- ifelse(any(temp==1), "TRUE", "FALSE")
}

# without an expposed class, the end of the epidemic is when "any.infx" first becomes false

### number of individuals infected ###
# theoretically, any individual infected should have a sum greater than 0 because only susceptible individuals will stay at 0 throughout
positives <- colSums(epidemic)
length(subset(positives, positives>0))
