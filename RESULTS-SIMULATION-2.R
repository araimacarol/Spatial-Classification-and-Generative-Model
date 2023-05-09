# ########################------------ Data Collection and analysis-----------------------------------------######################

source("RunEpicSimandMeasures.R")

source("different-average-degrees.R")

source("Different-degree-distribution.R")

source("FINAL-PIPELINE.R")

###---Packages----###
library(showtext)

library(png)

library(grid)

library(gridExtra)

#############-----------DATA--COLLECTION-----##########################################
########--- 5 NODES-------------############
set.seed(3456)
nsim=10 # number of simulated graphs for each graph structure
nreps=10 # number of epidemic simulations repeated on each graph
num_vertices_5=5 
nticks=100 
G_5_avgdeg_4 <-makeGraphs_5_avgdeg_4(nSamples=nsim,order = num_vertices_5)# 10 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) # transmission rate beta
gammaVals=c(0.033,0.143,0.25,0.5) # recovery rate gamma

observed_networks_5_avrgdeg_4=ParalleEpicSimOnGraphs(G_5_avgdeg_4, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observed_networks_5_avrgdeg_4.csv",report="i")

########--- 10 NODES-------------############
set.seed(3456)
nsim=10 
nreps=10
num_vertices_10=10
nticks=100 
G_10_avgdeg_4 <-makeGraphs_10_avgdeg_4(nSamples=nsim,order = num_vertices_10)

betaVals=c(0.01,0.1,0.33,0.5)
gammaVals=c(0.033,0.143,0.25,0.5)

observed_networks_10_avrgdeg_4=ParalleEpicSimOnGraphs(G_10_avgdeg_4, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observed_networks_10_avrgdeg_4.csv",report="i")


############--- 20 NODES-------------############
set.seed(3456)
nsim=10 
nreps=10
num_vertices_20=20 
nticks=100 
G_20_avgdeg_4 <-makeGraphs_20_avgdeg_4(nSamples=nsim,order = num_vertices_20)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)

observed_networks_20_avrgdeg_4=ParalleEpicSimOnGraphs(G_20_avgdeg_4, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observed_networks_20_avrgdeg_4.csv",report="i")

########--- 30 NODES-------------############
set.seed(3456)
nsim=10 
nreps=10
num_vertices_30=30 
nticks=100 
G_30_avgdeg_4 <-makeGraphs_30_avgdeg_4(nSamples=nsim,order = num_vertices_30)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)

observed_networks_30_avrgdeg_4=ParalleEpicSimOnGraphs(G_30_avgdeg_4, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observed_networks_30_avrgdeg_4.csv",report="i")


########--- 40 NODES-------------############
set.seed(3456)
nsim=10 
nreps=10
num_vertices_40=40 
nticks=100 
G_40_avgdeg_4 <-makeGraphs_40_avgdeg_4(nSamples=nsim,order = num_vertices_40)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)

observed_networks_40_avrgdeg_4=ParalleEpicSimOnGraphs(G_40_avgdeg_4, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observed_networks_40_avrgdeg_4.csv",report="i")



########--- 50 NODES-------------############
set.seed(3456)
nsim=10 
nreps=10
num_vertices_50=50 
nticks=100 
G_50_avgdeg_4 <-makeGraphs_50_avgdeg_4(nSamples=nsim,order = num_vertices_50)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5)
gammaVals=c(0.033,0.143,0.25,0.5)

observed_networks_50_avrgdeg_4=ParalleEpicSimOnGraphs(G_50_avgdeg_4, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observed_networks_50_avrgdeg_4.csv",report="i")


########--- 150 NODES-------------############
nsim=10 
nreps=10
num_vertices_150=150 
nticks=100 
G_150_avgdeg_4 <-makeGraphs_150_avgdeg_4(nSamples=nsim,order = num_vertices_150)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)

observed_networks_150_avrgdeg_4=ParalleEpicSimOnGraphs(G_150_avgdeg_4, 
                                                       nticks=nticks, beta=betaVals,gamma=gammaVals,
                                                       nreps=nreps,output_file="observed_networks_150_avrgdeg_4.csv ",report="i")

########--- 250 NODES-------------############
set.seed(3456)
nsim=10 
nreps=10
num_vertices_250=250 
nticks=100 
G_250_avgdeg_4 <-makeGraphs_250_avgdeg_4(nSamples=nsim,order = num_vertices_250)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)

observed_networks_250_avrgdeg_4=ParalleEpicSimOnGraphs(G_250_avgdeg_4, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observed_networks_250_avrgdeg_4.csv",report="i")

########--- 500 NODES-------------############
set.seed(3456)
nsim=10 
nreps=10
num_vertices_500=500 
nticks=100 
G_500_avgdeg_4 <-makeGraphs_500_avgdeg_4(nSamples=nsim,order = num_vertices_500)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)

observed_networks_500_avrgdeg_4=ParalleEpicSimOnGraphs(G_500_avgdeg_4, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observed_networks_500_avrgdeg_4.csv",report="i")

###############----1000 nodes----------####
set.seed(3456)
nsim=10 
nreps=10
num_vertices_1000=1000 
nticks=100 
G_1000_avgdeg_4 <-makeGraphs_1000_avgdeg_4(nSamples=nsim,order = num_vertices_1000)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)

observed_networks_1000_avrgdeg_4=ParalleEpicSimOnGraphs(G_1000_avgdeg_4, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observed_networks_1000_avrgdeg_4.csv",report="i")


##################----Observed Networks with average degree 6----------#############################

##### 50 nodes ###############  
set.seed(3456)
nsim=10 
nreps=10
num_vertices_50=50 
nticks=100 
G_50_avgdeg_6 <-makeGraphs_50_avgdeg_6(nSamples=nsim,order = num_vertices_50)# 5 samples for each graph
betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)
observed_networks_50_avrgdeg_6=ParalleEpicSimOnGraphs(G_50_avgdeg_6, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observed_networks_50_avrgdeg_6.csv",report="i")

##########--250 node--############
set.seed(3456)
nsim=10 
nreps=10
num_vertices_250=250 
#size=((4* num_vertices_50)/2) #target number of edges in each graph
nticks=100 
G_250_avgdeg_6 <-makeGraphs_250_avgdeg_6(nSamples=nsim,order = num_vertices_250)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)
observed_networks_250_avrgdeg_6=ParalleEpicSimOnGraphs(G_250_avgdeg_6, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observed_networks_250_avrgdeg_6.csv",report="i")


####----500 nodes----------####
set.seed(3456)
nsim=10 
nreps=10
num_vertices_500=500 
#size=((4* num_vertices_50)/2) #target number of edges in each graph
nticks=100 
G_500_avgdeg_6 <-makeGraphs_500_avgdeg_6(nSamples=nsim,order = num_vertices_500)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)
observed_networks_500_avrgdeg_6=ParalleEpicSimOnGraphs(G_500_avgdeg_6, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observed_networks_500_avrgdeg_6.csv",report="i")


######################----1000 nodes----------####
set.seed(3456)
nsim=10 
nreps=10
num_vertices_1000=1000 
#size=((4* num_vertices_50)/2) #target number of edges in each graph
nticks=100 
G_1000_avgdeg_6 <-makeGraphs_1000_avgdeg_6_1(nSamples=nsim,order = num_vertices_1000)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)

observed_networks_1000_avrgdeg_6=ParalleEpicSimOnGraphs(G_1000_avgdeg_6, 
                                                        nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observenetwork6B.csv",report="i")




##################----Observed Networks with average degree 8----------############################# 

#####----50 nodes----------####
set.seed(3456)
nsim=10
nreps=10
num_vertices_50=50 
nticks=100 
G_50_avgdeg_8 <-makeGraphs_50_avgdeg_8(nSamples=nsim,order = num_vertices_50)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)#  

observed_networks_50_avrgdeg_8=ParalleEpicSimOnGraphs(G_50_avgdeg_8, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observed_networks_50_avrgdeg_8.csv",report="i")


#####----250 nodes----------####
set.seed(3456)
nsim=10
nreps=10
num_vertices_250=250 
#size=((4* num_vertices_50)/2) #target number of edges in each graph
nticks=100 
G_250_avgdeg_8 <-makeGraphs_250_avgdeg_8(nSamples=nsim,order = num_vertices_250)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)#  

observed_networks_250_avrgdeg_8=ParalleEpicSimOnGraphs(G_250_avgdeg_8, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observed_networks_250_avrgdeg_8.csv",report="i")


#####----500 nodes----------####
set.seed(3456)
nsim=10
nreps=10
num_vertices_500=500 
#size=((4* num_vertices_50)/2) #target number of edges in each graph
nticks=100 
G_500_avgdeg_8 <-makeGraphs_500_avgdeg_8(nSamples=nsim,order = num_vertices_500)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)#  

observed_networks_500_avrgdeg_8=ParalleEpicSimOnGraphs(G_500_avgdeg_8, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observed_networks_500_avrgdeg_8.csv",report="i")


######################----1000 nodes----------################
set.seed(3456)
nsim=10 
nreps=10
num_vertices_1000=1000 
#size=((4* num_vertices_50)/2) #target number of edges in each graph
nticks=100 
G_1000_avgdeg_8 <-makeGraphs_1000_avgdeg_8(nSamples=nsim,order = num_vertices_1000)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)#  

observed_networks_1000_avrgdeg_8=ParalleEpicSimOnGraphs(G_1000_avgdeg_8, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observed_networks_1000_avrgdeg_8.csv",report="i")




#####################----Observed Networks with average degree 10----------#############################

#####----50 nodes----------####
set.seed(3456)
nsim=10
nreps=10
num_vertices_50=50 
nticks=100 
G_50_avgdeg_10 <-makeGraphs_50_avgdeg_10(nSamples=nsim,order = num_vertices_50)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)#  

observed_networks_50_avrgdeg_10=ParalleEpicSimOnGraphs(G_50_avgdeg_10, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observed_networks_50_avrgdeg_10.csv",report="i")



###########--- 150 NODES-------------############
nsim=10 
nreps=10
num_vertices_150=150 
nticks=100 
G_150_avgdeg_10 <-makeGraphs_150_avgdeg_10(nSamples=nsim,order = num_vertices_150)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)

observed_networks_150_avrgdeg_10=ParalleEpicSimOnGraphs(G_150_avgdeg_10, 
                                                        nticks=nticks, beta=betaVals,gamma=gammaVals,
                                                        nreps=nreps,output_file="observed_networks_150_avrgdeg_10.csv ",report="i")



#####----250 nodes----------####
set.seed(3456)
nsim=10
nreps=10
num_vertices_250=250 
#size=((4* num_vertices_50)/2) #target number of edges in each graph
nticks=100 
G_250_avgdeg_10 <-makeGraphs_250_avgdeg_10(nSamples=nsim,order = num_vertices_250)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)#  

observed_networks_250_avrgdeg_10=ParalleEpicSimOnGraphs(G_250_avgdeg_10, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observed_networks_250_avrgdeg_10.csv",report="i")




#####----500 nodes----------####
set.seed(3456)
nsim=10 # number of simulated graphs
nreps=10 # number of epidemic simulations to run on each graph

num_vertices_500=500 
#size=((4* num_vertices_50)/2) #target number of edges in each graph
nticks=100 
G_500_avgdeg_10 <-makeGraphs_500_avgdeg_10(nSamples=nsim,order = num_vertices_500)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)#  

observed_networks_500_avrgdeg_10=ParalleEpicSimOnGraphs(G_500_avgdeg_10, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observed_networks_500_avrgdeg_10.csv",report="i")


##################----1000 nodes----------####
set.seed(3456)
nsim=10 
nreps=10
num_vertices_1000=1000 
#size=((4* num_vertices_50)/2) #target number of edges in each graph
nticks=100 
G_1000_avgdeg_10 <-makeGraphs_1000_avgdeg_10(nSamples=nsim,order = num_vertices_1000)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)#  

observed_networks_1000_avrgdeg_10=ParalleEpicSimOnGraphs(G_1000_avgdeg_10, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observed_networks_1000_avrgdeg_10.csv",report="i")



##################----Observed Networks with average degree 12----------#############################


#######--- 20 NODES-------------############
set.seed(3456)
nsim=10 
nreps=10
num_vertices_20=20 
nticks=100 
G_20_avgdeg_12 <-makeGraphs_20_avgdeg_12(nSamples=nsim,order = num_vertices_20)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)

observed_networks_20_avrgdeg_12=ParalleEpicSimOnGraphs(G_20_avgdeg_12, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observed_networks_20_avrgdeg_12.csv",report="i")


########--- 30 NODES-------------############
set.seed(3456)
nsim=10 
nreps=10
num_vertices_30=30 
nticks=100 
G_30_avgdeg_12 <-makeGraphs_30_avgdeg_12(nSamples=nsim,order = num_vertices_30)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)

observed_networks_30_avrgdeg_12=ParalleEpicSimOnGraphs(G_30_avgdeg_12, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observed_networks_30_avrgdeg_12.csv",report="i")


########--- 40 NODES-------------############
set.seed(3456)
nsim=10 
nreps=10
num_vertices_40=40 
nticks=100 
G_40_avgdeg_12<-makeGraphs_40_avgdeg_12(nSamples=nsim,order = num_vertices_40)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)

observed_networks_40_avrgdeg_12=ParalleEpicSimOnGraphs(G_40_avgdeg_12, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observed_networks_40_avrgdeg_12.csv",report="i")

########--- 50 NODES-------------############
set.seed(3456)
nsim=10 
nreps=10
num_vertices_50=50 
nticks=100 
G_50_avgdeg_12 <-makeGraphs_50_avgdeg_12(nSamples=nsim,order = num_vertices_50)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)

observed_networks_50_avrgdeg_12=ParalleEpicSimOnGraphs(G_50_avgdeg_12, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observed_networks_50_avrgdeg_12.csv",report="i")



########--- 150 NODES-------------############
nsim=10 
nreps=10
num_vertices_150=150 
nticks=100 
G_150_avgdeg_12<-makeGraphs_150_avgdeg_12(nSamples=nsim,order = num_vertices_150)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)

observed_networks_150_avrgdeg_12=ParalleEpicSimOnGraphs(G_150_avgdeg_12, 
                                                        nticks=nticks, beta=betaVals,gamma=gammaVals,
                                                        nreps=nreps,output_file="observed_networks_150_avrgdeg_12.csv ",report="i")




######################----250 NODES-------------############        
set.seed(3456)
nsim=10 
nreps=10
num_vertices_250=250 
#size=((4* num_vertices_50)/2) #target number of edges in each graph
nticks=100 
G_250_avgdeg_12 <-makeGraphs_250_avgdeg_12(nSamples=nsim,order = num_vertices_250)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)

observed_networks_250_avrgdeg_12=ParalleEpicSimOnGraphs(G_250_avgdeg_12, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observed_networks_250_avrgdeg_12.csv",report="i")

###############-----500 NODES-------------############    
set.seed(3456)
nsim=10 
nreps=10
num_vertices_500=500 
#size=((4* num_vertices_50)/2) #target number of edges in each graph
nticks=100 
G_500_avgdeg_12 <-makeGraphs_500_avgdeg_12(nSamples=nsim,order = num_vertices_500)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)

observed_networks_500_avrgdeg_12=ParalleEpicSimOnGraphs(G_500_avgdeg_12, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="observed_networks_500_avrgdeg_12.csv",report="i")



################----1000 nodes----------##############
set.seed(3456)
nsim=10 
nreps=10
num_vertices_1000=1000 
nticks=100 
G_1000_avgdeg_12<-makeGraphs_1000_avgdeg_12_1(nSamples=nsim,order = num_vertices_1000)# 5 samples for each graph

betaVals=c(0.01,0.1,0.33,0.5) 
gammaVals=c(0.033,0.143,0.25,0.5)#

observed_networks_1000_avrgdeg_12=ParalleEpicSimOnGraphs(G_1000_avgdeg_12, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="A1_1000_avrgdeg_12.csv",report="i")



######################-------------DATA--ANALYSIS--------------######################################           
##--Epidemic--Parameters--##
betaVals=c(0.01,0.1,0.33,0.5) ##transmission rate
gammaVals=c(0.033,0.143,0.25,0.5) ##recovery rate


#########--------RESULTS FOR AVERAGE DEGREE  OF 4-----------##################

################--------------PROPORTION OF INFECTED INDIVIDUALS---------###################

#########----5 nodes-------------------#########
set.seed((234567))
nsim=10
ER_5_avgd4=plotfunc_obs(Name = "ER", ID=1,betaVal = .5,gammaVal=.143,Data = "observed_networks_5_avrgdeg_4.csv", nreps = nsim,plot_title="ER",net_size = 5)
SW_5_avgd4=plotfunc_obs(Name = "SW", ID=1,betaVal = .5,gammaVal=.143,Ylabel = "",Data = "observed_networks_5_avrgdeg_4.csv", nreps = nsim,plot_title="SW",net_size = 5)
SF_5_avgd4=plotfunc_obs(Name = "SF", ID=1,betaVal = .5,gammaVal=.143,Ylabel = "",Data = "observed_networks_5_avrgdeg_4.csv", nreps = nsim,plot_title="SF",net_size = 5)
SP_5_avgd4=plotfunc_obs(Name = "SP", ID=1,betaVal = .5,gammaVal=.143,Ylabel = "",Data = "observed_networks_5_avrgdeg_4.csv", nreps = nsim,plot_title="SP",net_size = 5)
Lat_5_avgd4=plotfunc_obs(Name = "Lat", ID=1,betaVal = .5,gammaVal=.143,Data = "observed_networks_5_avrgdeg_4.csv", nreps = nsim,plot_title="Lat",net_size = 5)
CG_5_avgd4=Complete_Graph(beta = .5,gamma=0.143,net_size = 5,sim = 10,ntime = 100,avrdg = 4,Ylabel = "")




#########----10 nodes-------------------#########
set.seed((234567))
nsim=10
ER_10_avgd4=plotfunc_obs(Name = "ER", ID=1,betaVal = .5,gammaVal=.143,Data = "observed_networks_10_avrgdeg_4.csv", nreps = nsim,plot_title="ER",net_size = 10)
SW_10_avgd4=plotfunc_obs(Name = "SW", ID=1,betaVal = .5,gammaVal=.143,Ylabel = "",Data = "observed_networks_10_avrgdeg_4.csv", nreps = nsim,plot_title="SW",net_size = 10)
SF_10_avgd4=plotfunc_obs(Name = "SF", ID=1,betaVal = .5,gammaVal=.143,Ylabel = "",Data = "observed_networks_10_avrgdeg_4.csv", nreps = nsim,plot_title="SF",net_size = 10)
SP_10_avgd4=plotfunc_obs(Name = "SP", ID=1,betaVal = .5,gammaVal=.143,Ylabel = "",Data = "observed_networks_10_avrgdeg_4.csv", nreps = nsim,plot_title="SP",net_size = 10)
Lat_10_avgd4=plotfunc_obs(Name = "Lat", ID=1,betaVal = .5,gammaVal=.143,Data = "observed_networks_10_avrgdeg_4.csv", nreps = nsim,plot_title="Lat",net_size = 10)
CG_10_avgd4=Complete_Graph(beta = 0.5,gamma=0.143,net_size = 10,sim = 10,ntime = 100,avrdg = 4,Ylabel = "")




#########----20 nodes-------------------#########
set.seed((234567))
nsim=10
ER_20_avgd4=plotfunc_obs(Name = "ER", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Data = "observed_networks_20_avrgdeg_4.csv", nreps = nsim,plot_title="ER",net_size = 20)
SW_20_avgd4=plotfunc_obs(Name = "SW", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Ylabel = "",Data = "observed_networks_20_avrgdeg_4.csv", nreps = nsim,plot_title="SW",net_size = 20)
SF_20_avgd4=plotfunc_obs(Name = "SF", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Ylabel = "",Data = "observed_networks_20_avrgdeg_4.csv", nreps = nsim,plot_title="SF",net_size = 20)
SP_20_avgd4=plotfunc_obs(Name = "SP", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Ylabel = "",Data = "observed_networks_20_avrgdeg_4.csv", nreps = nsim,plot_title="SP",net_size = 20)
Lat_20_avgd4=plotfunc_obs(Name = "Lat", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Data = "observed_networks_20_avrgdeg_4.csv", nreps = nsim,plot_title="Lat",net_size = 20)
CG_20_avgd4=Complete_Graph(beta = 0.5,gamma=0.143,net_size = 20,sim = 10,ntime = 100,avrdg = 4,Xlabel = "",Ylabel = "")




#########----30 nodes-------------------#########
set.seed((234567))
nsim=10
ER_30_avgd4=plotfunc_obs(Name = "ER", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Data = "observed_networks_30_avrgdeg_4.csv", nreps = nsim,plot_title="ER",net_size = 30)

SW_30_avgd4=plotfunc_obs(Name = "SW", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Ylabel = "",Data = "observed_networks_30_avrgdeg_4.csv", nreps = nsim,plot_title="SW",net_size = 30)
SF_30_avgd4=plotfunc_obs(Name = "SF", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Ylabel = "",Data = "observed_networks_30_avrgdeg_4.csv", nreps = nsim,plot_title="SF",net_size = 30)
SP_30_avgd4=plotfunc_obs(Name = "SP", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Ylabel = "",Data = "observed_networks_30_avrgdeg_4.csv", nreps = nsim,plot_title="SP",net_size = 30)
Lat_30_avgd4=plotfunc_obs(Name = "Lat", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Data = "observed_networks_30_avrgdeg_4.csv", nreps = nsim,plot_title="Lat",net_size = 30)
CG_30_avgd4=Complete_Graph(beta = 0.5,gamma=0.143,net_size = 30,sim = 10,ntime = 100,avrdg = 4,Xlabel = "",Ylabel = "")


#########----40 nodes-------------------#########
set.seed((234567))
nsim=10
ER_40_avgd4=plotfunc_obs(Name = "ER", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Data = "observed_networks_40_avrgdeg_4.csv", nreps = nsim,plot_title="ER",net_size = 40)
SW_40_avgd4=plotfunc_obs(Name = "SW", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Ylabel = "",Data = "observed_networks_40_avrgdeg_4.csv", nreps = nsim,plot_title="SW",net_size = 40)
SF_40_avgd4=plotfunc_obs(Name = "SF", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Ylabel = "",Data = "observed_networks_40_avrgdeg_4.csv", nreps = nsim,plot_title="SF",net_size = 40)
SP_40_avgd4=plotfunc_obs(Name = "SP", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Ylabel = "",Data = "observed_networks_40_avrgdeg_4.csv", nreps = nsim,plot_title="SP",net_size = 40)
Lat_40_avgd4=plotfunc_obs(Name = "Lat", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Data = "observed_networks_40_avrgdeg_4.csv", nreps = nsim,plot_title="Lat",net_size = 40)
CG_40_avgd4=Complete_Graph(beta = 0.5,gamma=0.143,net_size = 40,sim = 10,ntime = 100,avrdg = 4,Xlabel = "",Ylabel = "")


#########----50 nodes-------------------#########
set.seed((234567))
nsim=10
ER_50_avgd4=plotfunc_obs(Name = "ER", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Data = "observed_networks_50_avrgdeg_4.csv", nreps = nsim,plot_title="ER",net_size = 50)
SW_50_avgd4=plotfunc_obs(Name = "SW", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Ylabel = "",Data = "observed_networks_50_avrgdeg_4.csv", nreps = nsim,plot_title="SW",net_size = 50)
SF_50_avgd4=plotfunc_obs(Name = "SF", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Ylabel = "",Data = "observed_networks_50_avrgdeg_4.csv", nreps = nsim,plot_title="SF",net_size = 50)
SP_50_avgd4=plotfunc_obs(Name = "SP", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Ylabel = "",Data = "observed_networks_50_avrgdeg_4.csv", nreps = nsim,plot_title="SP",net_size = 50)
Lat_50_avgd4=plotfunc_obs(Name = "Lat", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Data = "observed_networks_50_avrgdeg_4.csv", nreps = nsim,plot_title="Lat",net_size = 50)
CG_50_avgd4=Complete_Graph(beta = 0.5,gamma=0.143,net_size = 50,sim = 10,ntime = 100,avrdg = 4,Xlabel = "",Ylabel = "")




#########----150 nodes-------------------#########
set.seed((234567))
nsim=10
ER_150_avgd4=plotfunc_obs(Name = "ER", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Data = "observed_networks_150_avrgdeg_4.csv", nreps = nsim,plot_title="ER",net_size = 150)
SW_150_avgd4=plotfunc_obs(Name = "SW", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Ylabel = "",Data = "observed_networks_150_avrgdeg_4.csv", nreps = nsim,plot_title="SW",net_size = 150)
SF_150_avgd4=plotfunc_obs(Name = "SF", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Ylabel = "",Data = "observed_networks_150_avrgdeg_4.csv", nreps = nsim,plot_title="SF",net_size = 150)
SP_150_avgd4=plotfunc_obs(Name = "SP", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Ylabel = "",Data = "observed_networks_150_avrgdeg_4.csv", nreps = nsim,plot_title="SP",net_size = 150)
Lat_150_avgd4=plotfunc_obs(Name = "Lat", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Data = "observed_networks_150_avrgdeg_4.csv", nreps = nsim,plot_title="Lat",net_size = 150)
CG_150_avgd4=Complete_Graph(beta = 0.5,gamma=0.143,net_size = 150,sim = 10,ntime = 100,avrdg = 4,Xlabel = "",Ylabel = "")


#########----250 nodes-------------------#########
set.seed((234567))
nsim=10
ER_250_avgd4=plotfunc_obs(Name = "ER", ID=1,betaVal = 0.5,gammaVal=.143,Xlabel = "",Data = "observed_networks_250_avrgdeg_4.csv", nreps = nsim,plot_title="ER",net_size = 250)
SW_250_avgd4=plotfunc_obs(Name = "SW", ID=1,betaVal = 0.5,gammaVal=.143,Xlabel = "",Ylabel = "",Data = "observed_networks_250_avrgdeg_4.csv", nreps = nsim,plot_title="SW",net_size = 250)
SF_250_avgd4=plotfunc_obs(Name = "SF", ID=1,betaVal = 0.5,gammaVal=.143,Xlabel = "",Ylabel = "",Data = "observed_networks_250_avrgdeg_4.csv", nreps = nsim,plot_title="SF",net_size = 250)
SP_250_avgd4=plotfunc_obs(Name = "SP", ID=1,betaVal = 0.5,gammaVal=.143,Xlabel = "",Ylabel = "",Data = "observed_networks_250_avrgdeg_4.csv", nreps = nsim,plot_title="SP",net_size = 250)
Lat_250_avgd4=plotfunc_obs(Name = "Lat", ID=1,betaVal = 0.5,gammaVal=.143,Xlabel = "",Data = "observed_networks_250_avrgdeg_4.csv", nreps = nsim,plot_title="Lat",net_size = 250)
CG_250_avgd4=Complete_Graph(beta = 0.5,gamma=0.143,net_size = 250,sim = 10,ntime = 100,avrdg = 4,Xlabel = "",Ylabel = "")

#########----500 nodes-------------------#########
set.seed((234567))
nsim=10
ER_500_avgd4=plotfunc_obs(Name = "ER", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Data = "observed_networks_500_avrgdeg_4.csv", nreps = nsim,plot_title="ER",net_size = 500)
SW_500_avgd4=plotfunc_obs(Name = "SW", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Ylabel = "",Data = "observed_networks_500_avrgdeg_4.csv", nreps = nsim,plot_title="SW",net_size = 500)
SF_500_avgd4=plotfunc_obs(Name = "SF", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Ylabel = "",Data = "observed_networks_500_avrgdeg_4.csv", nreps = nsim,plot_title="SF",net_size = 500)
SP_500_avgd4=plotfunc_obs(Name = "SP", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Ylabel = "",Data = "observed_networks_500_avrgdeg_4.csv", nreps = nsim,plot_title="SP",net_size = 500)
Lat_500_avgd4=plotfunc_obs(Name = "Lat", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Data = "observed_networks_500_avrgdeg_4.csv", nreps = nsim,plot_title="Lat",net_size = 500)
CG_500_avgd4=Complete_Graph(beta = 0.5,gamma=0.143,net_size = 500,sim = 10,ntime = 100,avrdg = 4,Xlabel = "",Ylabel = "")


#########----1000 nodes------------##############
set.seed((234567))
nsim=10
ER_1000_avgd4=plotfunc_obs(Name = "ER", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Data = "observed_networks_1000_avrgdeg_4.csv", nreps = nsim,plot_title="ER",net_size = 1000)
SW_1000_avgd4=plotfunc_obs(Name = "SW", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Ylabel = "",Data = "observed_networks_1000_avrgdeg_4.csv", nreps = nsim,plot_title="SW",net_size = 1000)
SF_1000_avgd4=plotfunc_obs(Name = "SF", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Ylabel = "",Data = "observed_networks_1000_avrgdeg_4.csv", nreps = nsim,plot_title="SF",net_size = 1000)
SP_1000_avgd4=plotfunc_obs(Name = "SP", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Ylabel = "",Data = "observed_networks_1000_avrgdeg_4.csv", nreps = nsim,plot_title="SP",net_size = 1000)
Lat_1000_avgd4=plotfunc_obs(Name = "Lat", ID=1,betaVal = .5,gammaVal=.143,Xlabel = "",Data = "observed_networks_1000_avrgdeg_4.csv", nreps = nsim,plot_title="Lat",net_size = 1000)
CG_1000_avgd4=Complete_Graph(beta = 0.5,gamma=0.143,net_size = 1000,sim = 10,ntime = 100,avrdg = 4,Xlabel = "",Ylabel = "")






#####################------------------NETWORK OF ALL SIZES AVERAGE DEGREE 4-----------################
#########-----PLOT_A------################
plotA_5=ggarrange(ER_5_avgd4,SW_5_avgd4,SF_5_avgd4,ncol=3,nrow = 1)
#plotA_5_annotate=annotate_figure(plotA_5, right = text_grob(" 5 nodes", color = "black", face = "bold", size = 26))

plotA_10=ggarrange(ER_10_avgd4,SW_10_avgd4,SF_10_avgd4,ncol=3,nrow = 1)

#plotA_10_annotate=annotate_figure(plotA_10, right = text_grob(" 10 nodes", color = "black", face = "bold", size = 26))
plotA_20=ggarrange(ER_20_avgd4,SW_20_avgd4,SF_20_avgd4,ncol=3,nrow = 1)
#plotA_20_annotate=annotate_figure(plotA_20, right = text_grob(" 20 nodes",color = "black", face = "bold", size = 26))

plotA_30=ggarrange(ER_30_avgd4,SW_30_avgd4,SF_30_avgd4,ncol=3,nrow = 1)
#plotA_30_annotate=annotate_figure(plotA_30, right = text_grob(" 30 nodes",color = "black", face = "bold", size = 26))

plotA_40=ggarrange(ER_40_avgd4,SW_40_avgd4,SF_40_avgd4,ncol=3,nrow = 1)
#plotA_40_annotate=annotate_figure(plotA_40, right = text_grob(" 40 nodes",color = "black", face = "bold", size = 26))

plotA_50=ggarrange(ER_50_avgd4,SW_50_avgd4,SF_50_avgd4,ncol=3,nrow = 1)
plotA_50_annotate=annotate_figure(plotA_50, right = text_grob(" 50 nodes", 
                                                              color = "black", face = "bold", size = 26))

plotA_150=ggarrange(ER_150_avgd4,SW_150_avgd4,SF_150_avgd4,ncol=3,nrow = 1)

#plotA_150_annotate=annotate_figure(plotA_150, right = text_grob(" 150 nodes", color = "black", face = "bold", size = 26))

plotA_250=ggarrange(ER_250_avgd4,SW_250_avgd4,SF_250_avgd4,ncol=3,nrow = 1)

#plotA_250_annotate=annotate_figure(plotA_250, right = text_grob(" 250 nodes",color = "black", face = "bold", size = 26))
plotA_500=ggarrange(ER_500_avgd4,SW_500_avgd4,SF_500_avgd4,ncol=3,nrow = 1)

#plotA_500_annotate=annotate_figure(plotA_500, right = text_grob(" 500 nodes",color = "black", face = "bold", size = 26))

plotA_1000=ggarrange(ER_1000_avgd4,SW_1000_avgd4,SF_1000_avgd4,ncol=3,nrow = 1)

#plotA_1000_annotate=annotate_figure(plotA_1000, right = text_grob(" 1000 nodes",color = "black", face = "bold", size = 26))


ALLPLOTA1=ggarrange(plotA_5,plotA_50,
                    plotA_500, nrow=3)

ggsave("ALLPLOTA1.png", width = 26, height = 28)
ggsave("ALLPLOTA1.pdf", width = 26, height = 28)

ALLPLOTA2=ggarrange(plotA_10,plotA_20,
                             plotA_30,plotA_40,plotA_150,
                             plotA_250,plotA_1000, nrow=7)

ggsave("ALLPLOTA2.png", width = 26, height = 40)
ggsave("ALLPLOTA2.pdf", width = 26, height = 40)



###--PLOTB----####
plotB_5=ggarrange(Lat_5_avgd4,SP_5_avgd4,CG_5_avgd4,ncol=3,nrow = 1)

#plotB_5_annotate=annotate_figure(plotB_5, right = text_grob(" 5 nodes", color = "black", face = "bold", size = 26))

plotB_10=ggarrange(Lat_10_avgd4,SP_10_avgd4,CG_10_avgd4,ncol=3,nrow = 1)

#plotB_10_annotate=annotate_figure(plotB_10, right = text_grob(" 10 nodes",color = "black", face = "bold", size = 26))

plotB_20=ggarrange(Lat_20_avgd4,SP_20_avgd4,CG_20_avgd4,ncol=3,nrow = 1)

#plotB_20_annotate=annotate_figure(plotB_20, right = text_grob(" 20 nodes", color = "black", face = "bold", size = 26))

plotB_30=ggarrange(Lat_30_avgd4,SP_30_avgd4,CG_30_avgd4,ncol=3,nrow = 1)
#plotB_30_annotate=annotate_figure(plotB_30, right = text_grob(" 30 nodes", color = "black", face = "bold", size = 26))

plotB_40=ggarrange(Lat_40_avgd4,SP_40_avgd4,CG_40_avgd4,ncol=3,nrow = 1)
#plotB_40_annotate=annotate_figure(plotB_40, right = text_grob(" 40 nodes", color = "black", face = "bold", size = 26))

plotB_50=ggarrange(Lat_50_avgd4,SP_50_avgd4,CG_50_avgd4,ncol=3,nrow = 1)

#plotB_50_annotate=annotate_figure(plotB_50, right = text_grob(" 50 nodes", color = "black", face = "bold", size = 26))

plotB_150=ggarrange(Lat_150_avgd4,SP_150_avgd4,CG_150_avgd4,ncol=3,nrow = 1)

#plotB_150_annotate=annotate_figure(plotB_150, right = text_grob(" 150 nodes", color = "black", face = "bold", size = 26))


plotB_250=ggarrange(Lat_250_avgd4,SP_250_avgd4,CG_250_avgd4,ncol=3,nrow = 1)

#plotB_250_annotate=annotate_figure(plotB_250, right = text_grob(" 250 nodes",  color = "black", face = "bold", size = 26))

plotB_500=ggarrange(Lat_500_avgd4,SP_500_avgd4,CG_500_avgd4,ncol=3,nrow = 1)

#plotB_500_annotate=annotate_figure(plotB_500, right = text_grob(" 500 nodes", color = "black", face = "bold", size = 26))

plotB_1000=ggarrange(Lat_1000_avgd4,SP_1000_avgd4,CG_1000_avgd4,ncol=3,nrow = 1)

#plotB_1000_annotate=annotate_figure(plotB_1000, right = text_grob(" 1000 nodes", color = "black", face = "bold", size = 26))


ALLPLOTB1=ggarrange(plotB_5,plotB_50,
                    plotB_500, nrow=3)

ggsave("ALLPLOTB1.png", width = 26, height = 28)

ggsave("ALLPLOTB1.pdf", width = 26, height = 28)


ALLPLOTB2=ggarrange(plotB_10,plotB_20,
                    plotB_30,plotB_40,plotB_150,
                    plotB_250,plotB_1000, nrow=7)

ggsave("ALLPLOTB2.png", width = 26, height = 40)
ggsave("ALLPLOTB2.pdf", width = 26, height = 40)

set.seed(82626)
nsim=10

############----------AVERAGE PROPORTION OF INFECTED FUNCTION------###############
AvrgPropOf<-function(Name="ER",ID=1,betaVal=.01,gammaVal=.2,nreps=nsim,Data="observed_networks_50_avrgdeg_4.csv",plot_title="",nticks=100,net_size=10){
  df_obs=data.frame(read.csv(Data,header = T, sep = ","))
  df_obs=df_obs%>% filter(GraphName==Name,GraphID==ID, beta==betaVal,gamma==gammaVal)
  df_obs=df_obs%>%select(c(GraphName,t0:paste("t",nticks,sep = "")))
  #y=f%>%group_by(GraphName)%>%summarize(colMeans(f[sapply(f, is.numeric)]))
  df=data.frame(colMeans(df_obs[sapply(df_obs, is.numeric)]))
  colnames(df)=c("value")
  df$Timesteps=1:(nticks+1)
  df$GraphName=Name
  df=df%>%select(GraphName, Timesteps,value)
  
  return(df)   
}

#########---------------------AVERAGE PROPORTION OF INFECTED (average degree of 4)------####################
###---Node 5---###
n_5 = 5
ER5=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_5_avrgdeg_4.csv")
SW5=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_5_avrgdeg_4.csv")
SF5=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_5_avrgdeg_4.csv")
SP5=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_5_avrgdeg_4.csv")
Lat5=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_5_avrgdeg_4.csv")
CG5=AvrgPropOf_CG(beta = 0.5,gamma=0.143,net_size = n_5,sim = 10,ntime = 100,avrdg = 4)

df=rbind(ER5,SW5,SF5,SP5,Lat5,CG5)
plot_obs_5_av4<- ggplot(df,aes(x=Timesteps, y=value/n_5, group=GraphName))+
  geom_line(show.legend = T,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title=" ")+
  xlab("Timesteps (days)")+ylab("Average-Prop-Infected")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_5_av4
ggsave("obs_5_av4.png", width = 24, height = 28)


###---Node 10---###
n_10= 10
ER10=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_10_avrgdeg_4.csv")
SW10=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_10_avrgdeg_4.csv")
SF10=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_10_avrgdeg_4.csv")
SP10=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_10_avrgdeg_4.csv")
Lat10=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_10_avrgdeg_4.csv")
CG10=AvrgPropOf_CG(beta = 0.5,gamma=0.143,net_size = n_10,sim = 10,ntime = 100,avrdg = 4)

df=rbind(ER10,SW10,SF10,SP10,Lat10,CG10)
plot_obs_10_av4<- ggplot(df,aes(x=Timesteps, y=value/n_10, group=GraphName))+
  geom_line(show.legend = T,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title=" ")+
  xlab("Timesteps (days)")+ylab(" ")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))
plot_obs_10_av4

combine_5_10_av4=ggarrange(plot_obs_5_av4,plot_obs_10_av4,ncol=2,nrow=1)


###---Node 20---###
n_20=20
ER20=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_20_avrgdeg_4.csv")
SW20=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_20_avrgdeg_4.csv")
SF20=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_20_avrgdeg_4.csv")
SP20=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_20_avrgdeg_4.csv")
Lat20=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_20_avrgdeg_4.csv")
CG20=AvrgPropOf_CG(beta = 0.5,gamma=0.143,net_size = n_20,sim = 10,ntime = 100,avrdg = 4)

df=rbind(ER20,SW20,SF20,SP20,Lat20,CG20)
plot_obs_20_av4<- ggplot(df,aes(x=Timesteps, y=value/n_20, group=GraphName))+
  geom_line(show.legend = T,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title=" ")+
  xlab(" ")+ylab("Average-Prop-Infected")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_20_av4

###---Node 30---###
n_30= 30
ER30=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_30_avrgdeg_4.csv")
SW30=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_30_avrgdeg_4.csv")
SF30=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_30_avrgdeg_4.csv")
SP30=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_30_avrgdeg_4.csv")
Lat30=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_30_avrgdeg_4.csv")
CG30=AvrgPropOf_CG(beta = 0.5,gamma=0.143,net_size = n_30,sim = 10,ntime = 100,avrdg = 4)

df=rbind(ER30,SW30,SF30,SP30,Lat30,CG30)
plot_obs_30_av4<- ggplot(df,aes(x=Timesteps, y=value/n_30, group=GraphName))+
  geom_line(show.legend = T,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title=" ")+
  xlab(" ")+ylab(" ")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_30_av4

combine_20_30_av4=ggarrange(plot_obs_20_av4,plot_obs_30_av4,ncol=2,nrow=1, common.legend = TRUE, legend="bottom")

###---Node 40---###
n_40 = 40
ER40=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_40_avrgdeg_4.csv")
SW40=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_40_avrgdeg_4.csv")
SF40=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_40_avrgdeg_4.csv")
SP40=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_40_avrgdeg_4.csv")
Lat40=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_40_avrgdeg_4.csv")
CG40=AvrgPropOf_CG(beta = 0.5,gamma=0.143,net_size = n_40,sim = 10,ntime = 100,avrdg = 4)

df=rbind(ER40,SW40,SF40,SP40,Lat40,CG40)
plot_obs_40_av4<- ggplot(df,aes(x=Timesteps, y=value/n_40, group=GraphName))+
  geom_line(show.legend = T,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title=" ")+
  xlab("Timesteps (days)")+ylab("Average-Prop-Infected")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))
plot_obs_40_av4

###---Node 50---###
n_50= 50
ER50=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_4.csv")
SW50=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_4.csv")
SF50=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_4.csv")
SP50=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_4.csv")
Lat50=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_4.csv")
CG50=AvrgPropOf_CG(beta = 0.5,gamma=0.143,net_size = n_50,sim = 10,ntime = 100,avrdg = 4)

df=rbind(ER50,SW50,SF50,SP50,Lat50,CG50)
plot_obs_50_av4<- ggplot(df,aes(x=Timesteps, y=value/n_50, group=GraphName))+
  geom_line(show.legend = T,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title=" ")+
  xlab("Timesteps (days)")+ylab("")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_50_av4

combine_40_50_av4=ggarrange(plot_obs_40_av4,plot_obs_50_av4,ncol=2,nrow=1)


###---Node 150---###
n_150= 150
ER150=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_4.csv")
SW150=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_4.csv")
SF150=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_4.csv")
SP150=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_4.csv")
Lat150=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_4.csv")
CG150=AvrgPropOf_CG(beta = 0.5,gamma=0.143,net_size = n_150,sim = 10,ntime = 100,avrdg = 4)

df=rbind(ER150,SW150,SF150,SP150,Lat150,CG150)
plot_obs_150_av4<- ggplot(df,aes(x=Timesteps, y=value/n_150, group=GraphName))+
  geom_line(show.legend = T,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title=" ")+
  xlab(" ")+ylab("Average-Prop-Infected")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_150_av4

###---Node 250---###
n_250 = 250
ER250=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_4.csv")
SW250=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_4.csv")
SF250=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_4.csv")
SP250=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_4.csv")
Lat250=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_4.csv")
CG250=AvrgPropOf_CG(beta = 0.5,gamma=0.143,net_size = n_250,sim = 10,ntime = 100,avrdg = 4)

df=rbind(ER250,SW250,SF250,SP250,Lat250,CG250)
plot_obs_250_av4<- ggplot(df,aes(x=Timesteps, y=value/n_250, group=GraphName))+
  geom_line(show.legend = T,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title=" ")+
  xlab(" ")+ylab(" ")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_250_av4

combine_150_250_av4=ggarrange(plot_obs_150_av4,plot_obs_250_av4,ncol=2,nrow=1)

###---Node 500---###
n_500 = 500
ER500=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_4.csv")
SW500=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_4.csv")
SF500=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_4.csv")
SP500=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_4.csv")
Lat500=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_4.csv")
CG500=AvrgPropOf_CG(beta = 0.5,gamma=0.143,net_size = n_500,sim = 10,ntime = 100,avrdg = 4)

df=rbind(ER500,SW500,SF500,SP500,Lat500,CG500)
plot_obs_500_av4<- ggplot(df,aes(x=Timesteps, y=value/n_500, group=GraphName))+
  geom_line(show.legend = T,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title=" ")+
  xlab(" ")+ylab("Average-Prop-Infected")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_500_av4

###---Node 1000---###
n_1000 = 1000
ER1000=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_4.csv")
SW1000=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_4.csv")
SF1000=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_4.csv")
SP1000=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_4.csv")
Lat1000=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_4.csv")
CG1000=AvrgPropOf_CG(beta = 0.5,gamma=0.143,net_size = n_1000,sim = 10,ntime = 100,avrdg = 4)

df=rbind(ER1000,SW1000,SF1000,SP1000,Lat1000,CG1000)
plot_obs_1000_av4<- ggplot(df,aes(x=Timesteps, y=value/n_1000, group=GraphName))+
  geom_line(show.legend = T,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title=" ")+
  xlab("")+ylab(" ")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_1000_av4

combine_500_1000_av4=ggarrange(plot_obs_500_av4,plot_obs_1000_av4,ncol=2,nrow=1,
                               common.legend = TRUE, legend="bottom")

NodesA1_av4_beta0.5_gamma0.143=ggarrange(combine_5_10_av4,combine_20_30_av4,
                                         nrow=2)
ggsave("NodesA1_av4_beta0.5_gamma0.143.png", width = 26, height = 28)
ggsave("AVRGDEGPLOT1.pdf", width = 26, height = 28)

NodesA2_av4_beta0.5_gamma0.143=ggarrange(combine_40_50_av4,combine_150_250_av4,
                                         combine_500_1000_av4,
                                         nrow=3)

ggsave("NodesA2_av4_beta0.5_gamma0.143.png", width = 26, height = 28)
ggsave("AVRGDEGPLOT2.pdf", width = 26, height = 28)


#############----Different infections (different betavalues), thus, beta=0.1,0.33,0.5, gamma=0.143-------------##########################

########------------DEGREE FOUR (4)---------------#######################
###---Node 50---###
n_50 = 50
##beta=0.1, gamma=0.143
ER50=AvrgPropOf(Name="ER",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_4.csv")
SW50=AvrgPropOf(Name="SW",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_4.csv")
SF50=AvrgPropOf(Name="SF",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_4.csv")
SP50=AvrgPropOf(Name="SP",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_4.csv")
Lat50=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_4.csv")
CG50=AvrgPropOf_CG(beta = 0.1,gamma=0.143,net_size = n_50,sim = 10,ntime = 100,avrdg = 4)

df=rbind(ER50,SW50,SF50,SP50,Lat50,CG50)
plot_obs_50_beta0.1<- ggplot(df,aes(x=Timesteps, y=value/n_50, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title=" ",subtitle = "beta:0.1")+
  xlab("Timesteps (days)")+ylab("Average-Prop-Infected")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_50_beta0.1

##beta=0.33,gamma=0.143
ER50=AvrgPropOf(Name="ER",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_4.csv")
SW50=AvrgPropOf(Name="SW",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_4.csv")
SF50=AvrgPropOf(Name="SF",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_4.csv")
SP50=AvrgPropOf(Name="SP",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_4.csv")
Lat50=AvrgPropOf(Name="Lat",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_4.csv")
CG50=AvrgPropOf_CG(beta =  0.33,gamma=0.143,net_size = n_50,sim = 10,ntime = 100,avrdg = 4)

df=rbind(ER50,SW50,SF50,SP50,Lat50,CG50)
plot_obs_50_beta0.33<- ggplot(df,aes(x=Timesteps, y=value/n_50, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="",subtitle = "beta:0.33")+
  xlab("Timesteps (days)")+ylab(" ")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_50_beta0.33


#beta=0.5,gamma=0.143
ER50=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_4.csv")
SW50=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_4.csv")
SF50=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_4.csv")
SP50=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_4.csv")
Lat50=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_4.csv")
CG50=AvrgPropOf_CG(beta = 0.5,gamma=0.143,net_size = n_50,sim = 10,ntime = 100,avrdg = 4)

df=rbind(ER50,SW50,SF50,SP50,Lat50,CG50)
plot_obs_50_beta0.5<- ggplot(df,aes(x=Timesteps, y=value/n_50, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="",subtitle = "beta:0.5")+
  xlab("Timesteps (days)")+ylab(" ")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))


plot_obs_50_beta0.5


### All Beta Plots for 50 nodes
All_betaplot_50_4=ggarrange(plot_obs_50_beta0.1,plot_obs_50_beta0.33,plot_obs_50_beta0.5,
                            nrow=1,ncol=3)


ggsave("50_all_beta_4.png", width = 30, height = 15)


###---Node 150---###
n_150= 150
##beta=0.1, gamma=0.143
ER150=AvrgPropOf(Name="ER",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_4.csv")
SW150=AvrgPropOf(Name="SW",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_4.csv")
SF150=AvrgPropOf(Name="SF",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_4.csv")
SP150=AvrgPropOf(Name="SP",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_4.csv")
Lat150=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_4.csv")
CG150=AvrgPropOf_CG(beta = 0.1,gamma=0.143,net_size = n_150,sim = 10,ntime = 100,avrdg = 4)

df=rbind(ER150,SW150,SF150,SP150,Lat150,CG150)
plot_obs_150_beta0.1<- ggplot(df,aes(x=Timesteps, y=value/n_150, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="",subtitle = "beta:0.1")+
  xlab("Timesteps (days)")+ylab("Average-Prop-Infected")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_150_beta0.1

##beta=0.33,gamma=0.143
ER150=AvrgPropOf(Name="ER",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_4.csv")
SW150=AvrgPropOf(Name="SW",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_4.csv")
SF150=AvrgPropOf(Name="SF",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_4.csv")
SP150=AvrgPropOf(Name="SP",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_4.csv")
Lat150=AvrgPropOf(Name="Lat",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_4.csv")
CG150=AvrgPropOf_CG(beta =  0.33,gamma=0.143,net_size = n_150,sim = 10,ntime = 100,avrdg = 4)

df=rbind(ER150,SW150,SF150,SP150,Lat150,CG150)
plot_obs_150_beta0.33<- ggplot(df,aes(x=Timesteps, y=value/n_150, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="",subtitle = "beta:0.33")+
  xlab("Timesteps (days)")+ylab(" ")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_150_beta0.33

#beta=0.5,gamma=0.143
ER150=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_4.csv")
SW150=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_4.csv")
SF150=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_4.csv")
SP150=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_4.csv")
Lat150=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_4.csv")
CG150=AvrgPropOf_CG(beta = 0.5,gamma=0.143,net_size = n_150,sim = 10,ntime = 100,avrdg = 4)

df=rbind(ER150,SW150,SF150,SP150,Lat150,CG150)
plot_obs_150_beta0.5<- ggplot(df,aes(x=Timesteps, y=value/n_150, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="",subtitle = "beta:0.5")+
  xlab("Timesteps (days)")+ylab("")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_150_beta0.5


### All Beta Plots for 150 nodes
All_betaplot_150_4=ggarrange(plot_obs_150_beta0.1,plot_obs_150_beta0.33,plot_obs_150_beta0.5,nrow=1,ncol=3)

###---Node 250---###
n_250= 250
##beta=0.1, gamma=0.143
ER250=AvrgPropOf(Name="ER",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_4.csv")
SW250=AvrgPropOf(Name="SW",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_4.csv")
SF250=AvrgPropOf(Name="SF",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_4.csv")
SP250=AvrgPropOf(Name="SP",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_4.csv")
Lat250=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_4.csv")
CG250=AvrgPropOf_CG(beta = 0.1,gamma=0.143,net_size = n_250,sim = 10,ntime = 100,avrdg = 4)

df=rbind(ER250,SW250,SF250,SP250,Lat250,CG250)
plot_obs_250_beta0.1<- ggplot(df,aes(x=Timesteps, y=value/n_250, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="",subtitle = "beta:0.1")+
  xlab(" ")+ylab("Average-Prop-Infected")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_250_beta0.1

##beta=0.33,gamma=0.143
ER250=AvrgPropOf(Name="ER",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_4.csv")
SW250=AvrgPropOf(Name="SW",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_4.csv")
SF250=AvrgPropOf(Name="SF",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_4.csv")
SP250=AvrgPropOf(Name="SP",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_4.csv")
Lat250=AvrgPropOf(Name="Lat",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_4.csv")
CG250=AvrgPropOf_CG(beta =  0.33,gamma=0.143,net_size = n_250,sim = 10,ntime = 100,avrdg = 4)

df=rbind(ER250,SW250,SF250,SP250,Lat250,CG250)
plot_obs_250_beta0.33<- ggplot(df,aes(x=Timesteps, y=value/n_250, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="",subtitle = "beta:0.33")+
  xlab(" ")+ylab(" ")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_250_beta0.33

#beta=0.5,gamma=0.143
ER250=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_4.csv")
SW250=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_4.csv")
SF250=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_4.csv")
SP250=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_4.csv")
Lat250=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_4.csv")
CG250=AvrgPropOf_CG(beta = 0.5,gamma=0.143,net_size = n_250,sim = 10,ntime = 100,avrdg = 4)

df=rbind(ER250,SW250,SF250,SP250,Lat250,CG250)
plot_obs_250_beta0.5<- ggplot(df,aes(x=Timesteps, y=value/n_250, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="",subtitle = "beta:0.5")+
  xlab(" ")+ylab("")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_250_beta0.5


### All Beta Plots for 250 nodes
All_betaplot_250_4=ggarrange(plot_obs_250_beta0.1,plot_obs_250_beta0.33,plot_obs_250_beta0.5,nrow=1,ncol=3)

###---Node 500---###
n_500= 500
##beta=0.1, gamma=0.143
ER500=AvrgPropOf(Name="ER",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_4.csv")
SW500=AvrgPropOf(Name="SW",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_4.csv")
SF500=AvrgPropOf(Name="SF",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_4.csv")
SP500=AvrgPropOf(Name="SP",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_4.csv")
Lat500=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_4.csv")
CG500=AvrgPropOf_CG(beta = 0.1,gamma=0.143,net_size = n_500,sim = 10,ntime = 100,avrdg = 4)

df=rbind(ER500,SW500,SF500,SP500,Lat500,CG500)
plot_obs_500_beta0.1<- ggplot(df,aes(x=Timesteps, y=value/n_500, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="",subtitle = "beta:0.1")+
  xlab(" ")+ylab("Average-Prop-Infected")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_500_beta0.1

##beta=0.33,gamma=0.143
ER500=AvrgPropOf(Name="ER",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_4.csv")
SW500=AvrgPropOf(Name="SW",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_4.csv")
SF500=AvrgPropOf(Name="SF",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_4.csv")
SP500=AvrgPropOf(Name="SP",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_4.csv")
Lat500=AvrgPropOf(Name="Lat",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_4.csv")
CG500=AvrgPropOf_CG(beta =  0.33,gamma=0.143,net_size = n_500,sim = 10,ntime = 100,avrdg = 4)

df=rbind(ER500,SW500,SF500,SP500,Lat500,CG500)
plot_obs_500_beta0.33<- ggplot(df,aes(x=Timesteps, y=value/n_500, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="",subtitle = "beta:0.33")+
  xlab(" ")+ylab(" ")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_500_beta0.33

#beta=0.5,gamma=0.143
ER500=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_4.csv")
SW500=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_4.csv")
SF500=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_4.csv")
SP500=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_4.csv")
Lat500=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_4.csv")
CG500=AvrgPropOf_CG(beta = 0.5,gamma=0.143,net_size = n_500,sim = 10,ntime = 100,avrdg = 4)

df=rbind(ER500,SW500,SF500,SP500,Lat500,CG500)
plot_obs_500_beta0.5<- ggplot(df,aes(x=Timesteps, y=value/n_500, group=GraphName))+
  geom_line(show.legend = T,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="",subtitle = "beta:0.5")+
  xlab(" ")+ylab("")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_500_beta0.5

### All Beta Plots for 500 nodes
All_betaplot_500_4=ggarrange(plot_obs_500_beta0.1,plot_obs_500_beta0.33,plot_obs_500_beta0.5,nrow=1,ncol=3,common.legend = TRUE, legend="bottom")

###---Node 1000---###
n_1000= 1000
##beta=0.1, gamma=0.143
ER1000=AvrgPropOf(Name="ER",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_4.csv")
SW1000=AvrgPropOf(Name="SW",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_4.csv")
SF1000=AvrgPropOf(Name="SF",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_4.csv")
SP1000=AvrgPropOf(Name="SP",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_4.csv")
Lat1000=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_4.csv")
CG1000=AvrgPropOf_CG(beta = 0.1,gamma=0.143,net_size = n_1000,sim = 10,ntime = 100,avrdg = 4)

df=rbind(ER1000,SW1000,SF1000,SP1000,Lat1000,CG1000)
plot_obs_1000_beta0.1<- ggplot(df,aes(x=Timesteps, y=value/n_1000, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="",subtitle = "beta:0.1")+
  xlab(" ")+ylab("Average-Prop-Infected")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_1000_beta0.1

##beta=0.33,gamma=0.143
ER1000=AvrgPropOf(Name="ER",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_4.csv")
SW1000=AvrgPropOf(Name="SW",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_4.csv")
SF1000=AvrgPropOf(Name="SF",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_4.csv")
SP1000=AvrgPropOf(Name="SP",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_4.csv")
Lat1000=AvrgPropOf(Name="Lat",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_4.csv")
CG1000=AvrgPropOf_CG(beta =  0.33,gamma=0.143,net_size = n_1000,sim = 10,ntime = 100,avrdg = 4)

df=rbind(ER1000,SW1000,SF1000,SP1000,Lat1000,CG1000)
plot_obs_1000_beta0.33<- ggplot(df,aes(x=Timesteps, y=value/n_1000, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="",subtitle = "beta:0.33")+
  xlab(" ")+ylab(" ")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_1000_beta0.33

#beta=0.5,gamma=0.143
ER1000=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_4.csv")
SW1000=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_4.csv")
SF1000=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_4.csv")
SP1000=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_4.csv")
Lat1000=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_4.csv")
CG1000=AvrgPropOf_CG(beta = 0.5,gamma=0.143,net_size = n_1000,sim = 10,ntime = 100,avrdg = 4)

df=rbind(ER1000,SW1000,SF1000,SP1000,Lat1000,CG1000)
plot_obs_1000_beta0.5<- ggplot(df,aes(x=Timesteps, y=value/n_1000, group=GraphName))+
  geom_line(show.legend = T,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="",subtitle = "beta:0.5")+
  xlab(" ")+ylab("")+scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "show",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_1000_beta0.5

All_betaplot_1000_4=ggarrange(plot_obs_1000_beta0.1,plot_obs_1000_beta0.33,
                              plot_obs_1000_beta0.5,nrow=1,ncol=3,common.legend = TRUE, legend="bottom")

###################---------COMBINED NETWORK ORDER PLOT FOR AVERAGE DEGREE 4--###############
combined_all4=ggarrange(All_betaplot_50_4,All_betaplot_150_4,All_betaplot_250_4,
                        All_betaplot_1000_4, All_betaplot_500_4,nrow = 5)

ggsave("AVRGDEGPLOT4.png", width = 20, height = 35)
ggsave("AVRGDEGPLOT4.pdf", width = 20, height = 35)

########------------DEGREE twelve (12)---------------#######################
###---Node 50---###
n = 50
##beta=0.1, gamma=0.143
ER50=AvrgPropOf(Name="ER",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_12.csv")
SW50=AvrgPropOf(Name="SW",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_12.csv")
SF50=AvrgPropOf(Name="SF",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_12.csv")
SP50=AvrgPropOf(Name="SP",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_12.csv")
Lat50=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_12.csv")
CG50=AvrgPropOf_CG(beta = 0.1,gamma=0.143,net_size = n,sim = 10,ntime = 100,avrdg = 12)

df=rbind(ER50,SW50,SF50,SP50,Lat50,CG50)
plot_obs_50_beta0.1<- ggplot(df,aes(x=Timesteps, y=value/n, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title=" ",subtitle = "beta:0.1")+
  xlab("Timesteps")+ylab("Average-Proportion")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x=element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot")

plot_obs_50_beta0.1

##beta=0.33,gamma=0.143
ER50=AvrgPropOf(Name="ER",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_12.csv")
SW50=AvrgPropOf(Name="SW",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_12.csv")
SF50=AvrgPropOf(Name="SF",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_12.csv")
SP50=AvrgPropOf(Name="SP",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_12.csv")
Lat50=AvrgPropOf(Name="Lat",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_12.csv")
CG50=AvrgPropOf_CG(beta =  0.33,gamma=0.143,net_size = n,sim = 10,ntime = 100,avrdg = 12)

df=rbind(ER50,SW50,SF50,SP50,Lat50,CG50)
plot_obs_50_beta0.33<- ggplot(df,aes(x=Timesteps, y=value/n, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="50-nodes",subtitle = "beta:0.33")+
  xlab("Timesteps")+ylab("Average-Proportion")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x=element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",plot.title = element_text(hjust = 0.6,vjust = -0.4))

plot_obs_50_beta0.33

#beta=0.5,gamma=0.143
ER50=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_12.csv")
SW50=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_12.csv")
SF50=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_12.csv")
SP50=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_12.csv")
Lat50=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_12.csv")
CG50=AvrgPropOf_CG(beta = 0.5,gamma=0.143,net_size = n,sim = 10,ntime = 100,avrdg = 12)

df=rbind(ER50,SW50,SF50,SP50,Lat50,CG50)
plot_obs_50_beta0.5<- ggplot(df,aes(x=Timesteps, y=value/n, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="",subtitle = "beta:0.5")+
  xlab("Timesteps")+ylab("Average-Proportion")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x=element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot")

plot_obs_50_beta0.5


### All Beta Plots for 50 nodes
All_betaplot_50_12=ggarrange(plot_obs_50_beta0.1,plot_obs_50_beta0.33,plot_obs_50_beta0.5,nrow=1,ncol=3)


All_betaplot_50_annotate_12=annotate_figure(All_betaplot_50_12,
                                            top = text_grob(" 50 nodes",
                                                            color = "black", face = "bold", size = 18))

ggsave("50_all_beta.png", width = 30, height = 15)

###---Node 150---###
n = 150
##beta=0.1, gamma=0.143
ER150=AvrgPropOf(Name="ER",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_12.csv")
SW150=AvrgPropOf(Name="SW",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_12.csv")
SF150=AvrgPropOf(Name="SF",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_12.csv")
SP150=AvrgPropOf(Name="SP",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_12.csv")
Lat150=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_12.csv")
CG150=AvrgPropOf_CG(beta = 0.1,gamma=0.143,net_size = n,sim = 10,ntime = 100,avrdg = 12)

df=rbind(ER150,SW150,SF150,SP150,Lat150,CG150)
plot_obs_150_beta0.1<- ggplot(df,aes(x=Timesteps, y=value/n, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title=" ")+
  xlab("Timesteps")+ylab("Average-Proportion")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x=element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot")

plot_obs_150_beta0.1

##beta=0.33,gamma=0.143
ER150=AvrgPropOf(Name="ER",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_12.csv")
SW150=AvrgPropOf(Name="SW",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_12.csv")
SF150=AvrgPropOf(Name="SF",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_12.csv")
SP150=AvrgPropOf(Name="SP",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_12.csv")
Lat150=AvrgPropOf(Name="Lat",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_12.csv")
CG150=AvrgPropOf_CG(beta =  0.33,gamma=0.143,net_size = n,sim = 10,ntime = 100,avrdg = 12)

df=rbind(ER150,SW150,SF150,SP150,Lat150,CG150)
plot_obs_150_beta0.33<- ggplot(df,aes(x=Timesteps, y=value/n, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="150-nodes")+
  xlab("Timesteps")+ylab("Average-Proportion")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x=element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",plot.title = element_text(hjust = 0.6,vjust = -0.4))

plot_obs_150_beta0.33

#beta=0.5,gamma=0.143
ER150=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_12.csv")
SW150=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_12.csv")
SF150=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_12.csv")
SP150=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_12.csv")
Lat150=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_12.csv")
CG150=AvrgPropOf_CG(beta = 0.5,gamma=0.143,net_size = n,sim = 10,ntime = 100,avrdg = 12)

df=rbind(ER150,SW150,SF150,SP150,Lat150,CG150)
plot_obs_150_beta0.5<- ggplot(df,aes(x=Timesteps, y=value/n, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="")+
  xlab("Timesteps")+ylab("Average-Proportion")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x=element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot")


plot_obs_150_beta0.5

### All Beta Plots for 150 nodes
All_betaplot_150_12=ggarrange(plot_obs_150_beta0.1,plot_obs_150_beta0.33,plot_obs_150_beta0.5,nrow=1,ncol=3)

All_betaplot_150_annotate_12=annotate_figure(All_betaplot_150_12,
                                             top = text_grob(" 150 nodes",
                                                             color = "black", face = "bold", size = 18))
ggsave("150_all_beta.png", width = 30, height = 15)


###---Node 250---###
n = 250
##beta=0.1, gamma=0.143
ER250=AvrgPropOf(Name="ER",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_12.csv")
SW250=AvrgPropOf(Name="SW",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_12.csv")
SF250=AvrgPropOf(Name="SF",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_12.csv")
SP250=AvrgPropOf(Name="SP",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_12.csv")
Lat250=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_12.csv")
CG250=AvrgPropOf_CG(beta = 0.1,gamma=0.143,net_size = n,sim = 10,ntime = 100,avrdg = 12)

df=rbind(ER250,SW250,SF250,SP250,Lat250,CG250)
plot_obs_250_beta0.1<- ggplot(df,aes(x=Timesteps, y=value/n, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="")+
  xlab("Timesteps")+ylab("Average-Proportion")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x=element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot")


plot_obs_250_beta0.1

n=250
ER250=AvrgPropOf(Name="ER",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_12.csv")
SW250=AvrgPropOf(Name="SW",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_12.csv")
SF250=AvrgPropOf(Name="SF",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_12.csv")
SP250=AvrgPropOf(Name="SP",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_12.csv")
Lat250=AvrgPropOf(Name="Lat",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_12.csv")
CG250=AvrgPropOf_CG(beta =  0.33,gamma=0.143,net_size = n,sim = 10,ntime = 100,avrdg = 12)

df=rbind(ER250,SW250,SF250,SP250,Lat250,CG250)
plot_obs_250_beta0.33<- ggplot(df,aes(x=Timesteps, y=value/n, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="250-nodes")+
  xlab("Timesteps")+ylab("Average-Proportion")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x=element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",plot.title = element_text(hjust = 0.6,vjust = -0.4))

plot_obs_250_beta0.33

#beta=0.5,gamma=0.143
n=250
ER250=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_12.csv")
SW250=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_12.csv")
SF250=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_12.csv")
SP250=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_12.csv")
Lat250=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_12.csv")
CG2250=AvrgPropOf_CG(beta = 0.5,gamma=0.143,net_size = n,sim = 10,ntime = 100,avrdg = 12)

df=rbind(ER250,SW250,SF250,SP250,Lat250,CG250)
plot_obs_250_beta0.5<- ggplot(df,aes(x=Timesteps, y=value/n, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="")+
  xlab("Timesteps")+ylab("Average-Proportion")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x=element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot")

plot_obs_250_beta0.5

### All Beta Plots for 250 nodes
All_betaplot_250_12=ggarrange(plot_obs_250_beta0.1,plot_obs_250_beta0.33,plot_obs_250_beta0.5,nrow=1,ncol=3)

All_betaplot_250_annotate_12=annotate_figure(All_betaplot_250_12,
                                             top = text_grob(" 250 nodes",
                                                             color = "black", face = "bold", size = 18))
ggsave("250_all_beta.png", width = 30, height = 15)

###############---Node 500---################
n = 500
##beta=0.1, gamma=0.143
ER500=AvrgPropOf(Name="ER",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_12.csv")
SW500=AvrgPropOf(Name="SW",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_12.csv")
SF500=AvrgPropOf(Name="SF",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_12.csv")
SP500=AvrgPropOf(Name="SP",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_12.csv")
Lat500=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_12.csv")
CG500=AvrgPropOf_CG(beta = 0.1,gamma=0.143,net_size = n,sim = 10,ntime = 100,avrdg = 12)

df=rbind(ER500,SW500,SF500,SP500,Lat500,CG500)
plot_obs_500_beta0.1<- ggplot(df,aes(x=Timesteps, y=value/n, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="")+
  xlab("Timesteps")+ylab("Average-Proportion")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x=element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot")


plot_obs_500_beta0.1

##beta=0.33,gamma=0.143
n=500
ER500=AvrgPropOf(Name="ER",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_12.csv")
SW500=AvrgPropOf(Name="SW",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_12.csv")
SF500=AvrgPropOf(Name="SF",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_12.csv")
SP500=AvrgPropOf(Name="SP",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_12.csv")
Lat500=AvrgPropOf(Name="Lat",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_12.csv")
CG500=AvrgPropOf_CG(beta =  0.33,gamma=0.143,net_size = n,sim = 10,ntime = 100,avrdg = 12)

df=rbind(ER500,SW500,SF500,SP500,Lat500,CG500)
plot_obs_500_beta0.33<- ggplot(df,aes(x=Timesteps, y=value/n, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="500-nodes")+
  xlab("Timesteps")+ylab("Average-Proportion")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x=element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",plot.title = element_text(hjust = 0.6,vjust = -0.4))

plot_obs_500_beta0.33

#beta=0.5,gamma=0.143
ER500=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_12.csv")
SW500=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_12.csv")
SF500=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_12.csv")
SP500=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_12.csv")
Lat500=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_12.csv")
CG500=AvrgPropOf_CG(beta = 0.5,gamma=0.143,net_size = n,sim = 10,ntime = 100,avrdg = 12)

df=rbind(ER500,SW500,SF500,SP500,Lat500,CG500)
plot_obs_500_beta0.5<- ggplot(df,aes(x=Timesteps, y=value/n, group=GraphName))+
  geom_line(show.legend = T,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="")+
  xlab("Timesteps")+ylab("Average-Proportion")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x=element_text(size = 26),
        legend.position = "none", 
        plot.title.position = "plot")

plot_obs_500_beta0.5

### All Beta Plots for 50 nodes
All_betaplot_500_12=ggarrange(plot_obs_500_beta0.1,
                              plot_obs_500_beta0.33,plot_obs_500_beta0.5,nrow=1,ncol=3)

All_betaplot_500_annotate_12=annotate_figure(All_betaplot_500_12,
                                             top = text_grob(" 500 nodes",
                                                             color = "black", face = "bold", size = 18))

ggsave("500_all_beta.png", width = 30, height = 15)

###############---Node 1000---################
n = 1000
##beta=0.1, gamma=0.143
ER1000=AvrgPropOf(Name="ER",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_12.csv")

SW1000=AvrgPropOf(Name="SW",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_12.csv")
SF1000=AvrgPropOf(Name="SF",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_12.csv")
SP1000=AvrgPropOf(Name="SP",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_12.csv")
Lat1000=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_12.csv")
CG1000=AvrgPropOf_CG(beta = 0.1,gamma=0.143,net_size = n,sim = 10,ntime = 100,avrdg = 12)

df=rbind(ER1000,SW1000,SF1000,SP1000,Lat1000,CG1000)
plot_obs_1000_beta0.1<- ggplot(df,aes(x=Timesteps, y=value/n, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="")+
  xlab("Timesteps")+ylab("Average-Proportion")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x=element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot")

plot_obs_1000_beta0.1

##beta=0.33,gamma=0.143
n=1000
ER1000=AvrgPropOf(Name="ER",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_12.csv")
SW1000=AvrgPropOf(Name="SW",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_12.csv")
SF1000=AvrgPropOf(Name="SF",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_12.csv")
SP1000=AvrgPropOf(Name="SP",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_12.csv")
Lat1000=AvrgPropOf(Name="Lat",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_12.csv")
CG1000=AvrgPropOf_CG(beta =  0.33,gamma=0.143,net_size = n,sim = 10,ntime = 100,avrdg = 12)

df=rbind(ER1000,SW1000,SF1000,SP1000,Lat1000,CG1000)
plot_obs_1000_beta0.33<- ggplot(df,aes(x=Timesteps, y=value/n, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="1000-nodes")+
  xlab("Timesteps")+ylab("Average-Proportion")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x=element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",plot.title = element_text(hjust = 0.6,vjust = -0.4))

plot_obs_1000_beta0.33

#beta=0.5,gamma=0.143
ER1000=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_12.csv")
SW1000=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_12.csv")
SF1000=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_12.csv")
SP1000=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_12.csv")
Lat1000=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_12.csv")
CG1000=AvrgPropOf_CG(beta = 0.5,gamma=0.143,net_size = n,sim = 10,ntime = 100,avrdg = 12)

df=rbind(ER1000,SW1000,SF1000,SP1000,Lat1000,CG1000)
plot_obs_1000_beta0.5<- ggplot(df,aes(x=Timesteps, y=value/n, group=GraphName))+
  geom_line(show.legend = T,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="")+
  xlab("Timesteps")+ylab("Average-Proportion")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x=element_text(size = 26),
        legend.position = "none", 
        plot.title.position = "plot")

plot_obs_1000_beta0.5

### All Beta Plots for 50 nodes
All_betaplot_1000_12=ggarrange(plot_obs_1000_beta0.1,
                               plot_obs_1000_beta0.33,plot_obs_1000_beta0.5,nrow=1,ncol=3,common.legend = TRUE, legend="bottom")

All_betaplot_1000_annotate_12=annotate_figure(All_betaplot_1000_12,
                                              top = text_grob(" 1000 nodes",
                                                              color = "black", face = "bold", size = 18))

ggsave("1000_all_beta.png", width = 30, height = 15)

###################---------COMBINED NETWORK ORDER PLOT FOR AVERAGE DEGREE 12--###############
combined_all12=ggarrange(All_betaplot_50_12,All_betaplot_150_12,All_betaplot_250_12,
                         All_betaplot_500_12,All_betaplot_1000_12,nrow = 5)
combined_all12

ggsave("combined_all_beta_avrg_12.png", width = 20, height = 35)

ggsave("AVRGDEGPLOT12.pdf", width = 20, height = 35)


########------------DEGREE twelve (24)---------------#######################
###---Node 50---###
n_50 = 50
##beta=0.1, gamma=0.143
ER50=AvrgPropOf(Name="ER",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_24.csv")
SW50=AvrgPropOf(Name="SW",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_24.csv")
SF50=AvrgPropOf(Name="SF",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_24.csv")
SP50=AvrgPropOf(Name="SP",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_24.csv")
Lat50=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_24.csv")

CG50=AvrgPropOf_CG(beta = 0.1,gamma= 0.143,net_size = n_50,sim = 10,ntime = 100,avrdg = 24)

df=rbind(ER50,SW50,SF50,SP50,Lat50,CG50)
plot_obs_50_beta0.1<- ggplot(df,aes(x=Timesteps, y=value/n_50, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title=" ",subtitle = "beta:0.1")+
  xlab("Timesteps (days)")+ylab("Average-Prop-Infected")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_50_beta0.1

##beta=0.33,gamma= 0.143
ER50=AvrgPropOf(Name="ER",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_24.csv")
SW50=AvrgPropOf(Name="SW",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_24.csv")
SF50=AvrgPropOf(Name="SF",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_24.csv")
SP50=AvrgPropOf(Name="SP",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_24.csv")
Lat50=AvrgPropOf(Name="Lat",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_24.csv")
CG50=AvrgPropOf_CG(beta =  0.33,gamma= 0.143,net_size = n_50,sim = 10,ntime = 100,avrdg = 24)

df=rbind(ER50,SW50,SF50,SP50,Lat50,CG50)
plot_obs_50_beta0.33<- ggplot(df,aes(x=Timesteps, y=value/n_50, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="",subtitle = "beta:0.33")+
  xlab("Timesteps (days)")+ylab(" ")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_50_beta0.33

#beta=0.5,gamma= 0.143
ER50=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_24.csv")
SW50=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_24.csv")
SF50=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_24.csv")
SP50=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_24.csv")
Lat50=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_50_avrgdeg_24.csv")
CG50=AvrgPropOf_CG(beta = 0.5,gamma= 0.143,net_size = n_50,sim = 10,ntime = 100,avrdg = 24)

df=rbind(ER50,SW50,SF50,SP50,Lat50,CG50)
plot_obs_50_beta0.5<- ggplot(df,aes(x=Timesteps, y=value/n_50, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="",subtitle = "beta:0.5")+
  xlab("Timesteps (days)")+ylab(" ")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))


plot_obs_50_beta0.5


### All Beta Plots for 50 nodes
All_betaplot_50_24=ggarrange(plot_obs_50_beta0.1,plot_obs_50_beta0.33,plot_obs_50_beta0.5,
                             nrow=1,ncol=3)

ggsave("50_all_beta_24.png", width = 30, height = 15)


###---Node 150---###
n_150= 150
##beta=0.1, gamma=0.143
ER150=AvrgPropOf(Name="ER",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_24.csv")
SW150=AvrgPropOf(Name="SW",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_24.csv")
SF150=AvrgPropOf(Name="SF",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_24.csv")
SP150=AvrgPropOf(Name="SP",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_24.csv")
Lat150=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_24.csv")
CG150=AvrgPropOf_CG(beta = 0.1,gamma= 0.143,net_size = n_150,sim = 10,ntime = 100,avrdg = 24)

df=rbind(ER150,SW150,SF150,SP150,Lat150,CG150)
plot_obs_150_beta0.1<- ggplot(df,aes(x=Timesteps, y=value/n_150, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="",subtitle = "beta:0.1")+
  xlab("Timesteps (days)")+ylab("Average-Prop-Infected")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_150_beta0.1

##beta=0.33,gamma= 0.143
ER150=AvrgPropOf(Name="ER",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_24.csv")
SW150=AvrgPropOf(Name="SW",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_24.csv")
SF150=AvrgPropOf(Name="SF",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_24.csv")
SP150=AvrgPropOf(Name="SP",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_24.csv")
Lat150=AvrgPropOf(Name="Lat",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_24.csv")
CG150=AvrgPropOf_CG(beta =  0.33,gamma= 0.143,net_size = n_150,sim = 10,ntime = 100,avrdg = 24)

df=rbind(ER150,SW150,SF150,SP150,Lat150,CG150)
plot_obs_150_beta0.33<- ggplot(df,aes(x=Timesteps, y=value/n_150, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="",subtitle = "beta:0.33")+
  xlab("Timesteps (days)")+ylab(" ")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_150_beta0.33

#beta=0.5,gamma= 0.143
ER150=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_24.csv")
SW150=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_24.csv")
SF150=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_24.csv")
SP150=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_24.csv")
Lat150=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_150_avrgdeg_24.csv")
CG150=AvrgPropOf_CG(beta = 0.5,gamma= 0.143,net_size = n_150,sim = 10,ntime = 100,avrdg = 24)

df=rbind(ER150,SW150,SF150,SP150,Lat150,CG150)
plot_obs_150_beta0.5<- ggplot(df,aes(x=Timesteps, y=value/n_150, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="",subtitle = "beta:0.5")+
  xlab("Timesteps (days)")+ylab("")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_150_beta0.5


### All Beta Plots for 150 nodes
All_betaplot_150_24=ggarrange(plot_obs_150_beta0.1,plot_obs_150_beta0.33,plot_obs_150_beta0.5,nrow=1,ncol=3)

###---Node 250---###
n_250= 250
##beta=0.1, gamma=0.143
ER250=AvrgPropOf(Name="ER",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_24.csv")
SW250=AvrgPropOf(Name="SW",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_24.csv")
SF250=AvrgPropOf(Name="SF",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_24.csv")
SP250=AvrgPropOf(Name="SP",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_24.csv")
Lat250=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_24.csv")
CG250=AvrgPropOf_CG(beta = 0.1,gamma= 0.143,net_size = n_250,sim = 10,ntime = 100,avrdg = 24)

df=rbind(ER250,SW250,SF250,SP250,Lat250,CG250)
plot_obs_250_beta0.1<- ggplot(df,aes(x=Timesteps, y=value/n_250, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="",subtitle = "beta:0.1")+
  xlab(" ")+ylab("Average-Prop-Infected")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_250_beta0.1

##beta=0.33,gamma= 0.143
ER250=AvrgPropOf(Name="ER",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_24.csv")
SW250=AvrgPropOf(Name="SW",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_24.csv")
SF250=AvrgPropOf(Name="SF",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_24.csv")
SP250=AvrgPropOf(Name="SP",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_24.csv")
Lat250=AvrgPropOf(Name="Lat",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_24.csv")
CG250=AvrgPropOf_CG(beta =  0.33,gamma= 0.143,net_size = n_250,sim = 10,ntime = 100,avrdg = 24)

df=rbind(ER250,SW250,SF250,SP250,Lat250,CG250)
plot_obs_250_beta0.33<- ggplot(df,aes(x=Timesteps, y=value/n_250, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="",subtitle = "beta:0.33")+
  xlab(" ")+ylab(" ")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_250_beta0.33

#beta=0.5,gamma= 0.143
ER250=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_24.csv")
SW250=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_24.csv")
SF250=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_24.csv")
SP250=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_24.csv")
Lat250=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_250_avrgdeg_24.csv")
CG250=AvrgPropOf_CG(beta = 0.5,gamma= 0.143,net_size = n_250,sim = 10,ntime = 100,avrdg = 24)

df=rbind(ER250,SW250,SF250,SP250,Lat250,CG250)
plot_obs_250_beta0.5<- ggplot(df,aes(x=Timesteps, y=value/n_250, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="",subtitle = "beta:0.5")+
  xlab(" ")+ylab("")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_250_beta0.5


### All Beta Plots for 250 nodes
All_betaplot_250_24=ggarrange(plot_obs_250_beta0.1,plot_obs_250_beta0.33,plot_obs_250_beta0.5,nrow=1,ncol=3)

###---Node 500---###
n_500= 500
##beta=0.1, gamma=0.143
ER500=AvrgPropOf(Name="ER",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_24.csv")
SW500=AvrgPropOf(Name="SW",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_24.csv")
SF500=AvrgPropOf(Name="SF",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_24.csv")
SP500=AvrgPropOf(Name="SP",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_24.csv")
Lat500=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_24.csv")
CG500=AvrgPropOf_CG(beta = 0.1,gamma= 0.143,net_size = n_500,sim = 10,ntime = 100,avrdg = 24)

df=rbind(ER500,SW500,SF500,SP500,Lat500,CG500)
plot_obs_500_beta0.1<- ggplot(df,aes(x=Timesteps, y=value/n_500, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="",subtitle = "beta:0.1")+
  xlab(" ")+ylab("Average-Prop-Infected")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_500_beta0.1

##beta=0.33,gamma= 0.143
ER500=AvrgPropOf(Name="ER",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_24.csv")
SW500=AvrgPropOf(Name="SW",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_24.csv")
SF500=AvrgPropOf(Name="SF",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_24.csv")
SP500=AvrgPropOf(Name="SP",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_24.csv")
Lat500=AvrgPropOf(Name="Lat",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_24.csv")
CG500=AvrgPropOf_CG(beta =  0.33,gamma= 0.143,net_size = n_500,sim = 10,ntime = 100,avrdg = 24)

df=rbind(ER500,SW500,SF500,SP500,Lat500,CG500)
plot_obs_500_beta0.33<- ggplot(df,aes(x=Timesteps, y=value/n_500, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="",subtitle = "beta:0.33")+
  xlab(" ")+ylab(" ")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_500_beta0.33

#beta=0.5,gamma= 0.143
ER500=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_24.csv")
SW500=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_24.csv")
SF500=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_24.csv")
SP500=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_24.csv")
Lat500=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_500_avrgdeg_24.csv")
CG500=AvrgPropOf_CG(beta = 0.5,gamma= 0.143,net_size = n_500,sim = 10,ntime = 100,avrdg = 24)

df=rbind(ER500,SW500,SF500,SP500,Lat500,CG500)
plot_obs_500_beta0.5<- ggplot(df,aes(x=Timesteps, y=value/n_500, group=GraphName))+
  geom_line(show.legend = T,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="",subtitle = "beta:0.5")+
  xlab(" ")+ylab("")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26,face = "italic"),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",
        plot.tag.position = c(0.1, 0.98))

plot_obs_500_beta0.5


### All Beta Plots for 500 nodes
All_betaplot_500_24=ggarrange(plot_obs_500_beta0.1,plot_obs_500_beta0.33,plot_obs_500_beta0.5,nrow=1,ncol=3,common.legend = TRUE, legend="bottom")

###############---Node 1000---################
n_1000= 1000
##beta=0.1, gamma=0.143
ER1000=AvrgPropOf(Name="ER",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_24.csv")
SW1000=AvrgPropOf(Name="SW",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_24.csv")
SF1000=AvrgPropOf(Name="SF",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_24.csv")
SP1000=AvrgPropOf(Name="SP",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_24.csv")
Lat1000=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_24.csv")
CG1000=AvrgPropOf_CG(beta = 0.1,gamma=0.143,net_size = n_1000,sim = 10,ntime = 100,avrdg = 24)

df=rbind(ER1000,SW1000,SF1000,SP1000,Lat1000,CG1000)
plot_obs_1000_beta0.1<- ggplot(df,aes(x=Timesteps, y=value/n_1000, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="")+
  xlab("Timesteps")+ylab("Average-Proportion")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x=element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot")


plot_obs_1000_beta0.1

##beta=0.33,gamma=0.143
ER1000=AvrgPropOf(Name="ER",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_24.csv")
SW1000=AvrgPropOf(Name="SW",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_24.csv")
SF1000=AvrgPropOf(Name="SF",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_24.csv")
SP1000=AvrgPropOf(Name="SP",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_24.csv")
Lat1000=AvrgPropOf(Name="Lat",ID=1,betaVal =  0.33,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_24.csv")
CG1000=AvrgPropOf_CG(beta =  0.33,gamma=0.143,net_size = n_1000,sim = 10,ntime = 100,avrdg = 24)

df=rbind(ER1000,SW1000,SF1000,SP1000,Lat1000,CG1000)
plot_obs_1000_beta0.33<- ggplot(df,aes(x=Timesteps, y=value/n_1000, group=GraphName))+
  geom_line(show.legend = F,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="1000-nodes")+
  xlab("Timesteps")+ylab("Average-Proportion")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x=element_text(size = 26),
        legend.position = "none",
        plot.title.position = "plot",plot.title = element_text(hjust = 0.6,vjust = -0.4))

plot_obs_1000_beta0.33

#beta=0.5,gamma=0.143
ER1000=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_24.csv")
SW1000=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_24.csv")
SF1000=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_24.csv")
SP1000=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_24.csv")
Lat1000=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,Data = "observed_networks_1000_avrgdeg_24.csv")
CG1000=AvrgPropOf_CG(beta = 0.5,gamma=0.143,net_size = n_1000,sim = 10,ntime = 100,avrdg = 24)

df=rbind(ER1000,SW1000,SF1000,SP1000,Lat1000,CG1000)
plot_obs_1000_beta0.5<- ggplot(df,aes(x=Timesteps, y=value/n_1000, group=GraphName))+
  geom_line(show.legend = T,aes(linetype=GraphName, color=GraphName),size = 2)+
  theme_classic()+labs(title="")+
  xlab("Timesteps")+ylab("Average-Proportion")+
  scale_y_continuous(limits = c(0,100))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x=element_text(size = 26),
        legend.position = "none", 
        plot.title.position = "plot")

plot_obs_1000_beta0.5

### All Beta Plots for 50 nodes
All_betaplot_1000_24=ggarrange(plot_obs_1000_beta0.1,
                               plot_obs_1000_beta0.33,plot_obs_1000_beta0.5,nrow=1,ncol=3,common.legend = TRUE, legend="bottom")

All_betaplot_1000_annotate_24=annotate_figure(All_betaplot_1000_24,
                                              top = text_grob(" 1000 nodes",
                                                              color = "black", face = "bold", size = 18))

ggsave("1000_all_beta.png", width = 20, height = 22)

###################---------COMBINED NETWORK ORDER PLOT FOR AVERAGE DEGREE 24--###############
combined_all24=ggarrange(All_betaplot_50_24,All_betaplot_150_24,All_betaplot_250_24,
                         All_betaplot_500_24,nrow = 4)
combined_all24
ggsave("AVRGDEGPLOT24.png", width = 20, height = 35)
ggsave("AVRGDEGPLOT24_.pdf", width = 20, height = 35)



##############---------EPIDEMIC MEASURES----------------------###############################
#############----Outbreak duration, Time to Infection peak and Maximum infected Individuals---##############

###########----------Box Plot of Average degree of 4------------#############
data4_1=as.data.frame(read.csv("observed_networks_50_avrgdeg_4_newdata.csv"))
data4_1$net_size=50

data4_2=as.data.frame(read.csv("observed_networks_150_avrgdeg_4_newdata.csv"))
data4_2$net_size=150

data4_3=as.data.frame(read.csv("observed_networks_250_avrgdeg_4_newdata.csv"))
data4_3$net_size=250

data4_4=as.data.frame(read.csv("observed_networks_500_avrgdeg_4_newdata.csv"))
data4_4$net_size=500


df_4=rbind(data4_1,data4_2,data4_3,data4_4)

write.csv(df_4,'C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi/network_size_50_150_250_500_avg4_newdata.csv')

D4="network_size_50_150_250_500_avg4_newdata.csv"
Data4=read.csv(D4,header = T, sep = ",")
as.factor(Data4$net_size)
#df=Data%>%filter(GraphID==1)
# m1<- subset(m1, beta != 0.1) 
# m1<- subset(m1, beta != 1.2) 
m4=filter(Data4,beta %in% c(0.01,0.1,0.33,0.5), gamma %in% c(0.143,0.25))

m4_1=m4%>%select(GraphName,beta,gamma,TimetoMax.i,net_size,OutbreakDuration,MaxNumOf.i)
# 
# str1=paste("\U03B2")
# str2=paste("\u0393")
colnames(m4_1)[2:7]=c("B","G","TimetoMaxInfected","size","OutbreakDuration","MaxInfected")



TimetoMaxInf_4=m4_1%>%gather(key, value = val,TimetoMaxInfected)%>%
  ggplot(aes(x = key, y = val, group_by=GraphName, color = GraphName)) +
  geom_boxplot(show.legend = T) +ylab("Timesteps (days)")+labs(title = "Average-degree-4")+
  facet_grid(B~G~size,labeller = label_both)+
  scale_y_continuous(limits = c(0,50))+
  theme(text = element_text(size =34),
        #strip.text.x = element_blank(),
        # axis.text.x = element_blank(),
        # axis.ticks.x = element_blank(),
        axis.title = element_text(size = 34),
        axis.text.y = element_text(size = 34),
        axis.text.x=element_text(size = 34),
        legend.position = "bottom",
        plot.title.position = "plot",plot.title = element_text(hjust = 0.6,vjust = -0.6))
TimetoMaxInf_4+theme_bw()

ggsave("TimetoMaxInfected_4_newdata.png", width = 24, height = 38)
ggsave("TimetoMaxInfected_4_newdata.pdf", width =24, height = 38)

OutDuratn_4=m4_1%>%gather(key, value = val,OutbreakDuration)%>%
  ggplot(aes(x = key, y = val, group_by=GraphName, color = GraphName)) +
  geom_boxplot() +ylab("Timesteps (days)")+xlab(" ")+labs(title = "Average-degree-4")+
  facet_grid(B~G~size,labeller = label_both)+
  scale_y_continuous(limits = c(0,50))+
  theme(text = element_text(size = 34),
        axis.title = element_text(size = 34),
        axis.text.y = element_text(size = 34),
        axis.text.x=element_text(size = 34),
        legend.position = "none",
        plot.title.position = "plot",plot.title = element_text(hjust = 0.6,vjust = -0.6))
OutDuratn_4+theme_bw()


ggsave("OutbreakDuration_4_newdata.png", width = 24, height = 38)
ggsave("OutbreakDuration_4_newdata.pdf", width = 24, height = 38)


MaxInf_4=m4_1%>%gather(key, value = val,MaxInfected)%>%
  ggplot(aes(x = key, y = val, group_by=GraphName, color = GraphName)) +
  geom_boxplot() +ylab(" Number-of-Infected-Individuals")+labs(title = "Average-degree-4")+
  facet_grid(B~G~size,labeller = label_both)+
  #scale_y_continuous(labels = scales::percent)+
  scale_y_continuous(limits=c(0,100))+
  theme(text = element_text(size = 34),
        axis.title = element_text(size = 34),
        axis.text.y = element_text(size = 34),
        axis.text.x=element_text(size = 34),
        legend.position = "none",
        plot.title.position = "plot",plot.title = element_text(hjust = 0.6,vjust = -0.6))

MaxInf_4+theme_bw()


ggsave("MaxInfected_4_newdata.png", width = 24, height = 38)

ggsave("MaxInfected_4_newdata.pdf", width = 24, height = 38)

##############################----------Box Plot of Average degree of 12------------#################
data12_1=as.data.frame(read.csv("observed_networks_50_avrgdeg_12_newdata.csv"))
data12_1$net_size=50

data12_2=as.data.frame(read.csv("observed_networks_150_avrgdeg_12_newdata.csv"))
data12_2$net_size=150

data12_3=as.data.frame(read.csv("observed_networks_250_avrgdeg_12_newdata.csv"))
data12_3$net_size=250

data12_4=as.data.frame(read.csv("observed_networks_500_avrgdeg_12_newdata.csv"))
data12_4$net_size=500

df_12=rbind(data12_1,data12_2,data12_3,data12_4)
write.csv(df_12,'C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi/network_size_50_150_250_500_avg12_newdata.csv')

D12="network_size_50_150_250_500_avg12_newdata.csv"
Data12=read.csv(D12,header = T, sep = ",")
as.factor(Data12$net_size)

#df=Data%>%filter(GraphID==1)
# m1<- subset(m1, beta != 0.1) 
# m1<- subset(m1, beta != 1.2) 
m12=filter(Data12,beta %in% c(0.01,0.1,0.33,0.5), gamma %in% c(0.143,0.25))

m12_1=m12%>%select(GraphName,beta,gamma,TimetoMax.i,net_size,OutbreakDuration,MaxNumOf.i)

colnames(m12_1)[2:7]=c("B","G","TimetoMaxInfected","size","OutbreakDuration","MaxInfected")

TimetoMaxInf_12=m12_1%>%gather(key, value = val,TimetoMaxInfected)%>%
  ggplot(aes(x = key, y = val, group_by=GraphName, color = GraphName)) +
  geom_boxplot(show.legend = T) +ylab("Timesteps (days)")+labs(title = "Average-degree-12")+
  facet_grid(B~G~size,labeller = label_both)+
  scale_y_continuous(limits = c(0,50))+
  theme(text = element_text(size = 34),
        axis.title = element_text(size = 34),
        axis.text.y = element_text(size = 34),
        axis.text.x=element_text(size = 34),
        legend.position = "bottom",
        plot.title.position = "plot",plot.title = element_text(hjust = 0.6,vjust = -0.6))
TimetoMaxInf_12+theme_bw()

ggsave("TimetoMaxInfected_12_newdata.png", width = 24, height = 38)

ggsave("TimetoMaxInfected_12_newdata.pdf", width = 24, height = 38)


OutDuratn_12=m12_1%>%gather(key, value = val,OutbreakDuration)%>%
  ggplot(aes(x = key, y = val, group_by=GraphName, color = GraphName)) +
  geom_boxplot(show.legend = T) +ylab("Timesteps (days)")+labs(title = "Average-degree-12")+
  facet_grid(B~G~size,labeller = label_both)+
  scale_y_continuous(limits = c(0,50))+
  theme(text = element_text(size = 34),
        axis.title = element_text(size = 34),
        axis.text.y = element_text(size = 34),
        axis.text.x=element_text(size = 34),
        legend.position = "bottom",
        plot.title.position = "plot",plot.title = element_text(hjust = 0.6,vjust = -0.6))
OutDuratn_12+theme_bw()

ggsave("OutbreakDuration_12_newdata.png", width = 24, height = 38)
ggsave("OutbreakDuration_12_newdata.pdf", width = 24, height = 38)


MaxInf_12=m12_1%>%gather(key, value = val,MaxInfected)%>%
  ggplot(aes(x = key, y = val, group_by=GraphName, color = GraphName)) +
  geom_boxplot(show.legend = T) +ylab("Number-of-Infected-Individuals")+labs(title = "Average-degree-12")+
  facet_grid(B~G~size,labeller = label_both)+
 # scale_y_continuous(labels = scales::percent)+
  scale_y_continuous(limits = c(0,100))+
  theme(text = element_text(size = 34),
        axis.title = element_text(size = 34),
        axis.text.y = element_text(size = 34),
        axis.text.x=element_text(size = 34),
        legend.position = "bottom",
        plot.title.position = "plot",plot.title = element_text(hjust = 0.6,vjust = -0.6))
MaxInf_12+theme_bw()


ggsave("MaxInfected_12_newdata.png", width = 24, height = 38 )
ggsave("MaxInfected_12_newdata.pdf", width = 24, height = 38)

##############################----------Box Plot of Average degree of 24------------#################
data24_1=as.data.frame(read.csv("observed_networks_50_avrgdeg_24_newdata.csv"))
data24_1$net_size=50

data24_2=as.data.frame(read.csv("observed_networks_150_avrgdeg_24_newdata.csv"))
data24_2$net_size=150

data24_3=as.data.frame(read.csv("observed_networks_250_avrgdeg_24_newdata.csv"))
data24_3$net_size=250

data24_4=as.data.frame(read.csv("observed_networks_500_avrgdeg_24_newdata.csv"))
data24_4$net_size=500

df_24=rbind(data24_1,data24_2,data24_3,data24_4)
write.csv(df_24,'C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi/network_size_50_150_250_500_avg24_newdata.csv')

D24="network_size_50_150_250_500_avg24_newdata.csv"
Data24=read.csv(D24,header = T, sep = ",")
as.factor(Data24$net_size)
#df=Data%>%filter(GraphID==1)
# m1<- subset(m1, beta != 0.1) 
# m1<- subset(m1, beta != 1.2) 
m24=filter(Data24,beta %in% c(0.01,0.1,0.33,0.5), gamma %in% c(0.143,0.25))

m1_24=m24%>%select(GraphName,beta,gamma,TimetoMax.i,net_size,OutbreakDuration,MaxNumOf.i)

colnames(m1_24)[2:7]=c("B","G","TimetoMaxInfected","size","OutbreakDuration","MaxInfected")

TimetoMaxInf_24=m1_24%>%gather(key, value = val,TimetoMaxInfected)%>%
  ggplot(aes(x = key, y = val, group_by=GraphName, color = GraphName)) +
  geom_boxplot(show.legend = T) +ylab("Timesteps (days)")+labs(title = "Average-degree-24")+
  facet_grid(B~G~size,labeller = label_both)+
  scale_y_continuous(limits = c(0,50))+
  theme(text = element_text(size = 34),
        axis.title = element_text(size = 34),
        axis.text.y = element_text(size = 34),
        axis.text.x=element_text(size = 34),
        legend.position = "none",
        plot.title.position = "plot",plot.title = element_text(hjust = 0.6,vjust = -0.6))

TimetoMaxInf_24+theme_bw()

ggsave("TimetoMaxInfected_24_newdata.png", width = 24, height = 28)
ggsave("TimetoMaxInfected_24_newdata.pdf", width = 24, height = 38)

OutDuratn_24=m1_24%>%gather(key, value = val,OutbreakDuration)%>%
  ggplot(aes(x = key, y = val, group_by=GraphName, color = GraphName)) +
  geom_boxplot(show.legend = T) +ylab("Timesteps (days)")+labs(title = "Average-degree-24")+
  facet_grid(B~G~size,labeller = label_both)+
  scale_y_continuous(limits = c(0,50))+
  theme(text = element_text(size = 34),
        axis.title = element_text(size = 34),
        axis.text.y = element_text(size = 34),
        axis.text.x=element_text(size = 34),
        legend.position = "none",
        plot.title.position = "plot",plot.title = element_text(hjust = 0.6,vjust = -0.6))

OutDuratn_24+theme_bw()
ggsave("OutbreakDuration_24_newdata.png", width = 24, height = 38)
ggsave("OutbreakDuration_24_newdata.pdf", width = 24, height = 38)


MaxInf_24=m1_24%>%gather(key, value = val,MaxInfected)%>%
  ggplot(aes(x = key, y = val, group_by=GraphName, color = GraphName)) +
  geom_boxplot(show.legend = T) +ylab("Number-of-Infected-Individuals)")+labs(title = "Average-degree-24")+
  facet_grid(B~G~size,labeller = label_both)+
  scale_y_continuous(limits=c(0,100))+
#  scale_y_continuous(labels=percent)+
  theme(text = element_text(size = 34),
        axis.title = element_text(size =34),
        axis.text.y = element_text(size = 34),
        axis.text.x=element_text(size = 34),
        legend.position = "none",
        plot.title.position = "plot",plot.title = element_text(hjust = 0.6,vjust = -0.6))
MaxInf_24+theme_bw()

ggsave("MaxInfected_24_newdata.png", width = 24, height = 38)
ggsave("MaxInfected_24_newdata.pdf", width = 24, height = 38)

###---Analysis--of--real-networks-and--their--equivalent--###
############------Dolphin----------######
kg <- read.table("mammalia-dolphin-social.edges") ## the dolphin network
kg1=graph_from_data_frame(as.matrix(kg),directed=FALSE)
ggk1=simplify(kg1,remove.multiple = T,remove.loops = T)
ggk1$type="Dolphin"
ggk1$id="1"
G=list(ggk1)

#ecount(ggk1)
#vcount(ggk1)

#######-----------------------Hyena----------------##############
dd <- read.table("mammalia-hyena-networkc.edges") ## the hyena network
gg=graph_from_data_frame(as.matrix(dd),directed=FALSE)
gg1=simplify(gg,remove.multiple = T,remove.loops = T)
gg1$type="Hyena"
gg1$id="1"
png("Hyena-Network.png", 600, 600)
plot(gg1,vertex.label.font=0.1,edge.color="black",
     vertex.label = NA,vertex.color="gray")
title("Hyena-Network", sub = "nodes:35,edges:509",
      cex.main = 2,   font.main= 4, col.main= "black",
      cex.sub = 2, font.sub = 4, col.sub = "black")

#title("Hyena-Network",cex.main=1,col.main="black")
dev.off()

#ggsave("Hyena-Network.png", width = 40, height = 40)

#ecount(gg1)
#vcount(gg1)
#plot(gg1)



######################------------GRAPH FEATURES OF DOLPHIN-----------#######################
dolph1=makeSimGraphs_dolph(nSamples = 100,order = vcount(ggk1))
dolph1$"real"=ggk1
d1=RunSimOnGraphFeatures(dolph1,nreps = 1)

d11=d1%>%select(1,2,4,8,9,10,12,14,16,15,17)
#save(t11, file = "tableofgraphfeatures.Rda")  
d2=write.csv(d11,'C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi/Dolphin_GraphFeat.csv')


######################------------GRAPH FEATURES OF HYENA-----------#######################
h=makeSimGraphs_hyena(nSamples = 100,order = vcount(gg1))
h$"real"=gg1
t1=RunSimOnGraphFeatures(h,nreps = 1)

t11=t1%>%select(1,2,4,8,9,10,12,14,16,15,17)
#save(t11, file = "tableofgraphfeatures.Rda")  
t2=write.csv(t11,'C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi/Hyena_GraphFeat.csv')


##########--------------PLOTS OF REAL AND SYNTHETIC EQUIVALENTS---------------------###################

### DOLPHIN NETWORK AND SYNTHETIC
png("Dolph.png", 400, 400)
#par(mfrow=c(1,2))
DOLPH=plot(dolph1[[402]],vertex.label = NA,vertex.size=5,vertex.color="red")

title(main = list("Dolphin:n=62,e=159", cex = 2, font= 2,
                  col = "black"))
dev.off()

### ERDOS NETWORK
png("ER-dolph.png", 400, 400)
#par(mfrow=c(1,2))
ERGRAPH_dolph=plot(dolph1[[1]],vertex.label = NA,vertex.size=5,vertex.color="red")

title(main = list("Erdos Renyi:n=62,e=159", cex = 2, font= 2,
                  col = "black"))
dev.off()

### SMALL WORLD NETWORK
png("Small-World-dolph.png", 400, 400)
#par(mfrow=c(1,2))
SWGRAPH_dolph=plot(dolph1[[101]],vertex.label = NA,vertex.size=5,vertex.color="red")


title(main = list("Small world:n=64,e=128", cex = 2, font= 2,
                  col = "black"))
dev.off()

### SCALE FREE NETWORK
png("Scale-Free-dolph.png", 400, 400)
#par(mfrow=c(1,2))
SFGRAPH_dolph=plot(dolph1[[205]],vertex.label = NA,vertex.size=5,vertex.color="red")


title(main = list("Scale-Free:n=62,e=121", cex = 2, font= 2,
                  col = "black"))
dev.off()

### SPATIAL NETWORK
png("Spatial-dolph.png", 400, 400)
#par(mfrow=c(1,2))
SPGRAPH_dolph=plot(dolph1[[388]],vertex.label = NA,vertex.size=5,vertex.color="red")


title(main = list("Spatial:n=62,e=159", cex = 2, font= 2,
                  col = "black"))
dev.off()

### LATTICE NETWORK
png("Lattice-dolph.png", 400, 400)
#par(mfrow=c(1,2))
LATGRAPH_dolph=plot(dolph1[[401]],vertex.label = NA,vertex.size=5,vertex.color="red")


title(main = list("Lattice:n=64,e=112", cex = 2, font= 2,
                  col = "black"))
dev.off()



####--PLOTS
plot1 <- readPNG('Dolph.png')
plot2 <- readPNG('Spatial-dolph.png')
plot3 <- readPNG('Lattice-dolph.png')
plot4 <- readPNG('Small-World-dolph.png')
plot5 <- readPNG('ER-dolph.png')
plot6 <- readPNG('Scale-Free-dolph.png')

plotsalldolph1<- arrangeGrob(rasterGrob(plot1),rasterGrob(plot2),
                             rasterGrob(plot3),rasterGrob(plot4),
                             rasterGrob(plot5),rasterGrob(plot6),nrow=2)
ggsave('all_dolph_synthetic1.png',plotsalldolph1,width=15,height=10)


####-----------------HYENA NETWORK AND SYNTHETIC----######
png("Hyena.png", 400, 400)

HYENA=plot(h[[402]],vertex.label = NA,vertex.size=5,vertex.color="red")
title(main = list("Hyena:n=35,e=509", cex = 2, font= 2,
                  col = "black"))
dev.off()


#######-----------------ERDOS RENYI------############
png("ER-Network.png", 400, 400)

ERGRAPH=plot(h[[1]],vertex.label = NA,vertex.size=5,vertex.color="red")

title(main = list("Erdos Renyi:n=35,e=509",  cex = 2, font= 2,
                  col = "black"))
dev.off()

####-----SMALL WORLD NETWORK----######
png("Small-World.png", 400, 400)

SWGRAPH=plot(h[[101]],vertex.label = NA,vertex.size=5,vertex.color="red")

title(main = list("Small world:n=36,e=540",  cex = 2, font= 2,
                  col = "black"))
dev.off()

### SCALE FREE NETWORK
png("Scale-Free.png", 400, 400)

SFGRAPH=plot(h[[205]],vertex.label = NA,vertex.size=5,vertex.color="red")
title(main = list("Scale-Free:n=35,e=517",  cex = 2, font= 2,
                  col = "black"))
dev.off()

### SPATIAL NETWORK
png("Spatial.png", 400, 400)

SPGRAPH=plot(h[[355]],vertex.label = NA,vertex.size=5,vertex.color="red")
title(main = list("Spatial:n=35,e=508",  cex = 2, font= 2,
                  col = "black"))
dev.off()

### LATTICE NETWORK
png("Lattice.png", 400, 400)

LATGRAPH=plot(h[[401]],vertex.label = NA,vertex.size=5,vertex.color="red")
title(main = list("Lattice:n=36,e=560", cex = 2, font= 2,
                  col = "black"))
dev.off()



plot1 <- readPNG('Hyena.png')
plot2 <- readPNG('Spatial.png')
plot3 <- readPNG('Lattice.png')
plot4 <- readPNG('Small-World.png')
plot5 <- readPNG('ER-Network.png')
plot6 <- readPNG('Scale-Free.png')

plotsall1<- arrangeGrob(rasterGrob(plot1),rasterGrob(plot2),
                        rasterGrob(plot3),rasterGrob(plot4),
                        rasterGrob(plot5),rasterGrob(plot6),nrow=2)

ggsave('all_hyena_synthetic1.png',plotsall1,width=15,height=10)


grid.arrange(rasterGrob(plot1),rasterGrob(plot2),
             rasterGrob(plot3),rasterGrob(plot4),
             rasterGrob(plot5),rasterGrob(plot6),nrow=2)

ggsave('all_hyena_synthetic3.png',plotsall1,width=25,height=30)

######---------------DATA SIMULATION of Dolphin and the synthetic networks------------#####
dolph1=makeSimGraphs_dolph(nSamples = 10,order = vcount(ggk1))
dolph1$"real"=ggk1
set.seed(3456)
#nsim=5 
nreps=10 # number of simulations on with each combination of beta and gamma
nticks=100 # number of timesteps

betaVals=c(0.1,0.28,0.33,0.5,0.76) 
gammaVals=c(0.033,0.143,0.25,0.43,0.8) 

Dolphin_vs_synthetic_networks=ParalleEpicSimOnGraphs(dolph1, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="Dolphine_vs_synthetic_network.csv",report="i")


#####----Some graph features on the real networks and their equivalents----#####

##Dolphin
x=Graphfeatures(Name="Dolphin",data="Dolphin_GraphFeat.csv")
min(x$minDegree)
max(x$maxDegree)
quantile(x$FiedlerValue,probs = c(0.25, 0.75))
quantile(x$modularity,probs = c(0.25, 0.75))
quantile(x$betweenness,probs = c(0.25, 0.75))
quantile(x$transitivity,probs = c(0.25, 0.75))


##ER
x=Graphfeatures(Name="ER",data="Dolphin_GraphFeat.csv")
min(x$minDegree)
max(x$maxDegree)
quantile(x$FiedlerValue,probs = c(0.25, 0.75))
quantile(x$modularity,probs = c(0.25, 0.75))
quantile(x$betweenness,probs = c(0.25, 0.75))
quantile(x$transitivity,probs = c(0.25, 0.75))

##SW
x=Graphfeatures(Name="SW",data="Dolphin_GraphFeat.csv")
min(x$minDegree)
max(x$maxDegree)
quantile(x$FiedlerValue,probs = c(0.25, 0.75))
quantile(x$modularity,probs = c(0.25, 0.75))
quantile(x$betweenness,probs = c(0.25, 0.75))
quantile(x$transitivity,probs = c(0.25, 0.75))

##SF
x=Graphfeatures(Name="SF",data="Dolphin_GraphFeat.csv")
min(x$minDegree)
max(x$maxDegree)
quantile(x$FiedlerValue,probs = c(0.25, 0.75))
quantile(x$modularity,probs = c(0.25, 0.75))
quantile(x$betweenness,probs = c(0.25, 0.75))
quantile(x$transitivity,probs = c(0.25, 0.75))


##SP
x=Graphfeatures(Name="SP",data="Dolphin_GraphFeat.csv")
min(x$minDegree)
max(x$maxDegree)
quantile(x$FiedlerValue,probs = c(0.25, 0.75))
quantile(x$modularity,probs = c(0.25, 0.75))
quantile(x$betweenness,probs = c(0.25, 0.75))
quantile(x$transitivity,probs = c(0.25, 0.75))


##Lat
x=Graphfeatures(Name="Lat",data="Dolphin_GraphFeat.csv")
min(x$minDegree)
max(x$maxDegree)
quantile(x$FiedlerValue,probs = c(0.25, 0.75))
quantile(x$modularity,probs = c(0.25, 0.75))
quantile(x$betweenness,probs = c(0.25, 0.75))
quantile(x$transitivity,probs = c(0.25, 0.75))



#########--------Epidemic dynamics on dolphin and its equivalent------#######

#####------Infected dynamics for dolphin----##########
set.seed((234567))
nsim=10
n_dolph=62 # number of nodes

#Beta1
ER_beta1=plotfunc_obs(Name = "ER", ID=1,betaVal = .5,gammaVal=.143,Data = "Dolphine_vs_synthetic_network.csv", nreps = nsim,plot_title="ER",net_size = 35)
SW_beta1=plotfunc_obs(Name = "SW", ID=1,betaVal = .5,gammaVal=.143,Data = "Dolphine_vs_synthetic_network.csv", nreps = nsim,plot_title="SW",net_size = 36)
SF_beta1=plotfunc_obs(Name = "SF", ID=1,betaVal = .5,gammaVal=.143,Data = "Dolphine_vs_synthetic_network.csv", nreps = nsim,plot_title="SF",net_size = 35)
SP_beta1=plotfunc_obs(Name = "SP", ID=1,betaVal = .5,gammaVal=.143,Data = "Dolphine_vs_synthetic_network.csv", nreps = nsim,plot_title="SP",net_size = 35)
Lat_beta1=plotfunc_obs(Name = "Lat", ID=1,betaVal = .5,gammaVal=.143,Data = "Dolphine_vs_synthetic_network.csv", nreps = nsim,plot_title="Lat",net_size = 36)
Dolph_beta1=plotfunc_obs(Name = "Dolphin", ID=1,betaVal = .5,gammaVal=.143,Data = "Dolphine_vs_synthetic_network.csv", nreps = nsim,plot_title="Hyena",net_size = 35)



plot1=ggarrange(ER_beta1,SW_beta1,SF_beta1,nrow = 1)
plot2=ggarrange(Lat_beta1,SP_beta1,Dolph_beta1,nrow=1)


AllPlotDolph=ggarrange(plot1,plot2,nrow = 2)

AllPlotDolph_annotate=annotate_figure(AllPlotDolph, 
                                      top = text_grob("Dolphin vs synthetic-networks-1", 
                                                      color = "black", face = "bold", size = 26),
                                      fig.lab = "a)",fig.lab.size = 40)

ggsave("allplot_Dolphin_synthetic.png", width = 28, height = 24)
ggsave("allplot_Dolphin_synthetic.pdf", width = 30, height = 38)


#####----AVerage proportion of infected for dolphin------------############
n_dolph=62
ER=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = .143,nreps = 10,Data = "Dolphine_vs_synthetic_network.csv")
SW=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = .143,nreps = 10,Data = "Dolphine_vs_synthetic_network.csv")
SF=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = .143,nreps = 10,Data = "Dolphine_vs_synthetic_network.csv")
SP=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = .143,nreps = 10,Data = "Dolphine_vs_synthetic_network.csv")
Lat=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = .143,nreps = 10,Data = "Dolphine_vs_synthetic_network.csv")
Dolphin=AvrgPropOf(Name="Dolphin",ID=1,betaVal = 0.5,gammaVal = .143,nreps = 10,Data = "Dolphine_vs_synthetic_network.csv")

df=rbind(ER,SW,SF,SP,Lat,Dolphin)
plotdolph<- ggplot(df,aes(x=Timesteps, y=value/n_dolph, group=GraphName))+
  geom_line(show.legend = T,aes(linetype=GraphName, color=GraphName),size = 1.5)+
  theme_classic()+labs(title="Dolphin vs synthetic-networks-1")+
  xlab("Timesteps (days)")+ylab("Average-proportion-of-infected")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 34),
        axis.title = element_text(size = 34),
        axis.text.y = element_text(size = 34),
        axis.text.x=element_text(size = 34),
        legend.position = "bottom",plot.title = element_text(size = 34,hjust = 0.3, face = "bold"))

plotdolph

plotdolph_annotate=annotate_figure(plotdolph,
                                   top = text_grob(" ",
                                                   color = "black", face = "bold",
                                                   size = 26),fig.lab = "a)",fig.lab.size = 40)



#####----AVERAGE PROPORTION FOR HYENA------------############
n_hyena=35
ER=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = .143,nreps = 10,Data = "hyena_vs_synthetic_netwok.csv")
SW=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = .143,nreps = 10,Data = "hyena_vs_synthetic_netwok.csv")
SF=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = .143,nreps = 10,Data = "hyena_vs_synthetic_netwok.csv")
SP=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = .143,nreps = 10,Data = "hyena_vs_synthetic_netwok.csv")
Lat=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = .143,nreps = 10,Data = "hyena_vs_synthetic_netwok.csv")
Hyena=AvrgPropOf(Name="Hyena",ID=1,betaVal = 0.5,gammaVal = .143,nreps = 10,Data = "hyena_vs_synthetic_netwok.csv")


df=rbind(ER,SW,SF,SP,Lat,Hyena)
plothyena<- ggplot(df,aes(x=Timesteps, y=value/n_hyena, group=GraphName))+
  geom_line(show.legend = T,aes(linetype=GraphName, color=GraphName),size = 1.5)+
  theme_classic()+labs(title="Hyena vs synthetic-network-2")+
  xlab("Timesteps (days)")+ylab("Average-proportion-of-infected")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 34),
        axis.title = element_text(size = 34),
        axis.text.y = element_text(size = 34),
        axis.text.x=element_text(size = 34),
        legend.position = "bottom",plot.title = element_text(size = 34,hjust = 0.3, face = "bold"))

plothyena

plothyena_annotate=annotate_figure(plothyena,
                                   top = text_grob(" ",
                                                   color = "black", face = "bold",
                                                   size = 26),fig.lab = "b)",fig.lab.size = 40)






All1=ggarrange(plotdolph_annotate,plothyena_annotate,nrow=1,labels="AUTO")
ggsave("allplot_hyena&dolph_synthetic.png", width = 18, height = 10)
ggsave("allplot_hyena&dolph_synthetic.pdf", width = 18, height = 10)

All2=ggarrange(plotdolph_annotate,plothyena_annotate,ncol=2)
ggsave("avd_allplot_hyena&dolph_synthetic.png", width = 15, height = 15)



############---------Epidemic dynamics With different beta for dolphin------------########################### 

#########Dolphin####
nticks=100
nreps=10

data_obs1=read.csv("Dolphine_vs_synthetic_network.csv",header = T, sep = ",")
df_obs1=data_obs1%>% filter(GraphName=="ER",GraphID==1)
df_obs1=df_obs1%>% select(c(GraphName,beta,gamma,t0:paste("t",nticks,sep = "")))
# View(df_obs1)


p_obs1=data.frame(0:(length(df_obs1[grep("t", names(df_obs1))])-1),t(df_obs1[grep("t", names(df_obs1))]))
colnames(p_obs1)= c("Timesteps (days)",paste("Infecteds",1:nreps,sep = "_")) #long format
df_obs_long_format1 <- melt(p_obs1, id="Timesteps (days)")  # convert to long format




###########----ERDOS-------##################
Dolph_ERbeta1=AvrgPropOf(Name="ER",ID=1,betaVal = 0.01,gammaVal = 0.143,nreps = 10,
                         Data = "Dolphine_vs_synthetic_network.csv")
Dolph_ERbeta1$Betavalues="beta:0.01"

Dolph_ERbeta2=AvrgPropOf(Name="ER",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,
                         Data = "Dolphine_vs_synthetic_network.csv")
Dolph_ERbeta2$Betavalues="beta:0.1"

Dolph_ERbeta3=AvrgPropOf(Name="ER",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,
                         Data = "Dolphine_vs_synthetic_network.csv")
Dolph_ERbeta3$Betavalues="beta:0.33"

Dolph_ERbeta4=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,
                         Data = "Dolphine_vs_synthetic_network.csv")
Dolph_ERbeta4$Betavalues="beta:0.5"


df=rbind(Dolph_ERbeta1,Dolph_ERbeta2,Dolph_ERbeta3,Dolph_ERbeta4)
n_dolph=62
plotER<- ggplot(df,aes(x=Timesteps, y=value/n_dolph, group=Betavalues))+
  geom_line(show.legend = T,
            aes(linetype=Betavalues, color=Betavalues),size=2)+
  #geom_point(colour = alpha("blue", 0.5))
  geom_point(aes(color=Betavalues))+
  theme_classic()+labs(title="ER")+
  xlab("Timesteps (days)")+ylab("Prop-infected")+
  scale_y_continuous(limits = c(0,1))+
  scale_linetype_manual(values=c("longdash", "solid", "dotted", "solid"))+
  scale_color_manual(values=c('#999999','#E69F00','red','black'))+
  scale_size_manual(values=c(2,1,2,1.5))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title = element_text(size = 34,hjust = 0.2,vjust = -0.6, face = "bold"))



plotER


###########----SMALL WORLD-------##################

Dolph_SWbeta1=AvrgPropOf(Name="SW",ID=1,betaVal = 0.01,gammaVal = 0.143,nreps = 10,
                         Data = "Dolphine_vs_synthetic_network.csv")
Dolph_SWbeta1$Betavalues="beta:0.01"

Dolph_SWbeta2=AvrgPropOf(Name="SW",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,
                         Data = "Dolphine_vs_synthetic_network.csv")
Dolph_SWbeta2$Betavalues="beta:0.1"

Dolph_SWbeta3=AvrgPropOf(Name="SW",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,
                         Data = "Dolphine_vs_synthetic_network.csv")
Dolph_SWbeta3$Betavalues="beta:0.33"

Dolph_SWbeta4=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,
                         Data = "Dolphine_vs_synthetic_network.csv")
Dolph_SWbeta4$Betavalues="beta:0.5"


df=rbind(Dolph_SWbeta1,Dolph_SWbeta2,Dolph_SWbeta3,Dolph_SWbeta4)
n_dolph=62
plotSW<- ggplot(df,aes(x=Timesteps, y=value/n_dolph, group=Betavalues))+
  geom_line(show.legend = T,
            aes(linetype=Betavalues, color=Betavalues),size=2)+
  #geom_point(colour = alpha("blue", 0.5))
  geom_point(aes(color=Betavalues))+
  theme_classic()+labs(title="SW")+
  xlab("Timesteps (days)")+ylab("Prop-infected")+
  scale_y_continuous(limits = c(0,1))+
  scale_linetype_manual(values=c("longdash", "solid", "dotted", "solid"))+
  scale_color_manual(values=c('#999999','#E69F00','red','black'))+
  scale_size_manual(values=c(2,1,2,1.5))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title = element_text(size = 34,hjust = 0.2,vjust = -0.6, face = "bold"))



plotSW

###########----SCALE FREE-------##################

Dolph_SFbeta1=AvrgPropOf(Name="SF",ID=1,betaVal = 0.01,gammaVal = 0.143,nreps = 10,
                         Data = "Dolphine_vs_synthetic_network.csv")
Dolph_SFbeta1$Betavalues="beta:0.01"

Dolph_SFbeta2=AvrgPropOf(Name="SF",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,
                         Data = "Dolphine_vs_synthetic_network.csv")
Dolph_SFbeta2$Betavalues="beta:0.1"

Dolph_SFbeta3=AvrgPropOf(Name="SF",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,
                         Data = "Dolphine_vs_synthetic_network.csv")
Dolph_SFbeta3$Betavalues="beta:0.33"

Dolph_SFbeta4=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,
                         Data = "Dolphine_vs_synthetic_network.csv")
Dolph_SFbeta4$Betavalues="beta:0.5"


df=rbind(Dolph_SFbeta1,Dolph_SFbeta2,Dolph_SFbeta3,Dolph_SFbeta4)
n_dolph=62
plotSF<- ggplot(df,aes(x=Timesteps, y=value/n_dolph, group=Betavalues))+
  geom_line(show.legend = T,
            aes(linetype=Betavalues, color=Betavalues),size=2)+
  #geom_point(colour = alpha("blue", 0.5))
  geom_point(aes(color=Betavalues))+
  theme_classic()+labs(title="SF")+
  xlab("Timesteps (days)")+ylab("Prop-infected")+
  scale_y_continuous(limits = c(0,1))+
  scale_linetype_manual(values=c("longdash", "solid", "dotted", "solid"))+
  scale_color_manual(values=c('#999999','#E69F00','red','black'))+
  scale_size_manual(values=c(2,1,2,1.5))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title = element_text(size = 34,hjust = 0.2,vjust = -0.6, face = "bold"))

plotSF

###########----SPATIAL-------##################

Dolph_SPbeta1=AvrgPropOf(Name="SP",ID=1,betaVal = 0.01,gammaVal = 0.143,nreps = 10,
                         Data = "Dolphine_vs_synthetic_network.csv")
Dolph_SPbeta1$Betavalues="beta:0.01"

Dolph_SPbeta2=AvrgPropOf(Name="SP",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,
                         Data = "Dolphine_vs_synthetic_network.csv")
Dolph_SPbeta2$Betavalues="beta:0.1"

Dolph_SPbeta3=AvrgPropOf(Name="SP",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,
                         Data = "Dolphine_vs_synthetic_network.csv")
Dolph_SPbeta3$Betavalues="beta:0.33"

Dolph_SPbeta4=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,
                         Data = "Dolphine_vs_synthetic_network.csv")
Dolph_SPbeta4$Betavalues="beta:0.5"


df=rbind(Dolph_SPbeta1,Dolph_SPbeta2,Dolph_SPbeta3,Dolph_SPbeta4)
n_dolph=62
plotSP<- ggplot(df,aes(x=Timesteps, y=value/n_dolph, group=Betavalues))+
  geom_line(show.legend = T,
            aes(linetype=Betavalues, color=Betavalues),size=2)+
  #geom_point(colour = alpha("blue", 0.5))
  geom_point(aes(color=Betavalues))+
  theme_classic()+labs(title="SP")+
  xlab("Timesteps (days)")+ylab("Prop-infected")+
  scale_y_continuous(limits = c(0,1))+
  scale_linetype_manual(values=c("longdash", "solid", "dotted", "solid"))+
  scale_color_manual(values=c('#999999','#E69F00','red','black'))+
  scale_size_manual(values=c(2,1,2,1.5))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title = element_text(size = 34,hjust = 0.2,vjust = -0.6, face = "bold"))

plotSP

###########----LATTICE-------##################

Dolph_Latbeta1=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.01,gammaVal = 0.143,nreps = 10,
                          Data = "Dolphine_vs_synthetic_network.csv")
Dolph_Latbeta1$Betavalues="beta:0.01"

Dolph_Latbeta2=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,
                          Data = "Dolphine_vs_synthetic_network.csv")
Dolph_Latbeta2$Betavalues="beta:0.1"

Dolph_Latbeta3=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,
                          Data = "Dolphine_vs_synthetic_network.csv")
Dolph_Latbeta3$Betavalues="beta:0.33"

Dolph_Latbeta4=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,
                          Data = "Dolphine_vs_synthetic_network.csv")
Dolph_Latbeta4$Betavalues="beta:0.5"


df=rbind(Dolph_Latbeta1,Dolph_Latbeta2,Dolph_Latbeta3,Dolph_Latbeta4)
n_dolph=62
plotLat<- ggplot(df,aes(x=Timesteps, y=value/n_dolph, group=Betavalues))+
  geom_line(show.legend = T,
            aes(linetype=Betavalues, color=Betavalues),size=2)+
  #geom_point(colour = alpha("blue", 0.5))
  geom_point(aes(color=Betavalues))+
  theme_classic()+labs(title="Lat")+
  xlab("Timesteps (days)")+ylab("Prop-infected")+
  scale_y_continuous(limits = c(0,1))+
  scale_linetype_manual(values=c("longdash", "solid", "dotted", "solid"))+
  scale_color_manual(values=c('#999999','#E69F00','red','black'))+
  scale_size_manual(values=c(2,1,2,1.5))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title = element_text(size = 34,hjust = 0.2,vjust = -0.6, face = "bold"))

plotLat

###########----Dolphin-------##################

Dolph_graphbeta1=AvrgPropOf(Name="Dolphin",ID=1,betaVal = 0.01,gammaVal = 0.143,nreps = 10,
                            Data = "Dolphine_vs_synthetic_network.csv")
Dolph_graphbeta1$Betavalues="beta:0.01"

Dolph_graphbeta2=AvrgPropOf(Name="Dolphin",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,
                            Data = "Dolphine_vs_synthetic_network.csv")
Dolph_graphbeta2$Betavalues="beta:0.1"

Dolph_graphbeta3=AvrgPropOf(Name="Dolphin",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,
                            Data = "Dolphine_vs_synthetic_network.csv")
Dolph_graphbeta3$Betavalues="beta:0.33"

Dolph_graphbeta4=AvrgPropOf(Name="Dolphin",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,
                            Data = "Dolphine_vs_synthetic_network.csv")
Dolph_graphbeta4$Betavalues="beta:0.5"


df=rbind(Dolph_graphbeta1,Dolph_graphbeta2,Dolph_graphbeta3,Dolph_graphbeta4)
n_dolph=62
plotDolphgraph<- ggplot(df,aes(x=Timesteps, y=value/n_dolph, group=Betavalues))+
  geom_line(show.legend = T,
            aes(linetype=Betavalues, color=Betavalues),size=2)+
  #geom_point(colour = alpha("blue", 0.5))
  geom_point(aes(color=Betavalues))+
  theme_classic()+labs(title="Dolph graph")+
  xlab("Timesteps (days)")+ylab("Prop-infected")+
  scale_y_continuous(limits = c(0,1))+
  scale_linetype_manual(values=c("longdash", "solid", "dotted", "solid"))+
  scale_color_manual(values=c('#999999','#E69F00','red','black'))+
  scale_size_manual(values=c(2,1,2,1.5))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title = element_text(size = 34,hjust = 0.2,vjust = -0.6, face = "bold"))

plotDolphgraph

plot1=ggarrange(plotER,plotSW,nrow=1)
plot2=ggarrange(plotLat,plotSF,nrow=1)
plot3=ggarrange(plotSP,plotDolphgraph,common.legend = TRUE, legend="bottom",nrow=1)

plotsalldolph=ggarrange(plot1,plot2,plot3,nrow=3,
                        common.legend = TRUE, legend="bottom")

Allannotatedolph=annotate_figure(
  plotsalldolph, top = text_grob("Dolphin vs synthetic-networks-1 ",
                                 color = "black", face = "bold", size = 30),
  fig.lab = "a)",fig.lab.size = 40)

#########HYENA####
nticks=100
nreps=10

data_obs=read.csv("hyena_vs_synthetic_netwok.csv",header = T, sep = ",")
df_obs=data_obs%>% filter(GraphName=="ER",GraphID==1)
df_obs=df_obs%>% select(c(GraphName,beta,gamma,t0:paste("t",nticks,sep = "")))
View(df_obs)


p_obs=data.frame(0:(length(df_obs[grep("t", names(df_obs))])-1),t(df_obs[grep("t", names(df_obs))]))
colnames(p_obs)= c("Timesteps (days)",paste("Infecteds",1:nreps,sep = "_")) #long format
df_obs_long_format <- melt(p_obs, id="Timesteps (days)")  # convert to long format


###########----ERDOS-------##################
Hyena_ERbeta1=AvrgPropOf(Name="ER",ID=1,betaVal = 0.01,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_ERbeta1$Betavalues="beta:0.01"

Hyena_ERbeta2=AvrgPropOf(Name="ER",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_ERbeta2$Betavalues="beta:0.1"

Hyena_ERbeta3=AvrgPropOf(Name="ER",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_ERbeta3$Betavalues="beta:0.33"

Hyena_ERbeta4=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_ERbeta4$Betavalues="beta:0.5"


df=rbind(Hyena_ERbeta1,Hyena_ERbeta2,Hyena_ERbeta3,Hyena_ERbeta4)
n_hyena=35
plotER_hyena<- ggplot(df,aes(x=Timesteps, y=value/n_hyena, group=Betavalues))+
  geom_line(show.legend = T,
            aes(linetype=Betavalues, color=Betavalues),size=2)+
  #geom_point(colour = alpha("blue", 0.5))
  geom_point(aes(color=Betavalues))+
  theme_classic()+labs(title="ER")+
  xlab("Timesteps (days)")+ylab("Prop-infected")+
  scale_y_continuous(limits = c(0,1))+
  scale_linetype_manual(values=c("longdash", "solid", "dotted", "solid"))+
  scale_color_manual(values=c('#999999','#E69F00','red','black'))+
  scale_size_manual(values=c(2,1,2,1.5))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title = element_text(size = 34,hjust = 0.2,vjust = -0.6, face = "bold"))



plotER_hyena

###########----SMALL WORLD-------##################
Hyena_SWbeta1=AvrgPropOf(Name="SW",ID=1,betaVal = 0.01,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_SWbeta1$Betavalues="beta:0.01"

Hyena_SWbeta2=AvrgPropOf(Name="SW",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_SWbeta2$Betavalues="beta:0.1"

Hyena_SWbeta3=AvrgPropOf(Name="SW",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_SWbeta3$Betavalues="beta:0.33"

Hyena_SWbeta4=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_SWbeta4$Betavalues="beta:0.5"


df=rbind(Hyena_SWbeta1,Hyena_SWbeta2,Hyena_SWbeta3,Hyena_SWbeta4)
n_hyena=35
plotSW_hyena<- ggplot(df,aes(x=Timesteps, y=value/n_hyena, group=Betavalues))+
  geom_line(show.legend = T,
            aes(linetype=Betavalues, color=Betavalues),size=2)+
  #geom_point(colour = alpha("blue", 0.5))
  geom_point(aes(color=Betavalues))+
  theme_classic()+labs(title="SW")+
  xlab("Timesteps (days)")+ylab("Prop-infected")+
  scale_y_continuous(limits = c(0,1))+
  scale_linetype_manual(values=c("longdash", "solid", "dotted", "solid"))+
  scale_color_manual(values=c('#999999','#E69F00','red','black'))+
  scale_size_manual(values=c(2,1,2,1.5))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title = element_text(size = 34,hjust = 0.2,vjust = -0.6, face = "bold"))



plotSW_hyena


###########----SCALE FREE-------##################
Hyena_SFbeta1=AvrgPropOf(Name="SF",ID=1,betaVal = 0.01,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_SFbeta1$Betavalues="beta:0.01"

Hyena_SFbeta2=AvrgPropOf(Name="SF",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_SFbeta2$Betavalues="beta:0.1"

Hyena_SFbeta3=AvrgPropOf(Name="SF",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_SFbeta3$Betavalues="beta:0.33"

Hyena_SFbeta4=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_SFbeta4$Betavalues="beta:0.5"


df=rbind(Hyena_SFbeta1,Hyena_SFbeta2,Hyena_SFbeta3,Hyena_SFbeta4)
n_hyena=35
plotSF_hyena<- ggplot(df,aes(x=Timesteps, y=value/n_hyena, group=Betavalues))+
  geom_line(show.legend = T,
            aes(linetype=Betavalues, color=Betavalues),size=2)+
  #geom_point(colour = alpha("blue", 0.5))
  geom_point(aes(color=Betavalues))+
  theme_classic()+labs(title="SF")+
  xlab("Timesteps (days)")+ylab("Prop-infected")+
  scale_y_continuous(limits = c(0,1))+
  scale_linetype_manual(values=c("longdash", "solid", "dotted", "solid"))+
  scale_color_manual(values=c('#999999','#E69F00','red','black'))+
  scale_size_manual(values=c(2,1,2,1.5))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title = element_text(size = 34,hjust = 0.2,vjust = -0.6, face = "bold"))

plotSF_hyena

###########----SPATIAL-------##################
Hyena_SPbeta1=AvrgPropOf(Name="SP",ID=1,betaVal = 0.01,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_SPbeta1$Betavalues="beta:0.01"

Hyena_SPbeta2=AvrgPropOf(Name="SP",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_SPbeta2$Betavalues="beta:0.1"

Hyena_SPbeta3=AvrgPropOf(Name="SP",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_SPbeta3$Betavalues="beta:0.33"

Hyena_SPbeta4=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_SPbeta4$Betavalues="beta:0.5"


df=rbind(Hyena_SPbeta1,Hyena_SPbeta2,Hyena_SPbeta3,Hyena_SPbeta4)
n_hyena=35
plotSP_hyena<- ggplot(df,aes(x=Timesteps, y=value/n_hyena, group=Betavalues))+
  geom_line(show.legend = T,
            aes(linetype=Betavalues, color=Betavalues),size=2)+
  #geom_point(colour = alpha("blue", 0.5))
  geom_point(aes(color=Betavalues))+
  theme_classic()+labs(title="SP")+
  xlab("Timesteps (days)")+ylab("Prop-infected")+
  scale_y_continuous(limits = c(0,1))+
  scale_linetype_manual(values=c("longdash", "solid", "dotted", "solid"))+
  scale_color_manual(values=c('#999999','#E69F00','red','black'))+
  scale_size_manual(values=c(2,1,2,1.5))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title = element_text(size = 34,hjust = 0.2,vjust = -0.6, face = "bold"))

plotSP_hyena

###########----LATTICE-------##################
Hyena_Latbeta1=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.01,gammaVal = 0.143,nreps = 10,
                          Data = "hyena_vs_synthetic_netwok.csv")
Hyena_Latbeta1$Betavalues="beta:0.01"

Hyena_Latbeta2=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,
                          Data = "hyena_vs_synthetic_netwok.csv")
Hyena_Latbeta2$Betavalues="beta:0.1"

Hyena_Latbeta3=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,
                          Data = "hyena_vs_synthetic_netwok.csv")


Hyena_Latbeta3$Betavalues="beta:0.33"

Hyena_Latbeta4=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,
                          Data = "hyena_vs_synthetic_netwok.csv")
Hyena_Latbeta4$Betavalues="beta:0.5"


df=rbind(Hyena_Latbeta1,Hyena_Latbeta2,Hyena_Latbeta3,Hyena_Latbeta4)
n_hyena=35
plotLat_hyena<- ggplot(df,aes(x=Timesteps, y=value/n_hyena, group=Betavalues))+
  geom_line(show.legend = T,
            aes(linetype=Betavalues, color=Betavalues),size=2)+
  #geom_point(colour = alpha("blue", 0.5))
  geom_point(aes(color=Betavalues))+
  theme_classic()+labs(title="Lat")+
  xlab("Timesteps (days)")+ylab("Prop-infected")+
  scale_y_continuous(limits = c(0,1))+
  scale_linetype_manual(values=c("longdash", "solid", "dotted", "solid"))+
  scale_color_manual(values=c('#999999','#E69F00','red','black'))+
  scale_size_manual(values=c(2,1,2,1.5))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title = element_text(size = 34,hjust = 0.2,vjust = -0.6, face = "bold"))

plotLat_hyena

###########----HYENA-------##################
Hyenabeta1=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.01,gammaVal = 0.143,nreps = 10,
                      Data = "hyena_vs_synthetic_netwok.csv")
Hyenabeta1$Betavalues="beta:0.01"

Hyenabeta2=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,
                      Data = "hyena_vs_synthetic_netwok.csv")
Hyenabeta2$Betavalues="beta:0.1"

Hyenabeta3=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,
                      Data = "hyena_vs_synthetic_netwok.csv")
Hyenabeta3$Betavalues="beta:0.33"

Hyenabeta4=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,
                      Data = "hyena_vs_synthetic_netwok.csv")
Hyenabeta4$Betavalues="beta:0.5"


df=rbind(Hyenabeta1,Hyenabeta2,Hyenabeta3,Hyenabeta4)
n_hyena=35

plotHyena<- ggplot(df,aes(x=Timesteps, y=value/n_hyena, group=Betavalues))+
  geom_line(show.legend = T,
            aes(linetype=Betavalues, color=Betavalues),size=2)+
  #geom_point(colour = alpha("blue", 0.5))
  geom_point(aes(color=Betavalues))+
  theme_classic()+labs(title="Hyena")+
  xlab("Timesteps (days)")+ylab("Prop-infected")+
  scale_y_continuous(limits = c(0,1))+
  scale_linetype_manual(values=c("longdash", "solid", "dotted", "solid"))+
  scale_color_manual(values=c('#999999','#E69F00','red','black'))+
  scale_size_manual(values=c(2,1,2,1.5))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title = element_text(size = 34,hjust = 0.2,vjust = -0.6, face = "bold"))

plotHyena

plot1_hyena=ggarrange(plotER_hyena,plotSW_hyena,nrow=1)
plot2_hyena=ggarrange(plotLat_hyena,plotSF_hyena,nrow=1)
plot3_hyena=ggarrange(plotSP_hyena,plotHyena,common.legend = TRUE, legend="bottom",nrow=1)

plotsallhyena=ggarrange(plot1_hyena,plot2_hyena,plot3_hyena,nrow=3,
                        common.legend = TRUE, legend="bottom")

Allannotatehyena=annotate_figure(
  plotsallhyena,
  top = text_grob("Hyena vs synthetic network-2 ",
                  color = "black", face = "bold", size = 30),
  fig.lab = "b)",fig.lab.size = 40)

##############----All Plots for VARYING BETA FOR HYENA AND DOLPHIN----#############

Allplot=ggarrange(Allannotatedolph,Allannotatehyena,nrow=2)
ggsave("allbeta_hyena&dolph_network.png", width = 30, height = 38)

ggsave("allbeta_hyena&dolph_network.pdf", width = 30, height = 38)

###############------------HYENA NETWORK------------#############

###########----------- Some graph features on hyena and its equivalent--------------
##ER
data="Hyena_GraphFeat.csv"
x=Graphfeatures(Name="ER",data)
min(x$minDegree)
max(x$maxDegree)
quantile(x$FiedlerValue,probs = c(0.25, 0.75))
quantile(x$modularity,probs = c(0.25, 0.75))
quantile(x$betweenness,probs = c(0.25, 0.75))
quantile(x$transitivity,probs = c(0.25, 0.75))

##SW
x=Graphfeatures(Name="SW",data)
min(x$minDegree)
max(x$maxDegree)
quantile(x$FiedlerValue,probs = c(0.25, 0.75))
quantile(x$modularity,probs = c(0.25, 0.75))
quantile(x$betweenness,probs = c(0.25, 0.75))
quantile(x$transitivity,probs = c(0.25, 0.75))

##SF
x=Graphfeatures(Name="SF",data)
min(x$minDegree)
max(x$maxDegree)
quantile(x$FiedlerValue,probs = c(0.25, 0.75))
quantile(x$modularity,probs = c(0.25, 0.75))
quantile(x$betweenness,probs = c(0.25, 0.75))
quantile(x$transitivity,probs = c(0.25, 0.75))


##SP
x=Graphfeatures(Name="SP",data)
min(x$minDegree)
max(x$maxDegree)
quantile(x$FiedlerValue,probs = c(0.25, 0.75))
quantile(x$modularity,probs = c(0.25, 0.75))
quantile(x$betweenness,probs = c(0.25, 0.75))
quantile(x$transitivity,probs = c(0.25, 0.75))


##Lat
x=Graphfeatures(Name="Lat",data)
min(x$minDegree)
max(x$maxDegree)
quantile(x$FiedlerValue,probs = c(0.25, 0.75))
quantile(x$modularity,probs = c(0.25, 0.75))
quantile(x$betweenness,probs = c(0.25, 0.75))
quantile(x$transitivity,probs = c(0.25, 0.75))

##Hyena
x=Graphfeatures(Name="Hyena",data)
min(x$minDegree)
max(x$maxDegree)
quantile(x$FiedlerValue,probs = c(0.25, 0.75))
quantile(x$modularity,probs = c(0.25, 0.75))
quantile(x$betweenness,probs = c(0.25, 0.75))
quantile(x$transitivity,probs = c(0.25, 0.75))


######---------------DATA SIMULATION------------#####
set.seed(3456)
nsim=10 
nreps=10
nticks=100 

betaVals=c(0.1,0.28,0.33,0.5,0.76,1.2) 
gammaVals=c(0.033,0.143,0.25,0.43,0.8,1.2) 

Hyena_vs_synthetic_networks=ParalleEpicSimOnGraphs(h, nticks=nticks, beta=betaVals,gamma=gammaVals, nreps=nreps,output_file="hyena_vs_synthetic_netwok.csv",report="i")



############---------Different Mean degree (Degree heterogeneity) ------------########################### 
#########HYENA####
nticks=100
nreps=10

data_obs=read.csv("hyena_vs_synthetic_netwok.csv",header = T, sep = ",")
df_obs=data_obs%>% filter(GraphName=="ER",GraphID==1)
df_obs=df_obs%>% select(c(GraphName,beta,gamma,t0:paste("t",nticks,sep = "")))
View(df_obs)


p_obs=data.frame(0:(length(df_obs[grep("t", names(df_obs))])-1),t(df_obs[grep("t", names(df_obs))]))
colnames(p_obs)= c("Timesteps (days)",paste("Infecteds",1:nreps,sep = "_")) #long format
df_obs_long_format <- melt(p_obs, id="Timesteps (days)")  # convert to long format


###########----ERDOS-------##################
Hyena_ERbeta1=AvrgPropOf(Name="ER",ID=1,betaVal = 0.01,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_ERbeta1$Betavalues="beta:0.01"

Hyena_ERbeta2=AvrgPropOf(Name="ER",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_ERbeta2$Betavalues="beta:0.1"

Hyena_ERbeta3=AvrgPropOf(Name="ER",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_ERbeta3$Betavalues="beta:0.33"

Hyena_ERbeta4=AvrgPropOf(Name="ER",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_ERbeta4$Betavalues="beta:0.5"


df=rbind(Hyena_ERbeta1,Hyena_ERbeta2,Hyena_ERbeta3,Hyena_ERbeta4)
n_hyena=35
plotER_hyena<- ggplot(df,aes(x=Timesteps, y=value/n_hyena, group=Betavalues))+
  geom_line(show.legend = T,
            aes(linetype=Betavalues, color=Betavalues),size=2)+
  #geom_point(colour = alpha("blue", 0.5))
  geom_point(aes(color=Betavalues))+
  theme_classic()+labs(title="ER")+
  xlab("Timesteps (days)")+ylab("Prop-infected")+
  scale_y_continuous(limits = c(0,1))+
  scale_linetype_manual(values=c("longdash", "solid", "dotted", "solid"))+
  scale_color_manual(values=c('#999999','#E69F00','red','black'))+
  scale_size_manual(values=c(2,1,2,1.5))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title = element_text(size = 34,hjust = 0.2,vjust = -0.6, face = "bold"))

plotER_hyena

###########----SMALL WORLD-------##################
Hyena_SWbeta1=AvrgPropOf(Name="SW",ID=1,betaVal = 0.01,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_SWbeta1$Betavalues="beta:0.01"

Hyena_SWbeta2=AvrgPropOf(Name="SW",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_SWbeta2$Betavalues="beta:0.1"

Hyena_SWbeta3=AvrgPropOf(Name="SW",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_SWbeta3$Betavalues="beta:0.33"

Hyena_SWbeta4=AvrgPropOf(Name="SW",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_SWbeta4$Betavalues="beta:0.5"


df=rbind(Hyena_SWbeta1,Hyena_SWbeta2,Hyena_SWbeta3,Hyena_SWbeta4)
n_hyena=35
plotSW_hyena<- ggplot(df,aes(x=Timesteps, y=value/n_hyena, group=Betavalues))+
  geom_line(show.legend = T,
            aes(linetype=Betavalues, color=Betavalues),size=2)+
  #geom_point(colour = alpha("blue", 0.5))
  geom_point(aes(color=Betavalues))+
  theme_classic()+labs(title="SW")+
  xlab("Timesteps (days)")+ylab("Prop-infected")+
  scale_y_continuous(limits = c(0,1))+
  scale_linetype_manual(values=c("longdash", "solid", "dotted", "solid"))+
  scale_color_manual(values=c('#999999','#E69F00','red','black'))+
  scale_size_manual(values=c(2,1,2,1.5))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title = element_text(size = 34,hjust = 0.2,vjust = -0.6, face = "bold"))



plotSW_hyena


###########----SCALE FREE-------##################
Hyena_SFbeta1=AvrgPropOf(Name="SF",ID=1,betaVal = 0.01,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_SFbeta1$Betavalues="beta:0.01"

Hyena_SFbeta2=AvrgPropOf(Name="SF",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_SFbeta2$Betavalues="beta:0.1"

Hyena_SFbeta3=AvrgPropOf(Name="SF",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_SFbeta3$Betavalues="beta:0.33"

Hyena_SFbeta4=AvrgPropOf(Name="SF",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_SFbeta4$Betavalues="beta:0.5"


df=rbind(Hyena_SFbeta1,Hyena_SFbeta2,Hyena_SFbeta3,Hyena_SFbeta4)
n_hyena=35
plotSF_hyena<- ggplot(df,aes(x=Timesteps, y=value/n_hyena, group=Betavalues))+
  geom_line(show.legend = T,
            aes(linetype=Betavalues, color=Betavalues),size=2)+
  #geom_point(colour = alpha("blue", 0.5))
  geom_point(aes(color=Betavalues))+
  theme_classic()+labs(title="SF")+
  xlab("Timesteps (days)")+ylab("Prop-infected")+
  scale_y_continuous(limits = c(0,1))+
  scale_linetype_manual(values=c("longdash", "solid", "dotted", "solid"))+
  scale_color_manual(values=c('#999999','#E69F00','red','black'))+
  scale_size_manual(values=c(2,1,2,1.5))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title = element_text(size = 34,hjust = 0.2,vjust = -0.6, face = "bold"))

plotSF_hyena

###########----SPATIAL-------##################
Hyena_SPbeta1=AvrgPropOf(Name="SP",ID=1,betaVal = 0.01,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_SPbeta1$Betavalues="beta:0.01"

Hyena_SPbeta2=AvrgPropOf(Name="SP",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_SPbeta2$Betavalues="beta:0.1"

Hyena_SPbeta3=AvrgPropOf(Name="SP",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_SPbeta3$Betavalues="beta:0.33"

Hyena_SPbeta4=AvrgPropOf(Name="SP",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,
                         Data = "hyena_vs_synthetic_netwok.csv")
Hyena_SPbeta4$Betavalues="beta:0.5"


df=rbind(Hyena_SPbeta1,Hyena_SPbeta2,Hyena_SPbeta3,Hyena_SPbeta4)
n_hyena=35
plotSP_hyena<- ggplot(df,aes(x=Timesteps, y=value/n_hyena, group=Betavalues))+
  geom_line(show.legend = T,
            aes(linetype=Betavalues, color=Betavalues),size=2)+
  #geom_point(colour = alpha("blue", 0.5))
  geom_point(aes(color=Betavalues))+
  theme_classic()+labs(title="SP")+
  xlab("Timesteps (days)")+ylab("Prop-infected")+
  scale_y_continuous(limits = c(0,1))+
  scale_linetype_manual(values=c("longdash", "solid", "dotted", "solid"))+
  scale_color_manual(values=c('#999999','#E69F00','red','black'))+
  scale_size_manual(values=c(2,1,2,1.5))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title = element_text(size = 34,hjust = 0.2,vjust = -0.6, face = "bold"))

plotSP_hyena

###########----LATTICE-------##################
Hyena_Latbeta1=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.01,gammaVal = 0.143,nreps = 10,
                          Data = "hyena_vs_synthetic_netwok.csv")
Hyena_Latbeta1$Betavalues="beta:0.01"

Hyena_Latbeta2=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,
                          Data = "hyena_vs_synthetic_netwok.csv")
Hyena_Latbeta2$Betavalues="beta:0.1"

Hyena_Latbeta3=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,
                          Data = "hyena_vs_synthetic_netwok.csv")


Hyena_Latbeta3$Betavalues="beta:0.33"

Hyena_Latbeta4=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,
                          Data = "hyena_vs_synthetic_netwok.csv")
Hyena_Latbeta4$Betavalues="beta:0.5"


df=rbind(Hyena_Latbeta1,Hyena_Latbeta2,Hyena_Latbeta3,Hyena_Latbeta4)
n_hyena=35
plotLat_hyena<- ggplot(df,aes(x=Timesteps, y=value/n_hyena, group=Betavalues))+
  geom_line(show.legend = T,
            aes(linetype=Betavalues, color=Betavalues),size=2)+
  #geom_point(colour = alpha("blue", 0.5))
  geom_point(aes(color=Betavalues))+
  theme_classic()+labs(title="Lat")+
  xlab("Timesteps (days)")+ylab("Prop-infected")+
  scale_y_continuous(limits = c(0,1))+
  scale_linetype_manual(values=c("longdash", "solid", "dotted", "solid"))+
  scale_color_manual(values=c('#999999','#E69F00','red','black'))+
  scale_size_manual(values=c(2,1,2,1.5))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title = element_text(size = 34,hjust = 0.2,vjust = -0.6, face = "bold"))

plotLat_hyena

###########----HYENA-------##################
Hyenabeta1=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.01,gammaVal = 0.143,nreps = 10,
                      Data = "hyena_vs_synthetic_netwok.csv")
Hyenabeta1$Betavalues="beta:0.01"

Hyenabeta2=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.1,gammaVal = 0.143,nreps = 10,
                      Data = "hyena_vs_synthetic_netwok.csv")
Hyenabeta2$Betavalues="beta:0.1"

Hyenabeta3=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.33,gammaVal = 0.143,nreps = 10,
                      Data = "hyena_vs_synthetic_netwok.csv")
Hyenabeta3$Betavalues="beta:0.33"

Hyenabeta4=AvrgPropOf(Name="Lat",ID=1,betaVal = 0.5,gammaVal = 0.143,nreps = 10,
                      Data = "hyena_vs_synthetic_netwok.csv")
Hyenabeta4$Betavalues="beta:0.5"


df=rbind(Hyenabeta1,Hyenabeta2,Hyenabeta3,Hyenabeta4)
n_hyena=35

plotHyena<- ggplot(df,aes(x=Timesteps, y=value/n_hyena, group=Betavalues))+
  geom_line(show.legend = T,
            aes(linetype=Betavalues, color=Betavalues),size=2)+
  #geom_point(colour = alpha("blue", 0.5))
  geom_point(aes(color=Betavalues))+
  theme_classic()+labs(title="Hyena")+
  xlab("Timesteps (days)")+ylab("Prop-infected")+
  scale_y_continuous(limits = c(0,1))+
  scale_linetype_manual(values=c("longdash", "solid", "dotted", "solid"))+
  scale_color_manual(values=c('#999999','#E69F00','red','black'))+
  scale_size_manual(values=c(2,1,2,1.5))+
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        plot.title = element_text(size = 34,hjust = 0.2,vjust = -0.6, face = "bold"))

plotHyena

plot1_hyena=ggarrange(plotER_hyena,plotSW_hyena,nrow=1)
plot2_hyena=ggarrange(plotLat_hyena,plotSF_hyena,nrow=1)
plot3_hyena=ggarrange(plotSP_hyena,plotHyena,common.legend = TRUE, legend="bottom",nrow=1)

plotsallhyena=ggarrange(plot1_hyena,plot2_hyena,plot3_hyena,nrow=3,
                        common.legend = TRUE, legend="bottom")

Allannotatehyena=annotate_figure(
  plotsallhyena,
  top = text_grob("Hyena vs synthetic network-2 ",
                  color = "black", face = "bold", size = 30),
  fig.lab = "b)",fig.lab.size = 40)

ggsave("allbeta_allnetwork.png", width = 24, height = 22)

##############----VARYING BETA FOR HYENA AND DOLPHIN----#############
Allplot=ggarrange(Allannotatedolph,Allannotatehyena,nrow=2)
ggsave("allbeta_hyena&dolph_network.png", width = 24, height = 24)