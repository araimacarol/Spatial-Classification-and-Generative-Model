install.packages("devtools")
#devtools::install_github('nfj1380/mrIML')
library(mrIML)
library(tidyverse)
library(tidymodels)
library(vip)
library(flashlight)
library(doParallel)
library(janitor)
library(viridis)
library(kernlab)
library(DataExplorer)
library(GGally)
library(ggfortify)
library(future.apply)

Simdata <- read.csv('NetworkDataSim.csv')

#note that - were changed to _ to avoid downstream issues in names. 
colnames(Simdata)[2] <- c('Graph.Name')
str(Simdata)  
Graph_metaData <- read.csv('Network_repository_metaData.csv')
str(Graph_metaData)

#convert characters to factors
Graph_metaData[sapply(Graph_metaData, is.character)] <- lapply(Graph_metaData[sapply(Graph_metaData, is.character)], 
                                       as.factor)
#merge

all_data <- left_join(Simdata, Graph_metaData , by = c("Graph.Name"))

str(all_data)

#remove features that are of limited value

All_data_reducedSet <- subset(all_data, select=-c(X,Mammal.,Citation)) #mammal not needed as we are looking across all groups

 
#remove missing data. SOme missing responses (removes 8 networks)
All_data_reducedSetNoNA<- All_data_reducedSet[complete.cases(All_data_reducedSet  ), ] 

#remove all networks with < 10 nodes
dataNoUnder10<- All_data_reducedSetNoNA  %>% filter( network_size >= 10) #removes 80 graphs

#most_infected_node removed as not a meaningful predictor in these models

Xworking <- dataNoUnder10[-c(1:7)] #str(Y)

##############################################################################################
#PCA 
##############################################################################################

pca_data <- Xworking[1:9]
pca_res <- prcomp(pca_data, scale. = TRUE)

autoplot(pca_res)

#with labels to work out outliers

autoplot(pca_res, data = Y, colour = 'Species',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3,
         label = TRUE)

#not labels

autoplot(pca_res, data = Y, colour = 'Species',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)

#Outlier (459) = complete voles network
#Outlier (500) ~complete tortoise network.

#kmeans

autoplot(kmeans(pca_data, 3), data = pca_data, label = TRUE, label.size = 3)

##############################################################################################
#Correlations
#############################################################################################

#species has too many levels so remove for the time being
Xfiltered <- subset(Xworking, select=-c(Source, Species, Edge_weight, most_infected_node)) #these features aren't useful here
ggpairs(Xfiltered) #lots of strong correlations

dataNoCOr <- Xfiltered %>%  select(-c('transitivity.gcc_s..type....global..','Rnot', 'Adj_val','Diameter', 'Mean_degree' ))
ggpairs(dataNoCOr)
#modularity and mean degree strongly correlated for example

#check out our responses
Ys <- dataNoUnder10[c(2:7)] 
ggpairs(Ys)
#all highly correlated   - create one model

create_report(Yworking)
create_report(Xs)
#beta 0.5 and beta 0.1 estimates strongly correlated. Choose 0.05 


##############################################################################################
#ML models
##############################################################################################

model1 <- 
  rand_forest(trees = 1000, mode = "regression", mtry = tune(), min_n = tune()) %>% #100 trees are set for brevity
  set_engine("ranger", importance = "impurity") #random forest doesn't work well with dummy variables

model1 <- linear_reg() %>% 
  set_engine("lm") %>% 
  set_mode("regression")

library(kernlab)
model1<-
  svm_rbf(cost = tune(), rbf_sigma = tune()) %>%
  set_engine("kernlab") %>%
  set_mode("regression")

#
Y <- Ys %>% select(avg_prop_infected_beta0.05 )

#X <- dataNoUnder10 %>% select(avg_prop_infected_beta0.05)

# Define set the outcomes of interest

Xraw <- dataNoCOr #up to hear  - need to add species again (no -s)

#add species as a predictor

X <-  bind_cols(Xraw,species=dataNoUnder10$Species); str(X) 

#If needed - need to upadte MrIML
library(fastDummies)

Xdummy <- dummy_cols(X, remove_first_dummy = TRUE, remove_selected_columns = TRUE)

#X %>% 
 # rename(`Captive_Semi-ranging` = Captive_Semi_ranging) #this isnt working hyphen

#colnames(Ydummy[6]) <- 'Captive_Semi_ranging' #r doesn't like '-' sometimes

#set up multicore
cl <- parallel::makeCluster(5)
plan(cluster, workers=cl)

yhats <- mrIMLpredicts(X=X,Y=Y , model1=model1, balance_data='no', mod='regression', tune_grid_size = 10, see=123)



saveRDS(yhats , 'yhT_avg_invasion_time_beta0.05_Sept2021')

#save(yhats, file='logreg_model')
ModelPerf <- mrIMLperformance(yhats, model1, Y=Y) 
ModelPerf 

VI <- mrVIptest (yhats, X=X)

 p <- plot_vi(VI=VI,  X=X,Y=Y, modelPerf=ModelPerf, cutoff= 0)+theme_bw()
 p <- plot_vi(VI=VI,  X=X,Y=Ydummy, modelPerf=ModelPerf, cutoff= 0)+theme_bw() #checks which particular species etc
 
flashlightObj <- mrFlashlight(yhats, Y=Y,X=X, response = "single", model='regression')

#plot prediction scatter for all responses.

plot(light_profile(flashlightObj, v = "Modularity", type = "ale",  n_bins =25))+
  theme_bw()+
  geom_rug(data=Y, aes(x =Modularity), inherit.aes = F)

plot(light_profile(flashlightObj, v = "network_size", type = "ale", n_bins =100))+
  theme_bw()+
  geom_rug(data=Y, aes(x =network_size), inherit.aes = F)

plot(light_profile(flashlightObj, v = "Centrality", type = "ale",  n_bins =50))+
  theme_bw()+
  geom_rug(data=Y, aes(x =Centrality), inherit.aes = F)
  
plot(light_profile(flashlightObj, v = "Fiedler", type = "ale", n_bins =100))+
  theme_bw()+
  geom_rug(data=Y, aes(x=Fiedler), inherit.aes = F)

plot(light_ice(flashlightObj, v = "Fiedler", center = "first"))+
  theme_bw()+
  geom_rug(data=Y, aes(x =Fiedler), inherit.aes = F)

plot(light_ice(flashlightObj, v = "Centrality", center = "first"))+
  theme_bw()+
  geom_rug(data=Y, aes(x = Centrality), inherit.aes = F)

plot(light_profile2d(flashlightObj, v = c("network_size", "Centrality")))

plot(light_profile2d(flashlightObj, v = c("network_size", "Mean_degree")))
geom_rug(data=Y, aes(x = network_size, y=Mean_degree), inherit.aes = F)+
  scale_fill_continuous(high = "#56B1F7", low = "#132B43")

plot(light_profile2d(flashlightObj, v = c("network_size", "Fiedler")))+
  geom_rug(data=Y, aes(x = network_size, y=Fiedler), inherit.aes = F)+
  #scale_fill_gradientn(colours = colorspace::heat_hcl(7))
  scale_fill_continuous(high = "#132B43", low = "#56B1F7")

plot(light_global_surrogate(flashlightObj))

st <- light_interaction(flashlightObj, grid_size = 30, n_max = 50, seed = 42, pairwise=T)+
  theme_bw()
plot(st)
