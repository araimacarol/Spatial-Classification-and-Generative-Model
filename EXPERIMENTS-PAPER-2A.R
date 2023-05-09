#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Set working directory
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
setwd("C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi-New")
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Loading packages
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

suppressPackageStartupMessages({
  if(!require(Matrix)) {
    install.packages("Matrix")
    library(Matrix)
  }
  if(!require(here)) {
    install.packages("here")
    library(here)
  }
  if(!require(stats)) {
    install.packages("stats")
    library(stats)
  }
  if(!require(igraph)) {
    install.packages("igraph")
    library(igraph)
  }
  if(!require(tidymodels)) {
    install.packages("tidymodels")
    library(tidymodels)
  }
  if(!require(tidyverse)) {
    install.packages("tidymodels")
    library(tidyverse)
  }
  if(!require(ggplot2)) {
    install.packages("ggplot2")
    library(ggplot2)
  }
  if(!require(janitor)) {
    install.packages("janitor")
    library(janitor)
  }
  if(!require(vip)) {
    install.packages("vip")
    library(vip)
  }
  if(!require(readxl)) {
    install.packages("readxl")
    library(readxl)
  }
  if(!require(naniar)) {
    install.packages("naniar")
    library(naniar)
  }
  if(!require(furrr)) {
    install.packages("furrr")
    library(furrr)
  }
})

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Helper functions
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

##########################################################################################################
# Input an igraph object from a file, treating it as a static graph
##########################################################################################################
getGraphFromFile <- function(file, simplify=TRUE, useBiggestComponent=TRUE, asUndirected=TRUE) {
  
  dat <- read.table(file) # just read static graph: ignoring third column
  G <- graph_from_data_frame(dat)
  if (asUndirected==TRUE) {
    G <- as.undirected(G, "collapse")
  }
  
  #  g_names <- gsub(".edges","",networks[i]) # one edge for each pair of connect vertices (not sure what this is for)
  
  if (useBiggestComponent==TRUE) {
    netclust <- clusters(G) #look for subgraphs
    gcc <- V(G)[netclust$membership == which.max(netclust$csize)]#select vertices from the largest sub-graph
    G <- induced.subgraph(G, gcc) #make it a igraph object.
  }
  if (simplify==TRUE) {
    G <- igraph::simplify(G, remove.multiple = TRUE, remove.loops = TRUE)
  }
  return(G)
}

###--Normalized Laplacina function--##
normalized_laplacian=function(Graphs){
  laplacian_matrix(Graphs,normalized = T)
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Calculate graph features
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
calcGraphFeatures <- function(Graphs=NULL) {
  
  features <- c(
    "order",                    # number of vertices
    "edges",                     # number of edges
    "connected",                # True / False
    "max_component",            # maximum component size (=order iff the graph is connected)
    "minDegree",                # minimum degree of any vertex
    "maxDegree",                # maximum degree of any vertex
    "mean_degree",                # average degree of any vertex
    "minCut",                   # minimum cut weight of the graph (might take a while to compute)
    "FiedlerValue",             # second-highest eigenvalue of the Laplacian matrix
    "Normalized_FiedlerValue",   # second-highest eigenvalue of the Normaized Laplacian matrix
    "closeness_centr",                # average inverse of distance between any pair of vertices
    "modularity",               # DEFINITION REQUIRED
    "diameter",                 # maximum distance between any two vertices (NAN if not connected)
    "betw_centr",              # max_{v} proportion of shortest paths going through vertex v
    "transitivity",             # aka Clustering Coefficient, is proportion of connected triples that form triangles: e.g., (a--b--c--a) when (a--b--c) is present.
    "threshold",                 # 1/max(eigen value of A)
    "spectral_radius"         # max (eigen value of A)
    
  )
  
  df <- as.data.frame(matrix(ncol=length(features),nrow=length(Graphs)))
  colnames(df)=features
  
  # Stuff that is simple to apply and needs no interim components:
  
  df$order = base::as.numeric(lapply(Graphs, gorder))
  df$edges = base::as.numeric(lapply(Graphs, gsize))
  df$connected = base::as.numeric(lapply(Graphs, is_connected))
  df$minCut = base::as.numeric(lapply(Graphs, min_cut))
  df$diameter = base::as.numeric(lapply(Graphs, diameter))
  df$transitivity = base::as.numeric(lapply(Graphs, transitivity))
  
  # stuff that needs interim things:
  degrees = lapply(Graphs,igraph::degree )
  df$minDegree = base::as.numeric(lapply(degrees, min))
  df$maxDegree = base::as.numeric(lapply(degrees, max))
  df$mean_degree = base::as.numeric(lapply(degrees, mean))
  
  # stuff that lapply doesn't like so has to be done in a loop:
  communities <-lapply(Graphs, cluster_walktrap) #lapply(Graphs, cluster_leading_eigen)
  Adj <- lapply(Graphs, as_adjacency_matrix)
  L <- lapply(Graphs, laplacian_matrix)
  Norm_Lap<-lapply(Graphs, normalized_laplacian)
  Fiedler.value=NULL
  norm.fiedler.value=NULL
  
  for (i in 1:length(Graphs)) {
    if (is.null(Graphs[[i]]$type)) { Graphs[[i]]$type = "untyped" }
    df$modularity[i] <- modularity(communities[[i]])
    df$spectral_radius[i] <- eigen(Adj[[i]], symmetric=TRUE, only.values=TRUE)$values[1]
    
    Fiedler.value[[i]]=eigen(L[[i]], symmetric=TRUE, only.values=TRUE)$values
    
    df$FiedlerValue[i] <- Fiedler.value[[i]][length(Fiedler.value[[i]])-1]
    
    norm.fiedler.value[[i]]=eigen(Norm_Lap[[i]], symmetric=TRUE, only.values=TRUE)$values
    
    df$Normalized_FiedlerValue[i] <- norm.fiedler.value[[i]][length(norm.fiedler.value[[i]])-1]
    
    df$eigen_centr[i] <- centr_eigen(Graphs[[i]])$centralization
    df$deg_centr[i] <- centr_degree(Graphs[[i]])$centralization
    
    df$betw_centr[i] <- centr_betw(Graphs[[i]])$centralization
    
    df$max_component[i] <- max(components(Graphs[[i]])$csize)
    df$mean_eccentr[i]<-mean(eccentricity(Graphs[[i]]))
    df$radius[i]<-radius(Graphs[[i]])
    df$mean_path_length[i]<-average.path.length(Graphs[[i]])
    #df$trace[i]<-sum(diag(Adj[[i]]))
    df$graph_energy[i]<-sum(abs(eigen(Adj[[i]], symmetric=TRUE, only.values=TRUE)$values))
    df$min_triangle[i]= min(count_triangles(Graphs[[i]]))
    df$mean_triangle[i]= mean(count_triangles(Graphs[[i]]))
    df$sd_triangle[i]= sd(count_triangles(Graphs[[i]]))
    df$max_triangle[i]= max(count_triangles(Graphs[[i]]))
    df$num_triangle[i]= sum(count_triangles(Graphs[[i]]))
    df$deg_assort_coef[i]=assortativity_degree(Graphs[[i]])
    
    df$threshold[i] <- 1/(df$spectral_radius[i])
    
    if (df$connected[i]==TRUE) {
      df$closeness_centr[i] = mean(closeness(Graphs[[i]]))
    } else { # handle the case where G isn't connected
      df$closeness_centr[i] = -1
    }
  }
  return (df)
}

#calcGraphFeatures(x)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Multiple Calculation of graph features
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Used to perform runs on multiple simulated graphs on any given network
RunSimOnGraphFeatures<-function(Graphs, nreps=nreps,output_file=NULL, seed=-1) {
  set.seed(1)
  # ### Definition and initialization of parameters for graphfeatures
  graphProperties=list()
  # ### Definition and initialization of Graph Prefix
  graphid=list();graphreplicate=list(); graphname=list(); GraphPrefix=list(); analysis=list()
  for (g in 1:length(Graphs)){
    for (reps in 1:nreps) {
      ### Calculate the graph features for each simulated graph of all the synthetic networks
      print(paste("Calculating graph features on", Graphs[[g]]$name))
      #graphProperties[[reps]] <- calcGraphFeatures(Graphs[g])
      graphProperties[[reps]] <- calcGraphFeatures(Graphs[g])
    }
    graphname[[g]]=Graphs[[g]]$type
    graphid[[g]]=Graphs[[g]]$id
    graphreplicate[[g]]=c(1:nreps)
    GraphPrefix=cbind(graphname[[g]],graphid[[g]],graphreplicate[[g]])
    colnames(GraphPrefix)=c("GraphName","GraphID","GraphReplicate")
    analysis[[g]]=as.data.frame(cbind(GraphPrefix,graphProperties[[reps]]))
    row.names(analysis[[g]])=1:nreps
  }
  All_results=do.call(rbind,analysis)
  # write.csv(All_results, file=output_file)
  return( All_results)
}


# dplyr::select(c(GraphName,order,edges,mean_eccentr,
#                 radius,mean_path_length,graph_energy,min_triangle,mean_triangle,
#                 sd_triangle,modularity,diameter,betw_centr,transitivity,threshold,
#                 spectral_radius,eigen_centr,deg_centr,
#                 mean_degree,minCut,FiedlerValue,Normalized_FiedlerValue,
#                 closeness_centr),max_triangle,num_triangle,deg_assort_coef)%>%
#   mutate_if(is.character,factor)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Data processing and cleaning
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
df1=read.csv("df.lat.csv")%>%
filter(order<=1000)


df1=df1%>%dplyr::select(c(GraphName,order,edges,
                          mean_eccentr,mean_path_length,graph_energy,
                          modularity,diameter,betw_centr,transitivity,
                          spectral_radius,eigen_centr,deg_centr,
                          mean_degree,minCut,FiedlerValue,Normalized_FiedlerValue,
                          closeness_centr))%>%
  mutate_if(is.character,factor)

# df2=rbind(read.csv("sbm.50.csv",sep = ",", header = T),
#           read.csv("sbm.100.csv",sep = ",", header = T),
#           read.csv("sbm.150.csv",sep = ",", header = T),
#           read.csv("sbm.150.csv",sep = ",", header = T))

df2=read.csv("df.sbm.csv")%>%
  filter(order<=1000)


df2=df2%>%dplyr::select(c(GraphName,order,edges,
                          mean_eccentr,mean_path_length,graph_energy,
                          modularity,diameter,betw_centr,transitivity,
                          spectral_radius,eigen_centr,deg_centr,
                          mean_degree,minCut,FiedlerValue,Normalized_FiedlerValue,
                          closeness_centr))%>%
  mutate_if(is.character,factor)

df3=rbind(read.csv("df.er.50.csv",sep = ",", header = T),
          read.csv("df.er.100.csv",sep = ",", header = T),
          read.csv("df.er.150.csv",sep = ",", header = T),
          read.csv("df.er.200.csv",sep = ",", header = T),
          read.csv("df.er.250.csv",sep = ",", header = T),
          read.csv("df.er.350.csv",sep = ",", header = T),
          read.csv("df.er.400.csv",sep = ",", header = T),
          read.csv("df.er.500.csv",sep = ",", header = T),
          read.csv("df.er.750.csv",sep = ",", header = T))


df3=df3%>%dplyr::select(c(GraphName,order,edges,
                          mean_eccentr,mean_path_length,graph_energy,
                          modularity,diameter,betw_centr,transitivity,
                          spectral_radius,eigen_centr,deg_centr,
                          mean_degree,minCut,FiedlerValue,Normalized_FiedlerValue,
                          closeness_centr))%>%
  mutate_if(is.character,factor)


df4=rbind(read.csv("df.SP.50.csv",sep = ",", header = T),
          read.csv("df.SP.100.csv",sep = ",", header = T),
          read.csv("df.SP.150.csv",sep = ",", header = T),
          read.csv("df.SP.200.csv",sep = ",", header = T),
          read.csv("df.SP.250.csv",sep = ",", header = T),
          read.csv("df.SP.300.csv",sep = ",", header = T),
          read.csv("df.SP.350.csv",sep = ",", header = T))

df4=df4%>%dplyr::select(c(GraphName,order,edges,
                          mean_eccentr,mean_path_length,graph_energy,
                          modularity,diameter,betw_centr,transitivity,
                          spectral_radius,eigen_centr,deg_centr,
                          mean_degree,minCut,FiedlerValue,Normalized_FiedlerValue,
                          closeness_centr))%>%
  mutate_if(is.character,factor)

df5=rbind(read.csv("SF1.50.csv",sep = ",", header = T),
          read.csv("SF1.100.csv",sep = ",", header = T),
          read.csv("SF1.150.csv",sep = ",", header = T),
          read.csv("SF1.200.csv",sep = ",", header = T))

df5=df5%>%rename("power"="prewireVals","mvals"="neiVals")

df5=df5%>%
  filter(power %in% c("2") & GraphReplicate <= 250)

df5=df5%>%dplyr::select(c(GraphName,order,edges,
                          mean_eccentr,mean_path_length,graph_energy,
                          modularity,diameter,betw_centr,transitivity,
                          spectral_radius,eigen_centr,deg_centr,
                          mean_degree,minCut,FiedlerValue,Normalized_FiedlerValue,
                          closeness_centr))%>%
  mutate_if(is.character,factor)


df6=rbind(read.csv("SW.50.csv",sep = ",", header = T),
          read.csv("SW.100.csv",sep = ",", header = T),
          read.csv("SW.150.csv",sep = ",", header = T),
          read.csv("SW.200.csv",sep = ",", header = T),
          read.csv("SW.250.csv",sep = ",", header = T),
          read.csv("SW.300.csv",sep = ",", header = T))


df6=df6%>%
  filter(prewireVals %in% c("0.2") & GraphReplicate <= 250)
#df6[df6$GraphReplicate <= 250,]
df6=df6%>%dplyr::select(c(GraphName,order,edges,
                          mean_eccentr,mean_path_length,graph_energy,
                          modularity,diameter,betw_centr,transitivity,
                          spectral_radius,eigen_centr,deg_centr,
                          mean_degree,minCut,FiedlerValue,Normalized_FiedlerValue,
                          closeness_centr))%>%
  mutate_if(is.character,factor)




Data=rbind(df1,df2,df3,df4,df5,df6)

##--Shuffle--data
Data = Data[sample(1:nrow(Data)), ]##shuffle row indices and randomly re order 
Data[Data==0] <- NA
Data%>%na.omit()
#print(na.omit(Data))
#replace(is.na(.), 0)#change NA's in deg_assort_coef to 0

# Change data frame to tibble (similar to data-frame)
df=as_tibble(Data)
df%>%count(GraphName,sort=T)

# Show first n row of the data
df%>%
  slice_head(n=10)
#df%>%View()

###--------Load the janitor package to clean the data---------###
df<-df%>%
  #mutate(GraphName=factor(GraphName,levels=c("ER","sbm","LAT")))%>%
  clean_names()

df%>%
  glimpse()

###----Data Partitioning----###
set.seed(384)

##sample and split
df_split=df[sample(1:nrow(df)), ]%>%
  na.omit()%>%
  #  replace(is.na(.), 0)%>%
  initial_split(prop = 0.7)

df.train=training(df_split)  
df.test=testing(df_split)



###----Print number of training and testing set----###
cat("training set:", nrow(df.train), "\n",
    "testing set :", nrow(df.test), "\n")

#########################################################################################
# Building models
#########################################################################################

##cross validation
df.cv.splits <- vfold_cv(df.train, v = 10)

##Preparing data
df.rec <- recipe(graph_name ~ ., data = df) %>%
  #--convert categorical variables to factors
  #recipes::step_string2factor(all_nominal()) %>%
  #-combine low frequency factor levels
  #recipes::step_other(all_nominal(), threshold = 0.01) %>%
  # remove no variance predictors which provide no predictive information 
  recipes::step_nzv(all_nominal()) %>%
  recipes::step_zv(all_predictors()) %>%
  recipes::step_normalize(all_numeric_predictors())
  #step_dummy(all_predictors())%>%
  
df.prep=df.rec%>%prep()
  

#+++++++++++++++++++++++++++++++++++++++++++
#+ MODELS
#+ +++++++++++++++++++++++++++++++++++++++++

####----[1]Random forest----####
rf.tune.df <- rand_forest(mtry = tune(), trees = tune()) %>%
  set_engine("ranger") %>%
  set_mode("classification")
# Hyperparameter grid
rf.grid <- rf.tune.df %>%
  parameters() %>%
  finalize(select(df, -graph_name)) %>%  
  grid_max_entropy(size = 10)
# Workflow bundling every step 
rf.wflow <- workflow() %>%
  add_recipe(df.prep) %>%
  add_model(rf.tune.df)

####----[2]Boost Trees----####
xgb.tune.df <- boost_tree(mtry = tune(), tree = tune(),
                             learn_rate = tune(), tree_depth = tune()) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

xgb.grid <- xgb.tune.df %>%
  parameters() %>%
  finalize(select(df, -graph_name)) %>%  
  grid_max_entropy(size = 10)

xgb.wflow <- workflow() %>%
  add_recipe(df.prep) %>%
  add_model(xgb.tune.df)

####----[3] Neural Network----####
nnet.tune.df <- mlp(hidden_units = tune(), penalty = tune(), activation = "relu") %>%
  set_engine("keras") %>%
  set_mode("classification")

nnet.grid <- nnet.tune.df %>%
  parameters() %>%
  grid_max_entropy(size = 10)

nnet.wflow <- workflow() %>%
  add_recipe(df.prep) %>%
  add_model(nnet.tune.df)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Fitting all Model
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
wflow_list <- list(rf.wflow,xgb.wflow,nnet.wflow)
#save trained workflows for the models
saveRDS(wflow_list,"wflow_list.rds")
grid_list <- list(rf.grid,xgb.grid,nnet.grid)
#save trained grid list
saveRDS(grid_list,"grid_list.rds")
plan(multiprocess, workers = 3)

trained_models_list <- future_map2(.x = wflow_list,
        .y = grid_list,
        ~tune_grid(.x , resamples = df.cv.splits , grid = .y))
#save trained models
saveRDS(trained_models_list,"trained_models_list.rds")

##----Checking best performing models
TrainedModelslist=readRDS("trained_models_list.rds")
TrainedModelslist[1:2] %>%
  map(show_best, metric = "accuracy", n = 1)

#extracting the best model
xgb.specs <- TrainedModelslist[[2]]

#saving best model specification
best.xgb.spec <- show_best(xgb.specs, "accuracy", 1)

#retrain best model
best.xgb.model <- boost_tree(mode = "classification", mtry = best.xgb.spec$mtry,
                             trees = best.xgb.spec$trees) %>%
                                       set_engine("xgboost")

df.final.rec <- recipe(graph_name ~ ., data = df) %>%
  #step_dummy(all_predictors())
  recipes::step_nzv(all_nominal()) %>%
  recipes::step_zv(all_predictors()) %>%
  recipes::step_normalize(all_numeric_predictors())

df.best.wflow <- workflow() %>%
  add_recipe(df.final.rec) %>%
  add_model(best.xgb.model)

best.fitted.model <- fit(df.best.wflow, data = df.train)

#confusion matrix of best model
best.model.predictions <- predict(best.fitted.model, new_data = df.test) %>%
  bind_cols(df.test)

#accuracy of best model
best.model.predictions %>%
  #mutate(graph_name = as.factor(graph_name)) %>%  
  accuracy(graph_name, .pred_class)

#confuion matrix of best model
best.model.predictions %>%
  #mutate(graph_name = as.factor(graph_name)) %>%  
  conf_mat(graph_name, .pred_class)

#If the predicted class and thruth were really close, we would need to use
#SMOTE() function in tidymodel to balance data but were good for our model so no need
# To train to make the model more accurate in a predicted class, 
# We will resample the training set, but by downsampling class ER/SP/SF/SBM and 
#upsampling class LAT.
# This can be done with the function SMOTE() from the {DMwR} package. However, 
# the testing set should have the same distribution as the population, so we can 
#not apply SMOTE() to the testing set. I will resplit the data, but this time with a
#95/5 % percent split; this way I have
# 5% of the original dataset used for testing, I can use SMOTE() on the 95% remaining training set

#eg
#train_smote = DMwR::SMOTE(graph_name ~ ., df.train, perc.over = 100, perc.under=200)


#########################################################################################
# Date processing, exploration and feature selection
#########################################################################################
# We need a  workflow, which we have already done
# to deal with the data preprocessing, explanation and then the training of the model
# To get explanation of data, we need a model and a way to get prediction
# The fitted workflow is able to make predictions
# We will use iml package which have several functions for qualitative anaylyses of the data
# and explanability

# We have to first define an object that takes the fitted final model as input 
#the design matrix, the target variable and the prediction function

library("iml")

#feature set
features <- df.test %>%
  select(-graph_name)

#target variable
target <- df.test %>%
 # mutate(job_search = as.factor(job_search)) %>%  
  select(graph_name)

#prediction function
predict.wrapper <- function(model, newdata){
  workflows:::predict.workflow(object = model, new_data = newdata)
}

#Because a workflow is a bit special, 
#we need to define this wrapper function that wraps the workflows:::predict.workflow() function
#predict() works but had some issues

##The predictor function enables to interrogate or get explanation of the data
predictorfunc <- Predictor$new(
  model = best.fitted.model,
  data = features, 
  y = target,
  predict.fun = predict.wrapper
)

##feature importance
feature.importance <- FeatureImp$new(predictorfunc, loss = "ce")
plot(feature.importance)

##########################################################
#ALE ANALYSIS
##########################################################


##########################################################
#SHAPELY ANALYSIS
##########################################################







####----Notes:Fitting Models with tidymodels----####
#--[1] fit()<- is used to train the model on the training data but not perform cross validation. eg
#fitted_model <- fit(model_formula, data = data_train)

#--[2] fit_resample()<- fits a model specification 
#(without varying hyper-parameters) to all the analysis sets contained in 
#the cross validation fold object (which contains the resampled training data for cross-validation) 
# however fit_resample() performs crossvaldation fit but not hyper parametr tunning
#eg
# fitted_resamples <- fit_resamples(model_wflow,
#     resamples = my_cv_splits,
#   control = control_resamples(save_pred = TRUE))

#--[3] tune_grid()<-performs cross validation and hyperparameter tunning
#eg 
# tuned_model <- tune_grid(model_wflow,
# resamples = my_cv_splits,
# grid = my_grid,
# control = control_resamples(save_pred = TRUE))

#map2()<- maps the +() function to each element of both vectors successively; so can be use
#to tune grids for mor than one model
#eg
#map2(c(1, 1, 1), c(2,2,2), `+`)
## [[1]]
## [1] 3
## 
## [[2]]
## [1] 3
## 
## [[3]]
## [1] 3

#note<-furrr::future_map2() will run one model per core, 
#and the way to do it is to simply define how many cores
# so 3 cores for 3 model

#setwd(dir = "https://networkrepository.com/asn.php")

# download a .zip file of the repository
# from the "Clone or download - Download ZIP" button
# on the GitHub repository of interest
#url="https://nrvis.com/download/data/asn/aves-barn-swallow-contact-network.zip"

# urls <- h %>% 
#   html_nodes('area') %>%    # get all `area` nodes
#   html_attr('href') %>%    # get the link attribute of each node
#   sub('.htm$', '.zip', .) %>%    # change file suffix
#   paste0('https://networkrepository.com/asn.php/.zip', .)    # append to base URL
# 
# destfile = "C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi-New/ANIMAL-SOCIAL-NETWORS-REPOSITORY/out.edges"
# #download.file(url,destfile)
# 
# 
# Map(function(u, d) download.file(u, d, mode="zip"), urls, destfile)
# # check it's there
# list.files(dir)
# 
# 
# 
#  library(rvest)
# # # get page source
#  h <- read_html('https://networkrepository.com/asn.php')
#  
#  urls <- h %>% 
#    html_nodes('area') %>%    # get all `area` nodes
#    html_attr('href') %>%    # get the link attribute of each node
#    sub('.htm$', '.zip', .) %>%    # change file suffix
#    paste0('https://networkrepository.com/', .)    # append to base URL
# # 
# # create a directory for it all
# dir <- file.path(tempdir(), 'asn')
# dir.create(dir)
# 
# # iterate and download
# lapply(urls, function(url) download.file(url, file.path(dir, basename(url))))
# 
# # check it's there
# list.files(dir)
