setwd("C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi-New")
source("SPATIAL-PIPELINE-NEW-MODEL.R")
# mu0.50=Sim_SpatialGraphNonComm(mu=rep(0,300),N=50)
# mu11.50=Sim_SpatialGraphNonComm(mu=rep(0.1,300),N=50)
# mu12.50=Sim_SpatialGraphNonComm(mu=rep(0.4,300),N=50)
# mu13.50=Sim_SpatialGraphNonComm(mu=rep(0.7,300),N=50)
# mu14.50=Sim_SpatialGraphNonComm(mu=rep(1,300),N=50)
# 
# ######## Saving all 50 nodes data to csv
# df_mu_50.new=rbind(mu0.50,mu11.50,mu12.50,mu13.50,mu14.50)
# dim(df_mu_50.new)
# df_mu_50.new<-df_mu_50.new %>%dplyr::filter(order==50)
# dim(df_mu_50.new)
# write.csv(df_mu_50.new,"df_mu_50.new.csv")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [1] Spatial Expander Propagation Graph with weighted edges, 
# node attributes, scale free degree distribution, small world effect and community structures
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#library(quadtree)
# Function to generate a spatial scale-free expander graph
# with the given parameters
#
# Arguments:
#   N: number of nodes
#   r: cutoff distance for adjacency matrix
#   beta: power law parameter for distance-dependent probability
#   m: number of edges added for each new node
#   alpha: parameter controlling strength of degree effect
#   mu: parameter controlling preference when connecting nodes. Thus, preference to care about 
#       spatial distance between nodes when connecting them, vs caring only about vertex degree. 0:care aboutdegree
#1: care about distance
#   prewire: rewiring probability for small world effect
#   node_attrs: list of node attributes to add to graph
#   edge_weights: logical indicating whether to add edge weights
#
# Returns:
#   A graph object of class igraph

spatial_scale_free_expander_noncomm <- function(N=50, beta=2, m=2, prewire = 0.1, 
                                                mu=0.2,add_edge_weight = FALSE, 
                                                add_node_attr = FALSE,r=0.1,alpha=0.2) {
  #Notes:
  ##--(1) use of of N or lambda
  ##--(2) fix L or ignore
  ###--(3) mu > 0 for preferential attachment: mu captures what p does in the doi:10.1088/1367-2630/9/6/190
  ###--(4) use unit square instead of torus
  
  ##--Generate points on a unit square
  points <- matrix(runif(N*2), ncol=2)
  # calculate distance matrix between all pairs of points
  dist_mat <- as.matrix(dist(points))
  # initialize graph with a single node and no edges
  graph <- list()
  graph$adj_mat <- matrix(0, nrow = N, ncol = N)
  graph$adj_mat[1, 1] <- 1
  # initialize the edge_list
  graph$edge_list <- as.data.frame(matrix(ncol=2, nrow=0))
  
  # Grow the network with new nodes added one at a time
  for (i in 3:N) {
    # preferential attachment effect with scale free technique
    graph$degrees= rowSums(graph$adj_mat[1:(i-1), ])
    deg_probs<- (graph$degrees[1:(i-1)]^alpha)/(sum(graph$degrees[1:(i-1)]^alpha))
    
    #--spatial distance effect incorporating short spatial distances (cut-off distance)
    spatial_probs <- ifelse(dist_mat[i,] <= r, 1/(1+(dist_mat[i,])^beta), 1/N)
    
    # total probability of attachment 
    totalprobs=((mu)*spatial_probs[1:(i-1)]) + ((1-(mu))*deg_probs)
    #--normalise probabilities
    #probattach=Total_ij/sum(Total_ij)
    
    
    #Add m edges to existing nodes with probability proportional to their degree and distance
    node_to_attach <- sample(1:(i-1), size = m,replace = TRUE, prob = totalprobs)
    #node_to_attach <- sample(1:(i-1), size = m,replace = TRUE, prob = probattach)
    
    graph$adj_mat[i, node_to_attach] <- 1
    graph$adj_mat[node_to_attach, i] <- 1
    #graph$degrees[i] <- graph$degrees[i] + 1
    #graph$degrees[node_to_attach] <- graph$degrees[node_to_attach] + 1
    
    #Bind to the empty edge_list
    graph$edge_list <- rbind(graph$edge_list, c(i, node_to_attach))
    
    # Small-world rewiring with probability p
    if (runif(1) < prewire) {
      # Select random node to rewire within cutoff distance
      neighbors <- which(graph$adj_mat[i,] == 1)
      if (length(neighbors) > 0) {
        new_neighbor <- sample(neighbors, 1)
        graph$adj_mat[i, new_neighbor] <- 0
        graph$adj_mat[new_neighbor, i] <- 0
        available_nodes <- setdiff(1:N, c(i, neighbors))
        new_neighbor <- sample(available_nodes, 1)
        graph$adj_mat[i, new_neighbor] <- 1
        graph$adj_mat[new_neighbor, i] <- 1
      }
    }
  } 
  
  ### Graph object
  GraphModel <- graph.adjacency(as.matrix(graph$adj_mat), mode="undirected")
  GraphModel=igraph::simplify(GraphModel,remove.loops = TRUE,remove.multiple = TRUE)
  graph$GraphObject <- GraphModel
  
  # Compute graph Laplacian for the expansion property
  D <- diag(rowSums(graph$adj_mat))
  L <- D - graph$adj_mat
  # Compute eigenvalues and eigenvectors of Laplacian matrix
  eig <- eigen(Matrix(L))
  eig_values <- Re(eig$values)
  eig_vectors <- eig$vectors
  
  # Compute the expansion ratio
  lambda_2 <- eig_values[2]
  lambda_n <- eig_values[nrow(L)]
  expansion_ratio <- lambda_2 / lambda_n
  # Check for strong expansion properties
  #spectral_gap=eigen(Matrix::t(L) %*% L, symmetric=TRUE, only.values=TRUE)$values[n-1]
  if(!is.connected(GraphModel)||expansion_ratio <= 0){  
    warning("Graph may not exhibit strong expansion properties or Graph is disconnected")
  }
  
  graph$ExpansionRatio=expansion_ratio
  
  if (add_edge_weight) {
    graph$adj_mat[graph$adj_mat == 1] <- runif(sum(graph$adj_mat == 1))
    graph$Weights=graph$adj_mat
  }
  
  # Create node attributes if requested
  if (add_node_attr) {
    node_attrs <- data.frame(x = points[,1], y = points[,2])
    graph$NodeAttributes=node_attrs
  } 
  return(graph)
}      

g=spatial_scale_free_expander_noncomm(N=50, beta=2, m=2, prewire = 0.1, 
                                      add_edge_weight = FALSE,mu=0, 
                                      add_node_attr = T,r=0.4,alpha=2)

g$degrees

plot(g$GraphObject,vertex.label=NA,vertex.size=2)

lo=layout.norm(as.matrix(cbind(g$NodeAttributes$x,g$NodeAttributes$y)))

plot(g$GraphObject,vertex.label=NA,vertex.size=2,layout=lo)


### Testing
SF=sample_pa(500,power = 2, m=1,directed = F,algorithm = c("psumtree"))

plot(SF,vertex.label=NA,vertex.size=1)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2] Spatial Expander Propagation Graph with weighted edges, 
# node attributes, scale free degree distribution, small world effect and community structures
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++ Code in R to grow a spatial scale-free expander graph with the following properties. 
#+This code should be built from scratch and not contain any igraph packages etc.

# 1) Nodes should be generated on a unit square where the initial set of nodes are distributed spatially with an underlying spatial structure, e.g. by not permitting nodes to land too near each other, or not too far from each other using an in built 'quadtree' algorithm built from scratch
# 2) Use the distance matrix to create an adjacency matrix for the graph.o create an adjacency matrix based on the distance matrix, you can define a cutoff distance 'r' and set the adjacency matrix element A_{ij} to 1 if dist_{ij} <= 'r' and 0 otherwise . This parameter 'r' controls the strength of the spatial distance effect
# 3) This network should be grown and nodes should be organized into 
#non-overlapping communities ( or blocks) with nodes connecting with each 
#other within the same community or between communities with an attachment probability 
#favouring short spatial distances and higher degree nodes that follows scale-free sdegree distribution. 
#Include 'alpha': the parameter that controls the strength of the degree effect
# 4)Add a rewiring probability 'prewire' for small world effect to the resulting graph to enhance the connectivity to long range nodes. Thus you can randomly rewire each edge with probability 'prewire' to a random node within a certain distance
# 5) Add arguments to create edge weights, node attributes etc
# 6)Check the expansion properties of the graph by computing the Laplacian matrix. The Laplacian matrix represents the graph's connectivity and its ability to spread information. A good spatial graph should have a small eigenvalue gap, indicating strong expansion properties.
# 7)This graph should be a single function where we can tun', N, 'alpha' , 'prewire', 'm', to generate different graph models

# count nodes for each clusters
#NumOfNodesPerCluster <- table(graph$clusters)
# num_nodes_per_comm <- rep(1, numofComm)

# for (i in 1:numofComm) {
#   for (j in 1:numofComm) {
#     if (i == j) {
#       P_ij[i, j] <- probWithin
#     } else {
#       P_ij[i, j] <- probBetween
#     }
#   }
# }

spatial_scale_free_expander_comm <- function(N=50, beta=2, m=1, prewire = 0.1, numofComm=4, 
                                             mu=0.2,add_edge_weight = FALSE,
                                             Prob_Betw_comm=0.1,add_node_attr = FALSE,r=0.1,alpha=0.2) {
  
  
  #--Step 1:Generate initial node positions
  #--Note: Profs Quad tree implementation here
  node_pos <- matrix(runif(N*2), ncol = 2)
  
  # Calculate distance matrix
  dist_mat <- as.matrix(dist(node_pos))
  
  # Initialize adjacency matrix
  graph <- list()
  graph$adj_mat <- matrix(0, nrow = N, ncol = N)
  
  # Add first node to network
  graph$adj_mat[1, 1] <- 1
  
  # initialize the edge_list and community membership
  graph$edge_list <- as.data.frame(matrix(ncol=2, nrow=0))
  
  # create initial community for first node. Assin first node to the first community
  graph$clusters <- rep(1, N)
  
  # Step 2: Create initial community/block probability matrix
  P_ij <- matrix(0, nrow = numofComm, ncol = numofComm)
  
  Prob_Within_comm=1-Prob_Betw_comm
  
  if(is.numeric(Prob_Betw_comm)){
    diag(P_ij) <- Prob_Within_comm
    P_ij[lower.tri(P_ij)] <- Prob_Betw_comm
    P_ij[upper.tri(P_ij)] <- t(P_ij)[upper.tri(P_ij)]
    
  }else {
    P_ij <- matrix(runif(numofComm^2), numofComm, numofComm)
  }
  
  # Step 3: Add nodes to the network one by one
  # Grow the network by adding new nodes and edges
  for (i in 3:N) {
    # Step 3i: Create arbitrary number of communities/blocks
    comm_probs<-rep(1/numofComm, numofComm)#sample.int(numofComm, (i-1), replace = TRUE)
    graph$clusters[i] <- sample(numofComm, 1, prob = comm_probs)
    
    # Step 3ii: Define community probability matrix
    P_comm <- P_ij[1:numofComm, 1:numofComm]
    
    ##--Step 3iii: Define probability of attachment function for within and between-community connections
    # P_unique_within <- rep(0, (i-1))
    # P_unique_between <- rep(0, (i-1))
    
    spatial_probs <- ifelse(dist_mat[1:(i-1), i] <= r, 1/(1+(dist_mat[1:(i-1), i])^beta), 0)
    graph$degrees<- rowSums(graph$adj_mat[1:(i-1), ])
    deg_probs <- (graph$degrees^alpha) / sum( graph$degrees^alpha)
    
    # Step 3iv: Add edges within and between communities based on community and attachment probabilty matrix 
    P_noncomm <- (mu*spatial_probs)+((1-mu)*deg_probs) 
    P_noncomm<-P_noncomm/sum(P_noncomm)
    if (numofComm==1){
      P_join <-P_comm[graph$clusters[(i-1)]] * P_noncomm
    }else{
      P_join <-P_comm[graph$clusters[1:(i-1)], graph$clusters[i]] * P_noncomm
    }
    P_join[is.na(P_join)] <- 0
    
    #Connect nodes based on attachment probability
    neighbors_within <- sample(1:(i-1), m, replace = TRUE, prob = P_join)
    graph$adj_mat[i, neighbors_within] <- 1
    graph$adj_mat[neighbors_within, i] <- 1
    graph$degrees[i] <- graph$degrees[i] + 1
    degrees[neighbors_within] <- degrees[neighbors_within] + 1
    graph$degrees[i] <- sum(graph$adj_mat[i, ])
    # Bind to the empty edge_list for within communities
    graph$edge_list <- rbind(graph$edge_list, c(i, neighbors_within))
    
    # Step 4: Small-world rewiring with probability p
    if (runif(1) < prewire) {
      # Select random node to rewire within cutoff distance
      neighbors <- which(graph$adj_mat[i,] == 1)
      if (length(neighbors) > 0) {
        new_neighbor <- sample(neighbors, 1)
        graph$adj_mat[i, new_neighbor] <- 0
        graph$adj_mat[new_neighbor, i] <- 0
        available_nodes <- setdiff(1:N, c(i, neighbors))
        new_neighbor <- sample(available_nodes, 1)
        graph$adj_mat[i, new_neighbor] <- 1
        graph$adj_mat[new_neighbor, i] <- 1
      }
    }
  }
  
  ### Graph object
  GraphModel <- graph.adjacency(as.matrix(graph$adj_mat), mode="undirected")
  GraphModel=simplify(GraphModel,remove.loops = TRUE,remove.multiple = TRUE)
  graph$GraphObject <- GraphModel
  
  # Compute graph Laplacian for the expansion property
  D <- diag(rowSums(graph$adj_mat))
  L <- D - graph$adj_mat
  # Compute eigenvalues and eigenvectors of Laplacian matrix
  eig <- eigen(Matrix(L))
  eig_values <- Re(eig$values)
  eig_vectors <- eig$vectors
  
  # Compute the expansion ratio
  lambda_2 <- eig_values[2]
  lambda_n <- eig_values[nrow(L)]
  expansion_ratio <- lambda_2 / lambda_n
  # Check for strong expansion properties
  #spectral_gap=eigen(Matrix::t(L) %*% L, symmetric=TRUE, only.values=TRUE)$values[n-1]
  if(!is.connected(GraphModel)||expansion_ratio <= 0){  
    warning("Graph may not exhibit strong expansion properties or Graph is disconnected")
  }
  
  graph$ExpansionRatio=expansion_ratio
  
  if (add_edge_weight) {
    graph$adj_mat[graph$adj_mat == 1] <- runif(sum(graph$adj_mat == 1))
    graph$Weights=graph$adj_mat
  }
  
  # Create node attributes if requested
  if (add_node_attr) {
    node_attrs <- data.frame(x = node_pos[,1], y = node_pos[,2])
    graph$NodeAttributes=node_attrs
  } 
  return(graph)
}      

g=spatial_scale_free_expander_comm(N=200, beta=2, m=10, prewire = 0.01, numofComm=2, 
                                   mu=0.2,add_edge_weight = FALSE,
                                   Prob_Betw_comm=0.0001,add_node_attr = T,r=0.1,alpha=0.2) 

lo=layout.norm(as.matrix(cbind(g$NodeAttributes$x,g$NodeAttributes$y)))


plot(g$GraphObject,vertex.label=NA,vertex.size=2)

plot(g$GraphObject,vertex.label=NA,vertex.size=2,layout=lo)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Training a machine learning regression model with random forest
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(tidymodels)

data(ames)

set.seed(4595)
data_split <- initial_split(ames, strata = "Sale_Price", prop = 0.75)

ames_train <- training(data_split)
ames_test  <- testing(data_split)


##########--------------------------MODEL--1--(RANDOM FOREST)-----------#####################
rfmodel<- rand_forest(mode = "regression")
rfmodel

# You can create create regression models with the tidymodels package 'parsnip'
# to predict continuous or numeric quantities.
# The parsnip package provides two different interfaces to fit a model:
#   
# [1] the non-formula interface (fit_xy()).
# [2] the formula interface (fit()), and

##------------[1] the non-formula model (fit_xy())-----------
preds <- c("Longitude", "Latitude", "Lot_Area", "Neighborhood", "Year_Sold")

rf_nonform_fit <- 
  rfmodel %>%
  set_engine("ranger") %>%
  fit_xy(
    x = ames_train[, preds],
    y = log10(ames_train$Sale_Price)
  )

rf_nonform_fit    

# The non-formula interface doesn't 
# do anything to the predictors before passing them to the underlying model function

##--Predicting with non-formula model
test_results_non_form<- 
  ames_test %>%
  select(Sale_Price) %>%
  mutate(Sale_Price = log10(Sale_Price)) %>%
  bind_cols(
    predict(rf_nonform_fit, new_data = ames_test[, preds])
  )
test_results_non_form
# Note that:
# If the model required indicator variables, we would have to create them manually prior to using fit() (perhaps using the recipes package).
# We had to manually log the outcome prior to modeling.

##-------[2] the formula model (fit()) with tunning parameters-----------
rf_form_fit=
  rand_forest(mode = "regression", mtry = 3, trees = 1000) %>%
  set_engine("ranger") %>%
  fit(
    log10(Sale_Price) ~ Longitude + Latitude + Lot_Area + Neighborhood + Year_Sold,
    data = ames_train
  )

rf_form_fit

##--Predicting with formula model
test_results_form <- 
  ames_test %>%
  select(Sale_Price) %>%
  mutate(Sale_Price = log10(Sale_Price)) %>%
  bind_cols(
    predict(rf_form_fit, new_data = ames_test[, preds])
  )
test_results_form

#Note:
# Suppose that we would like to use the randomForest package instead of ranger. To do so, 
#the only part of the syntax that needs to change is the set_engine() argument:

rand_forest(mode = "regression", mtry = 3, trees = 1000) %>%
  set_engine("randomForest") %>%
  fit(
    log10(Sale_Price) ~ Longitude + Latitude + Lot_Area + Neighborhood + Year_Sold,
    data = ames_train
  )

# Now suppose that we want to modify the value of mtry 
# based on the number of predictors in the data. Usually,
# a good default value is floor(sqrt(num_predictors)) but a 
# pure bagging model requires an mtry value equal to the total number of parameters.
# There may be cases where you may not know how many predictors are going to be present 
# when the model will be fit (perhaps due to the generation of indicator variables or 


# When the model it being fit by parsnip, data descriptors are made available. 
# These attempt to let you know what you will have available when the model is fit. 
# When a model object is created (say using rand_forest()), 
# the values of the arguments that you give it are immediately evaluated unless you delay them. 
# To delay the evaluation of any argument, you can used rlang::expr() to make an expression.
# 
# Two relevant data descriptors for our example model are:
#   
# .preds(): the number of predictor variables in the data set that are associated with the predictors prior to dummy variable creation.
# .cols(): the number of predictor columns after dummy variables (or other encodings) are created.
# Since ranger won't create indicator values, .preds() would be appropriate for mtry for a bagging model.
#For example, let's use an expression with the .preds() descriptor to fit a bagging model:

rand_forest(mode = "regression", mtry = .preds(), trees = 1000) %>%
  set_engine("ranger") %>%
  fit(
    log10(Sale_Price) ~ Longitude + Latitude + Lot_Area + Neighborhood + Year_Sold,
    data = ames_train
  )

##########---------------MODEL--2--(Regularised Linear regression)-----------#####################
# A linear model might work for this data set as well. We can use the linear_reg() parsnip model. 
# There are two engines that can perform regularization/penalization, 
# the glmnet and sparklyr package. The glmnet package only implements a non-formula method,
# but parsnip will allow either one to be used

#Note:
# When regularization is used, the predictors should first be centered and scaled before
# being passed to the model. The formula method won't do that automatically so we will need to do this
# ourselves. We'll use the recipes package for these steps.

norm_recipe <- 
  recipe(
    Sale_Price ~ Longitude + Latitude + Lot_Area + Neighborhood + Year_Sold, 
    data = ames_train
  ) %>%
  step_other(Neighborhood) %>% 
  step_dummy(all_nominal()) %>%
  step_center(all_predictors()) %>%
  step_scale(all_predictors()) %>%
  step_log(Sale_Price, base = 10) %>% 
  # estimate the means and standard deviations
  prep(training = ames_train, retain = TRUE)

##--Now let's fit the model using the processed version of the data
glmn_fit <- 
  linear_reg(penalty = 0.001, mixture = 0.5) %>% 
  set_engine("glmnet") %>%
  fit(Sale_Price ~ ., data = bake(norm_recipe, new_data = NULL))

glmn_fit

#Note: If penalty were not specified, all of the lambda values would be computed.
#To get the predictions for this specific value of lambda (aka penalty):

# First, get the processed version of the test set predictors:
test_normalized <- bake(norm_recipe, new_data = ames_test, all_predictors())

# Predict with the fitted model
test_results <- 
  test_results %>%
  rename(`random forest` = .pred) %>%
  bind_cols(
    predict(glmn_fit, new_data = test_normalized) %>%
      rename(glmnet = .pred)
  )
test_results 

#####--Evaluation-Metrics--#####
test_results %>% metrics(truth = Sale_Price, estimate = glmnet) 

###--Predicted --Plots--with--models--
test_results %>% 
  gather(model, prediction, -Sale_Price) %>% 
  ggplot(aes(x = prediction, y = Sale_Price)) + 
  geom_abline(col = "green", lty = 2) + 
  geom_point(alpha = .4) + 
  facet_wrap(~model) + 
  coord_fixed()

##############################################################################################
setwd("C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi-New")
source("SPATIAL-PIPELINE-NEW-MODEL.R")
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [1] Spatial Expander Propagation Graph with weighted edges, 
# node attributes, scale free degree distribution, small world effect and community structures
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#library(quadtree)
# Function to generate a spatial scale-free expander graph
# with the given parameters
#
# Arguments:
#   N: number of nodes
#   r: cutoff distance for adjacency matrix
#   beta: power law parameter for distance-dependent probability
#   m: number of edges added for each new node
#   alpha: parameter controlling strength of degree effect
#   mu: parameter controlling preference when connecting nodes. Thus, preference to care about 
#       spatial distance between nodes when connecting them, vs caring only about vertex degree. 0:care aboutdegree
#1: care about distance
#   prewire: rewiring probability for small world effect
#   node_attrs: list of node attributes to add to graph
#   edge_weights: logical indicating whether to add edge weights
#
# Returns:
#   A graph object of class igraph

spatial_scale_free_expander_noncomm <- function(N=50, beta=2, m=2, prewire = 0.1, 
                                                mu=0.2,add_edge_weight = FALSE, 
                                                add_node_attr = FALSE,r=0.1,alpha=0.2) {
  #Notes:
  ##--(1) use of of N or lambda
  ##--(2) fix L or ignore
  ###--(3) mu > 0 for preferential attachment: mu captures what p does in the doi:10.1088/1367-2630/9/6/190
  ###--(4) use unit square instead of torus
  
  ##--Generate points on a unit square
  points <- matrix(runif(N*2), ncol=2)
  # calculate distance matrix between all pairs of points
  dist_mat <- as.matrix(dist(points))
  # initialize graph with a single node and no edges
  graph <- list()
  graph$adj_mat <- matrix(0, nrow = N, ncol = N)
  graph$adj_mat[1, 1] <- 1
  # initialize the edge_list
  graph$edge_list <- as.data.frame(matrix(ncol=2, nrow=0))
  
  # Grow the network with new nodes added one at a time
  # Grow the network with new nodes added one at a time
  for (i in 3:N) {
    # preferential attachment effect with scale free technique
    graph$degrees= rowSums(graph$adj_mat[1:(i-1), ])
    deg_probs<- (graph$degrees^alpha)/(sum(graph$degrees^alpha))
    
    #--spatial distance effect incorporating short spatial distances (cut-off distance)
    spatial_probs <- ifelse(dist_mat[1:(i-1), i] <= r, 1/(1+(dist_mat[1:(i-1), i])^beta), 0)#ifelse(dist_mat[i,] <= r, 1/(1+(dist_mat[i,])^beta), 0)
    
    # total probability of attachment 
    totalprobs=(mu*spatial_probs)+((1-mu)*deg_probs) 
    norm_totalprobs=totalprobs/sum(totalprobs)
    # normalize total probability of attachment
    norm_totalprobs[which(norm_totalprobs == Inf)] <- 0 # remove infinite probabilities
    norm_totalprobs[which(is.na(norm_totalprobs))] <- 0
    
    if (mu!=1){
      #Add m edges to existing nodes with probability proportional to their degree and distance
      node_to_attach <- sample(1:(i-1), size = m,replace = TRUE, prob = norm_totalprobs)
      graph$adj_mat[i, node_to_attach] <- 1
      graph$adj_mat[node_to_attach, i] <- 1
    }else{
      node_to_attach<- sample(1:(i-1), size = m,replace = TRUE, prob = NULL)
      graph$adj_mat[i, node_to_attach] <- ifelse(spatial_probs[(i-1)]<= r, 1, 0)
      graph$adj_mat[node_to_attach, i] <- graph$adj_mat[i, node_to_attach]
      
    }
    #Bind to the empty edge_list
    graph$edge_list <- rbind(graph$edge_list, c(i, node_to_attach))
    
    # Small-world rewiring with probability p
    if (runif(1) < prewire) {
      # Select random node to rewire within cutoff distance
      neighbors <- which(graph$adj_mat[i,] == 1)
      if (length(neighbors) > 0) {
        new_neighbor <- sample(neighbors, 1)
        graph$adj_mat[i, new_neighbor] <- 0
        graph$adj_mat[new_neighbor, i] <- 0
        available_nodes <- setdiff(1:N, c(i, neighbors))
        new_neighbor <- sample(available_nodes, 1)
        graph$adj_mat[i, new_neighbor] <- 1
        graph$adj_mat[new_neighbor, i] <- 1
      }
    }
  } 
  
  ### Graph object
  GraphModel <- graph.adjacency(as.matrix(graph$adj_mat), mode="undirected")
  GraphModel=igraph::simplify(GraphModel,remove.loops = TRUE,remove.multiple = TRUE)
  graph$GraphObject <- GraphModel
  
  # Compute graph Laplacian for the expansion property
  D <- diag(rowSums(graph$adj_mat))
  L <- D - graph$adj_mat
  # Compute eigenvalues and eigenvectors of Laplacian matrix
  eig <- eigen(Matrix(L))
  eig_values <- Re(eig$values)
  eig_vectors <- eig$vectors
  
  # Compute the expansion ratio
  lambda_2 <- eig_values[2]
  lambda_n <- eig_values[nrow(L)]
  expansion_ratio <- lambda_2 / lambda_n
  # Check for strong expansion properties
  #spectral_gap=eigen(Matrix::t(L) %*% L, symmetric=TRUE, only.values=TRUE)$values[n-1]
  if(!is.connected(GraphModel)||expansion_ratio <= 0){  
    warning("Graph may not exhibit strong expansion properties or Graph is disconnected")
  }
  
  graph$ExpansionRatio=expansion_ratio
  
  if (add_edge_weight) {
    graph$adj_mat[graph$adj_mat == 1] <- runif(sum(graph$adj_mat == 1))
    graph$Weights=graph$adj_mat
  }
  
  # Create node attributes if requested
  if (add_node_attr) {
    node_attrs <- data.frame(x = points[,1], y = points[,2])
    graph$NodeAttributes=node_attrs
  } 
  return(graph)
}      

g=spatial_scale_free_expander_noncomm(N=50, beta=0.05, m=4, prewire = 0.3, 
                                      add_edge_weight = FALSE,mu=0.4, 
                                      add_node_attr = T,r=Inf,alpha=6)


plot(g$GraphObject,vertex.label=NA,vertex.size=2)

lo=layout.norm(as.matrix(cbind(g$NodeAttributes$x,g$NodeAttributes$y)))

plot(g$GraphObject,vertex.label=NA,vertex.size=2,layout=lo)


### Testing
SF=sample_pa(500,power = 2, m=1,directed = F,algorithm = c("psumtree"))

plot(SF,vertex.label=NA,vertex.size=1)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2] Spatial Expander Propagation Graph with weighted edges, 
# node attributes, scale free degree distribution, small world effect and community structures
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++ Code in R to grow a spatial scale-free expander graph with the following properties. 
#+This code should be built from scratch and not contain any igraph packages etc.

# 1) Nodes should be generated on a unit square where the initial set of nodes are distributed spatially with an underlying spatial structure, e.g. by not permitting nodes to land too near each other, or not too far from each other using an in built 'quadtree' algorithm built from scratch
# 2) Use the distance matrix to create an adjacency matrix for the graph.o create an adjacency matrix based on the distance matrix, you can define a cutoff distance 'r' and set the adjacency matrix element A_{ij} to 1 if dist_{ij} <= 'r' and 0 otherwise . This parameter 'r' controls the strength of the spatial distance effect
# 3) This network should be grown and nodes should be organized into 
#non-overlapping communities ( or blocks) with nodes connecting with each 
#other within the same community or between communities with an attachment probability 
#favouring short spatial distances and higher degree nodes that follows scale-free sdegree distribution. 
#Include 'alpha': the parameter that controls the strength of the degree effect
# 4)Add a rewiring probability 'prewire' for small world effect to the resulting graph to enhance the connectivity to long range nodes. Thus you can randomly rewire each edge with probability 'prewire' to a random node within a certain distance
# 5) Add arguments to create edge weights, node attributes etc
# 6)Check the expansion properties of the graph by computing the Laplacian matrix. The Laplacian matrix represents the graph's connectivity and its ability to spread information. A good spatial graph should have a small eigenvalue gap, indicating strong expansion properties.
# 7)This graph should be a single function where we can tun', N, 'alpha' , 'prewire', 'm', to generate different graph models

# count nodes for each clusters
#NumOfNodesPerCluster <- table(graph$clusters)
# num_nodes_per_comm <- rep(1, numofComm)

# for (i in 1:numofComm) {
#   for (j in 1:numofComm) {
#     if (i == j) {
#       P_ij[i, j] <- probWithin
#     } else {
#       P_ij[i, j] <- probBetween
#     }
#   }
# }

spatial_scale_free_expander_comm <- function(N=50, beta=2, m=1, prewire = 0.1, numofComm=4, 
                                             mu=0.2,add_edge_weight = FALSE,
                                             Prob_Betw_comm=0.1,add_node_attr = FALSE,r=0.1,alpha=0.2) {
  
  
  #--Step 1:Generate initial node positions
  #--Note: Profs Quad tree implementation here
  node_pos <- matrix(runif(N*2), ncol = 2)
  
  # Calculate distance matrix
  dist_mat <- as.matrix(dist(node_pos))
  
  # Initialize adjacency matrix
  graph <- list()
  graph$adj_mat <- matrix(0, nrow = N, ncol = N)
  
  # Add first node to network
  graph$adj_mat[1, 1] <- 1
  
  # initialize the edge_list and community membership
  graph$edge_list <- as.data.frame(matrix(ncol=2, nrow=0))
  
  # create initial community for first node. Assin first node to the first community
  graph$clusters <- rep(1, N)
  
  # Step 2: Create initial community/block probability matrix
  P_ij <- matrix(0, nrow = numofComm, ncol = numofComm)
  
  Prob_Within_comm=1-Prob_Betw_comm
  
  if(is.numeric(Prob_Betw_comm)){
    diag(P_ij) <- Prob_Within_comm
    P_ij[lower.tri(P_ij)] <- Prob_Betw_comm
    P_ij[upper.tri(P_ij)] <- t(P_ij)[upper.tri(P_ij)]
    
  }else {
    P_ij <- matrix(runif(numofComm^2), numofComm, numofComm)
  }
  
  # Step 3: Add nodes to the network one by one
  # Grow the network by adding new nodes and edges
  for (i in 3:N) {
    # Step 3i: Create arbitrary number of communities/blocks
    comm_probs<-rep(1/numofComm, numofComm)#sample.int(numofComm, (i-1), replace = TRUE)
    graph$clusters[i] <- sample(numofComm, 1, prob = comm_probs)
    
    # Step 3ii: Define community probability matrix
    P_comm <- P_ij[1:numofComm, 1:numofComm]
    
    ##--Step 3iii: Define probability of attachment function for within and between-community connections
    # P_unique_within <- rep(0, (i-1))
    # P_unique_between <- rep(0, (i-1))
    
    spatial_probs <- ifelse(dist_mat[1:(i-1), i] <= r, 1/(1+(dist_mat[1:(i-1), i])^beta), 0)
    graph$degrees<- rowSums(graph$adj_mat[1:(i-1), ])
    deg_probs <- (graph$degrees^alpha) / sum( graph$degrees^alpha)
    
    # Step 3iv: Add edges within and between communities based on community and attachment probabilty matrix 
    P_noncomm <- (mu*spatial_probs)+((1-mu)*deg_probs) 
    P_noncomm<-P_noncomm/sum(P_noncomm)
    if (numofComm==1){
      P_join <-P_comm[graph$clusters[(i-1)]] * P_noncomm
      P_join[is.na(P_join)] <- 0
      if (mu!=1){
        #Add m edges to existing nodes with probability proportional to their degree and distance
        node_to_attach <- sample(1:(i-1), size = m,replace = TRUE, prob = P_join)
        graph$adj_mat[i, node_to_attach] <- 1
        graph$adj_mat[node_to_attach, i] <- 1
        graph$degrees[i] <- graph$degrees[i] + 1
        graph$degrees[node_to_attach] <- graph$degrees[node_to_attach] + 1
        graph$edge_list <- rbind(graph$edge_list, c(i, node_to_attach))
      }else{
        node_to_attach<- sample(1:(i-1), size = m,replace = TRUE, prob = NULL)
        graph$adj_mat[i, node_to_attach] <- ifelse(spatial_probs[(i-1)]<= r, 1, 0)
        graph$adj_mat[node_to_attach, i] <- graph$adj_mat[i, node_to_attach]
        graph$degrees[i] <- graph$degrees[i] + 1
        graph$degrees[node_to_attach] <- graph$degrees[node_to_attach] + 1
        graph$edge_list <- rbind(graph$edge_list, c(i, node_to_attach))
      }
    }else{
      P_join <-P_comm[graph$clusters[1:(i-1)], graph$clusters[i]] * P_noncomm
      P_join[is.na(P_join)] <- 0
      node_to_attach<- sample(1:(i-1), m, replace = TRUE, prob = P_join)
      graph$adj_mat[i, node_to_attach] <- 1
      graph$adj_mat[node_to_attach, i] <- 1
      graph$degrees[i] <- graph$degrees[i] + 1
      graph$degrees[node_to_attach] <- graph$degrees[node_to_attach] + 1
      graph$edge_list <- rbind(graph$edge_list, c(i, node_to_attach))
    }
    #P_join[is.na(P_join)] <- 0
    
    #Connect nodes based on attachment probability
    #node_to_attach<- sample(1:(i-1), m, replace = TRUE, prob = P_join)
    #graph$adj_mat[i, node_to_attach] <- 1
    #graph$adj_mat[node_to_attach, i] <- 1
    #graph$degrees[i] <- graph$degrees[i] + 1
    #degrees[node_to_attach] <- degrees[node_to_attach] + 1
    #graph$degrees[i] <- sum(graph$adj_mat[i, ])
    # Bind to the empty edge_list for within communities
    #graph$edge_list <- rbind(graph$edge_list, c(i, node_to_attach))
    
    # Step 4: Small-world rewiring with probability p
    if (runif(1) < prewire) {
      # Select random node to rewire within cutoff distance
      neighbors <- which(graph$adj_mat[i,] == 1)
      if (length(neighbors) > 0) {
        new_neighbor <- sample(neighbors, 1)
        graph$adj_mat[i, new_neighbor] <- 0
        graph$adj_mat[new_neighbor, i] <- 0
        available_nodes <- setdiff(1:N, c(i, neighbors))
        new_neighbor <- sample(available_nodes, 1)
        graph$adj_mat[i, new_neighbor] <- 1
        graph$adj_mat[new_neighbor, i] <- 1
      }
    }
  }
  
  ### Graph object
  GraphModel <- graph.adjacency(as.matrix(graph$adj_mat), mode="undirected")
  GraphModel=simplify(GraphModel,remove.loops = TRUE,remove.multiple = TRUE)
  graph$GraphObject <- GraphModel
  
  # Compute graph Laplacian for the expansion property
  D <- diag(rowSums(graph$adj_mat))
  L <- D - graph$adj_mat
  # Compute eigenvalues and eigenvectors of Laplacian matrix
  eig <- eigen(Matrix(L))
  eig_values <- Re(eig$values)
  eig_vectors <- eig$vectors
  
  # Compute the expansion ratio
  lambda_2 <- eig_values[2]
  lambda_n <- eig_values[nrow(L)]
  expansion_ratio <- lambda_2 / lambda_n
  # Check for strong expansion properties
  #spectral_gap=eigen(Matrix::t(L) %*% L, symmetric=TRUE, only.values=TRUE)$values[n-1]
  if(!is.connected(GraphModel)||expansion_ratio <= 0){  
    warning("Graph may not exhibit strong expansion properties or Graph is disconnected")
  }
  
  graph$ExpansionRatio=expansion_ratio
  
  if (add_edge_weight) {
    graph$adj_mat[graph$adj_mat == 1] <- runif(sum(graph$adj_mat == 1))
    graph$Weights=graph$adj_mat
  }
  
  # Create node attributes if requested
  if (add_node_attr) {
    node_attrs <- data.frame(x = node_pos[,1], y = node_pos[,2])
    graph$NodeAttributes=node_attrs
  } 
  return(graph)
}      

g=spatial_scale_free_expander_comm(N=200, beta=2, m=10, prewire = 0.01, numofComm=1, 
                                   mu=0.2,add_edge_weight = FALSE,
                                   Prob_Betw_comm=0.0001,add_node_attr = T,r=0.1,alpha=0.2) 

lo=layout.norm(as.matrix(cbind(g$NodeAttributes$x,g$NodeAttributes$y)))


plot(g$GraphObject,vertex.label=NA,vertex.size=2)

plot(g$GraphObject,vertex.label=NA,vertex.size=2,layout=lo)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Training a machine learning regression model with random forest
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(tidymodels)

data(ames)

set.seed(4595)
data_split <- initial_split(ames, strata = "Sale_Price", prop = 0.75)

ames_train <- training(data_split)
ames_test  <- testing(data_split)


##########--------------------------MODEL--1--(RANDOM FOREST)-----------#####################
rfmodel<- rand_forest(mode = "regression")
rfmodel

# You can create create regression models with the tidymodels package 'parsnip'
# to predict continuous or numeric quantities.
# The parsnip package provides two different interfaces to fit a model:
#   
# [1] the non-formula interface (fit_xy()).
# [2] the formula interface (fit()), and

##------------[1] the non-formula model (fit_xy())-----------
preds <- c("Longitude", "Latitude", "Lot_Area", "Neighborhood", "Year_Sold")

rf_nonform_fit <- 
  rfmodel %>%
  set_engine("ranger") %>%
  fit_xy(
    x = ames_train[, preds],
    y = log10(ames_train$Sale_Price)
  )

rf_nonform_fit    

# The non-formula interface doesn't 
# do anything to the predictors before passing them to the underlying model function

##--Predicting with non-formula model
test_results_non_form<- 
  ames_test %>%
  select(Sale_Price) %>%
  mutate(Sale_Price = log10(Sale_Price)) %>%
  bind_cols(
    predict(rf_nonform_fit, new_data = ames_test[, preds])
  )
test_results_non_form
# Note that:
# If the model required indicator variables, we would have to create them manually prior to using fit() (perhaps using the recipes package).
# We had to manually log the outcome prior to modeling.

##-------[2] the formula model (fit()) with tunning parameters-----------
rf_form_fit=
  rand_forest(mode = "regression", mtry = 3, trees = 1000) %>%
  set_engine("ranger") %>%
  fit(
    log10(Sale_Price) ~ Longitude + Latitude + Lot_Area + Neighborhood + Year_Sold,
    data = ames_train
  )

rf_form_fit

##--Predicting with formula model
test_results_form <- 
  ames_test %>%
  select(Sale_Price) %>%
  mutate(Sale_Price = log10(Sale_Price)) %>%
  bind_cols(
    predict(rf_form_fit, new_data = ames_test[, preds])
  )
test_results_form

#Note:
# Suppose that we would like to use the randomForest package instead of ranger. To do so, 
#the only part of the syntax that needs to change is the set_engine() argument:

rand_forest(mode = "regression", mtry = 3, trees = 1000) %>%
  set_engine("randomForest") %>%
  fit(
    log10(Sale_Price) ~ Longitude + Latitude + Lot_Area + Neighborhood + Year_Sold,
    data = ames_train
  )

# Now suppose that we want to modify the value of mtry 
# based on the number of predictors in the data. Usually,
# a good default value is floor(sqrt(num_predictors)) but a 
# pure bagging model requires an mtry value equal to the total number of parameters.
# There may be cases where you may not know how many predictors are going to be present 
# when the model will be fit (perhaps due to the generation of indicator variables or 


# When the model it being fit by parsnip, data descriptors are made available. 
# These attempt to let you know what you will have available when the model is fit. 
# When a model object is created (say using rand_forest()), 
# the values of the arguments that you give it are immediately evaluated unless you delay them. 
# To delay the evaluation of any argument, you can used rlang::expr() to make an expression.
# 
# Two relevant data descriptors for our example model are:
#   
# .preds(): the number of predictor variables in the data set that are associated with the predictors prior to dummy variable creation.
# .cols(): the number of predictor columns after dummy variables (or other encodings) are created.
# Since ranger won't create indicator values, .preds() would be appropriate for mtry for a bagging model.
#For example, let's use an expression with the .preds() descriptor to fit a bagging model:

rand_forest(mode = "regression", mtry = .preds(), trees = 1000) %>%
  set_engine("ranger") %>%
  fit(
    log10(Sale_Price) ~ Longitude + Latitude + Lot_Area + Neighborhood + Year_Sold,
    data = ames_train
  )

##########---------------MODEL--2--(Regularised Linear regression)-----------#####################
# A linear model might work for this data set as well. We can use the linear_reg() parsnip model. 
# There are two engines that can perform regularization/penalization, 
# the glmnet and sparklyr package. The glmnet package only implements a non-formula method,
# but parsnip will allow either one to be used

#Note:
# When regularization is used, the predictors should first be centered and scaled before
# being passed to the model. The formula method won't do that automatically so we will need to do this
# ourselves. We'll use the recipes package for these steps.

norm_recipe <- 
  recipe(
    Sale_Price ~ Longitude + Latitude + Lot_Area + Neighborhood + Year_Sold, 
    data = ames_train
  ) %>%
  step_other(Neighborhood) %>% 
  step_dummy(all_nominal()) %>%
  step_center(all_predictors()) %>%
  step_scale(all_predictors()) %>%
  step_log(Sale_Price, base = 10) %>% 
  # estimate the means and standard deviations
  prep(training = ames_train, retain = TRUE)

##--Now let's fit the model using the processed version of the data
glmn_fit <- 
  linear_reg(penalty = 0.001, mixture = 0.5) %>% 
  set_engine("glmnet") %>%
  fit(Sale_Price ~ ., data = bake(norm_recipe, new_data = NULL))

glmn_fit

#Note: If penalty were not specified, all of the lambda values would be computed.
#To get the predictions for this specific value of lambda (aka penalty):

# First, get the processed version of the test set predictors:
test_normalized <- bake(norm_recipe, new_data = ames_test, all_predictors())

# Predict with the fitted model
test_results <- 
  test_results %>%
  rename(`random forest` = .pred) %>%
  bind_cols(
    predict(glmn_fit, new_data = test_normalized) %>%
      rename(glmnet = .pred)
  )
test_results 

#####--Evaluation-Metrics--#####
test_results %>% metrics(truth = Sale_Price, estimate = glmnet) 

###--Predicted --Plots--with--models--
test_results %>% 
  gather(model, prediction, -Sale_Price) %>% 
  ggplot(aes(x = prediction, y = Sale_Price)) + 
  geom_abline(col = "green", lty = 2) + 
  geom_point(alpha = .4) + 
  facet_wrap(~model) + 
  coord_fixed()
