#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+[1] Multi class or Multinomial classification
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
setwd("C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi-New")
source("SPATIAL-PIPELINE-NEW.R")
library(tidymodels)
df=read.csv("TrainingGraphFeatures.csv",sep = ",", header = T)
df=as_tibble(df)
df%>%count(GraphName,sort=T)
df%>%View()

Data= df%>%dplyr::select(c(GraphName,order,edges,mean_degree,minCut,
                           FiedlerValue,Normalized_FiedlerValue,closeness,modularity,diameter,
                           betweenness,transitivity,spectral_radius,centrality_eigen))%>%
  mutate_if(is.character,factor)


##----Data Prepocessing---#######
#Create bootstrap sample and measure training and testing performance on subset of each
df.boot=bootstraps(Data)

#If You have class imbalance say value of one class is lower or very high, use "step_smote"
#to generate new sample from the minority class using nearest neighbour technique
# library(themis). eg

# df_recipe=recipe(GraphName~., data=Data)%>%
#   update_role(order,new_role = "numberOfNodes")%>%
#   step_other(....)%>%
#   step_dummy(....)%>%
#   step_zv(all_predictors())%>%
#   step_normalize(all_predictors())%>%
#   step_smote(GraphName)


df_recipe=recipe(GraphName~., data=Data)%>%
  update_role(order,new_role = "numberOfNodes")%>%
   step_zv(all_predictors())%>%
   step_normalize(all_predictors())

df_prep=prep(df_recipe)
juice(df_prep)


###---Model-----####
rfModel=rand_forest(trees=1000)%>%
  set_mode("classification") %>%
  set_engine("ranger")

Graph_wf=workflow()%>%
add_recipe(df_recipe)%>%
  add_model(rfModel)

##--Fitting--Model--
Graph_res=fit_resamples(
  Graph_wf,
  resamples = df.boot,
  control=control_resamples(save_pred = T,verbose = T)
)
  
##--Result Exploration--##
Graph_res%>%
  collect_metrics()


Graph_res%>%
  collect_predictions()%>%
  conf_mat(GraphName,.pred_class)

Graph_res%>%
  collect_predictions()%>%
  ppv(GraphName,.pred_class)

Graph_res%>%
  collect_predictions()%>%View()

Graph_res%>%
  collect_predictions()%>%
  mutate(correct=GraphName==.pred_class)%>%count(correct)

Data%>%mutate(.row=row_number())

Graph_pred=Graph_res%>%
  collect_predictions()%>%
  mutate(correct=GraphName==.pred_class)%>%
  left_join(Data%>%mutate(.row=row_number()))

# Graph_res %>%
#   group_by(edges)%>%
#   collect_predictions()%>%
#   ppv(GraphName,.pred_class)

##---Variable--Importance--
library(vip)

rfModel %>%
  set_engine("ranger", importance = "permutation") %>%
  fit(
    GraphName ~.,
    data=Data,
  ) %>%
  vip(geom = "point")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+[2] Multi class or Multinomial classification
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(visdat)

#Visualizing all data structure
vis_dat(Data)

#Checking levels of factor or character nd creating a table
library(gt)
Data %>% 
  count(GraphName,
        sort = TRUE)%>%
  gt()


##--check for missing data
vis_miss(Data, sort_miss = TRUE)


##--Plotting a histogram/pairs plot of each numeric features
library(GGally)

Data1=Data%>%select(c(order,edges,mean_degree,minCut,
         FiedlerValue,Normalized_FiedlerValue,closeness,modularity,diameter,
         betweenness,transitivity,spectral_radius,centrality_eigen))%>%
  ggscatmat(alpha = 0.2)

Data1

# Data1=Data%>%select(c(order,edges,mean_degree,minCut,
#                       FiedlerValue,Normalized_FiedlerValue,closeness,modularity,diameter,
#                       betweenness,transitivity,spectral_radius,centrality_eigen))%>%
#   ggpairs()
# 
# Data1

##--Data--Splitting------------------##
Data %>% 
  ggplot(aes(GraphName)) +
  geom_bar()

# To actually split the data, we can use the rsample package (included in tidymodels)
# to create an object that contains the information on how to split the data 
# (which we call data_split), and
#then two more rsample functions to create data frames for the training and testing sets:

# This enables the analysis to be reproducible 
set.seed(123)

# Put 3/4 of the data into the training set 
df.split <- initial_split(Data, 
                            prop = 3/4, 
                            strata = GraphName)

# Create dataframes for the two sets:
train_data <- training(df.split) 
test_data <- testing(df.split)

###---Numerical variables
#We can use boxplots to check, if we actually find differences in our numeric variables for the different levels of our dependent categorical variable:
  
  Data %>% 
  ggplot(aes(x = GraphName, y = FiedlerValue, 
             fill = GraphName, color = GraphName)) +
  geom_boxplot(alpha=0.4) 
  
  

  
  
print_boxplot <- function(.y_var){
#---convert strings to variable
y_var <- dplyr::sym(.y_var) 
    
#--unquote variables using {{}}
    Data %>% 
      ggplot(aes(x = GraphName, y = {{y_var}},
                 fill = GraphName, color = GraphName)) +
      geom_boxplot(alpha=0.4) 
    
  }

#--Obtain all of the names of the y-variables we want to use for our plots:
y_var <- Data %>% 
  select(where(is.numeric), -GraphName) %>% 
  variable.names() # obtain name    

#--The map function applys the function print_boxplot to each element of our atomic vector y_var and returns the according plot:
library(purrr)
map(y_var, print_boxplot)

#--Data--preparation
glimpse(Data)




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Data prepropecessing recipe
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# The type of data preprocessing is dependent on the data and the type of model being fit. The excellent book "Tidy Modeling with R" provides an appendix with recommendations for baseline levels of preprocessing that are needed for various model functions.
# Let's create a base recipe for all of our classification models. Note that the sequence of steps matter:

# The recipe() function has two arguments:
#   
#A formula. Any variable on the left-hand side of the tilde (~) is considered the model outcome (here, price_category). On the right-hand side of the tilde are the predictors. Variables may be listed by name (separated by a +), or you can use the dot (.) to indicate all other variables as predictors.
# The data. A recipe is associated with the data set used to create the model. This will typically be the training set, so data = train_data here.
# update_role(): This step of adding roles to a recipe is optional; the purpose of using it here is that those two variables can be retained in the data but not included in the model. This can be convenient when, after the model is fit, we want to investigate some poorly predicted value. These ID columns will be available and can be used to try to understand what went wrong.
# 
# step_naomit() removes observations (rows of data) if they contain NA or NaN values. We use skip = TRUE because we don't want to perform this part to new data so that the number of samples in the assessment set is the same as the number of predicted values (even if they are NA).
# 
# Note that instead of deleting missing values we could also easily substitute (i.e., impute) missing values of variables by one of the following methods (using the training set):
#   
# median,
# mean,
# mode,
# k-nearest neighbors,
# linear model,
# bagged tree models

# 
# Take a look at the recipes reference for an overview about all possible imputation methods.
# step_novel() converts all nominal variables to factors and takes care of other issues related to categorical variables.
# step_log() will log transform data (since some of our numerical variables are right-skewed). Note that this step can not be performed on negative numbers.
# step_normalize() normalizes (center and scales) the numeric variables to have a standard deviation of one and a mean of zero. (i.e., z-standardization).
# step_dummy() converts our factor column ocean_proximity into numeric binary (0 and 1) variables.
# Note that this step may cause problems if your categorical variable has too many levels - especially if some of the levels are very infrequent. In this case you should either drop the variable or pool infrequently occurring values into an "other" category with step_other. This steps has to be performed befor step_dummy.
# step_zv(): removes any numeric variables that have zero variance.
# step_corr(): will remove predictor variables that have large correlations with other predictor variables.
# Note that the package themis contains extra steps for the recipes package for dealing with imbalanced data. A classification data set with skewed class proportions is called imbalanced. Classes that make up a large proportion of the data set are called majority classes. Those that make up a smaller proportion are minority classes (see Google Developers for more details). Themis provides various methods for over-sampling (e.g. SMOTE) and under-sampling. However, we don't have to use this methods since our data is not imbalanced.

df_recipe=recipe(GraphName~., data=train_data)%>%
  update_role(order,new_role = "numberOfNodes")%>%
  step_zv(all_predictors())%>%
  step_normalize(all_predictors())%>%
 step_corr(all_predictors(), threshold = 0.7, method = "spearman")

summary(df_recipe)

#If we would like to check if all of our preprocessing steps from above actually worked,
#we can proceed as follows:
  
  df_prep <- 
    df_recipe %>% # use the recipe object
  prep() %>% # perform the recipe on training data
  juice() # extract only the preprocessed dataframe 

  #Take a look at the data structure:
  glimpse(df_prep)
  
  
#--Validation--set
  # Remember that we already partitioned our data set into a training set and test set. 
  # This lets us judge whether a given model will generalize well to new data. 
  # However, using only two partitions may be insufficient when doing many rounds of
  # hyperparameter tuning (which we don't perform in this tutorial but it is always 
  #                        recommended to use a validation set).
  # 
  # Therefore, it is usually a good idea to create a so called validation set. 
  # Watch this short video from Google's Machine Learning crash course to learn more about the value of a validation set.
  # 
  # We use k-fold crossvalidation to build a set of 5 validation folds with the function 
  # vfold_cv. We also use stratified sampling:

  set.seed(100)
  cv_folds <-
    vfold_cv(train_data, 
             v = 10, 
             strata = GraphName)
  
##-----Model building----
  #Specify models
  # The process of specifying our models is always as follows:
  ###Pick a model type
  ###set the engine
  ###Set the mode: regression or classification


##--[A]--Random--Forest------
library(ranger)
  
  rf.model <- 
    rand_forest() %>% 
    set_engine("ranger", importance = "impurity") %>% 
    set_mode("classification")

##--[B]--Boosted--tree--(XGBoost)
  library(xgboost)
  xgb.model <- 
    boost_tree() %>% 
    set_engine("xgboost") %>% 
    set_mode("classification") 

##--[C]--Neural--network
  # library(keras)
  # nnet.model <-
  #   mlp() %>%
  #   set_mode("classification") %>% 
  #   set_engine("keras", verbose = 0) 

  
##---Create--workflows
  # To combine the data preparation recipe with the model building, 
  # we use the package workflows. 
  #A workflow is an object that can bundle together your pre-processing recipe,
  #modeling, and even post-processing requests (like calculating the RMSE).   
  
  rf.wkflow <-
    workflow() %>%
    add_recipe(df_recipe) %>% 
    add_model(rf.model) 
 
    
  xgb.wkflow <-
    workflow() %>%
    add_recipe(df_recipe) %>% 
    add_model(xgb.model)
  
  # nnet.wkflow <-
  #   workflow() %>%
  #   add_recipe(df_recipe) %>% 
  #   add_model(nnet.model)
  # 

##--Evaluate models
  # Now we can use our validation 
  # set (cv_folds) to estimate the performance of our models using the fit_resamples() 
  # function to fit the models on each of the folds and store the results.
  # 
  # Note that fit_resamples() will fit our model to each resample and 
  # evaluate on the heldout set from each resample. 
  # The function is usually only used for computing performance metrics across some 
  # set of resamples to evaluate our models (like accuracy) - the models are
  # not even stored. However, in our example we save the predictions in order to
  # visualize the model fit and residuals with control_resamples(save_pred = TRUE).
  #Finally, we collect the performance metrics with collect_metrics() and 
  #pick the model that does best on the validation set.    

  
  rf.res <-
    rf.wkflow %>% 
    fit_resamples(
      resamples = cv_folds, 
      metrics = metric_set(
        recall, precision, f_meas, 
        accuracy, kap,
        roc_auc, sens, spec),
      control = control_resamples(save_pred = TRUE)
    ) 
  
  rf.res %>%  collect_metrics(summarize = TRUE)
  
  
  xgb.res <- 
    xgb.wkflow %>% 
    fit_resamples(
      resamples = cv_folds, 
      metrics = metric_set(
        recall, precision, f_meas, 
        accuracy, kap,
        roc_auc, sens, spec),
      control = control_resamples(save_pred = TRUE)
    ) 
  
  xgb.res %>% collect_metrics(summarize = TRUE)
  
  
  
  # nnet.res <- 
  #   nnet.wkflow %>% 
  #   fit_resamples(
  #     resamples = cv_folds, 
  #     metrics = metric_set(
  #       recall, precision, f_meas, 
  #       accuracy, kap,
  #       roc_auc, sens, spec),
  #     control = control_resamples(save_pred = TRUE)
  #   ) 
  

###-----Compare----models-----
# Extract metrics from our models to compare them:
    
  rf.metrics <- 
    rf.res %>% 
    collect_metrics(summarise = TRUE) %>%
    mutate(model = "Random Forest")
  
  xgb.metrics <- 
    xgb.res %>% 
    collect_metrics(summarise = TRUE) %>%
    mutate(model = "XGBoost")
  
   # nnet.metrics <- 
   #   nnet_res %>% 
   #   collect_metrics(summarise = TRUE) %>%
   #   mutate(model = "Neural Net")
   # 
##--create dataframe with all models
  model.compare <- bind_rows(
    #log_metrics,
    rf.metrics,
    xgb.metrics
    #knn_metrics,
    # nnet_metrics
  ) 
  
##--change data structure
  model.comp <- 
    model.compare %>% 
    select(model, .metric, mean, std_err) %>% 
    pivot_wider(names_from = .metric, values_from = c(mean, std_err)) 
  
##--show mean F1-Score for every model
  library(forcats)
  model.comp %>% 
    arrange(mean_f_meas) %>% 
    mutate(model = fct_reorder(model, mean_f_meas)) %>% # order results
    ggplot(aes(model, mean_f_meas, fill=model)) +
    geom_col() +
    coord_flip() +
    scale_fill_brewer(palette = "Blues") +
    geom_text(
      size = 3,
      aes(label = round(mean_f_meas, 2), y = mean_f_meas + 0.08),
      vjust = 1
    )    
  
##--show mean area under the curve (auc) per model
  model.comp %>% 
    arrange(mean_roc_auc) %>% 
    mutate(model = fct_reorder(model, mean_roc_auc)) %>%
    ggplot(aes(model, mean_roc_auc, fill=model)) +
    geom_col() +
    coord_flip() +
    scale_fill_brewer(palette = "Blues") + 
    geom_text(
      size = 3,
      aes(label = round(mean_roc_auc, 2), y = mean_roc_auc + 0.08),
      vjust = 1
    )
  
  

model.comp %>% slice_max(mean_f_meas)  

##---Last evaluation on test set
#Tidymodels provides the function last_fit() which fits a model to the whole training data and evaluates it on the test set. We just need to provide the workflow object of the best model as well as the data split object (not the training data).

best.model.fit2 <- last_fit(rf.wkflow,
                        split = df.split,
                        metrics = metric_set(
                          recall, precision, f_meas,
                          accuracy, kap,
                          roc_auc, sens, spec)
)

best.model.fit <- rf.wkflow %>%
                    fit(data=train_data)   
                       

#--Show performance metrics
best.model.fit2 %>% 
  collect_metrics()

 # best.model.fit %>% 
 #   extract_fit_parsnip()%>%
 # tidy()


##--Variabe--Imortance
library(vip)
best.model.fit2%>% 
  pluck(".workflow", 1) %>%   
  pull_workflow_fit() %>% 
  vip(num_features = 12)

##--confusion matrix:
best.model.fit2 %>%
  collect_predictions() %>% 
  conf_mat(GraphName, .pred_class) %>% 
  autoplot(type = "heatmap")


#pred=predict(best.model.fit,....)
 # best.model.fit %>% 
 #    collect_predictions() %>% 
 #    roc_curve(GraphName, .pred_ER) %>% 
 #    autoplot()
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # # tab_header(
#   title = "California median house prices",
#   subtitle = "Districts above and below 150.000$"
# ) %>%
#   cols_label(
#     price_category = "Price",
#     districts_total = "Districts",
#     percent = "Percent

##--convert to numeric if not already done
# Data  <- 
#   Data %>% 
#   mutate(
#     housing_median_age = as.numeric(housing_median_age),
#     median_house_value = as.numeric(median_house_value)
#   )

##--convert all remaining character variables to factors 
# Data  <- 
#   Data  %>% 
#   mutate(across(where(is.character), as.factor))


# Data %>% 
#   count(GraphName, # count observations
#         instances ="total")

#--Create new variables
# One very important thing you may want to do at the beginning of your data science project is to create new variable combinations. For example:
#   
#   the total number of rooms in a district is not very useful if you don't know how many households there are. What you really want is the number of rooms per household.
# 
# Similarly, the total number of bedrooms by itself is not very useful: you probably want to compare it to the number of rooms.
# 
# And the population per household also seems like an interesting attribute combination to look at.
# 
# Let's create these new attributes:
#   
#   housing_df <- 
#   housing_df %>% 
#   mutate(rooms_per_household = total_rooms/households,
#          bedrooms_per_room = total_bedrooms/total_rooms,
#          population_per_household = population/households)




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Predict labels for empirical networks
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

largest_comp=function(df){
  G=graph_from_data_frame(as.matrix(df),directed=FALSE)
  df.graph=igraph::simplify(G,remove.multiple = T,remove.loops = T)
  Isolated = which(igraph::degree(df.graph)==0)
  Net=igraph::delete.vertices(df.graph, Isolated)
  components = igraph::clusters(Net, mode="weak")
  biggest_cluster_id = which.max(components$csize)
  vert_ids = V(df.graph)[components$membership== biggest_cluster_id]
  graph=igraph::induced_subgraph(df.graph, vert_ids)
}
#create graph object
#FolderName = c('C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi/Networks')
#function(x)read.table(x,header=TRUE,sep="\t")
library(data.table)
Network.Summary <- function(FolderName){
  
  Net.lists = list.files(FolderName, full.names=TRUE)
  Net.Feat.Summary = data.frame()
  
  g = lapply(Net.lists,fread)
  GraphName = cbind(lapply(Net.lists,function(x) gsub(".edges","",x)))
  
  G = lapply(g, largest_comp)
  
  g.features=calcGraphFeatures(G)
  all_graphs=cbind(GraphName, g.features)
  results=list(G,all_graphs)
  return(results)
}

FolderName = c('Animal-Social-Networks')
empdata=Network.Summary(FolderName)
empdata[[1]]
empdata[[2]]



RealNet=empdata[[2]]
RealNet= RealNet%>%dplyr::select(c(GraphName,order,edges,mean_degree,minCut,
                                   FiedlerValue,Normalized_FiedlerValue,closeness,modularity,diameter,
                                   betweenness,transitivity,spectral_radius,centrality_eigen))

#pred=predict(best.model.fit,....)
new_data=as_tibble(RealNet)
new.data=new_data%>%dplyr::select(c(GraphName,order,edges,mean_degree,minCut,
                           FiedlerValue,Normalized_FiedlerValue,closeness,modularity,diameter,
                           betweenness,transitivity,spectral_radius,centrality_eigen))%>%
  mutate_if(is.character,factor)

predicted_label <- predict(best.model.fit, new.data)
predicted_label

#confusionMatrix(predicted_label, train_data$GraphName)

#####Avesaves-wildbird-net1 predicted as Erdos Renyi
er.aves1=AvesWbirdSynth.Net1 [[1]]
sw.aves1=AvesWbirdSynth.Net1[[110]]
sf.linear.aves1=AvesWbirdSynth.Net1[[210]]
sf.nonlinear.aves1=AvesWbirdSynth.Net1[[310]]
sp.aves1=AvesWbirdSynth.Net1[[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[9]],vertex.label=NA,vertex.size=2)
plot(er.aves1,vertex.label=NA,vertex.size=2)
plot(sw.aves1,vertex.label=NA,vertex.size=2)
plot(sf.linear.aves1,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.aves1,vertex.label=NA,vertex.size=2)
plot(sp.aves1,vertex.label=NA,vertex.size=2)


#####Aves wildbird net2 predicted as Spatial
er.aves2=AvesWbirdSynth.Net2 [[1]]
sw.aves2=AvesWbirdSynth.Net2[[110]]
sf.linear.aves2=AvesWbirdSynth.Net2[[210]]
sf.nonlinear.aves2=AvesWbirdSynth.Net2[[310]]
sp.aves2=AvesWbirdSynth.Net2[[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[10]],vertex.label=NA,vertex.size=2)
plot(er.aves2,vertex.label=NA,vertex.size=2)
plot(sw.aves2,vertex.label=NA,vertex.size=2)
plot(sf.linear.aves2,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.aves2,vertex.label=NA,vertex.size=2)
plot(sp.aves2,vertex.label=NA,vertex.size=2)

#####Tortoise predicted as scale-free with linear attachment
er.tortoise=Tortoise.Net[[1]]
sw.tortoise=Tortoise.Net[[110]]
sf.linear.tortoise=Tortoise.Net[[210]]
sf.nonlinear.tortoise=Tortoise.Net[[310]]
sp.tortoise=Tortoise.Net[[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[29]],vertex.label=NA,vertex.size=2)
plot(er.tortoise,vertex.label=NA,vertex.size=2)
plot(sw.tortoise,vertex.label=NA,vertex.size=2)
plot(sf.linear.tortoise,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.tortoise,vertex.label=NA,vertex.size=2)
plot(sp.tortoise,vertex.label=NA,vertex.size=2)

##### Bat roosting predicted as spatial
er.bats1=Bat.Net1 [[1]]
sw.bats1=Bat.Net1[[110]]
sf.linear.bats1=Bat.Net1[[210]]
sf.nonlinear.bats1=Bat.Net1[[310]]
sp.bats1=Bat.Net1[[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[19]],vertex.label=NA,vertex.size=2)
plot(er.bats1,vertex.label=NA,vertex.size=2)
plot(sw.bats1,vertex.label=NA,vertex.size=2)
plot(sf.linear.bats1,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.bats1,vertex.label=NA,vertex.size=2)
plot(sp.bats1,vertex.label=NA,vertex.size=2)

##### Vampire Bat foodsharing predicted as small world
er.bats2=Bat.Net2 [[1]]
sw.bats2=Bat.Net2[[110]]
sf.linear.bats2=Bat.Net2[[210]]
sf.nonlinear.bats2=Bat.Net2[[310]]
sp.bats2=Bat.Net2[[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[28]],vertex.label=NA,vertex.size=2)
plot(er.bats2,vertex.label=NA,vertex.size=2)
plot(sw.bats2,vertex.label=NA,vertex.size=2)
plot(sf.linear.bats2,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.bats2,vertex.label=NA,vertex.size=2)
plot(sp.bats2,vertex.label=NA,vertex.size=2)


#####Ant Colony net1 (day 1) predicted as spatial
er.ants1=Ant.Net1 [[1]]
sw.ants1=Ant.Net1[[110]]
sf.linear.ants1=Ant.Net1[[210]]
sf.nonlinear.ants1=Ant.Net1[[310]]
sp.ants1=Ant.Net1[[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[16]],vertex.label=NA,vertex.size=2)
plot(er.ants1,vertex.label=NA,vertex.size=2)
plot(sw.ants1,vertex.label=NA,vertex.size=2)
plot(sf.linear.ants1,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.ants1,vertex.label=NA,vertex.size=2)
plot(sp.ants1,vertex.label=NA,vertex.size=2)


#####Ant Colony Net2 (day 2) predicted as spatial
er.ants2=Ant.Net2 [[1]]
sw.ants2=Ant.Net2[[110]]
sf.linear.ants2=Ant.Net2[[210]]
sf.nonlinear.ants2=Ant.Net2[[310]]
sp.ants2=Ant.Net2[[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[17]],vertex.label=NA,vertex.size=2)
plot(er.ants2,vertex.label=NA,vertex.size=2)
plot(sw.ants2,vertex.label=NA,vertex.size=2)
plot(sf.linear.ants2,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.ants2,vertex.label=NA,vertex.size=2)
plot(sp.ants2,vertex.label=NA,vertex.size=2)

#####Hyena net1 predicted as Spatial
er.hyena1=Hyena.Net1 [[1]]
sw.hyena1=Hyena.Net1 [[110]]
sf.linear.hyena1=Hyena.Net1 [[210]]
sf.nonlinear.hyena1=Hyena.Net1 [[310]]
sp.hyena1=Hyena.Net1 [[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[22]],vertex.label=NA,vertex.size=2)
plot(er.hyena1,vertex.label=NA,vertex.size=2)
plot(sw.hyena1,vertex.label=NA,vertex.size=2)
plot(sf.linear.hyena1,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.hyena1,vertex.label=NA,vertex.size=2)
plot(sp.hyena1,vertex.label=NA,vertex.size=2)


#####Hyena net2 predicted as Spatial
er.hyena2=Hyena.Net2 [[1]]
sw.hyena2=Hyena.Net2 [[110]]
sf.linear.hyena2=Hyena.Net2 [[210]]
sf.nonlinear.hyena2=Hyena.Net2 [[310]]
sp.hyena2=Hyena.Net2 [[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[23]],vertex.label=NA,vertex.size=2)
plot(er.hyena2,vertex.label=NA,vertex.size=2)
plot(sw.hyena2,vertex.label=NA,vertex.size=2)
plot(sf.linear.hyena2,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.hyena2,vertex.label=NA,vertex.size=2)
plot(sp.hyena2,vertex.label=NA,vertex.size=2)


#####Hyena net3 predicted as Spatial
er.hyena3=Hyena.Net3 [[1]]
sw.hyena3=Hyena.Net3 [[110]]
sf.linear.hyena3=Hyena.Net3 [[210]]
sf.nonlinear.hyena3=Hyena.Net3 [[310]]
sp.hyena3=Hyena.Net3 [[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[24]],vertex.label=NA,vertex.size=2)
plot(er.hyena3,vertex.label=NA,vertex.size=2)
plot(sw.hyena3,vertex.label=NA,vertex.size=2)
plot(sf.linear.hyena3,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.hyena3,vertex.label=NA,vertex.size=2)
plot(sp.hyena3,vertex.label=NA,vertex.size=2)


#####Cattle dominance predicted as Scale-free with non linear attachment
er.cattle=CattleDom.Net [[1]]
sw.cattle=CattleDom.Net [[110]]
sf.linear.cattle=CattleDom.Net [[210]]
sf.nonlinear.cattle=CattleDom.Net [[310]]
sp.cattle=CattleDom.Net [[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[20]],vertex.label=NA,vertex.size=2)
plot(er.cattle,vertex.label=NA,vertex.size=2)
plot(sw.cattle,vertex.label=NA,vertex.size=2)
plot(sf.linear.cattle,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.cattle,vertex.label=NA,vertex.size=2)
plot(sp.cattle,vertex.label=NA,vertex.size=2)

#####Dolphine social predicted as Erdos Renyi
er.dolphine=DolphinSoc.Net [[1]]
sw.dolphine=DolphinSoc.Net [[110]]
sf.linear.dolphine=DolphinSoc.Net [[210]]
sf.nonlinear.dolphine=DolphinSoc.Net [[310]]
sp.dolphine=DolphinSoc.Net [[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[21]],vertex.label=NA,vertex.size=2)
plot(er.dolphine,vertex.label=NA,vertex.size=2)
plot(sw.dolphine,vertex.label=NA,vertex.size=2)
plot(sf.linear.dolphine,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.dolphine,vertex.label=NA,vertex.size=2)
plot(sp.dolphine,vertex.label=NA,vertex.size=2)

#####sparrowlyon net1 predicted spatial
er.sparrowlyon1=Sparrowlyon.Net1 [[1]]
sw.sparrowlyon1=Sparrowlyon.Net1 [[110]]
sf.linear.sparrowlyon1=Sparrowlyon.Net1 [[210]]
sf.nonlinear.sparrowlyon1=Sparrowlyon.Net1 [[310]]
sp.sparrowlyon1=Sparrowlyon.Net1 [[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[5]],vertex.label=NA,vertex.size=2)
plot(er.sparrowlyon1,vertex.label=NA,vertex.size=2)
plot(sw.sparrowlyon1,vertex.label=NA,vertex.size=2)
plot(sf.linear.sparrowlyon1,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.sparrowlyon1,vertex.label=NA,vertex.size=2)
plot(sp.sparrowlyon1,vertex.label=NA,vertex.size=2)

#####sparrowlyon net2 predicted spatial
er.sparrowlyon2=Sparrowlyon.Net2 [[1]]
sw.sparrowlyon2=Sparrowlyon.Net2 [[110]]
sf.linear.sparrowlyon2=Sparrowlyon.Net2 [[210]]
sf.nonlinear.sparrowlyon2=Sparrowlyon.Net2 [[310]]
sp.sparrowlyon2=Sparrowlyon.Net2 [[405]]
par(mfrow=c(2,3))
plot(empdata[[1]][[6]],vertex.label=NA,vertex.size=2)
plot(er.sparrowlyon2,vertex.label=NA,vertex.size=2)
plot(sw.sparrowlyon2,vertex.label=NA,vertex.size=2)
plot(sf.linear.sparrowlyon2,vertex.label=NA,vertex.size=2)
plot(sf.nonlinear.sparrowlyon2,vertex.label=NA,vertex.size=2)
plot(sp.sparrowlyon2,vertex.label=NA,vertex.size=2)
# spweaver=fastSpatialNetwork(n =vcount(empdata[[1]][[1]]), r = 0.14, makeConnected=T)
# spweaver
# sfweaver=sample_pa(n=vcount(empdata[[1]][[2]]),power=1,m=3,directed=FALSE)
# sfweaver
# par(mfrow=c(1,3))
# plot(empdata[[1]][[2]]);plot(spweaver);plot(sfweaver)#plot(Hyena.Net2[[410]])#plot(AvesWbirdSynth.Net5[[310]])

# ####--Compare--Graph--Features--for--empirical--and--the--synthetic--class
# emp.g=empdata[[1]][[7]]
# emp.g$name="wildbirdnet7"
# emp.g$type="wildbirdnet7"
# emp.g$id="1"
# 
# synth.g=AvesWbirdSynth.Net7[[401]]
# t1=list(emp.g,synth.g)
# RunSimOnGraphFeatures(t1,nreps = 1)

# 
# 
# #varImpPlot(RfClassificationmodel)
# 
# 
# 
# library(tidymodels)
# 
# data(ames)
# 
 #set.seed(4595)
# data_split <- initial_split(ames, strata = "Sale_Price", prop = 0.75)
# 
# ames_train <- training(data_split)
# ames_test  <- testing(data_split)
# 
# 
# ##########--------------------------MODEL--1--(RANDOM FOREST)-----------#####################
# rfmodel<- rand_forest(mode = "regression")
# rfmodel
# 
# # You can create create regression models with the tidymodels package 'parsnip'
# # to predict continuous or numeric quantities.
# # The parsnip package provides two different interfaces to fit a model:
# #   
# # [1] the non-formula interface (fit_xy()).
# # [2] the formula interface (fit()), and
# 
# ##------------[1] the non-formula model (fit_xy())-----------
# preds <- c("Longitude", "Latitude", "Lot_Area", "Neighborhood", "Year_Sold")
# 
# rf_nonform_fit <- 
#   rfmodel %>%
#   set_engine("ranger") %>%
#   fit_xy(
#     x = ames_train[, preds],
#     y = log10(ames_train$Sale_Price)
#   )
# 
# rf_nonform_fit    
# 
# # The non-formula interface doesn't 
# # do anything to the predictors before passing them to the underlying model function
# 
# ##--Predicting with non-formula model
# test_results_non_form<- 
#   ames_test %>%
#   select(Sale_Price) %>%
#   mutate(Sale_Price = log10(Sale_Price)) %>%
#   bind_cols(
#     predict(rf_nonform_fit, new_data = ames_test[, preds])
#   )
# test_results_non_form
# # Note that:
# # If the model required indicator variables, we would have to create them manually prior to using fit() (perhaps using the recipes package).
# # We had to manually log the outcome prior to modeling.
# 
# ##-------[2] the formula model (fit()) with tunning parameters-----------
rf_form_fit=
  rand_forest(mode = "regression", mtry = 3, trees = 1000) %>%
  set_engine("ranger") %>%
  fit(
    log10(Sale_Price) ~ Longitude + Latitude + Lot_Area + Neighborhood + Year_Sold,
    data = ames_train
  )

rf_form_fit

# ##--Predicting with formula model
test_results_form <-
  ames_test %>%
  select(Sale_Price) %>%
  mutate(Sale_Price = log10(Sale_Price)) %>%
  bind_cols(
    predict(rf_form_fit, new_data = ames_test[,])
  )
test_results_form

# #Note:
# # Suppose that we would like to use the randomForest package instead of ranger. To do so, 
# #the only part of the syntax that needs to change is the set_engine() argument:
# 
# rand_forest(mode = "regression", mtry = 3, trees = 1000) %>%
#   set_engine("randomForest") %>%
#   fit(
#     log10(Sale_Price) ~ Longitude + Latitude + Lot_Area + Neighborhood + Year_Sold,
#     data = ames_train
#   )
# 
# # Now suppose that we want to modify the value of mtry 
# # based on the number of predictors in the data. Usually,
# # a good default value is floor(sqrt(num_predictors)) but a 
# # pure bagging model requires an mtry value equal to the total number of parameters.
# # There may be cases where you may not know how many predictors are going to be present 
# # when the model will be fit (perhaps due to the generation of indicator variables or 
# 
# 
# # When the model it being fit by parsnip, data descriptors are made available. 
# # These attempt to let you know what you will have available when the model is fit. 
# # When a model object is created (say using rand_forest()), 
# # the values of the arguments that you give it are immediately evaluated unless you delay them. 
# # To delay the evaluation of any argument, you can used rlang::expr() to make an expression.
# # 
# # Two relevant data descriptors for our example model are:
# #   
# # .preds(): the number of predictor variables in the data set that are associated with the predictors prior to dummy variable creation.
# # .cols(): the number of predictor columns after dummy variables (or other encodings) are created.
# # Since ranger won't create indicator values, .preds() would be appropriate for mtry for a bagging model.
# #For example, let's use an expression with the .preds() descriptor to fit a bagging model:
# 
# rand_forest(mode = "regression", mtry = .preds(), trees = 1000) %>%
#   set_engine("ranger") %>%
#   fit(
#     log10(Sale_Price) ~ Longitude + Latitude + Lot_Area + Neighborhood + Year_Sold,
#     data = ames_train
#   )


## rf_form_fit=
#   rand_forest(mode = "regression", mtry = 3, trees = 1000) %>%
#   set_engine("ranger") %>%
#   fit(
#     log10(Sale_Price) ~ Longitude + Latitude + Lot_Area + Neighborhood + Year_Sold,
#     data = ames_train
#   )
# 
# rf_form_fit
# 
# ##--Predicting with formula model
# test_results_form <- 
#   ames_test %>%
#   select(Sale_Price) %>%
#   mutate(Sale_Price = log10(Sale_Price)) %>%
#   bind_cols(
#     predict(rf_form_fit, new_data = ames_test[, preds])
#   )
# test_results_form
# 
# ##########---------------MODEL--2--(Regularised Linear regression)-----------#####################
# # A linear model might work for this data set as well. We can use the linear_reg() parsnip model. 
# # There are two engines that can perform regularization/penalization, 
# # the glmnet and sparklyr package. The glmnet package only implements a non-formula method,
# # but parsnip will allow either one to be used
# 
# #Note:
# # When regularization is used, the predictors should first be centered and scaled before
# # being passed to the model. The formula method won't do that automatically so we will need to do this
# # ourselves. We'll use the recipes package for these steps.
# 
# norm_recipe <- 
#   recipe(
#     Sale_Price ~ Longitude + Latitude + Lot_Area + Neighborhood + Year_Sold, 
#     data = ames_train
#   ) %>%
#   step_other(Neighborhood) %>% 
#   step_dummy(all_nominal()) %>%
#   step_center(all_predictors()) %>%
#   step_scale(all_predictors()) %>%
#   step_log(Sale_Price, base = 10) %>% 
#   # estimate the means and standard deviations
#   prep(training = ames_train, retain = TRUE)
# 
# ##--Now let's fit the model using the processed version of the data
# glmn_fit <- 
#   linear_reg(penalty = 0.001, mixture = 0.5) %>% 
#   set_engine("glmnet") %>%
#   fit(Sale_Price ~ ., data = bake(norm_recipe, new_data = NULL))
# 
# glmn_fit
# 
# #Note: If penalty were not specified, all of the lambda values would be computed.
# #To get the predictions for this specific value of lambda (aka penalty):
# 
# # First, get the processed version of the test set predictors:
# test_normalized <- bake(norm_recipe, new_data = ames_test, all_predictors())
# 
# # Predict with the fitted model
# test_results <- 
#   test_results %>%
#   rename(`random forest` = .pred) %>%
#   bind_cols(
#     predict(glmn_fit, new_data = test_normalized) %>%
#       rename(glmnet = .pred)
#   )
# test_results 
# 
# #####--Evaluation-Metrics--#####
# test_results %>% metrics(truth = Sale_Price, estimate = glmnet) 
# 
# ###--Predicted --Plots--with--models--
# test_results %>% 
#   gather(model, prediction, -Sale_Price) %>% 
#   ggplot(aes(x = prediction, y = Sale_Price)) + 
#   geom_abline(col = "green", lty = 2) + 
#   geom_point(alpha = .4) + 
#   facet_wrap(~model) + 
#   coord_fixed()
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ PREDICTION MODEL
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# #+ Training a machine learning regression model with random forest
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# #Regression models
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
# ##--Data importation and processing
# library(dplyr)
# df50=read.csv("df_mu_50.csv",sep = ",", header = T)
# df100=read.csv("df_mu_100.csv",sep = ",", header = T)
# df150=read.csv("df_mu_150.csv",sep = ",", header = T)
# 
# Data=rbind(df50,df100,df150)
# Data=df50
# 
# df=as_tibble(Data)
# df%>%count(GraphName,sort=T)
# df%>%View()
# 
# Data= df%>%dplyr::select(-c(connected,max_component,minDegree,maxDegree,
#                             threshold,GraphReplicate,GraphID,X,GraphName))%>%
#   mutate_if(is.character,factor)
# ##--Shuffle--data
# Data<- Data[sample(1:nrow(Data)), ]##shuffle row indices and randomly re order 
# head(Data)
# 
# 
# #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # Data prepropecessing recipe
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# library(visdat)
# 
# #Visualizing all data structure
# vis_dat(Data)
# 
# #Checking levels of factor or character nd creating a table
# library(gt)
# Data %>% 
#   count(GraphName,
#         sort = TRUE)%>%
#   gt()
# 
# set.seed(4595)
# data_split <- rsample::initial_split(Data, strata = "mu")#, prop = 0.75)
# 
# train.data <- training(data_split)
# test.data  <- testing(data_split)
# 
# dim(train.data)
# 
# ##---Creating--Data--Recipe---and--Prep--
# df_recipe=recipe(mu~., data=train.data)%>%
#   update_role(order,new_role = "numberOfNodes")%>%
#   step_zv(all_predictors())%>%
#   step_normalize(all_predictors())%>%
#   step_corr(all_predictors(), threshold = 0.7, method = "spearman")
# 
# summary(df_recipe)
# 
# df_prep <- 
#   df_recipe %>% # use the recipe object
#   prep() %>% # perform the recipe on training data
#   juice() # extract only the preprocessed dataframe 
# 
# #Take a look at the data structure:
# glimpse(df_prep)
# 
# ####------Setting--up--Cross--validation--set--to--tunes--the--models--we--are--going--to--compare
# set.seed(100)
# mu_folds <-
#   bootstraps(train.data,strata = mu)
# 
# mu.folds <- 
#   recipes::bake(
#     xgb.recipe, 
#     new_data = training(data_split)
#   ) %>%  
#   rsample::vfold_cv(v = 5)
# # mu_folds <-
# #   vfold_cv(train.data, 
# #            v = 10, 
# #            strata = mu)
# 
# library(usemodels)
# #This package contains information for options on setting up common types of models
# #eg
# use_ranger(mu~.,data=train.data)
# ##########------------SECIFYING--MODEL:SETTING--MODEL--ENGINE--AND--WORKFLOW----------#####################
# 
# # The process of specifying our models is always as follows:
# ###Pick a model type
# ###set the engine
# ###Set the mode: regression or classification
# 
# #######################--Setting---up---model--and--engine--###############################
# 
# ##--[A]--Random--Forest---model-engine and workflow---
# library(ranger)
# 
# rf.recipe=recipe(formula=mu~., data=train.data)%>%
#   update_role(order,new_role = "numberOfNodes")
# # step_zv(all_predictors())%>%
# #  step_normalize(all_predictors())%>%
# # step_corr(all_predictors(), threshold = 0.7, method = "spearman")
# 
# summary(rf.recipe)
# 
# rf.prep <- 
#   rf.recipe %>% # use the recipe object
#   prep() %>% # perform the recipe on training data
#   juice() # extract only the preprocessed dataframe 
# 
# #Take a look at the data structure:
# glimpse(rf.prep)
# 
# ### Parallel running
# #cores <- parallel::detectCores()
# #cores
# 
# rf.model <- 
#   rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>% 
#   set_engine("ranger") %>% 
#   set_mode("regression")
# 
# ##---Create--workflows
# # To combine the data preparation recipe with the model building, 
# # we use the package workflows. 
# #A workflow is an object that can bundle together your pre-processing recipe,
# #modeling, and even post-processing requests (like calculating the RMSE).   
# 
# rf.wkflow <-
#   workflow() %>%
#   add_recipe(rf.recipe) %>% 
#   add_model(rf.model) 
# 
# ## shows num of params to be tuned
# extract_parameter_set_dials(rf.wkflow)
# 
# set.seed(90223)
# doParallel::registerDoParallel()
# 
# rf.tune.params <-
#   tune_grid(rf.wkflow, 
#             resamples = mu_folds, 
#             grid = 11)
# 
# 
# ###---Explore results from tuned model
# show_best(rf.tune.params, metric="rmse")
# show_best(rf.tune.params, metric="rsq")
# 
# autoplot(rf.tune.params )
# 
# rf.final.wkflow<-rf.wkflow%>%
#   finalize_workflow(select_best(rf.tune.params ))
# 
# ##--Final fitted model
# rf.model.fit=last_fit(rf.final.wkflow,data_split )
# ##--collect metrics
# collect_metrics(rf.model.fit)
# 
# collect_predictions(rf.model.fit)%>%
#   ggplot(aes(mu,.pred))+
#   geom_abline(lty=2,color="blue")+
#   geom_point(alpha=0.5,color="red")+
#   coord_fixed()
# 
# ##---Making Predictions
# predicted.vals=predict(rf.model.fit$.workflow[[1]],test.data[,])
# 
# ##--Variable--importance
# library(vip)
# varimp.model=rf.model%>%
#   finalize_model(select_best(rf.tune.params ))%>%
#   set_engine("ranger",importance="permutation")
# 
# workflow() %>%
#   add_recipe(rf.recipe) %>% 
#   add_model(varimp.model)%>% 
#   fit(train.data)%>% 
#   pull_workflow_fit()%>% 
#   vip(aesthetics=list(alpha=0.8,fill="midnightblue"))
# 
# ##### simple example
# rf_form_fit=
#   rand_forest(mode = "regression", mtry = 24, trees = 1000) %>%
#   set_engine("ranger") %>%
#   fit(mu~.,
#       data = train.data
#   )
# # 
# rf_form_fit
# 
# ##--Predicting with formula model
# train_results_form <-
#   train.data %>%
#   select(mu) %>%
#   bind_cols(
#     predict(rf_form_fit, new_data = train.data[,])
#   )
# 
# train_results_form
# 
# test_results_form <-
#   test.data %>%
#   select(mu) %>%
#   bind_cols(
#     predict(rf_form_fit, new_data = test.data[,])
#   )
# test_results_form
# 
# #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ##--[B]--Boosted--tree--(XGBoost)
# #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# library(xgboost)
# xgb.recipe <- 
#   recipe(formula = mu ~ ., data = train.data) %>% 
#   step_zv(all_predictors()) 
# 
# 
# xgb.model <- 
#   boost_tree(rees = tune(), min_n = tune(), tree_depth = tune(), learn_rate = tune(),
#              loss_reduction = tune(), sample_size = tune()) %>% 
#   set_engine("xgboost") %>% 
#   set_mode("regression")
# 
# 
# xgb.wkflow <-
#   workflow() %>%
#   add_recipe(df_recipe) %>% 
#   add_model(xgb.model)
# 
# set.seed(90223)
# doParallel::registerDoParallel()
# 
# xgb.tune.params <-
#   tune_grid(xgb.wkflow, 
#             resamples = mu_folds, 
#             grid = 11)
# 
# ##------General--example----###
# xgb.recipe <- 
#   recipes::recipe(mu ~ ., data = training(data_split)) %>%
#   # convert categorical variables to factors
#   recipes::step_string2factor(all_nominal()) %>%
#   # combine low frequency factor levels
#   recipes::step_other(all_nominal(), threshold = 0.01) %>%
#   # remove no variance predictors which provide no predictive information 
#   recipes::step_nzv(all_nominal()) %>%
#   recipes::step_zv(all_predictors()) %>%
#   prep()
# 
# 
# # XGBoost model specification
# xgb.model <- 
#   parsnip::boost_tree(
#     mode = "regression",
#     trees = 1000,
#     min_n = tune(),
#     tree_depth = tune(),
#     learn_rate = tune(),
#     loss_reduction = tune()
#   ) %>%
#   set_engine("xgboost", objective = "reg:squarederror")
# 
# 
# # grid specification
# xgb.params <- 
#   dials::parameters(
#     min_n(),
#     tree_depth(),
#     learn_rate(),
#     loss_reduction()
#   )
# 
# # grid space
# xgb.grid <- 
#   dials::grid_max_entropy(
#     xgb.params, 
#     size = 60
#   )
# knitr::kable(head(xgb.grid))
# 
# # Workflow
# xgb.wkflow <- 
#   workflows::workflow() %>%
#   add_model(xgb.model) %>% 
#   add_formula(mu ~ .)
# 
# # hyperparameter tuning
# xgb.tuned <- tune::tune_grid(
#   object = xgb.wkflow,
#   resamples = mu.folds,
#   grid = xgb.grid,
#   metrics = yardstick::metric_set(rmse, rsq, mae),
#   control = tune::control_grid(verbose = TRUE)
# )
# 
# # Best hyper parameter values which minimises rmse
# xgb.tuned %>%
#   tune::show_best(metric = "rmse") %>%
#   knitr::kable()
# 
# 
# #isolate the best performing hyperparameter values.
# xgb.best.params <- xgb.tuned %>%
#   tune::select_best("rmse")
# knitr::kable(xgb.best.params)
# 
# #Finalize the XGBoost model to use the best tuning parameters.
# xgb.final.model<- xgb.model %>% 
#   finalize_model(xgb.best.params)
# 
# ##--Evauating--Performance---of---model--on--train--and--test--data
# train_processed <- bake(xgb.recipe,  new_data = training(data_split))
# 
# train_prediction <- xgb.final.model %>%
#   # fit the model on all the training data
#   fit(
#     formula = mu ~ ., 
#     data    = train_processed
#   ) %>%
#   # predict mu for the training data
#   predict(new_data = train_processed) %>%
#   bind_cols(training(data_split))
# 
# xgb.score.train <- 
#   train_prediction %>%
#   yardstick::metrics(mu, .pred) %>%
#   mutate(.estimate = format(round(.estimate, 2), big.mark = ","))
# knitr::kable(xgb.score.train)
# 
# #cbind(head(train_prediction$mu),head(train.data$mu))
# 
# ##------Predicting on test set
# test_processed  <- bake(xgb.recipe, new_data = testing(data_split))
# test_prediction <- xgb.final.model %>%
#   # fit the model on all the training data
#   fit(
#     formula = mu ~ ., 
#     data    = train_processed
#   ) %>%
#   # use the training model fit to predict the test data
#   predict(new_data = test_processed) %>%
#   bind_cols(testing(data_split))
# 
# # measure the accuracy of our model using `yardstick`
# xgb.score.test <- 
#   test_prediction %>%
#   yardstick::metrics(mu, .pred) %>%
#   mutate(.estimate = format(round(.estimate, 2), big.mark = ","))
# 
# knitr::kable(xgb.score.test)
# 
# 
# 
# #To quickly check that there is not an obvious issue with our model's predictions, 
# #let's plot the test data residuals.
# 
# mu_prediction_residual <- test_prediction %>%
#   arrange(.pred) %>%
#   mutate(residual_pct = (mu - .pred) / .pred) %>%
#   select(.pred, residual_pct)
# ggplot(mu_prediction_residual, aes(x = .pred, y = residual_pct)) +
#   geom_point() +
#   xlab("Predicted mu") +
#   ylab("Residual (%)") +
#   scale_x_continuous(labels = scales::dollar_format()) +
#   scale_y_continuous(labels = scales::percent)
# 
# 
# 
# 
# 
# 
# ##--[C]--GLM--network
# glm.model <-
#   linear_reg(penalty = tune(),mixture = 1) %>% 
#   set_engine("glmnet") %>%
#   set_mode("regression")
# 
# glm.wkflow <-
#   workflow() %>%
#   add_recipe(df_recipe) %>% 
#   add_model(glm.model)
# lr_reg_grid <- tibble(penalty = 10^seq(-4, -1, length.out = 30))
# 
# 
# 
# ##############-----------------Evaluate------models--------------------############
# # Now we can use our validation 
# # set (cv_folds) to estimate the performance of our models using the fit_resamples() 
# # function to fit the models on each of the folds and store the results.
# # 
# # Note that fit_resamples() will fit our model to each resample and 
# # evaluate on the heldout set from each resample. 
# # The function is usually only used for computing performance metrics across some 
# # set of resamples to evaluate our models (like accuracy) - the models are
# # not even stored. However, in our example we save the predictions in order to
# # visualize the model fit and residuals with control_resamples(save_pred = TRUE).
# #Finally, we collect the performance metrics with collect_metrics() and 
# #pick the model that does best on the validation set.    
# 
# rf.fit <-
#   rf.wkflow %>% 
#   fit(mu~.,data = train.data)
# 
# 
# rf.res %>%  collect_metrics(summarize = TRUE)
# 
# 
# xgb.res <- 
#   xgb.wkflow %>% 
#   fit_resamples(
#     resamples = cv_folds, 
#     metrics = metric_set(
#       recall, precision, f_meas, 
#       accuracy, kap,
#       roc_auc, sens, spec),
#     control = control_resamples(save_pred = TRUE)
#   ) 
# 
# xgb.res %>% collect_metrics(summarize = TRUE)

