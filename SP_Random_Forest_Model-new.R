library(dplyr)
library(ggplot2)
library("mgcv")
library(elasticnet)
library("gam")
library("brnn")
library("nnet")
library("glmnet")
library('caret')
library("clusterSim")
library("neuralnet")
library("randomForest")
library("gbm")
library("plyr")
library("rpart")
library("mlbench")
library("MLmetrics")
library("corrplot")
library("readxl")
library("varImp")
#install.packages("caret", dependencies = c("Depends", "Suggests"))


#####----Random---Forest--(RF)---#########
# Random forest Aggregate of the results of multiple predictors gives a better prediction 
# than the best individual predictor. A group of predictors is called an ensemble.
# Thus, this technique is called Ensemble Learning
# RF trainS a group of Decision Tree classifiers, each on a different random subset of the train set.
# To make a prediction, we just obtain the predictions of all individuals trees, 
# then predict the class that gets the most votes

#####---[1]--Imprting--Data--#####
# data_train <- read.csv("https://raw.githubusercontent.com/guru99-edu/R-Programming/master/train.csv")
# data_train= na.omit(data_train)
# data_test <- read.csv("https://raw.githubusercontent.com/guru99-edu/R-Programming/master/test.csv") 
# data_test= na.omit(data_train)

df4=read.csv("spdta.csv",header = T, sep=",")

data_df <- df4[sample(1:nrow(df4)), ]

data_df$GraphID=1:100
head(data_df)
# Drop variables
clean_df_data<- data_df %>%
  select(-c(GraphID,order, GraphReplicate, order,max_component,minDegree,
            threshold,maxDegree,X,GraphName,connected)) 

glimpse(clean_df_data)

df.sample <- sample(c(TRUE, FALSE), nrow(clean_df_data), replace=TRUE, prob=c(0.8,0.2))
train  <- clean_df_data[df.sample, ]
test   <- clean_df_data[!df.sample, ]



# One way to evaluate the performance of a model is to train it on a number of 
# different smaller datasets and evaluate them over the other smaller testing set. 
# This is called the F-fold cross-validation feature.
# R has a function to randomly split number of datasets of almost the same size.
# For example, if k=9, the model is evaluated over the nine folder and 
# tested on the remaining test set. This process is repeated until all 
# the subsets have been evaluated. This technique is widely used for model selection, 
#especially when the model has parameters to tune.

# We can now choose the parameters that best generalized the data.
# Random forest chooses a random subset of features and builds many Decision Trees. 
# The model averages out all the predictions of the Decisions trees.


# Random forest has some parameters that can be changed to improve the generalization
# of the prediction. You will use the function RandomForest() to train the model.

#Syntax for Randon Forest is
#RandomForest(formula, ntree=n, mtry=FALSE, maxnodes = NULL)
# Arguments:
#   - Formula: Formula of the fitted model
# - ntree: number of trees in the forest
# - mtry: Number of candidates draw to feed the algorithm. By default, it is the square of the number of columns.
# - maxnodes: Set the maximum amount of terminal nodes in the forest
# - importance=TRUE: Whether independent variables importance in the random forest be assessed

# Note: Random forest can be trained multiple parameters.
# 
# Tuning a model is very tedious work.There are lot of combination possible between the parameters.
# we dont have the time to try all of them.
# A good alternative is to let the machine find the best combination for you.
# There are two methods available to do that in RF:
# (A) Grid Search
# (B) Random Search

#(A) Grid Search
# The grid search method is simple,
# the model will be evaluated over all the combination you pass in the function,
# using cross-validation.For instance, you want to try the model with 10, 20, 30 number of trees
# and each tree will be tested over a number of mtry equals to 1, 2, 3, 4, 5.
# Then the machine will test 15 different models:

# The algorithm will evaluate:
# RandomForest(formula, ntree=10, mtry=1)
# RandomForest(formula, ntree=10, mtry=2)
# ......... etc
# 
# The number of experimentations becomes very high with the cross validation approach in the
# Grid search when using RF. Random search overcomes this caveats posed by grid search

#(B) Random search
# Random search will not evaluate all the combination of hyperparameter in the
# searching space like the grid search. Instead, it will randomly choose combination 
# at every iteration. The advantage is it lower the computational cost.

#--------Set---the---control--parameter
# We will proceed as follow to construct and evaluate the model:
# .Evaluate the model with the default setting
# .Find the best number of mtry
# .Find the best number of maxnodes
# .Find the best number of ntrees
# .Evaluate the model on the test dataset


###----Default setting----###
# K-fold cross validation is controlled by the trainControl() function

#trainControl(method = "cv", number = n, search ="grid")
# arguments
# - method = "cv": The method used to resample the dataset. 
# - number = n: Number of folders to create
# - search = "grid": Use the search grid method. For randomized method, use "rand"
# Note: You can refer to the vignette to see the other arguments of the function

####---Define the control
#eg
split1=trainControl(method = "cv",number = 10, search = "grid")
split2=trainControl(method = "repeatedcv",number = 10,repeats = 4,savePredictions="final",
classProbs = TRUE) #classProbs used for only classification models
# We will use caret library to evaluate the model.
# The library has one function called train() to evaluate almost all machine learning algorithm.

#####--[2]--Training--model--####
###--Basic--Syntax--###
# train(formula, df, method = "rf", metric= "Accuracy", trControl = trainControl(), tuneGrid = NULL)
# argument
# - `formula`: Define the formula of the algorithm
# - `method`: Define which model to train. Note, at the end of the tutorial, there is a list of all the models that can be trained
# - `metric` = "Accuracy": Define how to select the optimal model
# - `trControl = trainControl()`: Define the control parameters
# - `tuneGrid = NULL`: Return a data frame with all the possible combination

#eg
set.seed(1234)
###---Run--the--model

my_model.1<- train(r~.,
                    data = train,
                    method = "rf",
                    metric = "RMSE",
                    trControl = split1)
# Print the results
print(my_model.1)
plot(my_model.1)

#####----Explanation-----####
# split=trainControl(method=”cv”, number=10, search=”grid”): Evaluate the model with a grid search of 10 folder
# train(…): Train a random forest model. Best model is chosen with the accuracy measure.

#####---[3]---Search--the--best--mtry(Number of candidates/features draw to feed the algorithm. By default, it is the square of the number of columns)
#You can test the model with values of mtry from 1 to 10
set.seed(1234)
tuneGrid <- expand.grid(.mtry = c(1: 10))#vector of possible mtry
my_model_mtry <- train(r~.,
                 data = train,
                 method = "rf",
                 metric = "RMSE",
                 tuneGrid = tuneGrid,
                 trControl = split1,
                 importance = TRUE,
                 nodesize = 14,
                 ntree = 300)
print(my_model_mtry)
plot(my_model_mtry)

#The best value of mtry is stored in the position: my_model.2$bestTune$mtry
#You can also store the value corresponding to this position: max(my_model.2$results$RMSE) and 
#use it when you need to tune the other parameters:

#####---[3]---Search--the--best--maxnodes
# We need to create a loop to evaluate the different values of maxnodes. Thus,
# .Create a list
# .Create a variable with the best value of the parameter mtry; Compulsory
# .Create the loop
# .Store the current value of maxnode
# .Summarize the results

store_maxnode <- list()
best_mtry=my_model_mtry$bestTune$mtry
tuneGrid <- expand.grid(.mtry = best_mtry)
for (maxnodes in c(15: 30)) {
  set.seed(1234)
  my_model_maxnode <- train(r~.,
                      data = train,
                      method = "rf",
                      metric = "RMSE",
                      tuneGrid = tuneGrid,
                      trControl = split1,
                      importance = TRUE,
                      nodesize = 14,
                      maxnodes = maxnodes,
                      ntree = 300)
  
  current_iteration <- toString(maxnodes)
  
  store_maxnode[[current_iteration]] <- my_model_maxnode
}
results_mtry <- resamples(store_maxnode)
summary(results_mtry)

#plot(results_mtry$values$`15~RMSE`, type="l")

# Code explanation:
# store_maxnode <- list(): The results of the model will be stored in this list
# expand.grid(.mtry=best_mtry): Use the best value of mtry
# for (maxnodes in c(15:25)) { … }: Compute the model with values of maxnodes starting from 15 to 25.
# maxnodes=maxnodes: For each iteration, maxnodes is equal to the current value of maxnodes. i.e 15, 16, 17, …
# key <- toString(maxnodes): Store as a string variable the value of maxnode.
# store_maxnode[[key]] <- my_model_maxnode: Save the result of the model in the list.
# resamples(store_maxnode): Arrange the results of the model
# summary(results_mtry): Print the summary of all the combination.
# We can try with  a higher score (higher maxnodes) to see if we can get a higher accuracy


#####---[4]---Search--the--best--ntrees----####
# Now that we have the best value of mtry and maxnode, 
# we can tune the number of trees. The method is exactly the same as maxnode.

store_maxtrees <- list()
for (ntree in c(250, 300, 350, 400, 450, 500, 550, 600, 800, 1000, 2000)) {
  set.seed(5678)
  my_model_maxtrees <- train(r~.,
                       data = train,
                       method = "rf",
                       metric = "RMSE",
                       tuneGrid = tuneGrid,
                       trControl = split1,
                       importance = TRUE,
                       nodesize = 14,
                       maxnodes = 24,
                       ntree = ntree)
  key <- toString(ntree)
  store_maxtrees[[key]] <- my_model_maxtrees
}
results_tree <- resamples(store_maxtrees)
summary(results_tree)

# We have our final model. You can train the random forest with the following parameters:
#   
# ntree =250:250 trees will be trained
# mtry=3: 3 features is chosen for each iteration
# maxnodes = 15: Maximum 15 nodes in the terminal nodes (leaves)

fit_my_model <- train(r~.,
                train,
                method = "rf",
                metric = "RMSE",
                tuneGrid = tuneGrid,
                trControl = split1,
                importance = TRUE,
                nodesize = 13,
                ntree = 250,
                maxnodes = 15)

fit_my_model
#####---[5]---Model--Evaluation----####
#The library caret has a function to make prediction.
#syntax: predict(model, newdata= df)
# argument
# - `model`: Define the model evaluated before. 
# - `newdata`: Define the dataset to make prediction
prediction <-predict(fit_my_model , test)
plot(prediction,type="l", lty=2)
# we can use the prediction to compute the confusion matrix (for classification) and see the accuracy score
# 
#confusionMatrix(prediction, test$r)


#####---[6]---Evaluating--feature--importance----####
#We canlook at the feature importance with the function varImp()
#the important features are likely to appear closer to the root of the tree, 
#while less important features will often appear closed to the leaves.

###--Variable--importance---####
v.imp1=caret::varImp(fit_my_model)
plot(v.imp1)


# cvarImpPlot(fit_my_model)
# partialPlot(fit_my_model)
# partialPlot(fit_my_model,  data.frame(train))



###############---------------Different---Method----------------------------------------##############

#RANDOM--FOREST--direct--fitting

# df <- sample(c(TRUE, FALSE), nrow(clean_df_data), replace=TRUE, prob=c(0.8,0.2))
# train  <- clean_df_data[df.sample, ]
# test   <- clean_df_data[!df.sample, ]

#Random Forest Modelling
set.seed(863)
df= data_df %>%
  select(-c(order, GraphReplicate, order,max_component,minDegree,
            threshold,maxDegree,X)) 
newdf=write.csv(df,"sp2.csv")



model.rand.1 <- randomForest::randomForest(r~ ., data = train, importance=TRUE) 
attributes(model.rand.1)

#Conditional=True, adjusts for correlations between predictors.
var_scores <- caret::varImp(model.rand.1, conditional=TRUE)
var_scores

#Gathering rownames in 'var'  and converting it to the factor
#to provide 'fill' parameter for the bar chart. 
var_scores <- var_scores %>% tibble::rownames_to_column("var") 
var_scores$var<- var_scores$var %>% as.factor()

#Plotting the bar and polar charts for comparing variables
var_bar <- ggplot(data = var_scores) + 
  geom_bar(
    stat = "identity",#it leaves the data without count and bin
    mapping = aes(x = var, y=Overall, fill = var), 
    show.legend = FALSE,
    width = 1
  ) + 
  labs(x = NULL, y = NULL)
var_bar + coord_polar() + theme_minimal()
var_bar + coord_flip() + theme_minimal()

####----Prediction---train------
pred_train=predict(model.rand.1,train)
head(pred_train)## train data
head(train$r)### response. The values for the predicted response and train response is close

####----Prediction---test------
pred_test=predict(model.rand.1,test)
head(pred_test)## train data
head(test$r)### response. The values for the predicted response and train response is close

####----Error--Rate----#########
plot(model.rand.1)

###---Tunning--the--random--forest--model------
tn=tuneRF(train[,-1],train[,1],stepFactor=0.5,plot = TRUE,ntreeTry = 300,
       trace = TRUE,improve = 0.03)


model.rand.2 <- randomForest::randomForest(r~ ., data = train,
                    mtry=4,ntree=300,proximity=TRUE,importance=TRUE)

print(model.rand.2)

####----Prediction---train------
pred_train.2=predict(model.rand.2,train)
head(pred_train.2)## train data
head(train$r)### response. The values for the predicted response and train response is close

####----Prediction---test------
pred_test.2=predict(model.rand.2,test)
head(pred_test.2)## train data
head(test$r)### response. The values for the predicted response and train response is close

###--Number--of--nodes--for--the--trees
hist(treesize(model.rand.2),main = "No. of nodes for the trees",col="red")

###---Variable---importance----
varImpPlot(model.rand.2, sort = T,n.var = 5,main = "Top five predictors")
importance(model.rand.2)
varUsed(model.rand.2)

####---Partial--dependency--plot---
#Gives graphical depiction of the marginal effect of a variable on the response (regression)
partialPlot(model.rand.2,train,FiedlerValue,"1")

#####---Extract--single--tree
getTree(model.rand.2,1,labelVar = T)

####--Multi--dimensional--scaling--plot--of--proximity--matrix---
#MDSplot(model.rand.2,train$r)----#for class mainly

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
# ```{r}
# #| echo: true
# #| include: true
# #| code-fold: true
# 
# ### Tuning the number of mtrys
# tn=tuneRF(train[,-1],train$R,stepFactor=0.5,plot = TRUE,ntreeTry = 500,
#           trace = TRUE,improve = 0.03) #works only if mtry
# #> the number of variables(features). This is because mtry is the number of randomly sampled variable as candidate at each split. Base case use mtry=2  when this happens
# 
# set.seed(6748)
# rfmodel <- randomForest::randomForest(R~ ., data = train,
#                                       trainControl=train.control,
#                                       mtry=4,ntree=500,proximity=TRUE,importance=TRUE)
# 
# print(rfmodel)
# 
# ####----Prediction---train------
# pred.train.rf=predict(rfmodel,train)
# # head(pred.train.rf)## train data
# # head(train$R)### response. The values for the predicted response and train response is close
# 
# ####----Prediction---test------
# pred.test.rf=predict(rfmodel,test)
# # head(pred.test.rf)## train data
# # head(test$R)### response. The values for the predicted response and train response is close
# s5=rbind(head(test$R),c(head(pred.test.rf)))
# row.names(s5)=c("ActualVals","PredictedVals")
# s5
# 
# ####################
# obs.rf=test$R
# pred.rf=pred.test.rf
# dist.R.rf=data.frame(y=obs.rf,x=pred.rf)
# 
# #################################
# ggplot(dist.R.rf,aes(x=x,y=y))+
#   xlab("predicted")+
#   ylab("actual")+
#   geom_point()+
#   geom_abline(color="darkblue")+
#   ggtitle("Plot of actual vs predicted radius for rfmodel")
# 
# ####--Test--error-MSE--test--###
# test.error.rf= pred.test.rf-test$R 
# ####--RMSE---for--test--set
# RMSE_test.rf=sqrt(mean((test.error.rf)^2))
# RMSE_test.rf
# caret::RMSE(test$R,pred.test.rf)
# 
# ####--RSQUARED---for--test--set
# actual.rf=test$R
# R_squared = 1 - (sum((actual.rf-pred.test.rf)^2))/(sum((actual.rf-mean(actual.rf))^2))
# R_squared
# ####----Error--Rate----#########
# plot(rfmodel)
# 
# 
# ###--Number--of--nodes--for--the--trees
# #hist(treesize(rfmodel2),main = "No. of nodes for the trees",col="red")
# 
# ###---Variable---importance----
# varImpPlot(rfmodel, sort = T,n.var = 10,main = "Top five predictors")
# importance(rfmodel)
# varUsed(rfmodel)
# 
# ####---Partial--dependency--plot---
# #Gives graphical depiction of the marginal effect of a variable on the response (regression)
# partialPlot(rfmodel2,train,FiedlerValue,"1")
# 
# #####---Extract--single--tree
# #getTree(rfmodel2,1,labelVar = T)
# 
# #Conditional=True, adjusts for correlations between predictors.
# var_scores <- caret::varImp(rfmodel, conditional=TRUE,scale=TRUE)
# var_scores
# #Gathering rownames in 'var'  and converting it to the factor
# #to provide 'fill' parameter for the bar chart. 
# var_scores <- var_scores %>% tibble::rownames_to_column("var") 
# var_scores$var<- var_scores$var %>% as.factor()
# 
# #Plotting the bar and polar charts for comparing variables
# var_bar <- ggplot(data = var_scores) + 
#   geom_bar(
#     stat = "identity",#it leaves the data without count and bin
#     mapping = aes(x = var, y=Overall, fill = var), 
#     show.legend = FALSE,
#     width = 1
#   ) + 
#   labs(x = NULL, y = NULL)
# var_bar + coord_polar() + theme_minimal()
# var_bar + coord_flip() + theme_minimal()
# ```
