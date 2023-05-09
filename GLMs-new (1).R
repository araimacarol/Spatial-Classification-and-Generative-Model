

setwd("C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi")
a1=read.csv("x1.csv")
a2=read.csv("x2.csv")
a3=read.csv("x3.csv")
a4=read.csv("x4.csv")
a5=read.csv("x5.csv")
a6=read.csv("x6.csv")
a7=read.csv("x7.csv")
a8=read.csv("x8.csv")

d=rbind(a1,a2,a3,a4,a5,a6,a7,a8)


data= d%>% dplyr::select(-c(GraphID,order, GraphReplicate, order,max_component,minDegree,
                            threshold,maxDegree,GraphName,connected,X))


##--Shuffle--data
df<- data[sample(1:nrow(data)), ]##shuffle row indices and randomly re order 
head(df)

####-Split-to-Train-and-Test
set.seed(123)
Train_and_Test <- function(Data, size = 0.8, train = TRUE) {
  n_row = nrow(Data)
  df_split=size*n_row
  train_sample = 1:df_split
  if (train == TRUE) {
    return (Data[train_sample, ])
  } else {
    return (Data[-train_sample, ])
  }
}


train=Train_and_Test(df,0.8,train = T)
test=Train_and_Test(df,0.8,train = F)

dim(train)

###--Lets--Plot--the distribution--of--the--radius-for-the-train-and-test-set
par(mar=c(1,1,1,1))
par(mfrow = c(1, 2))
radius.dist.train=hist(train$R,main = "R-trian-dist",xlim = range(train$R), xlab = "radius")
radius.dist.test=hist(test$R,main = "R-test-dist",xlim = range(train$R),xlab = "radius")


###--Check--mean-and-variance--of--both--learning--set
Mean=rbind(mean(train$R),mean(test$R))
Variance=rbind(var(train$R),var(test$R))
x=data.frame(Mean,Variance);row.names(x)=c("train","test")

set.seed(9347)

###++++++++++++++--Simple GLM--+++++++++++++++
lm.model.1 <- glm(
  R~ .,
  data = train,
  family = gaussian()
)

summary(lm.model.1)

#--Make--predictions
pred.lm1 <- lm.model.1 %>% predict(test)

s1=rbind(head(test$R),head(pred.lm1))
row.names(s1)=c("ActualVals","PredictedVals")
s1

###++++++++++++++--Caret GLM Model--+++++++++++++++
train.control <- trainControl(method = "repeatedcv", number = 10)

lm.model.2 <- caret::train(
  form = R ~ ., 
  data = train, 
  trControl = train.control,
  method = "glm", 
  family = gaussian()
)

summary(lm.model.2$finalModel)
#--Make--predictions
pred.lm2 <- lm.model.2 %>% predict(test)

s2=rbind(head(test$R),head(pred.lm2))
row.names(s2)=c("ActualVals","PredictedVals")
s2

## R-squared
caret::R2(pred = pred.lm2, obs = test$R)
## RMSE
caret::RMSE(pred = pred.lm2, obs = test$R)

####----Penalization--of-features---####
# We use a regularized GLM if your data is  wide, sparse, collinear or big data
library(glmnet)
#[1]--Ridge--Regression--(L2 regularization)
grid = 10^seq(10, -2, length = 100)

x.train=train[,-1]
x.test=as.matrix(test[,-1])


y.train= train$R

y.test=test$R

ridge.mod = glmnet(x.train, y.train, alpha=0, lambda = grid, thresh = 1e-12)
plot(ridge.mod)# Draw plot of coefficients

ridge.pred = predict(ridge.mod, s = 4, newx = x.test)

##MSE
mean((ridge.pred - y.test)^2)

s3=rbind(head(y.test),c(head(ridge.pred)))
row.names(s3)=c("ActualVals","PredictedVals")
s3

###--Tunning the lambda parameter with cv.glmnet built in cross validation
set.seed(1)
cv.out = cv.glmnet(as.matrix(x.train), y.train, alpha = 0) # Fit ridge regression model on training data
plot(cv.out) # Draw plot of training MSE as a function of lambda

bestlam = cv.out$lambda.min  # Select lamda that minimizes training MSE
bestlam

ridge.pred.bestLambda = predict(ridge.mod, s = bestlam, newx = x.test) # Use best lambda to predict test data

#MSE
mean((ridge.pred.bestLambda - y.test)^2) # Calculate test MSE

s4=rbind(head(y.test),c(head(ridge.pred.bestLambda)))
row.names(s4)=c("ActualVals","PredictedVals")
s4

## R-squared
caret::R2(pred =ridge.pred.bestLambda, obs = y.test)
## RMSE
caret::RMSE(pred = ridge.pred.bestLambda, obs = y.test)

##--Fit Final Ridge model on full data set to estimate the coefficient of variables
x=df[,-1]
y=df$R
FinalRidgeModel = glmnet(x, y, alpha = 0) # Fit ridge regression model on full dataset
FinalRidgeModelPred=predict(FinalRidgeModel, type = "coefficients", s = bestlam)[1:13,] # Display coefficients using lambda chosen by CV


#--MODEL1:"ridge.pred.bestLambda"===important

#[2]--Lasso (Least absolut shrinkage and selection operator) Regression--(L2 regularization)

# Lasso (least absolute shrinkage and selection operator) 
# (also Lasso or LASSO) is a regression analysis method that 
# performs both variable selection and regularization in order to enhance 
# the prediction accuracy and interpretability of the statistical model it produces
# The use of alpha=1 gives a Lasso model
set.seed(1)
lasso.mod = glmnet(x.train, y.train, alpha = 1, lambda = grid) # Fit lasso model on training data
plot(lasso.mod)

cv.out = cv.glmnet(as.matrix(x.train), y.train, alpha = 1) # Fit lasso model on training data
plot(cv.out) # Draw plot of training MSE as a function of lambda

bestlam = cv.out$lambda.min # Select lamda that minimizes training MSE
lasso.pred.bestLambda = predict(lasso.mod, s = bestlam, newx = x.test) # Use best lambda to predict test data


s5=rbind(head(y.test),c(head(lasso.pred.bestLambda)))
row.names(s5)=c("ActualVals","PredictedVals")
s5

##--Plot--of--preficted--vs--actual
plot(y=y.test,x=lasso.pred.bestLambda,type="p",xlab="predicted",ylab="actual")
## MSE
mean((lasso.pred - y.test)^2) # Calculate test MSE

## R-squared
caret::R2(pred =lasso.pred, obs = y.test)
## RMSE
caret::RMSE(pred = lasso.pred, obs = y.test)

##Fit Final Lasso model on full data set to estimate the coefficient of variables
##--Shows you which features are important to keep
out = glmnet(x, y, alpha = 1, lambda = grid) # Fit lasso model on full dataset
lasso_coef = predict(out, type = "coefficients", s = bestlam)[1:12,] # Display coefficients using lambda chosen by CV
lasso_coef

lasso_coef[lasso_coef!=0] # Display only non-zero coefficients

# Unless you want to know what features are relevant, 
# or there is a cost associated in collecting all of the attributes,
# rather than just some of them, then don't perform feature selection at all, 
# and just use regularisation (e.g. ridge-regression) to prevent over-fitting instead.
# This is essentially the advice given by Miller in his 
# monograph on subset selection in regression.
# Feature selection is tricky and often makes predictive performance worse, 
# rather than better. 
# The reason is that it is easy to over-fit the feature selection criterion, 
# as there is esentially one degree of freedom for each attribute, so you get
# a subset of attributes that works really well on one particular sample of data, 
# but not necessarily on any other. Regularisation is easier, as you generally only
# have one degree of freedom, hence you tend to get less over-fitting (although the problem doesn't 
# #go away completely).


# Summarize the model
# summary(model)
# # Make predictions
# predictions <- model %>% predict(test.data)
# # Model performance
# # (a) Prediction error, RMSE
# RMSE(predictions, test.data$sales)
# # (b) R-square
# R2(predictions, test.data$sales)




#| echo: true
#| include: true
#| code-fold: true

# train.scale=scale(train[,-1])
# train.scale=cbind(train[1],train.scale) 
# test.scale=scale(test[,-1])
# test.scale=cbind(test[1],test.scale)

set.seed(9347)
###[1]++++++++++++++--Simple GLM--+++++++++++++++
lm.model.1 <- glm(
  R~ .,
  data = train.scale,
  family = gaussian()
)

summary(lm.model.1)

#--Make--predictions
pred.lm1 <- lm.model.1 %>% predict(test.scale)

s1=rbind(head(test$R),head(pred.lm1))
row.names(s1)=c("ActualVals","PredictedVals")
s1

###[2]+++++++++++++--Caret GLM Model--+++++++++++++++
train.control <- trainControl(method = "cv", number = 10,savePredictions = "all")

set.seed(753)
glm.model.caret <- caret::train(
  form = R ~ ., 
  data = train, 
  preProcess=c("center","scale"),
  trControl = train.control,
  method = "glm", 
  family = gaussian()
)

summary(glm.model.caret$finalModel)
#--Make--predictions
pred.glm.caret <- glm.model.caret %>% predict(test)

s2=rbind(head(test$R),head(pred.glm.caret))
row.names(s2)=c("ActualVals","PredictedVals")
s2

## R-squared
caret::R2(pred = pred.glm.caret, obs = test$R)
## RMSE
caret::RMSE(pred = pred.glm.caret, obs = test$R)

###################
obs.glm=test$R
pred.glm=pred.glm.caret
dist.R.glm=data.frame(y=obs.glm,x=pred.glm)

#################################
ggplot(dist.R.glm,aes(x=x,y=y))+
  xlab("predicted")+
  ylab("actual")+
  geom_point()+
  geom_abline(color="darkblue")+
  ggtitle("Plot of actual vs predicted radius for glm")

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

####--------------------Penalization--of-features---####
# We use a regularized GLM if your data is  wide, sparse, collinear or big data
# Ridge and lasso regression deals with matrix inputs so we will convert our test and train data as 
#such

## Train and Test split as matrix
#df=as.data.frame(df)
# split.data=createDataPartition(df$R,p=0.8,list=FALSE,times=1)
#train.df= scale(train[,-1])


## K fold cross validation (10 fold cross validation as training method)
train.control <- trainControl(method = "cv", number = 10,savePredictions = "all")
## List of potential Lambda values
lambda.vect=10^seq(10, -2, length = 100)

x.train=as.matrix(train[,-1])#train and test sets excluding the response
x.test=as.matrix(test[,-1])

y.train= train$R #train and test sets of the response
y.test=test$R

#[3]--Ridge--Regression--(L2 regularization)
set.seed(9234)
ridge.mod = glmnet(x.train, y.train, alpha = 0, lambda = lambda.vect,thresh = 1e-12)

ridge.pred = predict(ridge.mod, s = 4, newx = x.test)

mean((ridge.pred - y.test)^2)

mean((mean(y.train) - y.test)^2)
# ## Best optimal tuning parameters
# ridge.mod$bestTune
# ## Best lambda
# bestRidgeLambda=ridge.mod$bestTune$lambda

par(mar=c(2,2,2,2))
plot(ridge.mod)# Draw plot of coefficients

###--Tunning the lambda parameter with cv.glmnet built in cross validation
set.seed(1)
cv.out.ridge = cv.glmnet(x.train, y.train, alpha = 0) # Fit ridge regression model on training data
plot(cv.out.ridge) # Draw plot of training MSE as a function of lambda

bestlambda.ridge = cv.out.ridge$lambda.min  # Select lamda that minimizes training MSE
bestlambda.ridge

ridge.pred.bestLambda = predict(ridge.mod, s = bestlambda.ridge, newx = x.test) # Use best lambda to predict test data
ridge.pred.bestLambda
##--Plot--of--preficted--vs--actual for ridge model
obs.ridge=y.test
pred.ridge=ridge.pred.bestLambda
dist.R.ridge=data.frame(y=obs.ridge,x=pred.ridge)

ggplot(dist.R.ridge,aes(x=s1,y=y))+
  xlab("predicted")+
  ylab("actual")+
  geom_point()+
  geom_abline(color="darkblue")+
  ggtitle("Plot of actual vs predicted radius for ridge")

#MSE of Tuned Ridge Model
mean((ridge.pred.bestLambda - y.test)^2) # Calculate test MSE

s4=rbind(head(y.test),c(head(ridge.pred.bestLambda)))
row.names(s4)=c("ActualVals","PredictedVals")
s4
## R-squared of Tuned Ridge Model
caret::R2(pred =ridge.pred.bestLambda, obs = y.test)
## RMSE of Tuned Ridge Model
caret::RMSE(pred = ridge.pred.bestLambda, obs = y.test)

##--Fit Final Ridge model on full data set to estimate the coefficient of variables
x.data=as.matrix(df[,-1])
y.data=c(df$R)
FinalRidgeModel = glmnet(x.data, y.data, alpha = 0) # Fit ridge regression model on full dataset
FinalRidgeModelPred=predict(FinalRidgeModel, type = "coefficients", s = bestlambda.ridge)[1:13,] # Display coefficients using lambda chosen by CV

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###[4]--Lasso (Least absolut shrinkage and selection operator) Regression--(L2 regularization)

# Lasso (least absolute shrinkage and selection operator) 
# (also Lasso or LASSO) is a regression analysis method that 
# performs both variable selection and regularization in order to enhance 
# the prediction accuracy and interpretability of the statistical model it produces
# The use of alpha=1 gives a Lasso model
set.seed(1839)
lasso.mod = glmnet(x.train, y.train, alpha = 1, lambda = lambda.vect)

#glmnet(x.train, y.train, alpha = 1, lambda = grid) # Fit lasso model on training data
plot(lasso.mod)

cv.out.lasso = cv.glmnet(x.train, y.train, alpha = 1) # Fit lasso model on training data
plot(cv.out.lasso) # Draw plot of training MSE as a function of lambda

bestlambda.lasso = cv.out.lasso$lambda.min # Select lamda that minimizes training MSE
bestlambda.lasso
lasso.pred.bestLambda = predict(lasso.mod, s = bestlambda.lasso, newx = x.test) # Use best lambda to predict test data

s5=rbind(head(y.test),c(head(lasso.pred.bestLambda)))
row.names(s5)=c("ActualVals","PredictedVals")
s5

##--Plot--of--predicted--vs--actual
obs.laso=y.test
pred.laso=lasso.pred.bestLambda
dist.R.lasso=data.frame(y=obs.laso,x=pred.laso)

ggplot(dist.R.lasso,aes(x=s1,y=y))+
  xlab("predicted")+
  ylab("actual")+
  geom_point()+
  geom_abline(color="darkblue")+
  ggtitle("Plot of actual vs predicted radius for lasso")

## MSE
mean((lasso.pred.bestLambda - y.test)^2) # Calculate test MSE

## R-squared
caret::R2(pred =lasso.pred.bestLambda, obs = y.test)
## RMSE
caret::RMSE(pred = lasso.pred.bestLambda, obs = y.test)

##Fit Final Lasso model on full data set to estimate the coefficient of variables
##--Shows you which features are important to keep
out.lasso = glmnet(x, y, alpha = 1, lambda = grid) # Fit lasso model on full dataset
lasso_coef = predict(out.lasso, type = "coefficients", s = bestlam.lasso)[1:12,] # Display coefficients using lambda chosen by CV
lasso_coef

lasso_coef[lasso_coef!=0] # Display only non-zero coefficients

#caret::varImp(lasso.mod)

######Function to get variable importance for lasso
varImp <- function(object, lambda = NULL, ...) {
  
  ## skipping a few lines
  
  beta <- predict(object, s = lambda, type = "coef")
  if(is.list(beta)) {
    out <- do.call("cbind", lapply(beta, function(x) x[,1]))
    out <- as.data.frame(out, stringsAsFactors = TRUE)
  } else out <- data.frame(Overall = beta[,1])
  out <- abs(out[rownames(out) != "(Intercept)",,drop = FALSE])
  out
}

#varImp(lasso.mod, lambda = lasso.mod$lambda.min)

# Unless you want to know what features are relevant, 
# or there is a cost associated in collecting all of the attributes,
# rather than just some of them, then don't perform feature selection at all, 
# and just use regularisation (e.g. ridge-regression) to prevent over-fitting instead.
# This is essentially the advice given by Miller in his 
# monograph on subset selection in regression.
# Feature selection is tricky and often makes predictive performance worse, 
# rather than better. 
# The reason is that it is easy to over-fit the feature selection criterion, 
# as there is esentially one degree of freedom for each attribute, so you get
# a subset of attributes that works really well on one particular sample of data, 
# but not necessarily on any other. Regularisation is easier, as you generally only
# have one degree of freedom, hence you tend to get less over-fitting (although the problem doesn't 
# #go away completely).

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ```{r}
# #| echo: true
# #| include: true
# #| code-fold: true
# 
# ####--------------------Penalization--of-features---####
# # We use a regularized GLM if your data is  wide, sparse, collinear or big data
# # Ridge and lasso regression deals with matrix inputs so we will convert our test and train data as 
# #such
# 
# ## Train and Test split as matrix
# #df=as.data.frame(df)
# # split.data=createDataPartition(df$R,p=0.8,list=FALSE,times=1)
# #train.df= scale(train[,-1])
# 
# 
# ## K fold cross validation (10 fold cross validation as training method)
# train.control <- trainControl(method = "cv", number = 10,savePredictions = "all")
# ## List of potential Lambda values
# lambda.vect=10^seq(10, -2, length = 500)
# 
# x.train=as.matrix(scale(train[,-1]))#train and test sets excluding the response
# x.test=as.matrix(scale(test[,-1]))
# 
# y.train= train$R #train and test sets of the response
# y.test=test$R
# 
# #[3]--Ridge--Regression--(L2 regularization)
# set.seed(9234)
# ridge.mod = glmnet(x.train, y.train, alpha = 0, lambda = lambda.vect,thresh = 1e-12,
#                    trControl=train.control)
# 
# ridge.pred = predict(ridge.mod, s = 4, newx = x.test)
# 
# mean((ridge.pred - y.test)^2)
# 
# mean((mean(y.train) - y.test)^2)
# # ## Best optimal tuning parameters
# # ridge.mod$bestTune
# # ## Best lambda
# # bestRidgeLambda=ridge.mod$bestTune$lambda
# 
# par(mar=c(2,2,2,2))
# plot(ridge.mod)# Draw plot of coefficients
# 
# 
# 
# ###--Tunning the lambda parameter with cv.glmnet built in cross validation
# set.seed(1)
# cv.out.ridge = cv.glmnet(x.train, y.train, alpha = 0) # Fit ridge regression model on training data
# plot(cv.out.ridge) # Draw plot of training MSE as a function of lambda
# 
# bestlambda.ridge = cv.out.ridge$lambda.min  # Select lamda that minimizes training MSE
# bestlambda.ridge
# 
# ridge.pred.bestLambda = predict(ridge.mod, s = bestlambda.ridge, newx = x.test) # Use best lambda to predict test data
# 
# ##--Plot--of--preficted--vs--actual for ridge model
# obs.ridge=y.test
# pred.ridge=ridge.pred.bestLambda
# dist.R.ridge=data.frame(y=obs.ridge,x=pred.ridge)
# 
# ggplot(dist.R.ridge,aes(x=s1,y=y))+
#   xlab("predicted")+
#   ylab("actual")+
#   geom_point()+
#   geom_abline(color="darkblue")+
#   ggtitle("Plot of actual vs predicted radius for ridge")
# 
# #MSE of Tuned Ridge Model
# mean((ridge.pred.bestLambda - y.test)^2) # Calculate test MSE
# 
# s4=rbind(head(y.test),c(head(ridge.pred.bestLambda)))
# row.names(s4)=c("ActualVals","PredictedVals")
# s4
# ## R-squared of Tuned Ridge Model
# caret::R2(pred =ridge.pred.bestLambda, obs = y.test)
# ## RMSE of Tuned Ridge Model
# caret::RMSE(pred = ridge.pred.bestLambda, obs = y.test)
# 
# ##--Fit Final Ridge model on full data set to estimate the coefficient of variables
# x=df[,-1]
# y=df$R
# FinalRidgeModel = glmnet(x, y, alpha = 0) # Fit ridge regression model on full dataset
# FinalRidgeModelPred=predict(FinalRidgeModel, type = "coefficients", s = bestlambda.ridge)[1:13,] # Display coefficients using lambda chosen by CV
# 
# ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ###[4]--Lasso (Least absolut shrinkage and selection operator) Regression--(L2 regularization)
# 
# # Lasso (least absolute shrinkage and selection operator) 
# # (also Lasso or LASSO) is a regression analysis method that 
# # performs both variable selection and regularization in order to enhance 
# # the prediction accuracy and interpretability of the statistical model it produces
# # The use of alpha=1 gives a Lasso model
# set.seed(1839)
# lasso.mod = glmnet(x.train, y.train, alpha = 1, lambda = lambda.vect,thresh = 1e-12,
#                    trControl=train.control)
# 
# #glmnet(x.train, y.train, alpha = 1, lambda = grid) # Fit lasso model on training data
# plot(lasso.mod)
# 
# cv.out.lasso = cv.glmnet(x.train, y.train, alpha = 1) # Fit lasso model on training data
# plot(cv.out.lasso) # Draw plot of training MSE as a function of lambda
# 
# bestlambda.lasso = cv.out.lasso$lambda.min # Select lamda that minimizes training MSE
# bestlambda.lasso
# lasso.pred.bestLambda = predict(lasso.mod, s = bestlambda.lasso, newx = x.test) # Use best lambda to predict test data
# 
# s5=rbind(head(y.test),c(head(lasso.pred.bestLambda)))
# row.names(s5)=c("ActualVals","PredictedVals")
# s5
# 
# ##--Plot--of--predicted--vs--actual
# obs.laso=y.test
# pred.laso=lasso.pred.bestLambda
# dist.R.lasso=data.frame(y=obs.laso,x=pred.laso)
# 
# ggplot(dist.R.lasso,aes(x=s1,y=y))+
#   xlab("predicted")+
#   ylab("actual")+
#   geom_point()+
#   geom_abline(color="darkblue")+
#   ggtitle("Plot of actual vs predicted radius for lasso")
# 
# ## MSE
# mean((lasso.pred.bestLambda - y.test)^2) # Calculate test MSE
# 
# ## R-squared
# caret::R2(pred =lasso.pred.bestLambda, obs = y.test)
# ## RMSE
# caret::RMSE(pred = lasso.pred.bestLambda, obs = y.test)
# 
# ##Fit Final Lasso model on full data set to estimate the coefficient of variables
# ##--Shows you which features are important to keep
# out.lasso = glmnet(x, y, alpha = 1, lambda = grid) # Fit lasso model on full dataset
# lasso_coef = predict(out.lasso, type = "coefficients", s = bestlam.lasso)[1:12,] # Display coefficients using lambda chosen by CV
# lasso_coef
# 
# lasso_coef[lasso_coef!=0] # Display only non-zero coefficients
# 
# caret::varImp(lasso.mod)
# 
# ######Function to get variable importance for lasso
# varImp <- function(object, lambda = NULL, ...) {
#   
#   ## skipping a few lines
#   
#   beta <- predict(object, s = lambda, type = "coef")
#   if(is.list(beta)) {
#     out <- do.call("cbind", lapply(beta, function(x) x[,1]))
#     out <- as.data.frame(out, stringsAsFactors = TRUE)
#   } else out <- data.frame(Overall = beta[,1])
#   out <- abs(out[rownames(out) != "(Intercept)",,drop = FALSE])
#   out
# }

#varImp(lasso.mod, lambda = lasso.mod$lambda.min)

# Unless you want to know what features are relevant, 
# or there is a cost associated in collecting all of the attributes,
# rather than just some of them, then don't perform feature selection at all, 
# and just use regularisation (e.g. ridge-regression) to prevent over-fitting instead.
# This is essentially the advice given by Miller in his 
# monograph on subset selection in regression.
# Feature selection is tricky and often makes predictive performance worse, 
# rather than better. 
# The reason is that it is easy to over-fit the feature selection criterion, 
# as there is essentially one degree of freedom for each attribute, so you get
# a subset of attributes that works really well on one particular sample of data, 
# but not necessarily on any other. Regularisation is easier, as you generally only
# have one degree of freedom, hence you tend to get less over-fitting (although the problem doesn't 
# #go away completely).
#```