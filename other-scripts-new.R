###+++++++++++++This contains some functions, scripts etc I wrote and dont use now++++++
#---(1) Some tools to import and clean data
# data="node117-r1000sample.csv"
# df.data=read.csv(data,header = T, sep=",")
# data<- df.data[sample(1:nrow(df.data)), ]##shuffle row indices and randomly re order 
# #--selecting same node size and droping non-relevant columns
# node.num=117
# df.data<-data %>%dplyr::filter(order==node.num)%>% 
#   dplyr::select(-c(GraphID,order, GraphReplicate, order,max_component,minDegree,
#                    threshold,maxDegree,GraphName,connected,X))
# 
# dplyr::glimpse(df.data)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#--(2) Calculating data variance (whether data has 0 varaiance, thus, if entire column is of the same data set, this could hurt the PCA)
#df.data.var=nearZeroVar(df.data[,-1], saveMetrics = T)#checking variance of predictors
#df_nzv <- df.data[c(rownames(df.data.var[df.data.var$percentUnique > 0.1,])) ]
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#---(3)Principal Component Analysis (PCA)
#--find correlations between variables

# corel.data=cor(df.data[,-1])
# corel.data

#--PCA-1
#set.seed(26836)
# PCA_model.1=princomp(df.data[,-1],cor = T)## select all features except the response
# #and find Correlations between them
# sum_pca=summary(PCA_model.1, loadings=T)
# sum_pca
#df.pca=write.csv(sum_pca$loadings,"PCA.csv")

#set.seed(26836)
##--PCA-2
#PCA_model.2=prcomp(df.data[,-1], scale=T)## select all features except the response and scale them

#summary(PCA_model.2)## variation in the data is explained mostly by only PC1 and a little by PC2

#plot(PCA_model.2, type="l")

###--PCA-2 pot

#biplot(PCA_model.2,scale=0)

#str(PCA_model.2)

#PCA_model.2$x

##binding data and PC1 and PC2
#df.new=cbind(df.data,PCA_model.2$x[,1:2])## 

###Plots of PCs
#ggplot(df.new, aes(PC1,PC2, col=R,fill=R))+
#  stat_ellipse(geom="polygon",col="black",alpha=0.5)+
#  geom_point(shape=21,col="black")

### correlation matrix of PC1 and PC2
#cor_matrix=cor(df.new[,-1],df.new[,14:15])

# In summary: A PCA biplot shows both PC scores of samples (dots) 
# and loadings of variables (vectors). The further away these vectors 
# are from a PC origin, the more influence they have on that PC. 
# Loading plots also hint at how variables correlate with one another: 
# a small angle implies positive correlation, a large one suggests negative 
# correlation, and a 90° angle indicates no correlation between two 
# characteristics. A scree plot displays how much variation each principal
# component captures from the data. If the first two or three PCs are
# sufficient to describe the essence of the data, the scree plot is a
# steep curve that bends quickly and flattens out.


# Train data
# x2 = iris[-s,-5]
# 
# #Standardize train data
# x2 = scale(x2)
# 
# #Principal components of train data
# pr2 = prcomp(x2, center = FALSE, scale = FALSE)
# 
# #Test data
# y2 = iris[s,-5]
# 
# #Standardize test data with the same parameters used on train data
# y2 = scale(y2, center = attr(x2,"scaled:center"), scale = attr(x2,"scaled:scale"))
# 
# #Rotate standardized test data to the same space as train data
# #You can also keep the first K columns in case you want to retain a X% ammount of variance
# y2 = y2 %*% pr$rotation
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#(3) Feature-selection
###--Feature--reduction---####

#--Whcih features are correlated
#df.corr=cor(df.data[,-1])

# df.data=df.data%>%dplyr::select(-c(modularity,diameter, betweenness, transitivity, centrality_eigen,Normalized_FiedlerValue,minCut)) 
# df.data=df.data%>%dplyr::select(c(R,FiedlerValue,Normalized_FiedlerValue,spectral_radius,transitivity,modularity))
# 
# #df.data=df.data%>%dplyr::select(R,edges)
# 
# head(df.data)
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#(4) Linear Regression



# train.scale=scale(train[,-1])
# train.scale=cbind(train[1],train.scale) 
# test.scale=scale(test[,-1])
# test.scale=cbind(test[1],test.scale)

# set.seed(9347)
# ###++++++++++++++--Simple GLM--+++++++++++++++
# lm.model.1 <- glm(
#   R~ .,
#   data = train.scale,
#   family = gaussian()
# )
# 
# summary(lm.model.1)
# 
# #--Make--predictions
# pred.lm1 <- lm.model.1 %>% predict(test.scale)
# 
# s1=rbind(head(test$R),head(pred.lm1))
# row.names(s1)=c("ActualVals","PredictedVals")
# s1
# 
# ###+++++++++++++--Caret GLM Model--+++++++++++++++
# train.control <- trainControl(method = "cv", number = 10,savePredictions = "all")
# 
# set.seed(753)
# glm.model.caret <- caret::train(
#   form = R ~ ., 
#   data = train, 
#   preProcess=c("center","scale"),
#   trControl = train.control,
#   method = "glm", 
#   family = gaussian()
# )
# 
# summary(glm.model.caret$finalModel)
# #--Make--predictions
# pred.glm.caret <- glm.model.caret %>% predict(test)
# 
# s2=rbind(head(test$R),head(pred.glm.caret))
# row.names(s2)=c("ActualVals","PredictedVals")
# s2
# 
# ## R-squared
# caret::R2(pred = pred.glm.caret, obs = test$R)
# ## RMSE
# caret::RMSE(pred = pred.glm.caret, obs = test$R)
# 
# ###################
# obs.glm=test$R
# pred.glm=pred.glm.caret
# dist.R.glm=data.frame(y=obs.glm,x=pred.glm)
# 
# #################################
# ggplot(dist.R.glm,aes(x=x,y=y))+
#   xlab("predicted")+
#   ylab("actual")+
#   geom_point()+
#   geom_abline(color="darkblue")+
#   ggtitle("Plot of actual vs predicted radius for glm")

# 
# # Feature Scaling
# # training_set = scale(training_set)
# # test_set = scale(test_set)
# 
# # Fitting Multiple Linear Regression to the Training set
# #++++++++--Linear Model 1
# lm.regressor1 = lm(formula = R ~ .,
#                    data = train)
# 
# print(lm.regressor1)
# 
# # Predicting the Test set results
# y_pred = predict(lm.regressor1 , newdata = test)
# 
# head(y_pred)
# head(test$R)
# 
# #Seems like we have correlated features so we do
# ##--Step wise feature selection
# ##--Stepwise regression model (select best predictors)
# step.model.1 <- stepAIC(lm.regressor1, direction = "both", 
#                         trace = FALSE)
# summary(step.model.1)
# 
# #stepAIC() [MASS package], which choose the best model by AIC. It has an option named direction, which can take the following values: i) "both" (for stepwise regression, both forward and backward selection); "backward" (for backward selection) and "forward" (for forward selection). It return the best final model.
# 
# #-Best model with stepAIC
# # lm(formula = R ~ edges + minCut + FiedlerValue + closeness + 
# #     modularity + diameter + betweenness + centrality_eigen, data = train)
# 
# #+++++++-----Linear-Model-2-with leaps packages (select best predictors)
# lm.regressor2 <- regsubsets(R~., data = train, nvmax = 5,
#                             method = "seqrep")
# summary(lm.regressor2)
# 
# # 
# # Note that, the train() function [caret package] provides an easy workflow to perform stepwise selections using the leaps and the MASS packages. It has an option named method, which can take the following values:
# # 
# # "leapBackward", to fit linear regression with backward selection
# # "leapForward", to fit linear regression with forward selection
# # "leapSeq", to fit linear regression with stepwise selection .
# 
# # You also need to specify the tuning parameter nvmax, which corresponds to the maximum number of predictors to be incorporated in the model
# 
# #+++++++-----Linear-Model-3-with--train--package
# set.seed(123)
# # Set up repeated k-fold cross-validation
# train.control <- trainControl(method = "repeatedcv", number = 10)
# # Train the model
# train.data=train
# lm.regressor3 <- train(R ~., data = train.data,
#                        method = "leapBackward", 
#                        tuneGrid = data.frame(nvmax = 1:12),
#                        trControl = train.control
# )
# lm.regressor3$results
# 
# # nvmax: the number of variable in the model. For example nvmax = 2, specify the best 2-variables model
# # RMSE and MAE are two different metrics measuring the prediction error of each model. The lower the RMSE and MAE, the better the model.
# # Rsquared indicates the correlation between the observed outcome values and the values predicted by the model. The higher the R squared, the better the model.
# 
# lm.regressor3$bestTune
# # 
# summary(lm.regressor3$finalModel)
# # An asterisk specifies that a given variable is included in the corresponding model. For example, it can be seen that the best 4-variables model contains Agriculture, Education, Catholic, Infant.Mortality (Fertility ~ Agriculture + Education + Catholic + Infant.Mortality).
# # The regression coefficients of the final model (id = 4) can be accessed as follow:
# 
# coef(lm.regressor3$finalModel, 10)
# 
# # Additionally, the caret package has method to compute stepwise regression using the MASS package (method = "lmStepAIC"):
# 
# ###+++++++++++++FINAL--MODELS----
# # Train the model
# 
# lm.regressor.final.2 <- lm(R ~., data = train)
# step <- stepAIC(lm.regressor.final.2, direction = "both", trace = FALSE)
# step
# 
# 
# lm.regressor.final.2 <- train(R ~., data = train.data,
#                               method = "lmStepAIC", 
#                               trControl = train.control,
#                               trace = FALSE
# )
# 
# # Best model from StepAIC and train with caret package from stepwise regression
# # lm(formula = R ~ edges + minCut + FiedlerValue + closeness + 
# #     modularity + diameter + betweenness + centrality_eigen, data = train)
# 
# # Model accuracy
# lm.regressor.final.2$results
# # Final model coefficients
# lm.regressor.final.2$finalModel
# # Summary of the model
# summary(lm.regressor.final.2$finalModel)
# 
# 
# pred=predict(lm.regressor.final.2,test)
# head(test$R)
# head(pred)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#(5) GLMs
#| echo: true
#| include: true
#| code-fold: true

# train.scale=scale(train[,-1])
# train.scale=cbind(train[1],train.scale) 
# test.scale=scale(test[,-1])
# test.scale=cbind(test[1],test.scale)

# set.seed(9347)
# ###[1]++++++++++++++--Simple GLM--+++++++++++++++
# lm.model.1 <- glm(
#   R~ .,
#   data = train.scale,
#   family = gaussian()
# )
# 
# summary(lm.model.1)
# 
# #--Make--predictions
# pred.lm1 <- lm.model.1 %>% predict(test.scale)
# 
# s1=rbind(head(test$R),head(pred.lm1))
# row.names(s1)=c("ActualVals","PredictedVals")
# s1
# 
# ###[2]+++++++++++++--Caret GLM Model--+++++++++++++++
# train.control <- trainControl(method = "cv", number = 10,savePredictions = "all")
# 
# set.seed(753)
# glm.model.caret <- caret::train(
#   form = R ~ ., 
#   data = train, 
#   preProcess=c("center","scale"),
#   trControl = train.control,
#   method = "glm", 
#   family = gaussian()
# )
# 
# summary(glm.model.caret$finalModel)
# #--Make--predictions
# pred.glm.caret <- glm.model.caret %>% predict(test)
# 
# s2=rbind(head(test$R),head(pred.glm.caret))
# row.names(s2)=c("ActualVals","PredictedVals")
# s2
# 
# ## R-squared
# caret::R2(pred = pred.glm.caret, obs = test$R)
# ## RMSE
# caret::RMSE(pred = pred.glm.caret, obs = test$R)
# 
# ###################
# obs.glm=test$R
# pred.glm=pred.glm.caret
# dist.R.glm=data.frame(y=obs.glm,x=pred.glm)
# 
# #################################
# ggplot(dist.R.glm,aes(x=x,y=y))+
#   xlab("predicted")+
#   ylab("actual")+
#   geom_point()+
#   geom_abline(color="darkblue")+
#   ggtitle("Plot of actual vs predicted radius for glm")
# 
# ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

####--------------------Penalization--of-features---####
# We use a regularized GLM if your data is  wide, sparse, collinear or big data
# Ridge and lasso regression deals with matrix inputs so we will convert our test and train data as 
#such

## Train and Test split as matrix
#df=as.data.frame(df)
# split.data=createDataPartition(df$R,p=0.8,list=FALSE,times=1)
#train.df= scale(train[,-1])


## K fold cross validation (10 fold cross validation as training method)
# train.control <- trainControl(method = "cv", number = 10,savePredictions = "all")
# ## List of potential Lambda values
# lambda.vect=10^seq(10, -2, length = 100)
# 
# x.train=as.matrix(train[,-1])#train and test sets excluding the response
# x.test=as.matrix(test[,-1])
# 
# y.train= train$R #train and test sets of the response
# y.test=test$R
# 
# #[3]--Ridge--Regression--(L2 regularization)
# set.seed(9234)
# ridge.mod = glmnet(x.train, y.train, alpha = 0, lambda = lambda.vect,thresh = 1e-12)
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
# ###--Tunning the lambda parameter with cv.glmnet built in cross validation
# set.seed(1)
# cv.out.ridge = cv.glmnet(x.train, y.train, alpha = 0) # Fit ridge regression model on training data
# plot(cv.out.ridge) # Draw plot of training MSE as a function of lambda
# 
# bestlambda.ridge = cv.out.ridge$lambda.min  # Select lamda that minimizes training MSE
# bestlambda.ridge
# 
# ridge.pred.bestLambda = predict(ridge.mod, s = bestlambda.ridge, newx = x.test) # Use best lambda to predict test data
# ridge.pred.bestLambda
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
# x.data=as.matrix(df[,-1])
# y.data=c(df$R)
# FinalRidgeModel = glmnet(x.data, y.data, alpha = 0) # Fit ridge regression model on full dataset
# FinalRidgeModelPred=predict(FinalRidgeModel, type = "coefficients", s = bestlambda.ridge)[1:13,] # Display coefficients using lambda chosen by CV

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###[4]--Lasso (Least absolut shrinkage and selection operator) Regression--(L2 regularization)

# Lasso (least absolute shrinkage and selection operator) 
# (also Lasso or LASSO) is a regression analysis method that 
# performs both variable selection and regularization in order to enhance 
# the prediction accuracy and interpretability of the statistical model it produces
# The use of alpha=1 gives a Lasso model
# set.seed(1839)
# lasso.mod = glmnet(x.train, y.train, alpha = 1, lambda = lambda.vect)
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
# #caret::varImp(lasso.mod)
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
# as there is esentially one degree of freedom for each attribute, so you get
# a subset of attributes that works really well on one particular sample of data, 
# but not necessarily on any other. Regularisation is easier, as you generally only
# have one degree of freedom, hence you tend to get less over-fitting (although the problem doesn't 
# #go away completely).
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=+++++
##(6)----Random--Forest
# set.seed(45678)
# model.rand.1 <- randomForest::randomForest(R~ ., data = train, importance=TRUE,scale=TRUE)
# print(model.rand.1)
# 
# attributes(model.rand.1)

#--Conditional=True, adjusts for correlations between predictors.
#var_scores <- caret::varImp(model.rand.1, conditional=TRUE,scale=TRUE)
#var_scores

#--Gathering rownames in 'var'  and converting it to the factor
#to provide 'fill' parameter for the bar chart. 

#var_scores <- var_scores %>% tibble::rownames_to_column("var") 
#var_scores$var<- var_scores$var %>% as.factor()

#--Plotting the bar and polar charts for comparing variables
#var_bar <- ggplot(data = var_scores) + 
# geom_bar(
#   stat = "identity",#it leaves the data without count and bin
#  mapping = aes(x = var, y=Overall, fill = var), 
# show.legend = FALSE,
#  width = 1
# ) + 
#  labs(x = NULL, y = NULL)
#var_bar + coord_polar() + theme_minimal()
#var_bar + coord_flip() + theme_minimal()

####----Prediction---train------
# pred_train=predict(model.rand.1,train)
# head(pred_train)## train data
# head(train$R)### response.

####----Prediction---test------
# pred_test=predict(model.rand.1,test)
# head(pred_test)## test data
# head(test$R)### response. The values for the predicted response and train response is close

# ####--Test--error-MSE--test--###
# test_error= pred_test-test$R 
# ####--RMSE---for--test--set
# RMSE_test=sqrt(mean((test_error)^2))
# caret::RMSE(test$R,pred_test)
# ####--RSQUARED---for--test--set
# actual=test$R
# R_squared = 1 - (sum((actual-pred_test)^2))/(sum((actual-mean(actual))^2))
# 
# ####----Error--Rate----#########
# plot(model.rand.1)

#-model.2
# trainCtrl <- trainControl(method = "cv", number = 10, savePredictions = TRUE)
# rF.Model.2 <- train(R ~., method = "rf", trControl = trainCtrl, preProcess = "pca", data = train, prox = TRUE,ntree=500,importance=TRUE)

#train$R = as.factor(train$R) # conver response to factor
###---Tunning--the--random--forest--model------

# ridge.mod = train(R~.,
#                   data=train.df,
#                   preProcess=c("center","scale"),
#                   method="glmnet",
#                   tuneGrid=expand.grid(alpha=1,lambda=lambda.vect),
#                   trControl=train.control,
#                   na.ction=na.omit)

# train.control <- trainControl(method = "cv", number = 10,savePredictions = "all")
# ## List of potential Lambda values
# lambda.vect=10^seq(5, -5, length = 500)

### In built cross validation
# cv.out.ridge = cv.glmnet(x.train, y.train, alpha = 0) # Fit ridge regression model on training data
# plot(cv.out.ridge) # Draw plot of training MSE as a function of lambda
# 
# bestRidgeLambda = cv.out.ridge$lambda.min  # Select lamda that minimizes training MSE
# bestRidgeLambda
# 
# ridge.pred = predict(ridge.mod, s = bestRidgeLambda, newx = x.test)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#--[7]--Data--stuffs
# r1=0.01;r2=0.05;r3=0.1;r4=0.15;r5=0.2;r6=.25;r7=.3;r8=.35
# r9=.4;r10=.45;r11=.5;r12=.55;r13=.6;r14=.65;r15=.7;r16=.75
# r17=.8;r18=.85;r19=.9;r20=.95
#R=c(0.15,0.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,.95)

# r.1=runif(300,0.15,0.25)
# r.2=runif(300,0.25,0.35)
# r.3=runif(300,0.35,0.45)
# r.4=runif(300,0.45,0.55)
# r.5=runif(300,0.55,0.65)
# r.6=runif(300,0.65,0.75)
# r.7=runif(300,0.75,0.85)
# r.8=runif(300,0.85,0.95)
# 
# node.num=100
# #  
# networks.fun=function(R=c(0.3,0.4,0.6,0.05)){
#   net=NULL;data=NULL
#   for (i in 1:length(R)){
#     net[i]=makeSpatialGraphs(node.size =node.num,Radius=R[i])
#     data=RunSimOnGraphFeatures(net,nreps = 1)
#   }
#   df=cbind(R,data)
#   return(df)
# }
# 
# x.1=networks.fun(R=r.1)
# write.csv(x.1,"x1.csv")
# networks.fun(R=r.2)
# x2=write.csv(x.2,"x2.csv")
# networks.fun(R=r.3)
# x3=write.csv(x.3,"x3.csv")
# x.4=networks.fun(R=r.4)
# x4=write.csv(x.4,"x4.csv")
# x.5=networks.fun(R=r.5)
# x5=write.csv(x.5,"x5.csv")
# x.6=networks.fun(R=r.6)
# x6=write.csv(x.6,"x6.csv")
# x.7=networks.fun(R=r.7)
# x7=write.csv(x.7,"x7.csv")
# x.8=networks.fun(R=r.8)
# x8=write.csv(x.8,"x8.csv")
# # 

# a1=read.csv("x1.csv")
# a2=read.csv("x2.csv")
# a3=read.csv("x3.csv")
# a4=read.csv("x4.csv")
# a5=read.csv("x5.csv")
# a6=read.csv("x6.csv")
# a7=read.csv("x7.csv")
# a8=read.csv("x8.csv")
# 
# df.plot=rbind(a1,a2,a3,a4,a5)
# 
# 
# df.plot= df.plot %>% dplyr::select(-c(GraphID,order, GraphReplicate, order,max_component,minDegree,
#                                       threshold,maxDegree,GraphName,connected,X))
# head(df.plot)
# 
# df.plot.long <- gather(df.plot, GraphFeatures,values, edges:centrality_eigen, factor_key=TRUE)

#df.plot$R=as.factor(df.plot$R)

#df.plot.long$R=as.factor(df.plot.long$R)
# head(df.plot.long

#---- Distribution of R
# par(mar=c(1,1,1,1))
# d=density(df.plot$R)
# plot(d,main="Distribution of the radius")
# polygon(d, col="blue", border="green")


#geom_dotplot(binaxis = "y",binwidth = 0.2,stackdir = "center",show.legend = F)

# p2=  df.plot.long %>% ggplot(aes(x=GraphFeatures,y=R))+
#   geom_violin(show.legend = F,alpha=0.5,adjust=0.5,draw_quantiles = c(0.25,0.5,0.75))
#    # geom_dotplot(binaxis = "y",binwidth = 0.5,show.legend = F)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#[8]---Prediction--on--Empirical--networks

# all.nets=function(model=rfmodel2,data="aves-barn-swallow-contact-network.edges",graphname="dolph"){
#   df <- read.table(data)  
#   df_new=graph_from_data_frame(as.matrix(df),directed=FALSE)
#   df.graph=igraph::simplify(df_new,remove.multiple = T,remove.loops = T)
#   components = igraph::clusters(df.graph , mode="weak")
#   biggest_cluster_id = which.max(components$csize)
#   vert_ids = V(df.graph)[components$membership== biggest_cluster_id]
#   graph=igraph::induced_subgraph(df.graph, vert_ids)
#   
#   graph$name=graphname
#   graph$type=graphname
#   graph$id="1"
#   G=list(graph)
#   
#   ###---Graph--Features--###
#   data=RunSimOnGraphFeatures(G,nreps = 1)
#   
#   newdata=data %>%dplyr::select(-c(GraphID,order, GraphReplicate,max_component,minDegree,threshold,maxDegree,GraphName,connected))
#   
#   #newdata=data %>%dplyr::select(c(FiedlerValue,Normalized_FiedlerValue,spectral_radius,transitivity,modularity))
#   
#   #  newdata=data %>%dplyr::select(c(edges,minCut,FiedlerValue,closeness, 
#   # modularity,diameter,betweenness,centrality_eigen))
#   
#   predicted_radius=abs(predict(model,newdata))
#   
#   spatial.net=makeSpatialGraphs(node.size=vcount(df.graph),Radius=predicted_radius)
#   #fastSpatialNetwork(n=vcount(graph),r=predicted_radius,makeConnected=TRUE, keepCellsSeparate=FALSE)
#   
#   ###--Graph--Features--for--empirical--and--spatial--net
#   # spatial.net$name="Spatial"
#   # spatial.net$type="Spatial"
#   # spatial.net$id="1"
#   net=c(spatial.net,G)
#   feat=RunSimOnGraphFeatures(net,nreps = 1)
#   
#   #new.list <- list(spatial.net[[1]],graph,feat)
#   x=c(G,spatial.net)
#   x$summary=feat
#   return(x)
# }
# 
# x.1=all.nets(model=lasso.pred.bestLambda,data="mammalia-hyena-networkb.edges",graphname="hyena")
# x.1$summary
# par(mfrow=c(1,2))
# plot(x.1[[1]]);plot(x.1[[2]])
# 
# x.2=all.nets(model=lm.regressor.final.2,data="mammalia-bat-roosting-indiana.edges",graphname="bat")
# x.2$summary
# par(mfrow=c(1,2))
# plot(x.2[[1]]);plot(x.2[[2]])
# 
# x.3=all.nets(model=lm.regressor.final.2,data="insecta-ant-colony1-day01.edges",graphname="ant")
# x.3$summary
# par(mfrow=c(1,2))
# plot(x.3[[1]]);plot(x.3[[2]])

setwd("C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi")
g.net="aves-weaver-social.edges"
g.data <- read.table(g.net)  
g.new=graph_from_data_frame(as.matrix(g.data),directed=FALSE)
g.graph=igraph::simplify(g.new,remove.multiple = T,remove.loops = T)
components = igraph::clusters(g.graph , mode="weak")
biggest_cluster_id = which.max(components$csize)
vert_ids = V(g.graph)[components$membership== biggest_cluster_id]
g=igraph::induced_subgraph(g.graph, vert_ids)


 simulate.spatial<-function(N=50,pred.rad=0.4,nsim=100){
   spatial.graph=NULL
   for (i in 1:nsim){
     spatial.graph[[i]]=makeSpatialGraphs(node.size=N,Radius=pred.radius)
   }
   return(spatial.graph)
 }

al= SPNET(g,0.169,100)
edgeCount(al)

#### edge counts of graphs
edgeCount=function(G){
   cnt=NULL;store=NULL
   for (i in 1:length(G)){
     store[i]=G[[i]]
     cnt[i]=ecount(store[[i]])
   }
   
 return(cnt)
   }



# set.seed(8765)
# lasso.net=SPNET(df.graph,predicted_radius.lasso,N=100)
# spatial.net.lasso=lasso.net[[39]]
# ridge.net=SPNET(df.graph,predicted_radius.ridge,N=100)
# spatial.net.ridge=ridge.net[[98]]
# rfmodel2.net=SPNET(df.graph,predicted_radius.rfmodel2,N=100)
# spatial.net.rfmodel2=rfmodel2.net[[18]]
# Glm.net=SPNET(df.graph,predicted_radius.glm,N=100)
# spatial.net.glm=Glm.net[[1]]

#edgeCount(Glm.net)

#### edge counts of graphs
# edgeCount=function(G){
#   cnt=NULL;store=NULL#;m=1
#   for (i in 1:length(G)){
#     store[i]=G[[i]]
#     cnt[i]=ecount(store[[i]])
#    # m=m+1
#   }
#   
#   return(cnt)
# }

# edgeCount(lasso.net)

# mz=lasso.net[1:100]
# mz[[1]]


#newdata=data %>%dplyr::select(c(FiedlerValue,Normalized_FiedlerValue,spectral_radius,transitivity,modularity))

#  newdata=data %>%dplyr::select(c(edges,minCut,FiedlerValue,closeness, 
# modularity,diameter,betweenness,centrality_eigen))
#plot(y=y.test,x=lasso.pred.bestLambda,type="p",xlab="predicted",ylab="actual")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#[9]PCA with random forest*

#++++PCA--on--training--set
#pca on only the independent variables (features)
# set.seed(384394)
# pc.train <- prcomp(train[,-1],
#                    center = TRUE,
#                    scale. = TRUE)
# attributes(pc.train)
# summary(pc.train)
# print(pc.train$x[,1:2])
# 
# 
# ####---Prediction with pca test and train
# x1=data.frame(train[1])
# train.pca.1=data.frame(predict(pc.train,train))
# train.pca.2=cbind(train.pca.1,x1)
# 
# x2=data.frame(test[1])
# test.pca.1=data.frame(predict(pc.train,test))
# test.pca.2=cbind(test.pca.1,x2)
# ########--random forest with pca
# 
# set.seed(6748)
# PCARfmodel <- randomForest::randomForest(R~PC1+PC2, data =train.pca.2,
#                                          mtry=2,ntree=500,proximity=TRUE,importance=TRUE)
# 
# print(PCARfmodel)
# 
# ####----Prediction---train------
# pred.test.rfpca=predict(PCARfmodel,test.pca.2)
# s6=rbind(head(test$R),c(head(pred.test.rfpca)))
# row.names(s)=c("ActualVals","PredictedVals")
# s6
# 
# ################
# obs.rfpca=test$R
# pred.rfpca=pred.test.rfpca
# dist.R.rfpca=data.frame(y=obs.rfpca,x=pred.rfpca)
# 
# #################################
# ggplot(dist.R.rfpca,aes(x=x,y=y))+
#   xlab("predicted")+
#   ylab("actual")+
#   geom_point()+
#   geom_abline(color="darkblue")+
#   ggtitle("Plot of actual vs predicted radius for rfpca")
# 
# ####----Prediction---test------
# ## response. The values for the predicted response and train response is close
# ####--Test--error-MSE--test--###
# test.error.rfpca= pred.test.rfpca-test.pca.2$R 
# ####--RMSE---for--test--set
# RMSE_test.rfpca=sqrt(mean((test.error.rfpca)^2))
# RMSE_test.rfpca
# caret::RMSE(test.pca.2$R,pred.test.rfpca)
# 
# ####--RSQUARED---for--test--set
# actual.rfpca=test.pca.2$R
# R_squared.rfpca = 1 - (sum((actual.rfpca-pred.test.rfpca)^2))/(sum((actual.rfpca-mean(actual.rfpca))^2))
# R_squared.rfpca
# ####----Error--Rate----#########
# plot(PCARfmodel)
# 

# Summarize the model
# summary(model)
# # Make predictions
# predictions <- model %>% predict(test.data)
# # Model performance
# # (a) Prediction error, RMSE
# RMSE(predictions, test.data$sales)
# # (b) R-square
# R2(predictions, test.data$sales)

# ggplot(marketing, aes(x = youtube, y = sales)) +
#   geom_point() +
#   stat_smooth()

##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#[9]PCA with random forest**
  

#++++PCA--on--training--set
#pca on only the independent variables (features)
# set.seed(384394)
# pc.train <- prcomp(train[,-1],
#                    center = TRUE,
#                    scale. = TRUE)
# attributes(pc.train)
# summary(pc.train)
# print(pc.train$x[,1:2])
# 
# 
# ####---Prediction with pca test and train
# x1=data.frame(train[1])
# train.pca.1=data.frame(predict(pc.train,train))
# train.pca.2=cbind(train.pca.1,x1)
# 
# x2=data.frame(test[1])
# test.pca.1=data.frame(predict(pc.train,test))
# test.pca.2=cbind(test.pca.1,x2)
# ########--random forest with pca
# 
# set.seed(6748)
# PCARfmodel <- randomForest::randomForest(R~PC1+PC2, data =train.pca.2,
#                                          mtry=2,ntree=500,proximity=TRUE,importance=TRUE)
# 
# print(PCARfmodel)
# 
# ####----Prediction---train------
# pred.test.rfpca=predict(PCARfmodel,test.pca.2)
# s6=rbind(head(test$R),c(head(pred.test.rfpca)))
# row.names(s)=c("ActualVals","PredictedVals")
# s6
# 
# ################
# obs.rfpca=test$R
# pred.rfpca=pred.test.rfpca
# dist.R.rfpca=data.frame(y=obs.rfpca,x=pred.rfpca)
# 
# #################################
# ggplot(dist.R.rfpca,aes(x=x,y=y))+
#   xlab("predicted")+
#   ylab("actual")+
#   geom_point()+
#   geom_abline(color="darkblue")+
#   ggtitle("Plot of actual vs predicted radius for rfpca")
# 
# ####----Prediction---test------
# ## response. The values for the predicted response and train response is close
# ####--Test--error-MSE--test--###
# test.error.rfpca= pred.test.rfpca-test.pca.2$R 
# ####--RMSE---for--test--set
# RMSE_test.rfpca=sqrt(mean((test.error.rfpca)^2))
# RMSE_test.rfpca
# caret::RMSE(test.pca.2$R,pred.test.rfpca)
# 
# ####--RSQUARED---for--test--set
# actual.rfpca=test.pca.2$R
# R_squared.rfpca = 1 - (sum((actual.rfpca-pred.test.rfpca)^2))/(sum((actual.rfpca-mean(actual.rfpca))^2))
# R_squared.rfpca
# ####----Error--Rate----#########
# plot(PCARfmodel)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##Empirical Network
#data="mammalia-dolphin-social.edges"
#graphname="dolphin"
# allnets=function(data="mammalia-dolphin-social.edges", graphname="dolphin",csvfile="allnets3.csv"){
#   df <- read.table(data)  
#   df_new=graph_from_data_frame(as.matrix(df),directed=FALSE)
#   df.graph=igraph::simplify(df_new,remove.multiple = T,remove.loops = T)
#   components = igraph::clusters(df.graph , mode="weak")
#   biggest_cluster_id = which.max(components$csize)
#   vert_ids = V(df.graph)[components$membership== biggest_cluster_id]
#   graph=igraph::induced_subgraph(df.graph, vert_ids)
#   
#   graph$name=graphname
#   graph$type=graphname
#   graph$id="1"
#   G=list(graph)
#   
#   ###---Graph--Features--###
#   data=RunSimOnGraphFeatures(G,nreps = 1)
#   
#   newdata=data %>%dplyr::select(-c(GraphID,order, GraphReplicate,max_component,minDegree,threshold,maxDegree,GraphName,connected))
#   
#   ###---Predicted--Radius---##
#   predicted.radius.lasso=predict(lasso.mod, s = bestlambda.lasso, newx = as.matrix(newdata))
#   predicted.radius.ridge=predict(ridge.mod , s = bestlambda.ridge, newx = as.matrix(newdata))
#   predicted.radius.rfmodel=predict(rfmodel,newdata)
#   #predicted.radius.rfmodelpca=predict(PCARfmodel,newdata)
#   predicted.radius.glm=predict(glm.model.caret,newdata)
#   
#   ###--Predictiv--models
#   spatial.net.lasso=makeSpatialGraphs(node.size=vcount(df.graph),Radius=predicted.radius.lasso)
#   spatial.net.ridge=makeSpatialGraphs(node.size=vcount(df.graph),Radius=predicted.radius.ridge)
#   spatial.net.rfmodel=makeSpatialGraphs(node.size=vcount(df.graph),Radius=predicted.radius.rfmodel)
#   #spatial.net.rfmodelpca=makeSpatialGraphs(node.size=vcount(df.graph),Radius=predicted_radius.rfmodelpca)
#   #spatial.net.glm=makeSpatialGraphs(node.size=vcount(df.graph),Radius=predicted.radius.glm)
#   
#   ###---Theoretical--models
#   z=makeTheoGraphs(nSamples=1, order=62,edges=159,dim.sw=2,nei.sw=1,
#                    p.sw=0.2,power.sf=4,m.sf=2,r.sp=0.19,dim.lat=2,nei.lat=1)
#   z
#   
#   th.models=z
#   
#   ##--ERGMs---
#   set.seed(569)
#   
#   ergm_net <- asNetwork(graph)
#   summary(ergm_net ~ edges)
#   ergm.model <- ergm(ergm_net ~ edges) #1st fitted ERGM
#   summary(ergm.model)
#   #gf1=gof(ergm.model)
#   ergm.sim.model <- simulate(ergm.model,nsim=100)
#   
#   ergm.sim.net=lapply(ergm.sim.model,asIgraph)
#   #edgeCount_ergm=c(lapply(ergm.sim.net, ecount))
#   
#   
#   ergm.graph=ergm.sim.net[10]#ergm.sim.net[69]
#   
#   #fastSpatialNetwork(n=vcount(graph),r=predicted_radius,makeConnected=TRUE, keepCellsSeparate=FALSE)
#   
#   ###--Graph--Features--for--empirical--and--spatial--net
#   # spatial.net$name="Spatial"
#   # spatial.net$type="Spatial"
#   # spatial.net$id="1"
#   set.seed(7349)
#   net=c(G,ergm.graph,spatial.net.lasso,spatial.net.ridge,spatial.net.rfmodel,spatial.net.glm,th.models)
#   feat=RunSimOnGraphFeatures(net,nreps = 1)
#   feat$GraphName=c(graphname,"ERGM","SPNet.Lasso","SPNet.Ridge","SPNet.RF","SPNet.GLM","ER","SW","SF","SP","LAT")
#   write.csv(feat,csvfile)
#   return(feat)
# }
# x1=allnets(data="mammalia-bat-roosting-indiana.edges", graphname="ant",csvfile="dolph&nets.csv")
# x1
#new.list <- list(spatial.net[[1]],graph,feat)
#x=c(G,spatial.net)
#x$summary=feat
#return(x)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+Violin Plot
#+# violin plots

#```{r}

## | fig-cap: "Violin plots of Graph Features."
## | fig-subcap:
## |   - "violin plot of Fiedler"
## |   - "violin plot of spectralRadius"
## |   - "violin plot of modularity"
## |   - "violin plot of NormFiedler"
## |   - "violin plot of MeanDeg"
## |   - "violin plot of Edges"
## |   - "violin plot of Closeness"
## |   - "violin plot of Diameter"
## |   - "violin plot of Transitivity"
## |   - "violin plot of Centrality"
## |   - "violin plot of Betweeness"
## |   - "violin plot of Mincut"
## | layout-nrow=6
## | layout-ncol: 2
## | column: page
## | code-fold: TRUE
## | fig-width: 4
## | fig-height: 6

# df.plot %>% ggplot(aes(x=R,y=FiedlerValue))+
#    geom_violin(show.legend = F,alpha=0.5,adjust=0.5,draw_quantiles = c(0.5))
# 
# df.plot %>% ggplot(aes(x=R,y=spectral_radius))+
#    geom_violin(show.legend = F,alpha=0.5,adjust=0.5,draw_quantiles = c(0.5))
# 
# df.plot %>% ggplot(aes(x=R,y=modularity))+
#    geom_violin(show.legend = F,alpha=0.5,adjust=0.5,draw_quantiles = c(0.5))
# 
# df.plot %>% ggplot(aes(x=R,y=Normalized_FiedlerValue))+
#    geom_violin(show.legend = F,alpha=0.5,adjust=0.5,draw_quantiles = c(0.5))
# 
# df.plot %>% ggplot(aes(x=R,y=mean_degree))+
#    geom_violin(show.legend = F,alpha=0.5,adjust=0.5,draw_quantiles = c(0.5))
# 
# df.plot %>% ggplot(aes(x=R,y=edges))+
#    geom_violin(show.legend = F,alpha=0.5,adjust=0.5,draw_quantiles = c(0.5))
# 
# df.plot %>% ggplot(aes(x=R,y=closeness))+
#    geom_violin(show.legend = F,alpha=0.5,adjust=0.5,draw_quantiles = c(0.5))
# 
# df.plot %>% ggplot(aes(x=R,y=diameter))+
#    geom_violin(show.legend = F,alpha=0.5,adjust=0.5,draw_quantiles = c(0.5))
# 
# df.plot %>% ggplot(aes(x=R,y=transitivity))+
#    geom_violin(show.legend = F,alpha=0.5,adjust=0.5,draw_quantiles = c(0.5))
# 
# df.plot %>% ggplot(aes(x=R,y=centrality_eigen))+
#    geom_violin(show.legend = F,alpha=0.5,adjust=0.5,draw_quantiles = c(0.5))
# 
# df.plot %>% ggplot(aes(x=R,y=betweenness))+
#    geom_violin(show.legend = F,alpha=0.5,adjust=0.5,draw_quantiles = c(0.5))
# 
# df.plot %>% ggplot(aes(x=R,y=minCut))+
#    geom_violin(show.legend = F,alpha=0.5,adjust=0.5,draw_quantiles = c(0.5))
#```
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+
# networks.hybrid=function(R=c(0.3,0.4,0.6,0.05),p=c(0.1,0.2,0.6,0.4),n=100){
#   net=NULL;data=NULL
#   k=1
#   for (i in 1:length(R)){
#     net[i]=makeSpatialHybrid(node.size =n,Radius = R[i],prob=p)
#     data=RunSimOnGraphFeatures(net,nreps = 1)
#     k=k+1
#   }
#   df=cbind(R,p,data)
#   return(df)
# }
# 
# m15.100=networks.hybrid(R=rep(0.15,100),p=rep(0.8,100),n=100)
# write.csv(m15.100,"m15-100.csv")
# 
# m25.100=networks.hybrid(R=rep(0.2,100),p=rep(0.8,100),n=100)
# write.csv(m25.100,"m25-100.csv")
# 
# m35.100=networks.hybrid(R=rep(0.3,100),p=rep(0.8,100),n=100)
# write.csv(m35.100,"m35-100.csv")
# 
# m45.100=networks.hybrid(R=rep(0.4,100),p=rep(0.8,100),n=100)
# write.csv(m45.100,"m45-100.csv")
# 
# m55.100=networks.hybrid(R=rep(0.8,100),p=rep(0.8,100),n=100)
# write.csv(m55.100,"m55-100.csv")
# 
# m65.100=networks.hybrid(R=rep(0.6,100),p=rep(0.8,100),n=100)
# write.csv(m65.100,"m65-100.csv")
# 
# m75.100=networks.hybrid(R=rep(0.7,100),p=rep(0.8,100),n=100)
# write.csv(m75.100,"m75-100.csv")
# 
# m85.100=networks.hybrid(R=rep(0.8,100),p=rep(0.8,100),n=100)
# write.csv(m85.100,"m85-100.csv")
# 
# m95.100=networks.hybrid(R=rep(0.9,100),p=rep(0.8,100),n=100)
# write.csv(m95.100,"m95-100.csv")
# 
# m10.5.100=networks.hybrid(R=rep(1,100),p=rep(0.8,100),n=100)
# write.csv(m10.5.100,"m10-5-100.csv")


##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++###
                    
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# **3.6 Test with empirical network for multiple outcome**
#   
#   ```{r}
# #| echo: true
# #| include: true
# #| code-fold: true
# library(ergm)
# #--empirical--network--test
# 
# allnets2=function(emp.data2="mammalia-dolphin-social.edges", graphname2="dolphin",csvfile2="allnets3.csv",sample.sim=4,nsim=20){
#   graph.data2<- read.table(emp.data2)  
#   graph.new2=graph_from_data_frame(as.matrix(graph.data2),directed=FALSE)
#   g.graph2=igraph::simplify(graph.new2,remove.multiple = T,remove.loops = T)
#   g=list(g.graph2)
#   inducedgraph2=induced.graph(g)
#   
#   inducedgraph2$name=graphname2
#   inducedgraph2$type=graphname2
#   inducedgraph2$id="1"
#   G2=list(inducedgraph2)
#   
#   ###---Graph--Features--for-empirical###
#   emp_graphfeat2=RunSimOnGraphFeatures(G2,nreps = 1)
#   
#   input.data2=emp_graphfeat2 %>%dplyr::select(c(order,edges,mean_degree,minCut,FiedlerValue,Normalized_FiedlerValue,closeness,
#                                                 modularity,diameter,betweenness,transitivity,spectral_radius,centrality_eigen))
#   ###---Predicted--Radius---##
#   predicted.radius.ridge2=predict(ridge2,s=ridge2$lambda.1se,newx = as.matrix(input.data2))
#   predicted.radius.ridge2=data.frame(predicted.radius.ridge2)
#   predicted.radius.lasso2=predict(las2,s=las2$lambda.1se,newx = as.matrix(input.data2))
#   predicted.radius.lasso2=data.frame(predicted.radius.lasso2)
#   predicted.radius.Elasticnet2=predict(Elasticnet2,s=Elasticnet2$lambda.1se,newx = as.matrix(input.data2))
#   predicted.radius.Elasticnet2=data.frame(predicted.radius.Elasticnet2)
#   #predicted.radius.elasticnet2=predict(Elasticnet,s=Elasticnet$lambda.1se,newx = as.matrix(Sp.Net))
#   
#   
#   
#   Hybrid.Ridge.sim=simulate.spatialhybrid(N=vcount(g.graph2),radius=predicted.radius.ridge2$y1.1,
#                                           p=predicted.radius.ridge2$y2.1,nsim=nsim)
#   
#   Hybrid.Lasso.sim=simulate.spatialhybrid(N=vcount(g.graph2),radius=predicted.radius.lasso2$y1.1,
#                                           p=predicted.radius.lasso2$y2.1,nsim=nsim)
#   
#   Hybrid.ElsticNet.sim=simulate.spatialhybrid(N=vcount(g.graph2),radius=predicted.radius.Elasticnet2$y1.1,p=predicted.radius.Elasticnet2$y2.1,nsim=nsim)
#   
#   
#   Hybrid.RidgeGraph=Hybrid.Ridge.sim[[sample.sim]]
#   Hybrid.LassoGraph=Hybrid.Lasso.sim[[sample.sim]]
#   Hybrid.ElsticNetGraph=Hybrid.ElsticNet.sim[[sample.sim]]
#   #SP.EnetGraph=makeSpatialGraphs(node.size=n,Radius=predicted.radius.elasticnet)
#   
#   # spatial.net.rfmodel=makeSpatialGraphs(node.size=vcount(inducedgraph),Radius=predicted.radius.rfmodel)
#   # spatial.hybrid.rfmodel=makeSpatialHybrid(node.size = vcount(inducedgraph),Radius = predicted.radius.rfmodel,p=p)
#   # spatial.net.rfmodel2=makeSpatialGraphs(node.size=vcount(inducedgraph),Radius=predicted.radius.rfmodel2)
#   # spatial.hybrid.rfmodel2=makeSpatialHybrid(node.size = vcount(inducedgraph),Radius = predicted.radius.rfmodel2,p=p)
#   ##--ERGMs---
#   set.seed(569)
#   ergm_net2 <- asNetwork(inducedgraph2)
#   #  summary(ergm_net ~ edges)
#   ergm.model2 <- ergm(ergm_net2 ~ edges) #1st fitted ERGM
#   #  summary(ergm.model)
#   #gf1=gof(ergm.model)
#   ergm.sim.model2 <- simulate(ergm.model2,nsim=nsim)
#   
#   ergm.sim.net2=lapply(ergm.sim.model2,asIgraph)
#   #edgeCount_ergm=c(lapply(ergm.sim.net, ecount))
#   
#   
#   ergm.graph2=ergm.sim.net2[sample.sim]#ergm.sim.net[69]
#   
#   set.seed(7349)
#   net2=c(G2,ergm.graph2,Hybrid.RidgeGraph, Hybrid.LassoGraph, Hybrid.ElsticNetGraph)
#   
#   ###---Graph--Features-for All-nets-###
#   predvals=NULL
#   graph.feat2=RunSimOnGraphFeatures(net2,nreps = 1)
#   predvals$R=c(predicted.radius.ridge2$y1.1,predicted.radius.lasso2$y1.1,predicted.radius.Elasticnet2$y1.1)
#   predvals$p=c(predicted.radius.ridge2$y2.1,predicted.radius.lasso2$y2.1,predicted.radius.Elasticnet2$y2.1)
#   graph.feat2$GraphName=c(graphname2,"ERGM","Ridge-Hybrid","Lasso-Hybrid","ElasticNet-Hybrid")
#   
#   write.csv(graph.feat2,csvfile2)
#   results=list(graph.feat2,predvals)
#   return(results)
# }
# #set.seed(7439)
# #macaque
# x1=allnets2(emp.data2="mammalia-macaque-contact-sits.edges",sample.sim=4,nsim=50, graphname2="macaque",csvfile2="1.csv")
# x1[[1]]
# #dolphin
# x2=allnets2(emp.data2="mammalia-dolphin-social.edges", sample.sim=2,nsim=50,
#             graphname2="dolphin.social",csvfile2="1.csv")
# x2[[1]]
# #hyena
# x3=allnets2(emp.data2="mammalia-hyena-networkb.edges", graphname2="hyena.b",csvfile2="1.csv")
# x3[[1]]
# #ant
# x4=allnets2(emp.data2="insecta-ant-colony1-day01.edges", graphname2="antd1",csvfile2="i1.csv")
# x4[[1]]
# #weaver
# x5=allnets2(emp.data2="aves-weaver-social.edges", graphname2="weaver.social",csvfile2="weaver.csv")
# x5[[1]]
# #kangaroo
# x6=allnets2(emp.data2="mammalia-kangaroo-interactions.edges",graphname2="kangaroo",csvfile2="kangaroo.csv")
# x6[[1]]
# #Asianelephant
# x7=allnets2(emp.data2="mammalia-asianelephant.edges", graphname2="asianeleph",csvfile2="elephant.csv")
# x7[[1]]
# #chesapeake
# x8=allnets(emp.data="road-chesapeake.edges", graphname="road.chesap",csvfile="chesap.csv",p=0.1)
# x8
# #bat
# x9=allnets(emp.data="mammalia-bat-roosting-indiana.edges", graphname="bat",csvfile="bat.csv",p=0.1)
# x9
# #tortoise
# x10=allnets(emp.data="reptilia-tortoise-network-fi.edges", graphname="tortoise",csvfile="tortoise.csv",p=0.1)
#x10
#```


##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##

# **3.3 Test with empirical network for one outcome**
#   
#   ```{r}
# #| echo: true
# #| include: true
# #| code-fold: true
# library(ergm)
# #--empirical--network--test
# allnets=function(emp.data="mammalia-dolphin-social.edges", graphname="dolphin",csvfile="allnets3.csv",p=0.4){
#   graph.data<- read.table(emp.data)  
#   graph.new=graph_from_data_frame(as.matrix(graph.data),directed=FALSE)
#   g.graph=igraph::simplify(graph.new,remove.multiple = T,remove.loops = T)
#   components = igraph::clusters(g.graph , mode="weak")
#   biggest_cluster_id = which.max(components$csize)
#   vert_ids = V(g.graph)[components$membership== biggest_cluster_id]
#   inducedgraph=igraph::induced_subgraph(g.graph, vert_ids)
#   
#   inducedgraph$name=graphname
#   inducedgraph$type=graphname
#   inducedgraph$id="1"
#   G=list(inducedgraph)
#   
#   ###---Graph--Features--for-empirical###
#   emp_graphfeat=RunSimOnGraphFeatures(G,nreps = 1)
#   
#   input.data=emp_graphfeat %>%dplyr::select(c(edges,mean_degree,minCut,FiedlerValue,Normalized_FiedlerValue,closeness,
#                                               modularity,diameter,betweenness,transitivity,spectral_radius,centrality_eigen))
#   ###---Predicted--Radius---##
#   predicted.radius.ridge=predict(ridge,s=ridge$lambda.1se,newx = as.matrix(input.data))
#   #predicted.radius.lasso=predict(las, s=las$lambda.1se, newx = as.matrix(input.data))
#   #predicted.radius.elastnet=predict(Elasticnet,s=Elasticnet$lambda.1se, newx = as.matrix(input.data))
#   predicted.radius.rfmodel=predict(rfmodel,input.data)
#   predicted.radius.rfmodel2=predict(rfmodel.2,input.data)
#   
#   
#   ###--Predictiv--models
#   #spatial.net.lasso=makeSpatialGraphs(node.size=vcount(inducedgraph),Radius=predicted.radius.lasso)
#   #spatial.net.elasnet=makeSpatialGraphs(node.size=vcount(inducedgraph),Radius=predicted.radius.elastnet)
#   spatial.net.ridge=makeSpatialGraphs(node.size=vcount(inducedgraph),Radius=predicted.radius.ridge)
#   spatial.hybrid.ridge=makeSpatialHybrid(node.size = vcount(inducedgraph),Radius = predicted.radius.ridge,p=p)
#   spatial.net.rfmodel=makeSpatialGraphs(node.size=vcount(inducedgraph),Radius=predicted.radius.rfmodel)
#   spatial.hybrid.rfmodel=makeSpatialHybrid(node.size = vcount(inducedgraph),Radius = predicted.radius.rfmodel,p=p)
#   spatial.net.rfmodel2=makeSpatialGraphs(node.size=vcount(inducedgraph),Radius=predicted.radius.rfmodel2)
#   spatial.hybrid.rfmodel2=makeSpatialHybrid(node.size = vcount(inducedgraph),Radius = predicted.radius.rfmodel2,p=p)
#   ##--ERGMs---
#   set.seed(569)
#   ergm_net <- asNetwork(inducedgraph)
#   #  summary(ergm_net ~ edges)
#   ergm.model <- ergm(ergm_net ~ edges) #1st fitted ERGM
#   #  summary(ergm.model)
#   #gf1=gof(ergm.model)
#   ergm.sim.model <- simulate(ergm.model,nsim=10)
#   
#   ergm.sim.net=lapply(ergm.sim.model,asIgraph)
#   #edgeCount_ergm=c(lapply(ergm.sim.net, ecount))
#   
#   
#   ergm.graph=ergm.sim.net[1]#ergm.sim.net[69]
#   
#   set.seed(7349)
#   net=c(G,ergm.graph,spatial.net.ridge,spatial.hybrid.ridge,spatial.net.rfmodel,spatial.hybrid.rfmodel,spatial.net.rfmodel2,spatial.hybrid.rfmodel2)
#   
#   ###---Graph--Features-for All-nets-###
#   graph.feat=RunSimOnGraphFeatures(net,nreps = 1)
#   graph.feat$GraphName=c(graphname,"ERGM","SPNet.Ridge","Hybrid.Ridge","SPNet.RF","Hybrid.RF","SPNet.RF2","Hybrid.RF2")
#   write.csv(graph.feat,csvfile)
#   return(graph.feat)
# }
# set.seed(7439)
# #macaque
# x1=allnets(emp.data="mammalia-macaque-contact-sits.edges", graphname="macaque",csvfile="1.csv",p=0.2)
# x1
# #dolphin
# x2=allnets(emp.data="mammalia-dolphin-social.edges", graphname="dolphin.social",csvfile="dolphin.csv",p=0.1)
# x2
# #hyena
# x3=allnets(emp.data="mammalia-hyena-networkb.edges", graphname="hyena.b",csvfile="hyena.csv",p=0.1)
# x3
# #ant
# x4=allnets(emp.data="insecta-ant-colony1-day01.edges", graphname="antd1",csvfile="insect-day1.csv",p=0.1)
# x4
# #weaver
# x5=allnets(emp.data="aves-weaver-social.edges", graphname="weaver.social",csvfile="weaver.csv",p=0.1)
# x5
# #kangaroo
# x6=allnets(emp.data="mammalia-kangaroo-interactions.edges", graphname="kangaroo",csvfile="kangaroo.csv",p=0.1)
# x6
# #Asianelephant
# x7=allnets(emp.data="mammalia-asianelephant.edges", graphname="asianeleph",csvfile="elephant.csv",p=0.1)
# x7
# #chesapeake
# x8=allnets(emp.data="road-chesapeake.edges", graphname="road.chesap",csvfile="chesap.csv",p=0.1)
# x8
# #bat
# x9=allnets(emp.data="mammalia-bat-roosting-indiana.edges", graphname="bat",csvfile="bat.csv",p=0.1)
# x9
# #tortoise
# x10=allnets(emp.data="reptilia-tortoise-network-fi.edges", graphname="tortoise",csvfile="tortoise.csv",p=0.1)
# x10
# ```

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# g <- igraph::sample_gnp(20, 1/20)
# 
# components <- igraph::clusters(g, mode="weak")
# biggest_cluster_id <- which.max(components$csize)
# 
# # ids
# vert_ids <- V(g)[components$membership == biggest_cluster_id]
# 
# # subgraph
# igraph::induced_subgraph(g, vert_ids)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# **3.1 Models: Allmodels for both 1 and 2 response variables**
#   
#   ```{r}
# #| echo: false
# #| include: false
# #| code-fold: true
# 
# 
# #########---------Turning traininga dn test set to matrix------####
# x.train=as.matrix(train[,-c(1,2)])#train and test sets excluding the response
# x.test=as.matrix(test[,-c(1,2)])
# 
# y.train.1= train$R #train and test sets of the response
# y.test.1=test$R
# 
# y.train.2= cbind(train$R,train$p) #train and test sets of the response
# y.test.2=cbind(test$R,test$p)
# 
# ##-------set controls
# train.control <- trainControl(method = "cv", number = 10,savePredictions = "all")
# ## List of potential Lambda values
# lambda.vect=10^seq(10, -2, length = 500)
# 
# set.seed(29394)
# ####---Ridge---regression--one--outcome
# ridge=cv.glmnet(x.train,y.train.1,type.measure = "mse",alpha=0,family="gaussian")
# ridgepredict=predict(ridge,s=ridge$lambda.1se,newx=x.test)
# mean((y.test.1-ridgepredict)^2)
# 
# ####---Ridge---regression--two--outcomes
# ridge2=cv.glmnet(x.train,y.train.2,type.measure = "mse",alpha=0,family="mgaussian")
# ridgepredict2=predict(ridge2,s=ridge2$lambda.1se,newx=x.test)
# 
# #### mse for Radius
# mean((as.data.frame(y.test.2)[,1]-as.data.frame(ridgepredict2)[,1])^2)
# #### mse for p
# mean((as.data.frame(y.test.2)[,2]-as.data.frame(ridgepredict2)[,2])^2)
# 
# #### RMSE for Radius
# caret::RMSE(pred = as.data.frame(ridgepredict2)[,1], obs =as.data.frame(y.test.2)[,1])
# #### RMSE for p
# caret::RMSE(pred = as.data.frame(ridgepredict2)[,2], obs =as.data.frame(y.test.2)[,2])
# #### R2for Radius
# caret::R2(pred = as.data.frame(ridgepredict2)[,1], obs =as.data.frame(y.test.2)[,1])
# #### R2 for p
# caret::R2(pred = as.data.frame(ridgepredict2)[,2], obs =as.data.frame(y.test.2)[,2])
# 
# ####---lasso---regression--for--one--outcome--###
# set.seed(29394)
# las=cv.glmnet(x.train,y.train.1,type.measure = "mse",alpha=1,family="gaussian")
# laspredict=predict(las,s=las$lambda.1se,newx=x.test)
# mean((y.test.1-laspredict)^2)
# 
# ####---lasso---regression--for--two--outcome--###
# las2=cv.glmnet(x.train,y.train.2,type.measure = "mse",alpha=1,family="mgaussian")
# laspredict2=predict(las2,s=las2$lambda.1se,newx=x.test)
# #### mse for Radius
# mean((as.data.frame(y.test.2)[,1]-as.data.frame(laspredict2)[,1])^2)
# #### mse for p
# mean((as.data.frame(y.test.2)[,2]-as.data.frame(laspredict2)[,2])^2)
# 
# #### RMSE for Radius
# caret::RMSE(pred = as.data.frame(laspredict2)[,1], obs =as.data.frame(y.test.2)[,1])
# #### RMSE for p
# caret::RMSE(pred = as.data.frame(laspredict2)[,2], obs =as.data.frame(y.test.2)[,2])
# #### R2for Radius
# caret::R2(pred = as.data.frame(laspredict2)[,1], obs =as.data.frame(y.test.2)[,1])
# #### R2 for p
# caret::R2(pred = as.data.frame(laspredict2)[,2], obs =as.data.frame(y.test.2)[,2])
# 
# 
# ####++++Elastic----net--for--one--and--two---outcomes-------#####
# ####
# set.seed(28443)
# Models=function(){
#   models.list=list();res=data.frame()
#   models.list2=list();res2=data.frame()
#   allresults=NULL
#   for (j in 0:10){
#     models.names=paste0("alpha",j/10)
#     ####### 1 outocme modes
#     models.list[[models.names]]=cv.glmnet(x.train,y.train.1,type.measure = "mse",alpha=j/10,
#                                           family="gaussian")
#     ####### 2 outocmes modes
#     models.names2=paste0("alpha",j/10)
#     models.list2[[models.names2]]=cv.glmnet(x.train,y.train.2,type.measure = "mse",alpha=j/10,
#                                             family="mgaussian")
#   }
#   
#   for (j in 0:10){
#     models.names=paste0("alpha",j/10)
#     models.names2=paste0("alpha",j/10)
#     pred.models=predict(models.list[[models.names]],s=models.list[[models.names]]$lambda.1se,newx=x.test)
#     
#     pred.models2=predict(models.list2[[models.names2]],s=models.list2[[models.names2]]$lambda.1se,newx=x.test)
#     
#     #### mse of Radius for 1 outcome
#     MSE.1=mean((y.test.1-pred.models)^2)
#     #### mse of Radius for 2 outcomes
#     MSE.2=mean((as.data.frame(y.test.2)[,1]-as.data.frame(pred.models2)[,1])^2)
#     #### RMSE of Radius for 1 outcome 
#     RMSE.1=caret::RMSE(pred = pred.models, obs =y.test.1)
#     #### RMSE of Radius for 2 outcomes
#     RMSE.2=caret::RMSE(pred = as.data.frame(pred.models2)[,1], obs =as.data.frame(y.test.2)[,1])
#     #### RSquared of Radius for 1 outcome
#     Rsquared.1=caret::R2(pred = pred.models, obs =y.test.1)
#     #### RSquared of Radius for 2 outcomes
#     Rsquared.2=caret::R2(pred = as.data.frame(pred.models2)[,1], obs =as.data.frame(y.test.2)[,1])
#     ### final model for 1 outcome
#     mod=data.frame(alpha=j/10,mse=MSE.1,RMSE=RMSE.1,Rsquared=Rsquared.1,models.names=models.names)
#     ### final model for 2 outcomes
#     mod2=data.frame(alpha=j/10,mse=MSE.2,RMSE=RMSE.2,Rsquared=Rsquared.2,models.names2=models.names2)
#     #### result for 1 outcome
#     res=rbind(res,mod)
#     ### result for two outcome
#     res2=rbind(res2,mod2)
#     
#     allresults$result1=res
#     allresults$result2=res2
#     # allres=list(res,res2)
#     
#   }
#   return(allresults)
# }
# 
# set.seed(8362)
# output=Models()
# ###--minimum--mse--for--1--outcome
# min(output$result1$mse)
# 
# ###--minimum--mse--for--2--outcomes
# min(output$result2$mse)
# ###--show--the--alpha--with--the--min--mse--for--1--outcome
# minalpha1=output$result1[output$result1$mse==min(output$result1$mse),]
# minalpha1
# 
# ###--show--the--alpha--with--the--min--mse--for--2--outcomes
# minalpha2=output$result2[output$result2$mse==min(output$result2$mse),]
# minalpha2
# 
# ###Ue minimum alphas to perform elastic net
# set.seed(29394)
# ####---elastic---regression--for--1--outcome
# Elasticnet=cv.glmnet(x.train,y.train.1,type.measure = "mse",alpha=minalpha1$alpha,family="gaussian")
# Elasticnet.predict=predict(Elasticnet,s=Elasticnet$lambda.1se,newx=x.test)
# mean((y.test.1-Elasticnet.predict)^2)
# ####---elastic---regression--for--2--outcomes
# Elasticnet2=cv.glmnet(x.train,y.train.2,type.measure = "mse",alpha=minalpha2$alpha,family="mgaussian")
# Elasticnet.predict2=predict(Elasticnet2,s=Elasticnet2$lambda.1se,newx=x.test)
# mean((y.test.2[,1]-Elasticnet.predict2)^2)
# 
# #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ###-----------------------Random--forest-----------#####
# ### Tuning the number of mtrys
# set.seed((6283))
# tn=tuneRF(train[,-1],train$R,stepFactor=0.5,plot = TRUE,ntreeTry = 500,
#           trace = TRUE,improve = 0.01) #works only if mtry
# #> the number of variables(features). This is because mtry is the number of randomly sampled variable as candidate at each split. Base case use mtry=2  when this happens
# 
# tn=as.data.frame(tn)
# tn.min=tn$mtry[tn$OOBError== min(tn$OOBError)] 
# 
# set.seed(6748)
# rfmodel <- randomForest::randomForest(R~ ., data = train[,-2],
#                                       trainControl=train.control,
#                                       mtry=tn.min,ntree=500,proximity=TRUE,importance=TRUE)
# print(rfmodel)
# ####----Prediction---train------
# pred.train.rf=predict(rfmodel,train)
# ####----Prediction---test------
# pred.test.rf=predict(rfmodel,test)
# 
# s5=rbind(head(test$R),c(head(pred.test.rf)))
# row.names(s5)=c("ActualVals","PredictedVals")
# s5
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
# caret::RMSE(test$R,pred.test.rf)
# ####--RSQUARED---for--test--set
# caret::R2(test$R,pred.test.rf)
# ####----Error--Rate----#########
# plot(rfmodel)
# 
# ##Random-forest-model-2
# set.seed(6748)
# rfmodel.2 <- randomForest::randomForest(R~ spectral_radius+centrality_eigen+closeness+minCut+transitivity+betweenness+FiedlerValue+modularity+Normalized_FiedlerValue+diameter, data = train[,-2],
#                                         trainControl=train.control,
#                                         mtry=tn.min,ntree=500,proximity=TRUE,importance=TRUE)
# print(rfmodel.2)
# # ####----Prediction---train------
# pred.train.rf2=predict(rfmodel.2,train)
# # ####----Prediction---test------
# pred.test.rf2=predict(rfmodel.2,test)
# # 
# # s6=rbind(head(test$R),c(head(pred.test.rf2)))
# # row.names(s6)=c("ActualVals","PredictedVals")
# # s6
# # 
# # ####--Test--error-MSE--test--###
# # test.error.rf2= pred.test.rf2-test$R 
# # ####--RMSE---for--test--set
# # caret::RMSE(test$R,pred.test.rf2)
# # ####--RSQUARED---for--test--set
# # caret::R2(test$R,pred.test.rf2)
# # ####----Error--Rate----#########
# # plot(rfmodel.2)
# ```
# 
# 
# 
# 
# 
# **Test with empirical networks for spatial and spatial hybrid models**
#   
#   ```{r}
# #| echo: true
# #| include: true
# #| code-fold: true
# library(ergm)
# #--empirical--network--test
# allnets=function(emp.data="mammalia-dolphin-social.edges", graphname="dolphin",
#                  csvfile="allnets3.csv",sample.sim=4,nsim=20){
#   graph.data<- read.table(emp.data)  
#   graph.new=graph_from_data_frame(as.matrix(graph.data),directed=FALSE)
#   g.graph=igraph::simplify(graph.new,remove.multiple = T,remove.loops = T)
#   g=list(g.graph)
#   inducedgraph=induced.graph(g)
#   
#   inducedgraph$name=graphname
#   inducedgraph$type=graphname
#   inducedgraph$id="1"
#   G=list(inducedgraph)
#   
#   #++++++++++++++++++++++++Graph--Features--for-empirical+++++++++++++++++++++++++++++++++++++++++++++++++++#
#   emp_graphfeat=RunSimOnGraphFeatures(G,nreps = 1)
#   input.data=emp_graphfeat %>%dplyr::select(c(order,edges,mean_degree,minCut,FiedlerValue,Normalized_FiedlerValue,closeness,                          modularity,diameter,betweenness,transitivity,spectral_radius,centrality_eigen))
#   
#   #+++++++++++++++++++++++++Predicted--Radius--for Spatial--models+++++++++++++++++++++++++++++++++++++#
#   predicted.radius.ridge=predict(ridge,s=ridge$lambda.1se,newx = as.matrix(input.data))
#   predicted.radius.lasso=predict(las, s=las$lambda.1se, newx = as.matrix(input.data))
#   predicted.radius.elastnet=predict(Elasticnet,s=Elasticnet$lambda.1se, newx = as.matrix(input.data))
#   #predicted.radius.rfmodel=predict(rfmodel,input.data)
#   #predicted.radius.rfmodel2=predict(rfmodel.2,input.data)
#   
#   #+++++++++++++++++++++++++Predicted--Radius--for Spatial-hybrid--models++++++++++++++++++++++++++++++++#
#   predicted.radius.ridge2=predict(ridge2,s=ridge2$lambda.1se,newx = as.matrix(input.data))
#   predicted.radius.ridge2=data.frame(predicted.radius.ridge2)
#   predicted.radius.lasso2=predict(las2,s=las2$lambda.1se,newx = as.matrix(input.data))
#   predicted.radius.lasso2=data.frame(predicted.radius.lasso2)
#   predicted.radius.Elasticnet2=predict(Elasticnet2,s=Elasticnet2$lambda.1se,newx = as.matrix(input.data))
#   predicted.radius.Elasticnet2=data.frame(predicted.radius.Elasticnet2)
#   
#   
#   #++++++++++++++++++++++++++++++++++++Simulate--Spatial--models++++++++++++++++++++++++++++++++++++#
#   ##--ridge
#   spatial.net.ridge=simulate.spatial(N=vcount(inducedgraph),radius=predicted.radius.ridge,nsim=nsim)
#   
#   ##--lasso
#   spatial.net.lasso=simulate.spatial(N=vcount(inducedgraph),radius=predicted.radius.lasso,nsim=nsim)
#   ##--elasticnet
#   spatial.net.elasnet=simulate.spatial(N=vcount(inducedgraph),radius=predicted.radius.elastnet,nsim=nsim)
#   
#   #++++++++++++++++++++++++++++++++++++Simulate--Spatial--hybrid--networks++++++++++++++++++++++++++++++++++#
#   ###--Ridge
#   Hybrid.Ridge.sim=simulate.spatialhybrid(N=vcount(g.graph),radius=predicted.radius.ridge2$y1.1,
#                                           p=predicted.radius.ridge2$y2.1,nsim=nsim)
#   ###--Lasso
#   Hybrid.Lasso.sim=simulate.spatialhybrid(N=vcount(g.graph),radius=predicted.radius.lasso2$y1.1,
#                                           p=predicted.radius.lasso2$y2.1,nsim=nsim)
#   ###--Elasticnetw
#   Hybrid.ElasticNet.sim=simulate.spatialhybrid(N=vcount(g.graph),radius=predicted.radius.Elasticnet2$y1.1,
#                                                p=predicted.radius.Elasticnet2$y2.1,nsim=nsim)
#   
#   #+++++++++++++++++++++++++++++++++++Get--smaple----graph--for--spatial--models++++++++++++++++++++++++++++#  
#   spatial.ridgegraph=spatial.net.ridge[[sample.sim]]
#   spatial.lassograph=spatial.net.lasso[[sample.sim]]
#   spatial.elasnetgraph=spatial.net.lasso[[sample.sim]]
#   
#   #+++++++++++++++++++++++++++++++++++Get--smaple----graph--for--spatial--hybrid--models++++++++++++++++++++++++++++#  
#   Hybrid.RidgeGraph=Hybrid.Ridge.sim[[sample.sim]]
#   Hybrid.LassoGraph=Hybrid.Lasso.sim[[sample.sim]]
#   Hybrid.ElasticNetGraph=Hybrid.ElasticNet.sim[[sample.sim]]
#   
#   #++++++++++++++++++++++++++++++++++++ERGMs+++++++++++++++++++++++++++++++++++++++++
#   set.seed(569)
#   ergm_net <- asNetwork(inducedgraph)
#   #  summary(ergm_net ~ edges)
#   ergm.model <- ergm(ergm_net ~ edges) #1st fitted ERGM
#   #  summary(ergm.model)
#   #gf1=gof(ergm.model)
#   ergm.sim.model <- simulate(ergm.model,nsim=nsim)
#   
#   ergm.sim.net=lapply(ergm.sim.model,asIgraph)
#   
#   ergm.graph=ergm.sim.net[sample.sim]
#   
#   #++++++++++++++++++++++++++++++++All--graph--models+++++++++++++++++++++++++++++++++++++++++++++++++++++++  
#   set.seed(7349)
#   net=c(G,ergm.graph,spatial.ridgegraph,Hybrid.RidgeGraph,spatial.lassograph,
#         Hybrid.LassoGraph,spatial.elasnetgraph,Hybrid.ElasticNetGraph)
#   
#   #++++++++++++++++++++++++++++++++Graph--Features-for-all-graph-models+++++++++++++++++++++++++++++++++++++#
#   graph.feat=RunSimOnGraphFeatures(net,nreps = 1)
#   
#   graph.feat$GraphName=c(graphname,"ERGM","SP.RidgeModel","Hybrid.RidgeModel",
#                          "SP.LassoModel","Hybrid.LassoModel","SP.ElasticNetModel","Hybrid.ElasticNetModel")
#   
#   
#   predvals=NULL
#   predvals$R=c(predicted.radius.ridge2$y1.1,predicted.radius.lasso2$y1.1,predicted.radius.Elasticnet2$y1.1)
#   predvals$p=c(predicted.radius.ridge2$y2.1,predicted.radius.lasso2$y2.1,predicted.radius.Elasticnet2$y2.1)
#   results=list(graph.feat,predvals)
#   
#   write.csv(graph.feat,csvfile)
#   return(results)
# }
# set.seed(7439)
# 
# x1=allnets(emp.data="mammalia-macaque-contact-sits.edges",sample.sim=4,nsim=50, graphname="macaque",csvfile="1.csv")
# x1[[1]]
# #dolphin
# x2=allnets(emp.data="mammalia-dolphin-social.edges", sample.sim=,nsim=50,
#            graphname="dolphin.social",csvfile="1.csv")
# x2[[1]]
# #hyena
# x3=allnets(emp.data="mammalia-hyena-networkb.edges", sample.sim=,nsim=50,
#            graphname="hyena.b",csvfile="1.csv")
# x3[[1]]
# #ant
# x4=allnets(emp.data="insecta-ant-colony1-day01.edges", sample.sim=,nsim=50,
#            graphname="antd1",csvfile="1.csv")
# x4[[1]]
# #weaver
# x5=allnets(emp.data="aves-weaver-social.edges", sample.sim=,nsim=50,
#            graphname="weaver.social",csvfile="weaver.csv")
# x5[[1]]
# #kangaroo
# x6=allnets(emp.data="mammalia-kangaroo-interactions.edges",sample.sim=,nsim=50,
#            graphname="kangaroo",csvfile="kangaroo.csv")
# x6[[1]]
# #Asianelephant
# x7=allnets(emp.data="mammalia-asianelephant.edges", sample.sim=,nsim=50,
#            graphname="asianeleph",csvfile="elephant.csv")
# x7[[1]]
# #chesapeake
# x8=allnets(emp.data="road-chesapeake.edges",sample.sim=,nsim=50,
#            graphname="road.chesap",csvfile="chesap.csv")
# x8[[1]]
# #bat
# x9=allnets(emp.data="mammalia-bat-roosting-indiana.edges", sample.sim=,nsim=50,
#            graphname="bat",csvfile="bat.csv")
# x9[[1]]
# #tortoise
# x10=allnets(emp.data="reptilia-tortoise-network-fi.edges",sample.sim=,nsim=50, graphname="tortoise",csvfile="tortoise.csv")
# x10[[1]]
# ```
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                       Empirical test
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

```{r}
#| echo: false
#| include: false
#| code-fold: true


#########---------Turning traininga dn test set to matrix------####
x.train=as.matrix(train[,-c(1,2)])#train and test sets excluding the response
x.test=as.matrix(test[,-c(1,2)])

y.train.1= train$R #train and test sets of the response
y.test.1=test$R

y.train.2= cbind(train$R,train$p) #train and test sets of the response
y.test.2=cbind(test$R,test$p)

##-------set controls
train.control <- trainControl(method = "cv", number = 10,savePredictions = "all")
## List of potential Lambda values
lambda.vect=10^seq(10, -2, length = 500)

set.seed(29394)
####---Ridge---regression--one--outcome
ridge=cv.glmnet(x.train,y.train.1,type.measure = "mse",alpha=0,family="gaussian")
ridgepredict=predict(ridge,s=ridge$lambda.1se,newx=x.test)
mean((y.test.1-ridgepredict)^2)

ridge1.dat=rbind(head(y.test.1),c(head(ridgepredict)))
row.names(ridge1.dat)=c("ActualVals","PredictedVals")
ridge1.dat

####---Ridge---regression--two--outcomes
ridge2=cv.glmnet(x.train,y.train.2,type.measure = "mse",alpha=0,family="mgaussian")
ridgepredict2=predict(ridge2,s=ridge2$lambda.1se,newx=x.test)

ridge2.dat=rbind(head(as.data.frame(y.test.2)[,1]),c(head(as.data.frame(ridgepredict2)[,1])))
row.names(ridge2.dat)=c("ActualVals","PredictedVals")
ridge2.dat

#### mse for Radius
mean((as.data.frame(y.test.2)[,1]-as.data.frame(ridgepredict2)[,1])^2)
#### mse for p
mean((as.data.frame(y.test.2)[,2]-as.data.frame(ridgepredict2)[,2])^2)

#### RMSE for Radius
caret::RMSE(pred = as.data.frame(ridgepredict2)[,1], obs =as.data.frame(y.test.2)[,1])
#### RMSE for p
caret::RMSE(pred = as.data.frame(ridgepredict2)[,2], obs =as.data.frame(y.test.2)[,2])
#### R2for Radius
caret::R2(pred = as.data.frame(ridgepredict2)[,1], obs =as.data.frame(y.test.2)[,1])
#### R2 for p
caret::R2(pred = as.data.frame(ridgepredict2)[,2], obs =as.data.frame(y.test.2)[,2])

####---lasso---regression--for--one--outcome--###
set.seed(29394)
las=cv.glmnet(x.train,y.train.1,type.measure = "mse",alpha=1,family="gaussian")
laspredict=predict(las,s=las$lambda.1se,newx=x.test)
mean((y.test.1-laspredict)^2)

####---lasso---regression--for--two--outcome--###
las2=cv.glmnet(x.train,y.train.2,type.measure = "mse",alpha=1,family="mgaussian")
laspredict2=predict(las2,s=las2$lambda.1se,newx=x.test)
#### mse for Radius
mean((as.data.frame(y.test.2)[,1]-as.data.frame(laspredict2)[,1])^2)
#### mse for p
mean((as.data.frame(y.test.2)[,2]-as.data.frame(laspredict2)[,2])^2)

#### RMSE for Radius
caret::RMSE(pred = as.data.frame(laspredict2)[,1], obs =as.data.frame(y.test.2)[,1])
#### RMSE for p
caret::RMSE(pred = as.data.frame(laspredict2)[,2], obs =as.data.frame(y.test.2)[,2])
#### R2for Radius
caret::R2(pred = as.data.frame(laspredict2)[,1], obs =as.data.frame(y.test.2)[,1])
#### R2 for p
caret::R2(pred = as.data.frame(laspredict2)[,2], obs =as.data.frame(y.test.2)[,2])


####++++Elastic----net--for--one--and--two---outcomes-------#####
####
set.seed(28443)
Models=function(){
  models.list=list();res=data.frame()
  models.list2=list();res2=data.frame()
  allresults=NULL
  for (j in 0:10){
    models.names=paste0("alpha",j/10)
    ####### 1 outocme modes
    models.list[[models.names]]=cv.glmnet(x.train,y.train.1,type.measure = "mse",alpha=j/10,
                                          family="gaussian")
    ####### 2 outocmes modes
    models.names2=paste0("alpha",j/10)
    models.list2[[models.names2]]=cv.glmnet(x.train,y.train.2,type.measure = "mse",alpha=j/10,
                                            family="mgaussian")
  }
  
  for (j in 0:10){
    models.names=paste0("alpha",j/10)
    models.names2=paste0("alpha",j/10)
    pred.models=predict(models.list[[models.names]],s=models.list[[models.names]]$lambda.1se,newx=x.test)
    
    pred.models2=predict(models.list2[[models.names2]],s=models.list2[[models.names2]]$lambda.1se,newx=x.test)
    
    #### mse of Radius for 1 outcome
    MSE.1=mean((y.test.1-pred.models)^2)
    #### mse of Radius for 2 outcomes
    MSE.2=mean((as.data.frame(y.test.2)[,1]-as.data.frame(pred.models2)[,1])^2)
    #### RMSE of Radius for 1 outcome 
    RMSE.1=caret::RMSE(pred = pred.models, obs =y.test.1)
    #### RMSE of Radius for 2 outcomes
    RMSE.2=caret::RMSE(pred = as.data.frame(pred.models2)[,1], obs =as.data.frame(y.test.2)[,1])
    #### RSquared of Radius for 1 outcome
    Rsquared.1=caret::R2(pred = pred.models, obs =y.test.1)
    #### RSquared of Radius for 2 outcomes
    Rsquared.2=caret::R2(pred = as.data.frame(pred.models2)[,1], obs =as.data.frame(y.test.2)[,1])
    ### final model for 1 outcome
    mod=data.frame(alpha=j/10,mse=MSE.1,RMSE=RMSE.1,Rsquared=Rsquared.1,models.names=models.names)
    ### final model for 2 outcomes
    mod2=data.frame(alpha=j/10,mse=MSE.2,RMSE=RMSE.2,Rsquared=Rsquared.2,models.names2=models.names2)
    #### result for 1 outcome
    res=rbind(res,mod)
    ### result for two outcome
    res2=rbind(res2,mod2)
    
    allresults$result1=res
    allresults$result2=res2
    # allres=list(res,res2)
    
  }
  return(allresults)
}

set.seed(8362)
output=Models()
###--minimum--mse--for--1--outcome
min(output$result1$mse)

###--minimum--mse--for--2--outcomes
min(output$result2$mse)
###--show--the--alpha--with--the--min--mse--for--1--outcome
minalpha1=output$result1[output$result1$mse==min(output$result1$mse),]
minalpha1

###--show--the--alpha--with--the--min--mse--for--2--outcomes
minalpha2=output$result2[output$result2$mse==min(output$result2$mse),]
minalpha2

###Ue minimum alphas to perform elastic net
set.seed(29394)
####---elastic---regression--for--1--outcome
Elasticnet=cv.glmnet(x.train,y.train.1,type.measure = "mse",alpha=minalpha1$alpha,family="gaussian")
Elasticnet.predict=predict(Elasticnet,s=Elasticnet$lambda.1se,newx=x.test)
mean((y.test.1-Elasticnet.predict)^2)
####---elastic---regression--for--2--outcomes
Elasticnet2=cv.glmnet(x.train,y.train.2,type.measure = "mse",alpha=minalpha2$alpha,family="mgaussian")
Elasticnet.predict2=predict(Elasticnet2,s=Elasticnet2$lambda.1se,newx=x.test)
mean((y.test.2[,1]-Elasticnet.predict2)^2)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###-----------------------Random--forest-----------#####
### Tuning the number of mtrys
set.seed((6283))
tn=tuneRF(train[,-1],train$R,stepFactor=0.5,plot = TRUE,ntreeTry = 500,
          trace = TRUE,improve = 0.01) #works only if mtry
#> the number of variables(features). This is because mtry is the number of randomly sampled variable as candidate at each split. Base case use mtry=2  when this happens

tn=as.data.frame(tn)
tn.min=tn$mtry[tn$OOBError== min(tn$OOBError)] 

set.seed(6748)
rfmodel <- randomForest::randomForest(R~ ., data = train[,-2],
                                      trainControl=train.control,
                                      mtry=tn.min,ntree=500,proximity=TRUE,importance=TRUE)
print(rfmodel)
####----Prediction---train------
pred.train.rf=predict(rfmodel,train)
####----Prediction---test------
pred.test.rf=predict(rfmodel,test)

s5=rbind(head(test$R),c(head(pred.test.rf)))
row.names(s5)=c("ActualVals","PredictedVals")
s5
####################
obs.rf=test$R
pred.rf=pred.test.rf
dist.R.rf=data.frame(y=obs.rf,x=pred.rf)

#################################
ggplot(dist.R.rf,aes(x=x,y=y))+
  xlab("predicted")+
  ylab("actual")+
  geom_point()+
  geom_abline(color="darkblue")+
  ggtitle("Plot of actual vs predicted radius for rfmodel")

####--Test--error-MSE--test--###
test.error.rf= pred.test.rf-test$R 
####--RMSE---for--test--set
caret::RMSE(test$R,pred.test.rf)
####--RSQUARED---for--test--set
caret::R2(test$R,pred.test.rf)
####----Error--Rate----#########
plot(rfmodel)

##Random-forest-model-2
set.seed(6748)
rfmodel.2 <- randomForest::randomForest(R~ spectral_radius+centrality_eigen+closeness+minCut+transitivity+betweenness+FiedlerValue+modularity+Normalized_FiedlerValue+diameter, data = train[,-2],
                                        trainControl=train.control,
                                        mtry=tn.min,ntree=500,proximity=TRUE,importance=TRUE)
print(rfmodel.2)
# ####----Prediction---train------
pred.train.rf2=predict(rfmodel.2,train)
# ####----Prediction---test------
pred.test.rf2=predict(rfmodel.2,test)
# 
# s6=rbind(head(test$R),c(head(pred.test.rf2)))
# row.names(s6)=c("ActualVals","PredictedVals")
# s6
# 
# ####--Test--error-MSE--test--###
# test.error.rf2= pred.test.rf2-test$R 
# ####--RMSE---for--test--set
# caret::RMSE(test$R,pred.test.rf2)
# ####--RSQUARED---for--test--set
# caret::R2(test$R,pred.test.rf2)
# ####----Error--Rate----#########
# plot(rfmodel.2)
```
