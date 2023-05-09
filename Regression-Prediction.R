
setwd("C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi")
source("SPATIAL-PIPELINE.R")

m_glm1 <- glm(R ~., family = Gamma(link = "log"),data=train)


m_glm2 <- glm(R ~FiedlerValue + Normalized_FiedlerValue + spectral_radius
              +centrality_eigen+
                minCut+diameter+betweenness,
              family = Gamma(link = "inverse"),data=train)


summary(m_glm2)

predeict.glm2=predict(m_glm2,test)

glm2.dataframe=rbind(tail(test$R),c(tail(predeict.glm2)))
row.names(glm2.dataframe)=c("ActualVals","PredictedVals")
glm2.dataframe


### Beta Regression
model.beta1 = betareg(R ~FiedlerValue + Normalized_FiedlerValue + spectral_radius
                      +centrality_eigen+transitivity+closeness+
                        minCut+diameter+betweenness, data=train)

summary(model.beta1)
##--Beta--prediction--with--Model--1
model.beta.predict1=predict(model.beta1,test)

mean((test$R-model.beta.predict1)^2)

##--Snipet--of--actual--and--predicted--values
betamodel1.dataframe=rbind(tail(test$R),c(tail(model.beta.predict1)))
row.names(betamodel1.dataframe)=c("ActualVals","PredictedVals")
betamodel1.dataframe



# Find the best lambda using cross-validation
x.train=as.matrix(train[,-1])#train and test sets excluding the response
x.test=as.matrix(test[,-1])

class(y.train)

y.train= train$R #train and test sets of the response
y.test=test$R


#### Lasso Model
set.seed(123) 
cv.lasso <- cv.glmnet(x.train, log(y.train), alpha = 1, family = "gaussian")
# Fit the final model on the training data
model.lasso <- glmnet(x.train, log(y.train), alpha = 1, family = "gaussian",
                lambda = cv.lasso$lambda.1se)
# Display regression coefficients
coef(model.lasso)
# Make predictions on the test data
lasso.predict <- model.lasso %>% predict( cv.lasso$lambda.1se,newx = x.test)
head(c(exp(lasso.predict)))
head(y.test)

#### Ridge Model
set.seed(123) 
cv.ridge <- cv.glmnet(x.train, log(y.train), alpha = 0, family = "gaussian")
# Fit the final model on the training data
model.ridge <- glmnet(x.train, log(y.train), alpha = 0, family = "gaussian",
                      lambda = cv.ridge$lambda.1se)
# Display regression coefficients
coef(model.ridge)
# Make predictions on the test data
ridge.predict <- model.ridge %>% predict(cv.ridge$lambda.min,newx = x.test)
head(c(exp(ridge.predict)))
head(y.test)

###Elastic Net
set.seed(28443)
ENET.Model=function(){
  models.list=list();res=data.frame()
  for (j in 0:10){
    models.names=paste0("alpha",j/10)
    models.list[[models.names]]=cv.glmnet(x.train,log(y.train),type.measure = "mse",alpha=j/10,
                                          family="gaussian")
    
  }
  for (j in 0:10){
    models.names=paste0("alpha",j/10)
    pred.logmodels=predict(models.list[[models.names]],s=models.list[[models.names]]$lambda.1se,newx=x.test)
    pred.models=exp(pred.logmodels)
    
    #### MSE of Radius for 1 outcome
    MSE=mean((y.test-pred.models)^2)
    
    #### RMSE of Radius for 1 outcome 
    RMSE=caret::RMSE(pred = pred.models, obs =y.test)
    
    #### RSquared of Radius for 1 outcome
    Rsquared=caret::R2(pred = pred.models, obs =y.test)
    
    ### final model for 1 outcome
    mod=data.frame(alpha=j/10,mse=MSE,RMSE=RMSE,Rsquared=Rsquared,models.names=models.names)
    
    #### result for 1 outcome
    res=rbind(res,mod)
  }
  return(res)
}

set.seed(8362)
output=ENET.Model()
###--minimum--mse--for--1--outcome
min(output$mse)
###--show--the--alpha--with--the--min--mse--for--1--outcome
minalpha=output[output$mse==min(output$mse),]
minalpha

###--Use minimum alphas to perform elastic net
set.seed(29394)
####---Elastic---Model--for--1--outcome
cv.elasticnet=cv.glmnet(x.train,log(y.train),type.measure = "mse",alpha=minalpha$alpha,family="gaussian")

# Fit the final model on the training data
model.enet <- glmnet(x.train, log(y.train),alpha=minalpha$alpha, family = "gaussian",
                      lambda = cv.ridge$lambda.1se)
# Make predictions on the test data
enet.predict <- model.enet  %>% predict(s=cv.elasticnet$lambda.1se,newx = x.test)
head(c(exp(enet.predict)))
head(y.test)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+             Prediction on empirical data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#"aves-weaver-social.edges",
#"insecta-ant-colony1-day01.edges",
#"mammalia-macaque-contact-sits.edges"
# mammalia-asianelephant.edges
# mammalia-dolphin-social.edges
# mammalia-hyena-networkb.edges
# mammalia-kangaroo-interactions.edges
# mammalia-bat-roosting-indiana.edges
# reptilia-tortoise-network-fi.edges
# road-chesapeake.edges
emp.data.1="insecta-ant-colony1-day01.edges"
graphname.1="emp"
graph.data.1<- read.table(emp.data.1)  
graph.new.1=graph_from_data_frame(as.matrix(graph.data.1),directed=FALSE)
g.graph.1=igraph::simplify(graph.new.1,remove.multiple = T,remove.loops = T)
components = igraph::clusters(g.graph.1, mode="weak")
biggest_cluster_id = which.max(components$csize)
# # ids
vert_ids = V(g.graph.1)[components$membership== biggest_cluster_id]
# # subgraph
G=igraph::induced_subgraph(g.graph.1, vert_ids)
inducedgraph.1=G
#g.1=list(g.graph.1)
#inducedgraph.1=induced.graph(g.1)

inducedgraph.1$name=graphname.1
inducedgraph.1$type=graphname.1
inducedgraph.1$id="1"
G.1=list(inducedgraph.1)



#++++++++++++++++++++++++Graph--Features--for-empirical-network++++++++++++++++++++

emp_graphfeat.1=RunSimOnGraphFeatures(G.1,nreps = 1)
input.data.1=emp_graphfeat.1 %>%dplyr::select(c(order,edges,mean_degree,minCut,
                                                FiedlerValue,Normalized_FiedlerValue,closeness,modularity,diameter,betweenness,
                                                transitivity,spectral_radius,centrality_eigen))

lasso.predict=predict(model.lasso,s=cv.lasso$lambda.1se,as.matrix(input.data.1))
transform.lasso.predict=exp(lasso.predict)

ridge.predict=predict(model.ridge,s=cv.ridge$lambda.1se,as.matrix(input.data.1))
transform.ridge.predict=exp(ridge.predict)

elasticnet.predict=predict(model.enet,s=cv.elasticnet$lambda.1se,as.matrix(input.data.1))
transform.elasticnet.predict=exp(elasticnet.predict)



spatial.lasso=simulate.spatial(N=vcount(inducedgraph.1),
                               radius=transform.lasso.predict,nsim=10)

spatial.ridge=simulate.spatial(N=vcount(inducedgraph.1),
                               radius=transform.ridge.predict,nsim=10)

spatial.elasticnet=simulate.spatial(N=vcount(inducedgraph.1),
                                    radius=transform.elasticnet.predict,nsim=10)

spatial.lassograph=spatial.lasso[[9]]
spatial.ridgegraph=spatial.ridge[[9]]
spatial.elasticnetgraph=spatial.elasticnet[[9]]

Nets=c(G.1,spatial.lassograph,spatial.ridgegraph,spatial.elasticnetgraph)
graphfeat=RunSimOnGraphFeatures(Nets,nreps = 1)
graphfeat$GraphName=c(graphname.1,"Sp.lasso","Sp.ridge","Sp.enet")
graphfeat

par(mfrow=c(1,3))
plot(inducedgraph.1,vertex.label=NA,vertex.size=2)
plot(spatial.lassograph[[1]],vertex.label=NA,vertex.size=2)
plot(spatial.ridgegraph[[1]],vertex.label=NA,vertex.size=2)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# betareg caret
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# library(caret)
# library(betareg)
# 
# # create a betareg caret model
# # type and library
# betaregression <- list(type='Regression',
#                        library='betareg',
#                        loop=NULL)
# 
# # parameters to tune
# prm <- data.frame(parameter=c("link", "type", "link.phi"),
#                   class= rep("character", 3))
# 
# # add to the model
# betaregression$parameters <- prm
# 
# # grid search, this is the default grid search, user can specify otherwise
# # creates 54 separate models, so if looking to speed up try fewer params in grid
# betaGrid  <- function(x, y, len=NULL, search="grid"){
#   if(search == "grid"){
#     out <- expand.grid(link=c("logit", "probit", "cloglog", "cauchit", "log", "loglog"),
#                        type=c("ML", "BC", "BR"),
#                        link.phi=c("identity", "log", "sqrt"), stringsAsFactors = F) # here force the strings as character,
#     # othewise get error that the model arguments
#     # were expecting 'chr' when fitting
#   }
#   out
# }
# 
# # add the grid search
# betaregression$grid <- betaGrid
# 
# # create the fit
# betaFit <- function(x, y, wts, param, lev, last, weights, classProbs, ...){
# 
#   dat <- if(is.data.frame(x)) x else as.data.frame(x)
#   dat$.outcome <- y
# 
#   theDots <- list(...)
# 
#   modelArgs <- c(list(formula = as.formula(".outcome ~ ."), data = dat, link=param$link, type=param$type), theDots)
# 
#   out <- do.call(betareg::betareg, modelArgs)
#   out$call <- NULL
#   out
# }
# 
# # betaregression fit
# betaregression$fit <- betaFit
# 
# # predict element
# betaPred <- function(modelFit, newdata, preProc=NULL, submodels=NULL){
#   if(!is.data.frame(newdata)) newdata <- as.data.frame(newdata)
#   betareg::predict(modelFit, newdata)
# }
# 
# # add the predict method
# betaregression$predict <- betaPred
# 
# # regression, no probabities calculated
# # just assigning NULL didnt work for some reason
# # wrapped in a function instead
# betaProb <- function(){
#   return(NULL)
# }
# betaregression$prob <- betaProb
# 
# 
# # test it on a dataset
# #data('GasolineYield', package = 'betareg')
# 
# # 10 fold cross validation
# fitControl <- trainControl(method='repeatedcv', number = 10, repeats = 5)
# 
# # betaregression, takes a min or so with full grid
# betareg <- train(R~FiedlerValue + Normalized_FiedlerValue + spectral_radius+
#                  centrality_eigen+transitivity+closeness+
#                    minCut+diameter+betweenness,
#                  data=train, method=betaregression,
#                  trControl = fitControl)
# 
# #look at output
# betareg

