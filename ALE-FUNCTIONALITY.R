library(ALEPlot)

N=500
x1 <- runif(N, min=0, max=1)
x2 <- runif(N, min=0, max=1)
x3 <- runif(N, min=0, max=1)

y = x1 + 2*x2^2 + rnorm(N, 0, 0.1)
DAT = data.frame(y, x1, x2, x3)
lm.DAT = lm(y ~ .^2 + I(x1^2) + I(x2^2) + I(x3^2), DAT)
## Define the predictive function (easy in this case, since \code{lm} has a built-in
## predict function that suffices)
yhat <- function(X.model, newdata) as.numeric(predict(X.model, newdata))
## Calculate and plot the ALE main and second-order interaction effects of x1, x2, x3
par(mfrow = c(2,3))
ALE.1=ALEPlot(DAT[,2:4], lm.DAT, pred.fun=yhat, J=1, K=50, NA.plot = TRUE)
ALE.2=ALEPlot(DAT[,2:4], lm.DAT, pred.fun=yhat, J=2, K=50, NA.plot = TRUE)
ALE.3=ALEPlot(DAT[,2:4], lm.DAT, pred.fun=yhat, J=3, K=50, NA.plot = TRUE)
ALE.12=ALEPlot(DAT[,2:4], lm.DAT, pred.fun=yhat, J=c(1,2), K=20, NA.plot = TRUE)
ALE.13=ALEPlot(DAT[,2:4], lm.DAT, pred.fun=yhat, J=c(1,3), K=20, NA.plot = TRUE)
ALE.23=ALEPlot(DAT[,2:4], lm.DAT, pred.fun=yhat, J=c(2,3), K=20, NA.plot = TRUE)



library(tidyverse)
library(tidymodels)
library(broom)
library(iml)

df=read.csv("data.csv",header = T,sep = ",")
df = df[sample(1:nrow(df)), ]##shuffle row indices and randomly re order 
df%>%dplyr::select(c(GraphName,order,edges,
                     mean_eccentr,mean_path_length,graph_energy,
                     modularity,diameter,betw_centr,transitivity,
                     spectral_radius,eigen_centr,deg_centr,
                     mean_degree,minCut,FiedlerValue,Normalized_FiedlerValue,
                     closeness_centr))%>%
  mutate_if(is.character,factor)

df=df%>%na.omit()%>%
  as_tibble()%>%
  clean_names()%>%
  select(-x)%>%
  mutate_if(is.character,factor)

# Start making some non-interpretable models 
rf.model <- rand_forest() %>% 
  set_mode("classification") %>% 
  set_engine("randomForest")

xgb.model <- boost_tree() %>% 
  set_mode("classification") %>% 
  set_engine("xgboost")

#Train models
set.seed(42)
rf.model <- rf.model %>% 
  fit(graph_name~., data = df)

xgb.model <- xgb.model %>% 
  fit(graph_name~., data = df)


x.data <- df %>% select(-graph_name) %>% as.data.frame()

#high_predictor <- Predictor$new(rfm, data = data, type = "prob", class="High")

rf.predictor <- Predictor$new(rf.model, data = x.data, y = df$graph_name)#, y = df$graph_name, prob=)

xgb.predictor <- Predictor$new(xgb.model, data = x.data, y = df$graph_name)

rf.feature.imp <- FeatureImp$new(rf.predictor, loss = "ce")

xgb.feature.imp <- FeatureImp$new(xgb.predictor, loss = "ce")

rf.feature.imp %>% plot()

xgb.feature.imp %>% plot()


###----ALE PLOTS----

rf.predictor.class <- Predictor$new(rf.model, data = x.data,
                              type = "prob")

xgb.predictor.class <- Predictor$new(xgb.model, data = df,y = df$graph_name,
                              type = "prob")

# ## modularity variable 
modularity.preds.rf <- FeatureEffect$new(rf.predictor.class,
                                   feature = "modularity",
                                   grid.size = 30)

modularity.preds.rf%>%
  plot()

modularity.preds.xgb <- FeatureEffect$new(xgb.predictor.class,
                                      feature = "modularity",
                                      grid.size = 30)

modularity.preds.xgb %>%
  plot()

### transitivity variable 
transitivity.preds.rf <- FeatureEffect$new(rf.predictor.class,
                                           feature = "transitivity",
                                           grid.size = 30)

transitivity.preds.rf%>%
  plot()

transitivity.preds.xgb <- FeatureEffect$new(xgb.predictor.class,
                                            feature = "transitivity",
                                            grid.size = 30)

transitivity.preds.xgb%>%
  plot()
 

### normalized_fiedler_value variable 
normalized_fiedler_value.preds.rf <- FeatureEffect$new(rf.predictor.class,
                                                       feature = "normalized_fiedler_value",
                                                       grid.size = 30)

normalized_fiedler_value.preds.rf%>%
  plot()

normalized_fiedler_value.preds.xgb <- FeatureEffect$new(xgb.predictor.class,
                                                        feature = "normalized_fiedler_value",
                                                        grid.size = 30)

normalized_fiedler_value.preds.xgb%>%
  plot()

### graph_energy variable 
graph_energy.preds.rf <- FeatureEffect$new(rf.predictor.class,
                                           feature = "graph_energy",
                                           grid.size = 30)

graph_energy.preds.rf%>%
  plot()

graph_energy.preds.xgb <- FeatureEffect$new(xgb.predictor.class,
                                            feature = "graph_energy",
                                            grid.size = 30)

graph_energy.preds.xgb%>%
  plot()


### deg_centr variable 
deg_centr.preds.rf <- FeatureEffect$new(rf.predictor.class,
                                        feature = "deg_centr",
                                        grid.size = 30)

deg_centr.preds.rf%>%
  plot()

deg_centr.preds.xgb <- FeatureEffect$new(rf.predictor.class,
                                         feature = "deg_centr",
                                         grid.size = 30)

deg_centr.preds.xgb%>%
  plot()


#############################################
#+ Interactions
#############################################
interactn1.preds<- FeatureEffect$new(xgb.predictor.class,
                            feature=c("modularity","transitivity"),
                            grid.size = 30)

interactn1.preds %>%
  plot()


#### All predictors
FeatureEffects$new(xgb.predictor.class)$plot(ncol = 2)

#rf.ale <- FeatureEffects$new(rf.predictor)
#xgb.ale <- FeatureEffects$new(xgb.predictor)

#rf.aleplot=rf.ale %>% plot()
#ggsave(rf.aleplot,width=30,height=32)

#xgb.aleplot=xgb.ale %>% plot()



# We train a random forest on the Boston dataset:
# #data("Boston", package = "MASS")
# library("rpart")
# rf <- rpart(medv ~ ., data = Boston)
# mod <- Predictor$new(rf, data = Boston)
# # Compute the accumulated local effects for the first feature
# eff <- FeatureEffect$new(mod, feature = "rm", grid.size = 30)
# eff$plot()
# 
# # FeatureEffect plots support up to two features:
# eff <- FeatureEffect$new(mod, feature = c("Sepal.Length", "Petal.Length"))
# eff$plot()
# # show where the actual data lies
# eff$plot(show.data = TRUE)
# 
# # For multiclass classification models, you can choose to only show one class:
# mod <- Predictor$new(rf, data = iris, type = "prob", class = 1)


# high_veg <- plot(FeatureEffect$new(high_predictor, feature = "veg", method = "ale")) 
# 


# # ALE plot
# low_predictor <- Predictor$new(rfm, data = data, type = "prob", class="Low")
# high_predictor <- Predictor$new(rfm, data = data, type = "prob", class="High")
# 
# ## Vegetation variable 
# low_veg <- plot(FeatureEffect$new(Low_predictor, feature = "veg", method = "ale")) 
# high_veg <- plot(FeatureEffect$new(high_predictor, feature = "veg", method = "ale")) 
# 
# plot(low_veg)
# plot(high_veg)
# 
# ## slope variable
# low_slope <- plot(FeatureEffect$new(Low_predictor, feature = "slope", method = "ale")) 
# high_slope <- plot(FeatureEffect$new(high_predictor, feature = "slope", method = "ale")) 
# 
# plot(low_slope)
# plot(high_slope)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Lime model 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
#For analyzing how models make individual predictions 
rf_lime <- LocalModel$new(rf.predictor, x.interest = x.data[1,])
rf_lime %>% plot()


xgb_lime <- LocalModel$new(xgb.predictor, x.interest = x.data[1,])
xgb_lime %>% plot()



##---The function ggscatmat from the package GGally creates a matrix with scatterplots, densities and correlations for numeric columns

df.new=df%>%select(-c(graph_name))

ggscatmat(df.new, 
          alpha=0.2,
          columns = 1:ncol(df.new),
          corMethod = "pearson")
##--[1]--Correlation plot
corell.plot=pairs.panels(df.new,
                         gap = 0,
                         bg = c("red", "yellow", "blue")[df$graph_name],
                         smooth = TRUE,
                         stars = TRUE, 
                         ci=TRUE,
                         pch=24)
corell.plot















