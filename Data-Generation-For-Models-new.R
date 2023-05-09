####Simulate networks-with-R-and-P
setwd("C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi-New")
source("SPATIAL-PIPELINE.R")


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                    Erdos Renyi  data

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MKERGraphs<- function(n=25,p=0.8) {
  Graphs = list()#;G=list() # set up the list of Graphs first
  initgraph= list()
  i= 1
  print("Creating ER Graphs")
  Graphs[[i]] = erdos.renyi.game(n,p,type = c("gnp"),directed = F,loops = F)

  Graphs[[i]]$type = "ER"
  Graphs[[i]]$id = "1"
  Graphs[[i]]$name="ER"
  i <- i+1
  return(Graphs)
}


sim.erdos.net=function(p=c(0.3,0.4,0.6,0.05),n=100){
  net=NULL;data=NULL
  k=1
  for (i in 1:length(p)){
    net[i]= MKERGraphs(n,p[i])
    data=RunSimOnGraphFeatures(net,nreps = 1)
    k=k+1
  }
  df=cbind(p,data)
  return(net)
}

### 50 nodes 
ER1.50=sim.erdos.net(p=rep(0.001,250),n=50)
write.csv(ER1.50,"ER1.50.csv")
ER2.50=sim.erdos.net(p=rep(0.005,250),n=50)
write.csv(ER2.50,"ER2.50.csv")
ER3.50=sim.erdos.net(p=rep(0.009,250),n=50)
write.csv(ER3.50,"ER3.50.csv")
ER4.50=sim.erdos.net(p=rep(0.01,250),n=50)
write.csv(ER4.50,"ER4.50.csv")
ER5.50=sim.erdos.net(p=rep(0.05,250),n=50)
write.csv(ER5.50,"ER5.50.csv")
ER6.50=sim.erdos.net(p=rep(0.09,250),n=50)
write.csv(ER6.50,"ER6.50.csv")
ER7.50=sim.erdos.net(p=rep(0.1,250),n=50)
write.csv(ER7.50,"ER7.50.csv")
ER8.50=sim.erdos.net(p=rep(0.5,250),n=50)
write.csv(ER8.50,"ER8.50.csv")
######## Saving all ER graphs 50 nodes data to csv
ER1=read.csv("ER1.50.csv");ER2=read.csv("ER2.50.csv");ER3=read.csv("ER3.50.csv")
ER4=read.csv("ER4.50.csv");ER5=read.csv("ER5.50.csv");ER6=read.csv("ER6.50.csv")
ER7=read.csv("ER7.50.csv");ER8=read.csv("ER8.50.csv")

df.er.50=rbind(ER1,ER2,ER3,ER4,ER5,ER6,ER7,ER8)
write.csv(df.er.50,"df.er.50.csv")

# dim(df.er.50)
# node.num=50
# df.er.50<-df.er.50 %>%dplyr::filter(order==node.num)
# dim(df.er.50)
# write.csv(df.er.50,"df.50.csv")


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

networks.spatial=function(R=c(0.3,0.4,0.6,0.05),n=100){
  net=NULL;data=NULL
  k=1
  for (i in 1:length(R)){
    net[i]= makeSpatialGraphs(node.size =n,Radius = R[i])
    data=RunSimOnGraphFeatures(net,nreps = 1)
    k=k+1
  }
  df=cbind(R,data)
  return(df)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#      Data for spatial network
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### 50 nodes 
a1.50=networks.spatial(R=rep(0.15,200),n=50)
write.csv(a1.50,"a1.50.csv")
a2.50=networks.spatial(R=rep(0.2,200),n=50)
write.csv(a2.50,"a2.50.csv")
a3.50=networks.spatial(R=rep(0.3,200),n=50)
write.csv(a3.50,"a3.50.csv")
a4.50=networks.spatial(R=rep(0.4,200),n=50)
write.csv(a4.50,"a4.50.csv")
a5.50=networks.spatial(R=rep(0.5,200),n=50)
write.csv(a5.50,"a5.50.csv")
a6.50=networks.spatial(R=rep(0.6,200),n=50)
write.csv(a6.50,"a6.50.csv")
a7.50=networks.spatial(R=rep(0.7,200),n=50)
write.csv(a7.50,"a7.50.csv")
a8.50=networks.spatial(R=rep(0.8,200),n=50)
write.csv(a8.50,"a8.50.csv")
a9.50=networks.spatial(R=rep(0.9,200),n=50)
write.csv(a9.50,"a9.50.csv")
a10.50=networks.spatial(R=rep(0.99,200),n=50)
write.csv(a10.50,"a10.50.csv")
a11.50=networks.spatial(R=rep(0.05,200),n=50)
write.csv(a11.50,"a11.50.csv")
a12.50=networks.spatial(R=rep(0.01,200),n=50)
write.csv(a12.50,"a12.50.csv")
a13.50=networks.spatial(R=rep(0.08,200),n=50)
write.csv(a13.50,"a13.50.csv")

######## Saving all 50 nodes data to csv
a1=read.csv("a1.50.csv");a2=read.csv("a2.50.csv");a3=read.csv("a3.50.csv")
a4=read.csv("a4.50.csv");a5=read.csv("a5.50.csv");a6=read.csv("a6.50.csv")
a7=read.csv("a7.50.csv");a8=read.csv("a8.50.csv");a9=read.csv("a9.50.csv")
a10=read.csv("a10.50.csv")

df.50=rbind(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
dim(df.50)
node.num=50
df.50<-df.50 %>%dplyr::filter(order==node.num)
dim(df.50)
write.csv(df.50,"df.50.csv")


##100 nodes
a1.100=networks.spatial(R=rep(0.15,200),n=100)
write.csv(a1.100,"a1.100.csv")
a2.100=networks.spatial(R=rep(0.2,200),n=100)
write.csv(a2.100,"a2.100.csv")
a3.100=networks.spatial(R=rep(0.3,200),n=100)
write.csv(a3.100,"a3.100.csv")
a4.100=networks.spatial(R=rep(0.4,200),n=100)
write.csv(a4.100,"a4.100.csv")
a5.100=networks.spatial(R=rep(0.5,200),n=100)
write.csv(a5.100,"a5.100.csv")
a6.100=networks.spatial(R=rep(0.6,200),n=100)
write.csv(a6.100,"a6.100.csv")
a7.100=networks.spatial(R=rep(0.7,200),n=100)
write.csv(a7.100,"a7.100.csv")
a8.100=networks.spatial(R=rep(0.8,200),n=100)
write.csv(a8.100,"a8.100.csv")
a9.100=networks.spatial(R=rep(0.9,200),n=100)
write.csv(a9.100,"a9.100.csv")
a10.100=networks.spatial(R=rep(0.99,200),n=100)
write.csv(a10.100,"a10.100.csv")
a11.100=networks.spatial(R=rep(0.05,200),n=100)
write.csv(a11.100,"a11.100.csv")
a12.100=networks.spatial(R=rep(0.01,200),n=100)
write.csv(a12.100,"a12.100.csv")
a13.100=networks.spatial(R=rep(0.08,200),n=100)
write.csv(a13.100,"a13.100.csv")

######## Saving all 100 nodes data to csv
a1=read.csv("a1.100.csv");a2=read.csv("a2.100.csv");a3=read.csv("a3.100.csv")
a4=read.csv("a4.100.csv");a5=read.csv("a5.100.csv");a6=read.csv("a6.100.csv")
a7=read.csv("a7.100.csv");a8=read.csv("a8.100.csv");a9=read.csv("a9.100.csv")
a10=read.csv("a10.100.csv")

df.100=rbind(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
dim(df.100)
node.num=100
df.100<-df.100 %>%dplyr::filter(order==node.num)
dim(df.100)
write.csv(df.100,"df.100.csv")

##150 nodes
a1.150=networks.spatial(R=rep(0.15,200),n=150)
write.csv(a1.150,"a1.150.csv")
a2.150=networks.spatial(R=rep(0.2,200),n=150)
write.csv(a2.150,"a2.150.csv")
a3.150=networks.spatial(R=rep(0.3,200),n=150)
write.csv(a3.150,"a3.150.csv")
a4.150=networks.spatial(R=rep(0.4,200),n=150)
write.csv(a4.150,"a4.150.csv")
a5.150=networks.spatial(R=rep(0.5,200),n=150)
write.csv(a5.150,"a5.150.csv")
a6.150=networks.spatial(R=rep(0.6,200),n=150)
write.csv(a6.150,"a6.150.csv")
a7.150=networks.spatial(R=rep(0.7,200),n=150)
write.csv(a7.150,"a7.150.csv")
a8.150=networks.spatial(R=rep(0.8,200),n=150)
write.csv(a8.150,"a8.150.csv")
a9.150=networks.spatial(R=rep(0.9,200),n=150)
write.csv(a9.150,"a9.150.csv")
a10.150=networks.spatial(R=rep(0.99,200),n=150)
write.csv(a10.150,"a10.150.csv")
a11.150=networks.spatial(R=rep(0.05,200),n=150)
write.csv(a11.150,"a11.150.csv")
a12.150=networks.spatial(R=rep(0.01,200),n=150)
write.csv(a12.150,"a12.150.csv")
a13.150=networks.spatial(R=rep(0.08,200),n=150)
write.csv(a13.150,"a13.150.csv")

######## Saving all 150 nodes data to csv
a1=read.csv("a1.150.csv");a2=read.csv("a2.150.csv");a3=read.csv("a3.150.csv")
a4=read.csv("a4.150.csv");a5=read.csv("a5.150.csv");a6=read.csv("a6.150.csv")
a7=read.csv("a7.150.csv");a8=read.csv("a8.150.csv");a9=read.csv("a9.150.csv")
a10=read.csv("a10.150.csv")

df.150=rbind(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
dim(df.150)
node.num=150
df.150<-df.150 %>%dplyr::filter(order==node.num)
dim(df.150)
write.csv(df.150,"df.150.csv")


##150 nodes
a1.150=networks.spatial(R=rep(0.15,200),n=150)
write.csv(a1.150,"a1.150.csv")
a2.150=networks.spatial(R=rep(0.2,200),n=150)
write.csv(a2.150,"a2.150.csv")
a3.150=networks.spatial(R=rep(0.3,200),n=150)
write.csv(a3.150,"a3.150.csv")
a4.150=networks.spatial(R=rep(0.4,200),n=150)
write.csv(a4.150,"a4.150.csv")
a5.150=networks.spatial(R=rep(0.5,200),n=150)
write.csv(a5.150,"a5.150.csv")
a6.150=networks.spatial(R=rep(0.6,200),n=150)
write.csv(a6.150,"a6.150.csv")
a7.150=networks.spatial(R=rep(0.7,200),n=150)
write.csv(a7.150,"a7.150.csv")
a8.150=networks.spatial(R=rep(0.8,200),n=150)
write.csv(a8.150,"a8.150.csv")
a9.150=networks.spatial(R=rep(0.9,200),n=150)
write.csv(a9.150,"a9.150.csv")
a10.150=networks.spatial(R=rep(0.99,200),n=150)
write.csv(a10.150,"a10.150.csv")
a11.150=networks.spatial(R=rep(0.05,200),n=150)
write.csv(a11.150,"a11.150.csv")
a12.150=networks.spatial(R=rep(0.01,200),n=150)
write.csv(a12.150,"a12.150.csv")
a13.150=networks.spatial(R=rep(0.08,200),n=150)
write.csv(a13.150,"a13.150.csv")

######## Saving all 150 nodes data to csv
a1=read.csv("a1.150.csv");a2=read.csv("a2.150.csv");a3=read.csv("a3.150.csv")
a4=read.csv("a4.150.csv");a5=read.csv("a5.150.csv");a6=read.csv("a6.150.csv")
a7=read.csv("a7.150.csv");a8=read.csv("a8.150.csv");a9=read.csv("a9.150.csv")
a10=read.csv("a10.150.csv")

df.150=rbind(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
dim(df.150)
node.num=150
df.150<-df.150 %>%dplyr::filter(order==node.num)
dim(df.150)
write.csv(df.150,"df.150.csv")


##200 nodes
a1.200=networks.spatial(R=rep(0.15,200),n=200)
write.csv(a1.200,"a1.200.csv")
a2.200=networks.spatial(R=rep(0.2,200),n=200)
write.csv(a2.200,"a2.200.csv")
a3.200=networks.spatial(R=rep(0.3,200),n=200)
write.csv(a3.200,"a3.200.csv")
a4.200=networks.spatial(R=rep(0.4,200),n=200)
write.csv(a4.200,"a4.200.csv")
a5.200=networks.spatial(R=rep(0.5,200),n=200)
write.csv(a5.200,"a5.200.csv")
a6.200=networks.spatial(R=rep(0.6,200),n=200)
write.csv(a6.200,"a6.200.csv")
a7.200=networks.spatial(R=rep(0.7,200),n=200)
write.csv(a7.200,"a7.200.csv")
a8.200=networks.spatial(R=rep(0.8,200),n=200)
write.csv(a8.200,"a8.200.csv")
a9.200=networks.spatial(R=rep(0.9,200),n=200)
write.csv(a9.200,"a9.200.csv")
a10.200=networks.spatial(R=rep(0.99,200),n=200)
write.csv(a10.200,"a10.200.csv")
a11.200=networks.spatial(R=rep(0.05,200),n=200)
write.csv(a11.200,"a11.200.csv")
a12.200=networks.spatial(R=rep(0.01,200),n=200)
write.csv(a12.200,"a12.200.csv")
a13.200=networks.spatial(R=rep(0.08,200),n=200)
write.csv(a13.200,"a13.200.csv")

######## Saving all 200 nodes data to csv
a1=read.csv("a1.200.csv");a2=read.csv("a2.200.csv");a3=read.csv("a3.200.csv")
a4=read.csv("a4.200.csv");a5=read.csv("a5.200.csv");a6=read.csv("a6.200.csv")
a7=read.csv("a7.200.csv");a8=read.csv("a8.200.csv");a9=read.csv("a9.200.csv")
a10=read.csv("a10.200.csv")

df.200=rbind(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
dim(df.200)
node.num=200
df.200<-df.200 %>%dplyr::filter(order==node.num)
dim(df.200)
write.csv(df.200,"df.200.csv")


##250 nodes
a1.250=networks.spatial(R=rep(0.15,250),n=250)
write.csv(a1.250,"a1.250.csv")
a2.250=networks.spatial(R=rep(0.2,250),n=250)
write.csv(a2.250,"a2.250.csv")
a3.250=networks.spatial(R=rep(0.3,250),n=250)
write.csv(a3.250,"a3.250.csv")
a4.250=networks.spatial(R=rep(0.4,250),n=250)
write.csv(a4.250,"a4.250.csv")
a5.250=networks.spatial(R=rep(0.5,250),n=250)
write.csv(a5.250,"a5.250.csv")
a6.250=networks.spatial(R=rep(0.6,250),n=250)
write.csv(a6.250,"a6.250.csv")
a7.250=networks.spatial(R=rep(0.7,250),n=250)
write.csv(a7.250,"a7.250.csv")
a8.250=networks.spatial(R=rep(0.8,250),n=250)
write.csv(a8.250,"a8.250.csv")
a9.250=networks.spatial(R=rep(0.9,250),n=250)
write.csv(a9.250,"a9.250.csv")
a10.250=networks.spatial(R=rep(0.99,250),n=250)
write.csv(a10.250,"a10.250.csv")
a11.250=networks.spatial(R=rep(0.05,250),n=250)
write.csv(a11.250,"a11.250.csv")
a12.250=networks.spatial(R=rep(0.01,250),n=250)
write.csv(a12.250,"a12.250.csv")
a13.250=networks.spatial(R=rep(0.08,250),n=250)
write.csv(a13.250,"a13.250.csv")

######## Saving all 250 nodes data to csv
a1=read.csv("a1.250.csv");a2=read.csv("a2.250.csv");a3=read.csv("a3.250.csv")
a4=read.csv("a4.250.csv");a5=read.csv("a5.250.csv");a6=read.csv("a6.250.csv")
a7=read.csv("a7.250.csv");a8=read.csv("a8.250.csv");a9=read.csv("a9.250.csv")
a10=read.csv("a10.250.csv")

df.250=rbind(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
dim(df.250)
node.num=250
df.250<-df.250 %>%dplyr::filter(order==node.num)
dim(df.250)
write.csv(df.250,"df.250.csv")

###300 nodes
a1.300=networks.spatial(R=rep(0.15,200),n=300)
write.csv(a1.300,"a1.300.csv")
a2.300=networks.spatial(R=rep(0.2,200),n=300)
write.csv(a2.300,"a2.300.csv")
a3.300=networks.spatial(R=rep(0.3,200),n=300)
write.csv(a3.300,"a3.300.csv")
a4.300=networks.spatial(R=rep(0.4,200),n=300)
write.csv(a4.300,"a4.300.csv")
a5.300=networks.spatial(R=rep(0.5,200),n=300)
write.csv(a5.300,"a5.300.csv")
a6.300=networks.spatial(R=rep(0.6,200),n=300)
write.csv(a6.300,"a6.300.csv")
a7.300=networks.spatial(R=rep(0.7,200),n=300)
write.csv(a7.300,"a7.300.csv")
a8.300=networks.spatial(R=rep(0.8,200),n=300)
write.csv(a8.300,"a8.300.csv")
a9.300=networks.spatial(R=rep(0.9,200),n=300)
write.csv(a9.300,"a9.300.csv")
a10.300=networks.spatial(R=rep(0.99,200),n=300)
write.csv(a10.300,"a10.300.csv")
a11.300=networks.spatial(R=rep(0.05,200),n=300)
write.csv(a11.300,"a11.300.csv")
a12.300=networks.spatial(R=rep(0.01,200),n=300)
write.csv(a12.300,"a12.300.csv")
a13.300=networks.spatial(R=rep(0.08,200),n=300)
write.csv(a13.300,"a13.300.csv")

######## Saving all 300 nodes data to csv
a1=read.csv("a1.300.csv");a2=read.csv("a2.300.csv");a3=read.csv("a3.300.csv")
a4=read.csv("a4.300.csv");a5=read.csv("a5.300.csv");a6=read.csv("a6.300.csv")
a7=read.csv("a7.300.csv");a8=read.csv("a8.300.csv");a9=read.csv("a9.300.csv")
a10=read.csv("a10.300.csv")

df.300=rbind(a1,a2,a3,a4,a5,a6,a7,a8,a9)#,a10)
dim(df.300)
node.num=300
df.300<-df.300 %>%dplyr::filter(order==node.num)
dim(df.300)
write.csv(df.300,"df.300.csv")

##350 nodes
a1.350=networks.spatial(R=rep(0.15,200),n=350)
write.csv(a1.350,"a1.350.csv")
a2.350=networks.spatial(R=rep(0.2,200),n=350)
write.csv(a2.350,"a2.350.csv")
a3.350=networks.spatial(R=rep(0.3,200),n=350)
write.csv(a3.350,"a3.350.csv")
a4.350=networks.spatial(R=rep(0.4,200),n=350)
write.csv(a4.350,"a4.350.csv")
a5.350=networks.spatial(R=rep(0.5,200),n=350)
write.csv(a5.350,"a5.350.csv")
a6.350=networks.spatial(R=rep(0.6,200),n=350)
write.csv(a6.350,"a6.350.csv")
a7.350=networks.spatial(R=rep(0.7,200),n=350)
write.csv(a7.350,"a7.350.csv")
a8.350=networks.spatial(R=rep(0.8,200),n=350)
write.csv(a8.350,"a8.350.csv")
a9.350=networks.spatial(R=rep(0.9,200),n=350)
write.csv(a9.350,"a9.350.csv")
a10.350=networks.spatial(R=rep(0.99,200),n=350)
write.csv(a10.350,"a10.350.csv")
a11.350=networks.spatial(R=rep(0.05,200),n=350)
write.csv(a11.350,"a11.350.csv")
a12.350=networks.spatial(R=rep(0.01,200),n=350)
write.csv(a12.350,"a12.350.csv")
a13.350=networks.spatial(R=rep(0.08,200),n=350)
write.csv(a13.350,"a13.350.csv")

######## Saving all 350 nodes data to csv
a1=read.csv("a1.350.csv");a2=read.csv("a2.350.csv");a3=read.csv("a3.350.csv")
a4=read.csv("a4.350.csv");a5=read.csv("a5.350.csv");a6=read.csv("a6.350.csv")
a7=read.csv("a7.350.csv");a8=read.csv("a8.350.csv")#;a9=read.csv("a9.350.csv")
#a10=read.csv("a10.350.csv")

df.350=rbind(a1,a2,a3,a4,a5,a6,a7,a8)#,a9,a10)
dim(df.350)
node.num=350
df.350<-df.350 %>%dplyr::filter(order==node.num)
dim(df.350)
write.csv(df.350,"df.350.csv")

### 400 nodes
a1.400=networks.spatial(R=rep(0.15,200),n=400)
write.csv(a1.400,"a1.400.csv")
a2.400=networks.spatial(R=rep(0.2,200),n=400)
write.csv(a2.400,"a2.400.csv")
a3.400=networks.spatial(R=rep(0.3,200),n=400)
write.csv(a3.400,"a3.400.csv")
a4.400=networks.spatial(R=rep(0.4,200),n=400)
write.csv(a4.400,"a4.400.csv")
a5.400=networks.spatial(R=rep(0.5,200),n=400)
write.csv(a5.400,"a5.400.csv")
a6.400=networks.spatial(R=rep(0.6,200),n=400)
write.csv(a6.400,"a6.400.csv")
a7.400=networks.spatial(R=rep(0.7,200),n=400)
write.csv(a7.400,"a7.400.csv")
a8.400=networks.spatial(R=rep(0.8,200),n=400)
write.csv(a8.400,"a8.400.csv")
a9.400=networks.spatial(R=rep(0.9,200),n=400)
write.csv(a9.400,"a9.400.csv")
a10.400=networks.spatial(R=rep(0.99,200),n=400)
write.csv(a10.400,"a10.400.csv")
a11.400=networks.spatial(R=rep(0.05,200),n=400)
write.csv(a11.400,"a11.400.csv")
a12.400=networks.spatial(R=rep(0.01,200),n=400)
write.csv(a12.400,"a12.400.csv")
a13.400=networks.spatial(R=rep(0.08,200),n=400)
write.csv(a13.400,"a13.400.csv")

######## Saving all 400 nodes data to csv
a1=read.csv("a1.400.csv");a2=read.csv("a2.400.csv");a3=read.csv("a3.400.csv")
a4=read.csv("a4.400.csv");a5=read.csv("a5.400.csv");a6=read.csv("a6.400.csv")
a7=read.csv("a7.400.csv");a8=read.csv("a8.400.csv");a9=read.csv("a9.400.csv")
a10=read.csv("a10.400.csv")

df.400=rbind(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
dim(df.400)
node.num=400
df.400<-df.400 %>%dplyr::filter(order==node.num)
dim(df.400)
write.csv(df.400,"df.400.csv")

a1.450=networks.spatial(R=rep(0.15,200),n=450)
write.csv(a1.450,"a1.450.csv")
a2.450=networks.spatial(R=rep(0.2,200),n=450)
write.csv(a2.450,"a2.450.csv")
a3.450=networks.spatial(R=rep(0.3,200),n=450)
write.csv(a3.450,"a3.450.csv")
a4.450=networks.spatial(R=rep(0.4,200),n=450)
write.csv(a4.450,"a4.450.csv")
a5.450=networks.spatial(R=rep(0.5,200),n=450)
write.csv(a5.450,"a5.450.csv")
a6.450=networks.spatial(R=rep(0.6,200),n=450)
write.csv(a6.450,"a6.450.csv")
a7.450=networks.spatial(R=rep(0.7,200),n=450)
write.csv(a7.450,"a7.450.csv")
a8.450=networks.spatial(R=rep(0.8,200),n=450)
write.csv(a8.450,"a8.450.csv")
a9.450=networks.spatial(R=rep(0.9,200),n=450)
write.csv(a9.450,"a9.450.csv")
a10.450=networks.spatial(R=rep(0.99,200),n=450)
write.csv(a10.450,"a10.450.csv")
a11.450=networks.spatial(R=rep(0.05,200),n=450)
write.csv(a11.450,"a11.450.csv")
a12.450=networks.spatial(R=rep(0.01,200),n=450)
write.csv(a12.450,"a12.450.csv")
a13.450=networks.spatial(R=rep(0.08,200),n=450)
write.csv(a13.450,"a13.450.csv")

######## Saving all 450 nodes data to csv
a1=read.csv("a1.450.csv");a2=read.csv("a2.450.csv");a3=read.csv("a3.450.csv")
a4=read.csv("a4.450.csv");a5=read.csv("a5.450.csv");a6=read.csv("a6.450.csv")
a7=read.csv("a7.450.csv");a8=read.csv("a8.450.csv");a9=read.csv("a9.450.csv")
a10=read.csv("a10.450.csv")

df.450=rbind(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
dim(df.450)
node.num=450
df.450<-df.450 %>%dplyr::filter(order==node.num)
dim(df.450)
write.csv(df.450,"df.450.csv")


###500 nodes
a1.500=networks.spatial(R=rep(0.15,200),n=500)
write.csv(a1.500,"a1.500.csv")
a2.500=networks.spatial(R=rep(0.2,200),n=500)
write.csv(a2.500,"a2.500.csv")
a3.500=networks.spatial(R=rep(0.3,200),n=500)
write.csv(a3.500,"a3.500.csv")
a4.500=networks.spatial(R=rep(0.4,200),n=500)
write.csv(a4.500,"a4.500.csv")
a5.500=networks.spatial(R=rep(0.5,200),n=500)
write.csv(a5.500,"a5.500.csv")
a6.500=networks.spatial(R=rep(0.6,200),n=500)
write.csv(a6.500,"a6.500.csv")
a7.500=networks.spatial(R=rep(0.7,200),n=500)
write.csv(a7.500,"a7.500.csv")
a8.500=networks.spatial(R=rep(0.8,200),n=500)
write.csv(a8.500,"a8.500.csv")
a9.500=networks.spatial(R=rep(0.9,200),n=500)
write.csv(a9.500,"a9.500.csv")
a10.500=networks.spatial(R=rep(0.99,200),n=500)
write.csv(a10.500,"a10.500.csv")
a11.500=networks.spatial(R=rep(0.05,200),n=500)
write.csv(a11.500,"a11.500.csv")
a12.500=networks.spatial(R=rep(0.01,200),n=500)
write.csv(a12.500,"a12.500.csv")
a13.500=networks.spatial(R=rep(0.08,200),n=500)
write.csv(a13.500,"a13.500.csv")

######## Saving all 500 nodes data to csv
a1=read.csv("a1.500.csv");a2=read.csv("a2.500.csv");a3=read.csv("a3.500.csv")
a4=read.csv("a4.500.csv");a5=read.csv("a5.500.csv");a6=read.csv("a6.500.csv")
a7=read.csv("a7.500.csv")#;a8=read.csv("a8.500.csv");a9=read.csv("a9.500.csv")
#a10=read.csv("a10.500.csv")

df.500=rbind(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
dim(df.500)
node.num=500
df.500<-df.500 %>%dplyr::filter(order==node.num)
dim(df.500)
write.csv(df.500,"df.500.csv")


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                    Spatial hybrid data

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++ 50 NODES +++++++++++++++++++++++++++++++++++++++++++++++++++++++
### 50 nodes with p=0.1
m1.50.p0.1=networks.hybrid(R=rep(0.15,100),p=0.1,n=50)
write.csv(m1.50.p0.1,"m1.50.p0.1.csv")
m2.50.p0.1=networks.hybrid(R=rep(0.2,100),p=0.1,n=50)
write.csv(m2.50.p0.1,"m2.50.p0.1.csv")
m3.50.p0.1=networks.hybrid(R=rep(0.3,100),p=0.1,n=50)
write.csv(m3.50.p0.1,"m3.50.p0.1.csv")
m4.50.p0.1=networks.hybrid(R=rep(0.4,100),p=0.1,n=50)
write.csv(m4.50.p0.1,"m4.50.p0.1.csv")
m5.50.p0.1=networks.hybrid(R=rep(0.5,100),p=0.1,n=50)
write.csv(m5.50.p0.1,"m5.50.p0.1.csv")
m6.50.p0.1=networks.hybrid(R=rep(0.6,100),p=0.1,n=50)
write.csv(m6.50.p0.1,"m6.50.p0.1.csv")
m7.50.p0.1=networks.hybrid(R=rep(0.7,100),p=0.1,n=50)
write.csv(m7.50.p0.1,"m7.50.p0.1.csv")
m8.50.p0.1=networks.hybrid(R=rep(0.8,100),p=0.1,n=50)
write.csv(m8.50.p0.1,"m8.50.p0.1.csv")
m9.50.p0.1=networks.hybrid(R=rep(0.9,100),p=0.1,n=50)
write.csv(m9.50.p0.1,"m9.50.p0.1.csv")
m10.50.p0.1=networks.hybrid(R=rep(0.99,100),p=0.1,n=50)
write.csv(m10.50.p0.1,"m10.50.p0.1.csv")







#Hybrid spatial network
networks.hybrid=function(R=c(0.3,0.4,0.6,0.05),p=c(0.5,0.4,0.1,0.05),n=100){
  net=NULL;data=NULL
  k=1
  for (i in 1:length(R)){
    net[i]= makeSpatialHybrid(node.size =n,Radius = R[i],prob=p[i])
    data=RunSimOnGraphFeatures(net,nreps = 1)
    k=k+1
  }
  df=cbind(R,p,data)
  return(df)
}
#+++++++++++++++++++++ 50 NODES +++++++++++++++++++++++++++++++++++++++++++++++++++++++
### 50 nodes with p=0.1
m1.50.p0.1=networks.hybrid(R=rep(0.15,100),p=0.1,n=50)
write.csv(m1.50.p0.1,"m1.50.p0.1.csv")
m2.50.p0.1=networks.hybrid(R=rep(0.2,100),p=0.1,n=50)
write.csv(m2.50.p0.1,"m2.50.p0.1.csv")
m3.50.p0.1=networks.hybrid(R=rep(0.3,100),p=0.1,n=50)
write.csv(m3.50.p0.1,"m3.50.p0.1.csv")
m4.50.p0.1=networks.hybrid(R=rep(0.4,100),p=0.1,n=50)
write.csv(m4.50.p0.1,"m4.50.p0.1.csv")
m5.50.p0.1=networks.hybrid(R=rep(0.5,100),p=0.1,n=50)
write.csv(m5.50.p0.1,"m5.50.p0.1.csv")
m6.50.p0.1=networks.hybrid(R=rep(0.6,100),p=0.1,n=50)
write.csv(m6.50.p0.1,"m6.50.p0.1.csv")
m7.50.p0.1=networks.hybrid(R=rep(0.7,100),p=0.1,n=50)
write.csv(m7.50.p0.1,"m7.50.p0.1.csv")
m8.50.p0.1=networks.hybrid(R=rep(0.8,100),p=0.1,n=50)
write.csv(m8.50.p0.1,"m8.50.p0.1.csv")
m9.50.p0.1=networks.hybrid(R=rep(0.9,100),p=0.1,n=50)
write.csv(m9.50.p0.1,"m9.50.p0.1.csv")
m10.50.p0.1=networks.hybrid(R=rep(0.99,100),p=0.1,n=50)
write.csv(m10.50.p0.1,"m10.50.p0.1.csv")

######## Saving all 50 nodes p:0.1 data to csv
m1=read.csv("m1.50.p0.1.csv");m2=read.csv("m2.50.p0.1.csv");m3=read.csv("m3.50.p0.1.csv")
m4=read.csv("m4.50.p0.1.csv");m5=read.csv("m5.50.p0.1.csv");m6=read.csv("m6.50.p0.1.csv")
m7=read.csv("m7.50.p0.1.csv");m8=read.csv("m8.50.p0.1.csv");m9=read.csv("m9.50.p0.1.csv")
m10=read.csv("m10.50.p0.1.csv")

data.50.p0.1=rbind(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)
write.csv(data.50.p0.1,"data.50.p0.1.csv")

### 50 nodes with p=0.2
m1.50.p0.2=networks.hybrid(R=rep(0.15,100),p=0.2,n=50)
write.csv(m1.50.p0.2,"m1.50.p0.2.csv")
m2.50.p0.2=networks.hybrid(R=rep(0.2,100),p=0.2,n=50)
write.csv(m2.50.p0.2,"m2.50.p0.2.csv")
m3.50.p0.2=networks.hybrid(R=rep(0.3,100),p=0.2,n=50)
write.csv(m3.50.p0.2,"m3.50.p0.2.csv")
m4.50.p0.2=networks.hybrid(R=rep(0.4,100),p=0.2,n=50)
write.csv(m4.50.p0.2,"m4.50.p0.2.csv")
m5.50.p0.2=networks.hybrid(R=rep(0.5,100),p=0.2,n=50)
write.csv(m5.50.p0.2,"m5.50.p0.2.csv")
m6.50.p0.2=networks.hybrid(R=rep(0.6,100),p=0.2,n=50)
write.csv(m6.50.p0.2,"m6.50.p0.2.csv")
m7.50.p0.2=networks.hybrid(R=rep(0.7,100),p=0.2,n=50)
write.csv(m7.50.p0.2,"m7.50.p0.2.csv")
m8.50.p0.2=networks.hybrid(R=rep(0.8,100),p=0.2,n=50)
write.csv(m8.50.p0.2,"m8.50.p0.2.csv")
m9.50.p0.2=networks.hybrid(R=rep(0.9,100),p=0.2,n=50)
write.csv(m9.50.p0.2,"m9.50.p0.2.csv")
m10.50.p0.2=networks.hybrid(R=rep(0.99,100),p=0.2,n=50)
write.csv(m10.50.p0.2,"m10.50.p0.2.csv")

######## Saving all 50 nodes p:0.2 data to csv
m1=read.csv("m1.50.p0.2.csv");m2=read.csv("m2.50.p0.2.csv");m3=read.csv("m3.50.p0.2.csv")
m4=read.csv("m4.50.p0.2.csv");m5=read.csv("m5.50.p0.2.csv");m6=read.csv("m6.50.p0.2.csv")
m7=read.csv("m7.50.p0.2.csv");m8=read.csv("m8.50.p0.2.csv");m9=read.csv("m9.50.p0.2.csv")
m10=read.csv("m10.50.p0.2.csv")

data.50.p0.2=rbind(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)
write.csv(data.50.p0.2,"data.50.p0.2.csv")


### 50 nodes with p=0.3
m1.50.p0.3=networks.hybrid(R=rep(0.15,100),p=0.3,n=50)
write.csv(m1.50.p0.3,"m1.50.p0.3.csv")
m2.50.p0.3=networks.hybrid(R=rep(0.2,100),p=0.3,n=50)
write.csv(m2.50.p0.3,"m2.50.p0.3.csv")
m3.50.p0.3=networks.hybrid(R=rep(0.3,100),p=0.3,n=50)
write.csv(m3.50.p0.3,"m3.50.p0.3.csv")
m4.50.p0.3=networks.hybrid(R=rep(0.4,100),p=0.3,n=50)
write.csv(m4.50.p0.3,"m4.50.p0.3.csv")
m5.50.p0.3=networks.hybrid(R=rep(0.5,100),p=0.3,n=50)
write.csv(m5.50.p0.3,"m5.50.p0.3.csv")
m6.50.p0.3=networks.hybrid(R=rep(0.6,100),p=0.3,n=50)
write.csv(m6.50.p0.3,"m6.50.p0.3.csv")
m7.50.p0.3=networks.hybrid(R=rep(0.7,100),p=0.3,n=50)
write.csv(m7.50.p0.3,"m7.50.p0.3.csv")
m8.50.p0.3=networks.hybrid(R=rep(0.8,100),p=0.3,n=50)
write.csv(m8.50.p0.3,"m8.50.p0.3.csv")
m9.50.p0.3=networks.hybrid(R=rep(0.9,100),p=0.3,n=50)
write.csv(m9.50.p0.3,"m9.50.p0.3.csv")
m10.50.p0.3=networks.hybrid(R=rep(0.99,100),p=0.3,n=50)
write.csv(m10.50.p0.3,"m10.50.p0.3.csv")


######## Saving all 50 nodes p:0.3 data to csv
m1=read.csv("m1.50.p0.3.csv");m2=read.csv("m2.50.p0.3.csv");m3=read.csv("m3.50.p0.3.csv")
m4=read.csv("m4.50.p0.3.csv");m5=read.csv("m5.50.p0.3.csv");m6=read.csv("m6.50.p0.3.csv")
m7=read.csv("m7.50.p0.3.csv");m8=read.csv("m8.50.p0.3.csv");m9=read.csv("m9.50.p0.3.csv")
m10=read.csv("m10.50.p0.3.csv")
data.50.p0.3=rbind(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)
write.csv(data.50.p0.3,"data.50.p0.3.csv")

######## Saving all 50 nodes and p data to csv
k1=read.csv("data.50.p0.1.csv")
k2=read.csv("data.50.p0.2.csv")
k3=read.csv("data.50.p0.3.csv")
k4=rbind(k1,k2,k3)
k_50_all_p=write.csv(k4,"data.50.all_p.csv")

k_50=read.csv(data.50.all_p.csv)

##+++++++++++++++ 100 NODES +++++++++++++++++++++++++++++++++++++++++++++++++++++++
### 100 nodes with p=0.1
m1.100.p0.1=networks.hybrid(R=rep(0.15,100),p=0.1,n=100)
write.csv(m1.100.p0.1,"m1.100.p0.1.csv")
m2.100.p0.1=networks.hybrid(R=rep(0.2,100),p=0.1,n=100)
write.csv(m2.100.p0.1,"m2.100.p0.1.csv")
m3.100.p0.1=networks.hybrid(R=rep(0.3,100),p=0.1,n=100)
write.csv(m3.100.p0.1,"m3.100.p0.1.csv")
m4.100.p0.1=networks.hybrid(R=rep(0.4,100),p=0.1,n=100)
write.csv(m4.100.p0.1,"m4.100.p0.1.csv")
m5.100.p0.1=networks.hybrid(R=rep(0.5,100),p=0.1,n=100)
write.csv(m5.100.p0.1,"m5.100.p0.1.csv")
m6.100.p0.1=networks.hybrid(R=rep(0.6,100),p=0.1,n=100)
write.csv(m6.100.p0.1,"m6.100.p0.1.csv")
m7.100.p0.1=networks.hybrid(R=rep(0.7,100),p=0.1,n=100)
write.csv(m7.100.p0.1,"m7.100.p0.1.csv")
m8.100.p0.1=networks.hybrid(R=rep(0.8,100),p=0.1,n=100)
write.csv(m8.100.p0.1,"m8.100.p0.1.csv")
m9.100.p0.1=networks.hybrid(R=rep(0.9,100),p=0.1,n=100)
write.csv(m9.100.p0.1,"m9.100.p0.1.csv")
m10.100.p0.1=networks.hybrid(R=rep(0.99,100),p=0.1,n=100)
write.csv(m10.100.p0.1,"m10.100.p0.1.csv")

######## Saving all 100 nodes p:0.1 data to csv
m1=read.csv("m1.100.p0.1.csv");m2=read.csv("m2.100.p0.1.csv");m3=read.csv("m3.100.p0.1.csv")
m4=read.csv("m4.100.p0.1.csv");m5=read.csv("m5.100.p0.1.csv");m6=read.csv("m6.100.p0.1.csv")
m7=read.csv("m7.100.p0.1.csv");m8=read.csv("m8.100.p0.1.csv");m9=read.csv("m9.100.p0.1.csv")
m10=read.csv("m10.100.p0.1.csv")
data.100.p0.1=rbind(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)
write.csv(data.100.p0.1,"data.100.p0.1.csv")

### 100 nodes with p=0.2
m1.100.p0.2=networks.hybrid(R=rep(0.15,100),p=0.2,n=100)
write.csv(m1.100.p0.2,"m1.100.p0.2.csv")
m2.100.p0.2=networks.hybrid(R=rep(0.2,100),p=0.2,n=100)
write.csv(m2.100.p0.2,"m2.100.p0.2.csv")
m3.100.p0.2=networks.hybrid(R=rep(0.3,100),p=0.2,n=100)
write.csv(m3.100.p0.2,"m3.100.p0.2.csv")
m4.100.p0.2=networks.hybrid(R=rep(0.4,100),p=0.2,n=100)
write.csv(m4.100.p0.2,"m4.100.p0.2.csv")
m5.100.p0.2=networks.hybrid(R=rep(0.5,100),p=0.2,n=100)
write.csv(m5.100.p0.2,"m5.100.p0.2.csv")
m6.100.p0.2=networks.hybrid(R=rep(0.6,100),p=0.2,n=100)
write.csv(m6.100.p0.2,"m6.100.p0.2.csv")
m7.100.p0.2=networks.hybrid(R=rep(0.7,100),p=0.2,n=100)
write.csv(m7.100.p0.2,"m7.100.p0.2.csv")
m8.100.p0.2=networks.hybrid(R=rep(0.8,100),p=0.2,n=100)
write.csv(m8.100.p0.2,"m8.100.p0.2.csv")
m9.100.p0.2=networks.hybrid(R=rep(0.9,100),p=0.2,n=100)
write.csv(m9.100.p0.2,"m9.100.p0.2.csv")
m10.100.p0.2=networks.hybrid(R=rep(0.99,100),p=0.2,n=100)
write.csv(m10.100.p0.2,"m10.100.p0.2.csv")


######## Saving all 100 nodes p:0.2 data to csv
m1=read.csv("m1.100.p0.2.csv");m2=read.csv("m2.100.p0.2.csv");m3=read.csv("m3.100.p0.2.csv")
m4=read.csv("m4.100.p0.2.csv");m5=read.csv("m5.100.p0.2.csv");m6=read.csv("m6.100.p0.2.csv")
m7=read.csv("m7.100.p0.2.csv");m8=read.csv("m8.100.p0.2.csv");m9=read.csv("m9.100.p0.2.csv")
m10=read.csv("m10.100.p0.2.csv")
data.100.p0.2=rbind(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)
write.csv(data.100.p0.2,"data.100.p0.2.csv")


### 100 nodes with p=0.3
m1.100.p0.3=networks.hybrid(R=rep(0.15,100),p=0.3,n=100)
write.csv(m1.100.p0.3,"m1.100.p0.3.csv")
m2.100.p0.3=networks.hybrid(R=rep(0.2,100),p=0.3,n=100)
write.csv(m2.100.p0.3,"m2.100.p0.3.csv")
m3.100.p0.3=networks.hybrid(R=rep(0.3,100),p=0.3,n=100)
write.csv(m3.100.p0.3,"m3.100.p0.3.csv")
m4.100.p0.3=networks.hybrid(R=rep(0.4,100),p=0.3,n=100)
write.csv(m4.100.p0.3,"m4.100.p0.3.csv")
m5.100.p0.3=networks.hybrid(R=rep(0.5,100),p=0.3,n=100)
write.csv(m5.100.p0.3,"m5.100.p0.3.csv")
m6.100.p0.3=networks.hybrid(R=rep(0.6,100),p=0.3,n=100)
write.csv(m6.100.p0.3,"m6.100.p0.3.csv")
m7.100.p0.3=networks.hybrid(R=rep(0.7,100),p=0.3,n=100)
write.csv(m7.100.p0.3,"m7.100.p0.3.csv")
m8.100.p0.3=networks.hybrid(R=rep(0.8,100),p=0.3,n=100)
write.csv(m8.100.p0.3,"m8.100.p0.3.csv")
m9.100.p0.3=networks.hybrid(R=rep(0.9,100),p=0.3,n=100)
write.csv(m9.100.p0.3,"m9.100.p0.3.csv")
m10.100.p0.3=networks.hybrid(R=rep(0.99,100),p=0.3,n=100)
write.csv(m10.100.p0.3,"m10.100.p0.3.csv")


######## Saving all 100 nodes p:0.3 data to csv
m1=read.csv("m1.100.p0.3.csv");m2=read.csv("m2.100.p0.3.csv");m3=read.csv("m3.100.p0.3.csv")
m4=read.csv("m4.100.p0.3.csv");m5=read.csv("m5.100.p0.3.csv");m6=read.csv("m6.100.p0.3.csv")
m7=read.csv("m7.100.p0.3.csv");m9=read.csv("m9.100.p0.3.csv")#;m8=read.csv("m8.100.p0.3.csv")
m10=read.csv("m10.100.p0.3.csv")
data.100.p0.3=rbind(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)
write.csv(data.100.p0.3,"data.100.p0.3.csv")

######## Saving all 100 nodes and p data to csv
k1=read.csv("data.100.p0.1.csv")
k2=read.csv("data.100.p0.2.csv")
k3=read.csv("data.100.p0.3.csv")
k4=rbind(k1,k2,k3)
k_100_all_p=write.csv(k4,"data.100.all_p.csv")

k_100=read.csv(data.100.all_p.csv)






##+++++++++++++++ 150 NODES +++++++++++++++++++++++++++++++++++++++++++++++++++++++
### 150 nodes with p=0.1
m1.150.p0.1=networks.hybrid(R=rep(0.15,100),p=0.1,n=150)
write.csv(m1.150.p0.1,"m1.150.p0.1.csv")
m2.150.p0.1=networks.hybrid(R=rep(0.2,100),p=0.1,n=150)
write.csv(m2.150.p0.1,"m2.150.p0.1.csv")
m3.150.p0.1=networks.hybrid(R=rep(0.3,100),p=0.1,n=150)
write.csv(m3.150.p0.1,"m3.150.p0.1.csv")
m4.150.p0.1=networks.hybrid(R=rep(0.4,100),p=0.1,n=150)
write.csv(m4.150.p0.1,"m4.150.p0.1.csv")
m5.150.p0.1=networks.hybrid(R=rep(0.5,100),p=0.1,n=150)
write.csv(m5.150.p0.1,"m5.150.p0.1.csv")
m6.150.p0.1=networks.hybrid(R=rep(0.6,100),p=0.1,n=150)
write.csv(m6.150.p0.1,"m6.150.p0.1.csv")
m7.150.p0.1=networks.hybrid(R=rep(0.7,150),p=0.1,n=150)
write.csv(m7.150.p0.1,"m7.150.p0.1.csv")
m8.150.p0.1=networks.hybrid(R=rep(0.8,150),p=0.1,n=150)
write.csv(m8.150.p0.1,"m8.150.p0.1.csv")
m9.150.p0.1=networks.hybrid(R=rep(0.9,100),p=0.1,n=150)
write.csv(m9.150.p0.1,"m9.150.p0.1.csv")
m10.150.p0.1=networks.hybrid(R=rep(0.99,100),p=0.1,n=150)
write.csv(m10.150.p0.1,"m10.150.p0.1.csv")

### 150 nodes with p=0.2
m1.150.p0.2=networks.hybrid(R=rep(0.15,100),p=0.2,n=150)
write.csv(m1.150.p0.2,"m1.150.p0.2.csv")
m2.150.p0.2=networks.hybrid(R=rep(0.2,100),p=0.2,n=150)
write.csv(m2.150.p0.2,"m2.150.p0.2.csv")
m3.150.p0.2=networks.hybrid(R=rep(0.3,100),p=0.2,n=150)
write.csv(m3.150.p0.2,"m3.150.p0.2.csv")
m4.150.p0.2=networks.hybrid(R=rep(0.4,100),p=0.2,n=150)
write.csv(m4.150.p0.2,"m4.150.p0.2.csv")
m5.150.p0.2=networks.hybrid(R=rep(0.5,100),p=0.2,n=150)
write.csv(m5.150.p0.2,"m5.150.p0.2.csv")
m6.150.p0.2=networks.hybrid(R=rep(0.6,100),p=0.2,n=150)
write.csv(m6.150.p0.2,"m6.150.p0.2.csv")
m7.150.p0.2=networks.hybrid(R=rep(0.7,100),p=0.2,n=150)
write.csv(m7.150.p0.2,"m7.150.p0.2.csv")
m8.150.p0.2=networks.hybrid(R=rep(0.8,100),p=0.2,n=150)
write.csv(m8.150.p0.2,"m8.150.p0.2.csv")
m9.150.p0.2=networks.hybrid(R=rep(0.9,100),p=0.2,n=150)
write.csv(m9.150.p0.2,"m9.150.p0.2.csv")
m10.150.p0.2=networks.hybrid(R=rep(0.99,100),p=0.2,n=150)
write.csv(m10.150.p0.2,"m10.150.p0.2.csv")

### 150 nodes with p=0.3
m1.150.p0.3=networks.hybrid(R=rep(0.15,100),p=0.3,n=150)
write.csv(m1.150.p0.3,"m1.150.p0.3.csv")
m2.150.p0.3=networks.hybrid(R=rep(0.2,100),p=0.3,n=150)
write.csv(m2.150.p0.3,"m2.150.p0.3.csv")
m3.150.p0.3=networks.hybrid(R=rep(0.3,100),p=0.3,n=150)
write.csv(m3.150.p0.3,"m3.150.p0.3.csv")
m4.150.p0.3=networks.hybrid(R=rep(0.4,100),p=0.3,n=150)
write.csv(m4.150.p0.3,"m4.150.p0.3.csv")
m5.150.p0.3=networks.hybrid(R=rep(0.5,100),p=0.3,n=150)
write.csv(m5.150.p0.3,"m5.150.p0.3.csv")
m6.150.p0.3=networks.hybrid(R=rep(0.6,100),p=0.3,n=150)
write.csv(m6.150.p0.3,"m6.150.p0.3.csv")
m7.150.p0.3=networks.hybrid(R=rep(0.7,100),p=0.3,n=150)
write.csv(m7.150.p0.3,"m7.150.p0.3.csv")
m8.150.p0.3=networks.hybrid(R=rep(0.8,100),p=0.3,n=150)
write.csv(m8.150.p0.3,"m8.150.p0.3.csv")
m9.150.p0.3=networks.hybrid(R=rep(0.9,100),p=0.3,n=150)
write.csv(m9.150.p0.3,"m9.150.p0.3.csv")
m10.150.p0.3=networks.hybrid(R=rep(0.99,100),p=0.3,n=150)
write.csv(m10.150.p0.3,"m10.150.p0.3.csv")

##+++++++++++++++ 200 NODES +++++++++++++++++++++++++++++++++++++++++++++++++++++++
### 200 nodes with p=0.1
m1.200.p0.1=networks.hybrid(R=rep(0.15,100),p=0.1,n=200)
write.csv(m1.200.p0.1,"m1.200.p0.1.csv")
m2.200.p0.1=networks.hybrid(R=rep(0.2,100),p=0.1,n=200)
write.csv(m2.200.p0.1,"m2.200.p0.1.csv")
m3.200.p0.1=networks.hybrid(R=rep(0.3,100),p=0.1,n=200)
write.csv(m3.200.p0.1,"m3.200.p0.1.csv")
m4.200.p0.1=networks.hybrid(R=rep(0.4,100),p=0.1,n=200)
write.csv(m4.200.p0.1,"m4.200.p0.1.csv")
m5.200.p0.1=networks.hybrid(R=rep(0.5,100),p=0.1,n=200)
write.csv(m5.200.p0.1,"m5.200.p0.1.csv")
m6.200.p0.1=networks.hybrid(R=rep(0.6,100),p=0.1,n=200)
write.csv(m6.200.p0.1,"m6.200.p0.1.csv")
m7.200.p0.1=networks.hybrid(R=rep(0.7,200),p=0.1,n=200)
write.csv(m7.200.p0.1,"m7.200.p0.1.csv")
m8.200.p0.1=networks.hybrid(R=rep(0.8,200),p=0.1,n=200)
write.csv(m8.200.p0.1,"m8.200.p0.1.csv")
m9.200.p0.1=networks.hybrid(R=rep(0.9,100),p=0.1,n=200)
write.csv(m9.200.p0.1,"m9.200.p0.1.csv")
m10.200.p0.1=networks.hybrid(R=rep(0.99,100),p=0.1,n=200)
write.csv(m10.200.p0.1,"m10.200.p0.1.csv")


######## Saving all 200 nodes p:0.1 data to csv
m1=read.csv("m1.200.p0.1.csv");m2=read.csv("m2.200.p0.1.csv");m3=read.csv("m3.200.p0.1.csv")
m4=read.csv("m4.200.p0.1.csv");m5=read.csv("m5.200.p0.1.csv");m6=read.csv("m6.200.p0.1.csv")
m7=read.csv("m7.200.p0.1.csv");m8=read.csv("m8.200.p0.1.csv");m9=read.csv("m9.200.p0.1.csv")
m10=read.csv("m10.200.p0.1.csv")
data.200.p0.1=rbind(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)
write.csv(data.200.p0.1,"data.200.p0.1.csv")


### 200 nodes with p=0.2
m1.200.p0.2=networks.hybrid(R=rep(0.15,100),p=0.2,n=200)
write.csv(m1.200.p0.2,"m1.200.p0.2.csv")
m2.200.p0.2=networks.hybrid(R=rep(0.2,100),p=0.2,n=200)
write.csv(m2.200.p0.2,"m2.200.p0.2.csv")
m3.200.p0.2=networks.hybrid(R=rep(0.3,100),p=0.2,n=200)
write.csv(m3.200.p0.2,"m3.200.p0.2.csv")
m4.200.p0.2=networks.hybrid(R=rep(0.4,100),p=0.2,n=200)
write.csv(m4.200.p0.2,"m4.200.p0.2.csv")
m5.200.p0.2=networks.hybrid(R=rep(0.5,100),p=0.2,n=200)
write.csv(m5.200.p0.2,"m5.200.p0.2.csv")
m6.200.p0.2=networks.hybrid(R=rep(0.6,100),p=0.2,n=200)
write.csv(m6.200.p0.2,"m6.200.p0.2.csv")
m7.200.p0.2=networks.hybrid(R=rep(0.7,100),p=0.2,n=200)
write.csv(m7.200.p0.2,"m7.200.p0.2.csv")
m8.200.p0.2=networks.hybrid(R=rep(0.8,100),p=0.2,n=200)
write.csv(m8.200.p0.2,"m8.200.p0.2.csv")
m9.200.p0.2=networks.hybrid(R=rep(0.9,100),p=0.2,n=200)
write.csv(m9.200.p0.2,"m9.200.p0.2.csv")
m10.200.p0.2=networks.hybrid(R=rep(0.99,100),p=0.2,n=200)
write.csv(m10.200.p0.2,"m10.200.p0.2.csv")

######## Saving all 200 nodes p:0.2 data to csv
m1=read.csv("m1.200.p0.2.csv");m2=read.csv("m2.200.p0.2.csv");m3=read.csv("m3.200.p0.2.csv")
m4=read.csv("m4.200.p0.2.csv");m5=read.csv("m5.200.p0.2.csv");m6=read.csv("m6.200.p0.2.csv")
m7=read.csv("m7.200.p0.2.csv");m8=read.csv("m8.200.p0.2.csv");m9=read.csv("m9.200.p0.2.csv")
m10=read.csv("m10.200.p0.2.csv")
data.200.p0.2=rbind(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)
write.csv(data.200.p0.2,"data.200.p0.2.csv")

### 200 nodes with p=0.3
m1.200.p0.3=networks.hybrid(R=rep(0.15,100),p=0.3,n=200)
write.csv(m1.200.p0.3,"m1.200.p0.3.csv")
m2.200.p0.3=networks.hybrid(R=rep(0.2,100),p=0.3,n=200)
write.csv(m2.200.p0.3,"m2.200.p0.3.csv")
m3.200.p0.3=networks.hybrid(R=rep(0.3,100),p=0.3,n=200)
write.csv(m3.200.p0.3,"m3.200.p0.3.csv")
m4.200.p0.3=networks.hybrid(R=rep(0.4,100),p=0.3,n=200)
write.csv(m4.200.p0.3,"m4.200.p0.3.csv")
m5.200.p0.3=networks.hybrid(R=rep(0.5,100),p=0.3,n=200)
write.csv(m5.200.p0.3,"m5.200.p0.3.csv")
m6.200.p0.3=networks.hybrid(R=rep(0.6,100),p=0.3,n=200)
write.csv(m6.200.p0.3,"m6.200.p0.3.csv")
m7.200.p0.3=networks.hybrid(R=rep(0.7,100),p=0.3,n=200)
write.csv(m7.200.p0.3,"m7.200.p0.3.csv")
m8.200.p0.3=networks.hybrid(R=rep(0.8,100),p=0.3,n=200)
write.csv(m8.200.p0.3,"m8.200.p0.3.csv")
m9.200.p0.3=networks.hybrid(R=rep(0.9,100),p=0.3,n=200)
write.csv(m9.200.p0.3,"m9.200.p0.3.csv")
m10.200.p0.3=networks.hybrid(R=rep(0.99,100),p=0.3,n=200)
write.csv(m10.200.p0.3,"m10.200.p0.3.csv")


######## Saving all 200 nodes p:0.3 data to csv
m1=read.csv("m1.200.p0.3.csv");m2=read.csv("m2.200.p0.3.csv");m3=read.csv("m3.200.p0.3.csv")
m4=read.csv("m4.200.p0.3.csv");m5=read.csv("m5.200.p0.3.csv");m6=read.csv("m6.200.p0.3.csv")
m7=read.csv("m7.200.p0.3.csv");m8=read.csv("m8.200.p0.3.csv");m9=read.csv("m9.200.p0.3.csv")
m10=read.csv("m10.200.p0.3.csv")
data.200.p0.3=rbind(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)
write.csv(data.200.p0.3,"data.200.p0.3.csv")


######## Saving all 200 nodes and p data to csv
k1=read.csv("data.200.p0.1.csv")
k2=read.csv("data.200.p0.2.csv")
k3=read.csv("data.200.p0.3.csv")
k4=rbind(k1,k2,k3)
k_200_all_p=write.csv(k4,"data.200.all_p.csv")



##+++++++++++++++ 300 NODES +++++++++++++++++++++++++++++++++++++++++++++++++++++++
### 300 nodes with p=0.1
m1.300.p0.1=networks.hybrid(R=rep(0.15,100),p=0.1,n=300)
write.csv(m1.300.p0.1,"m1.300.p0.1.csv")
m2.300.p0.1=networks.hybrid(R=rep(0.2,100),p=0.1,n=300)
write.csv(m2.300.p0.1,"m2.300.p0.1.csv")
m3.300.p0.1=networks.hybrid(R=rep(0.3,100),p=0.1,n=300)
write.csv(m3.300.p0.1,"m3.300.p0.1.csv")
m4.300.p0.1=networks.hybrid(R=rep(0.4,100),p=0.1,n=300)
write.csv(m4.300.p0.1,"m4.300.p0.1.csv")
m5.300.p0.1=networks.hybrid(R=rep(0.8,100),p=0.1,n=300)
write.csv(m5.300.p0.1,"m5.300.p0.1.csv")
m6.300.p0.1=networks.hybrid(R=rep(0.6,100),p=0.1,n=300)
write.csv(m6.300.p0.1,"m6.300.p0.1.csv")
m7.300.p0.1=networks.hybrid(R=rep(0.7,300),p=0.1,n=300)
write.csv(m7.300.p0.1,"m7.300.p0.1.csv")
m8.300.p0.1=networks.hybrid(R=rep(0.8,300),p=0.1,n=300)
write.csv(m8.300.p0.1,"m8.300.p0.1.csv")
m9.300.p0.1=networks.hybrid(R=rep(0.9,100),p=0.1,n=300)
write.csv(m9.300.p0.1,"m9.300.p0.1.csv")
m10.300.p0.1=networks.hybrid(R=rep(0.99,100),p=0.1,n=300)
write.csv(m10.300.p0.1,"m10.300.p0.1.csv")

### 300 nodes with p=0.2
m1.300.p0.2=networks.hybrid(R=rep(0.15,100),p=0.2,n=300)
write.csv(m1.300.p0.2,"m1.300.p0.2.csv")
m2.300.p0.2=networks.hybrid(R=rep(0.2,100),p=0.2,n=300)
write.csv(m2.300.p0.2,"m2.300.p0.2.csv")
m3.300.p0.2=networks.hybrid(R=rep(0.3,100),p=0.2,n=300)
write.csv(m3.300.p0.2,"m3.300.p0.2.csv")
m4.300.p0.2=networks.hybrid(R=rep(0.4,100),p=0.2,n=300)
write.csv(m4.300.p0.2,"m4.300.p0.2.csv")
m5.300.p0.2=networks.hybrid(R=rep(0.8,100),p=0.2,n=300)
write.csv(m5.300.p0.2,"m5.300.p0.2.csv")
m6.300.p0.2=networks.hybrid(R=rep(0.6,100),p=0.2,n=300)
write.csv(m6.300.p0.2,"m6.300.p0.2.csv")
m7.300.p0.2=networks.hybrid(R=rep(0.7,100),p=0.2,n=300)
write.csv(m7.300.p0.2,"m7.300.p0.2.csv")
m8.300.p0.2=networks.hybrid(R=rep(0.8,100),p=0.2,n=300)
write.csv(m8.300.p0.2,"m8.300.p0.2.csv")
m9.300.p0.2=networks.hybrid(R=rep(0.9,100),p=0.2,n=300)
write.csv(m9.300.p0.2,"m9.300.p0.2.csv")
m10.300.p0.2=networks.hybrid(R=rep(0.99,100),p=0.2,n=300)
write.csv(m10.300.p0.2,"m10.300.p0.2.csv")

### 300 nodes with p=0.3
m1.300.p0.3=networks.hybrid(R=rep(0.15,100),p=0.3,n=300)
write.csv(m1.300.p0.3,"m1.300.p0.3.csv")
m2.300.p0.3=networks.hybrid(R=rep(0.2,100),p=0.3,n=300)
write.csv(m2.300.p0.3,"m2.300.p0.3.csv")
m3.300.p0.3=networks.hybrid(R=rep(0.3,100),p=0.3,n=300)
write.csv(m3.300.p0.3,"m3.300.p0.3.csv")
m4.300.p0.3=networks.hybrid(R=rep(0.4,100),p=0.3,n=300)
write.csv(m4.300.p0.3,"m4.300.p0.3.csv")
m5.300.p0.3=networks.hybrid(R=rep(0.8,100),p=0.3,n=300)
write.csv(m5.300.p0.3,"m5.300.p0.3.csv")
m6.300.p0.3=networks.hybrid(R=rep(0.6,100),p=0.3,n=300)
write.csv(m6.300.p0.3,"m6.300.p0.3.csv")
m7.300.p0.3=networks.hybrid(R=rep(0.7,100),p=0.3,n=300)
write.csv(m7.300.p0.3,"m7.300.p0.3.csv")
m8.300.p0.3=networks.hybrid(R=rep(0.8,100),p=0.3,n=300)
write.csv(m8.300.p0.3,"m8.300.p0.3.csv")
m9.300.p0.3=networks.hybrid(R=rep(0.9,100),p=0.3,n=300)
write.csv(m9.300.p0.3,"m9.300.p0.3.csv")
m10.300.p0.3=networks.hybrid(R=rep(0.99,100),p=0.3,n=300)
write.csv(m10.300.p0.3,"m10.300.p0.3.csv")

##+++++++++++++++++++++++ 500 NODES +++++++++++++++++++++++++++++++++++++++++++++++++++++++
### 500 nodes with p=0.1
m1.500.p0.1=networks.hybrid(R=rep(0.15,100),p=0.1,n=500)
write.csv(m1.500.p0.1,"m1.500.p0.1.csv")
m2.500.p0.1=networks.hybrid(R=rep(0.2,100),p=0.1,n=500)
write.csv(m2.500.p0.1,"m2.500.p0.1.csv")
m3.500.p0.1=networks.hybrid(R=rep(0.3,100),p=0.1,n=500)
write.csv(m3.500.p0.1,"m3.500.p0.1.csv")
m4.500.p0.1=networks.hybrid(R=rep(0.4,100),p=0.1,n=500)
write.csv(m4.500.p0.1,"m4.500.p0.1.csv")
m5.500.p0.1=networks.hybrid(R=rep(0.8,100),p=0.1,n=500)
write.csv(m5.500.p0.1,"m5.500.p0.1.csv")
m6.500.p0.1=networks.hybrid(R=rep(0.6,100),p=0.1,n=500)
write.csv(m6.500.p0.1,"m6.500.p0.1.csv")
m7.500.p0.1=networks.hybrid(R=rep(0.7,500),p=0.1,n=500)
write.csv(m7.500.p0.1,"m7.500.p0.1.csv")
m8.500.p0.1=networks.hybrid(R=rep(0.8,500),p=0.1,n=500)
write.csv(m8.500.p0.1,"m8.500.p0.1.csv")
m9.500.p0.1=networks.hybrid(R=rep(0.9,100),p=0.1,n=500)
write.csv(m9.500.p0.1,"m9.500.p0.1.csv")
m10.500.p0.1=networks.hybrid(R=rep(0.99,100),p=0.1,n=500)
write.csv(m10.500.p0.1,"m10.500.p0.1.csv")

### 500 nodes with p=0.2
m1.500.p0.2=networks.hybrid(R=rep(0.15,100),p=0.2,n=500)
write.csv(m1.500.p0.2,"m1.500.p0.2.csv")
m2.500.p0.2=networks.hybrid(R=rep(0.2,100),p=0.2,n=500)
write.csv(m2.500.p0.2,"m2.500.p0.2.csv")
m3.500.p0.2=networks.hybrid(R=rep(0.3,100),p=0.2,n=500)
write.csv(m3.500.p0.2,"m3.500.p0.2.csv")
m4.500.p0.2=networks.hybrid(R=rep(0.4,100),p=0.2,n=500)
write.csv(m4.500.p0.2,"m4.500.p0.2.csv")
m5.500.p0.2=networks.hybrid(R=rep(0.8,100),p=0.2,n=500)
write.csv(m5.500.p0.2,"m5.500.p0.2.csv")
m6.500.p0.2=networks.hybrid(R=rep(0.6,100),p=0.2,n=500)
write.csv(m6.500.p0.2,"m6.500.p0.2.csv")
m7.500.p0.2=networks.hybrid(R=rep(0.7,100),p=0.2,n=500)
write.csv(m7.500.p0.2,"m7.500.p0.2.csv")
m8.500.p0.2=networks.hybrid(R=rep(0.8,100),p=0.2,n=500)
write.csv(m8.500.p0.2,"m8.500.p0.2.csv")
m9.500.p0.2=networks.hybrid(R=rep(0.9,100),p=0.2,n=500)
write.csv(m9.500.p0.2,"m9.500.p0.2.csv")
m10.500.p0.2=networks.hybrid(R=rep(0.99,100),p=0.2,n=500)
write.csv(m10.500.p0.2,"m10.500.p0.2.csv")

### 500 nodes with p=0.3
m1.500.p0.3=networks.hybrid(R=rep(0.15,100),p=0.3,n=500)
write.csv(m1.500.p0.3,"m1.500.p0.3.csv")
m2.500.p0.3=networks.hybrid(R=rep(0.2,100),p=0.3,n=500)
write.csv(m2.500.p0.3,"m2.500.p0.3.csv")
m3.500.p0.3=networks.hybrid(R=rep(0.3,100),p=0.3,n=500)
write.csv(m3.500.p0.3,"m3.500.p0.3.csv")
m4.500.p0.3=networks.hybrid(R=rep(0.4,100),p=0.3,n=500)
write.csv(m4.500.p0.3,"m4.500.p0.3.csv")
m5.500.p0.3=networks.hybrid(R=rep(0.8,100),p=0.3,n=500)
write.csv(m5.500.p0.3,"m5.500.p0.3.csv")
m6.500.p0.3=networks.hybrid(R=rep(0.6,100),p=0.3,n=500)
write.csv(m6.500.p0.3,"m6.500.p0.3.csv")
m7.500.p0.3=networks.hybrid(R=rep(0.7,100),p=0.3,n=500)
write.csv(m7.500.p0.3,"m7.500.p0.3.csv")
m8.500.p0.3=networks.hybrid(R=rep(0.8,100),p=0.3,n=500)
write.csv(m8.500.p0.3,"m8.500.p0.3.csv")
m9.500.p0.3=networks.hybrid(R=rep(0.9,100),p=0.3,n=500)
write.csv(m9.500.p0.3,"m9.500.p0.3.csv")
m10.500.p0.3=networks.hybrid(R=rep(0.99,100),p=0.3,n=500)
write.csv(m10.500.p0.3,"m10.500.p0.3.csv")

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++ 50 NODES +++++++++++++++++++++++++++++++++++++++++++++++++++++++
### 50 nodes with p=runif(100,0.1,0.35)
e1.50.p0.1=networks.hybrid(R=rep(0.15,100),p=runif(100,0.1,0.35),n=50)
write.csv(e1.50.p0.1,"e1.50.p0.1.csv")
e2.50.p0.1=networks.hybrid(R=rep(0.2,100),p=runif(100,0.1,0.35),n=50)
write.csv(e2.50.p0.1,"e2.50.p0.1.csv")
e3.50.p0.1=networks.hybrid(R=rep(0.3,100),p=runif(100,0.1,0.35),n=50)
write.csv(e3.50.p0.1,"e3.50.p0.1.csv")
e4.50.p0.1=networks.hybrid(R=rep(0.4,100),p=runif(100,0.1,0.35),n=50)
write.csv(e4.50.p0.1,"e4.50.p0.1.csv")
e5.50.p0.1=networks.hybrid(R=rep(0.5,100),p=runif(100,0.1,0.35),n=50)
write.csv(e5.50.p0.1,"e5.50.p0.1.csv")
e6.50.p0.1=networks.hybrid(R=rep(0.6,100),p=runif(100,0.1,0.35),n=50)
write.csv(e6.50.p0.1,"e6.50.p0.1.csv")
e7.50.p0.1=networks.hybrid(R=rep(0.7,100),p=runif(100,0.1,0.35),n=50)
write.csv(e7.50.p0.1,"e7.50.p0.1.csv")
e8.50.p0.1=networks.hybrid(R=rep(0.8,100),p=runif(100,0.1,0.35),n=50)
write.csv(e8.50.p0.1,"e8.50.p0.1.csv")
e9.50.p0.1=networks.hybrid(R=rep(0.9,100),p=runif(100,0.1,0.35),n=50)
write.csv(e9.50.p0.1,"e9.50.p0.1.csv")
e10.50.p0.1=networks.hybrid(R=rep(0.99,100),p=runif(100,0.1,0.35),n=50)
write.csv(e10.50.p0.1,"e10.50.p0.1.csv")

######## Saving all 50 nodes p:0.1 data to csv
e1=read.csv("e1.50.p0.1.csv");e2=read.csv("e2.50.p0.1.csv");e3=read.csv("e3.50.p0.1.csv")
e4=read.csv("e4.50.p0.1.csv");e5=read.csv("e5.50.p0.1.csv");e6=read.csv("e6.50.p0.1.csv")
e7=read.csv("e7.50.p0.1.csv");e8=read.csv("e8.50.p0.1.csv");e9=read.csv("e9.50.p0.1.csv")
e10=read.csv("e10.50.p0.1.csv")

data.50.p0.1=rbind(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)
write.csv(data.50.p0.1,"data.50.p0.1.csv")

#+++++++++++++++++++++ 100 NODES +++++++++++++++++++++++++++++++++++++++++++++++++++++++
### 100 nodes with p=runif(100,0.1,0.35)
e1.100.p0.1=networks.hybrid(R=rep(0.15,100),p=runif(100,0.1,0.35),n=100)
write.csv(e1.100.p0.1,"e1.100.p0.1.csv")
e2.100.p0.1=networks.hybrid(R=rep(0.2,100),p=runif(100,0.1,0.35),n=100)
write.csv(e2.100.p0.1,"e2.100.p0.1.csv")
e3.100.p0.1=networks.hybrid(R=rep(0.3,100),p=runif(100,0.1,0.35),n=100)
write.csv(e3.100.p0.1,"e3.100.p0.1.csv")
e4.100.p0.1=networks.hybrid(R=rep(0.4,100),p=runif(100,0.1,0.35),n=100)
write.csv(e4.100.p0.1,"e4.100.p0.1.csv")
e5.100.p0.1=networks.hybrid(R=rep(0.5,100),p=runif(100,0.1,0.35),n=100)
write.csv(e5.100.p0.1,"e5.100.p0.1.csv")
e6.100.p0.1=networks.hybrid(R=rep(0.6,100),p=runif(100,0.1,0.35),n=100)
write.csv(e6.100.p0.1,"e6.100.p0.1.csv")
e7.100.p0.1=networks.hybrid(R=rep(0.7,100),p=runif(100,0.1,0.35),n=100)
write.csv(e7.100.p0.1,"e7.100.p0.1.csv")
e8.100.p0.1=networks.hybrid(R=rep(0.8,100),p=runif(100,0.1,0.35),n=100)
write.csv(e8.100.p0.1,"e8.100.p0.1.csv")
e9.100.p0.1=networks.hybrid(R=rep(0.9,100),p=runif(100,0.1,0.35),n=100)
write.csv(e9.100.p0.1,"e9.100.p0.1.csv")
e10.100.p0.1=networks.hybrid(R=rep(0.99,100),p=runif(100,0.1,0.35),n=100)
write.csv(e10.100.p0.1,"e10.100.p0.1.csv")

######## Saving all 100 nodes p:0.1 data to csv
e1=read.csv("e1.100.p0.1.csv");e2=read.csv("e2.100.p0.1.csv");e3=read.csv("e3.100.p0.1.csv")
e4=read.csv("e4.100.p0.1.csv");e5=read.csv("e5.100.p0.1.csv");e6=read.csv("e6.100.p0.1.csv")
e7=read.csv("e7.100.p0.1.csv");e8=read.csv("e8.100.p0.1.csv");e9=read.csv("e9.100.p0.1.csv")
e10=read.csv("e10.100.p0.1.csv")

data.100.p0.1=rbind(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)
write.csv(data.100.p0.1,"data.100.p0.1.csv")



#+++++++++++++++++++++ 150 NODES +++++++++++++++++++++++++++++++++++++++++++++++++++++++
### 150 nodes with p=runif(150,0.1,0.35)
e1.150.p0.1=networks.hybrid(R=rep(0.15,100),p=runif(100,0.1,0.35),n=150)
write.csv(e1.150.p0.1,"e1.150.p0.1.csv")
e2.150.p0.1=networks.hybrid(R=rep(0.2,100),p=runif(100,0.1,0.35),n=150)
write.csv(e2.150.p0.1,"e2.150.p0.1.csv")
e3.150.p0.1=networks.hybrid(R=rep(0.3,100),p=runif(100,0.1,0.35),n=150)
write.csv(e3.150.p0.1,"e3.150.p0.1.csv")
e4.150.p0.1=networks.hybrid(R=rep(0.4,100),p=runif(100,0.1,0.35),n=150)
write.csv(e4.150.p0.1,"e4.150.p0.1.csv")
e5.150.p0.1=networks.hybrid(R=rep(0.5,100),p=runif(100,0.1,0.35),n=150)
write.csv(e5.150.p0.1,"e5.150.p0.1.csv")
e6.150.p0.1=networks.hybrid(R=rep(0.6,100),p=runif(100,0.1,0.35),n=150)
write.csv(e6.150.p0.1,"e6.150.p0.1.csv")
e7.150.p0.1=networks.hybrid(R=rep(0.7,100),p=runif(100,0.1,0.35),n=150)
write.csv(e7.150.p0.1,"e7.150.p0.1.csv")
e8.150.p0.1=networks.hybrid(R=rep(0.8,100),p=runif(100,0.1,0.35),n=150)
write.csv(e8.150.p0.1,"e8.150.p0.1.csv")
e9.150.p0.1=networks.hybrid(R=rep(0.9,100),p=runif(100,0.1,0.35),n=150)
write.csv(e9.150.p0.1,"e9.150.p0.1.csv")
e10.150.p0.1=networks.hybrid(R=rep(0.99,100),p=runif(100,0.1,0.35),n=150)
write.csv(e10.150.p0.1,"e10.150.p0.1.csv")

######## Saving all 150 nodes p:0.1 data to csv
e1=read.csv("e1.150.p0.1.csv");e2=read.csv("e2.150.p0.1.csv");e3=read.csv("e3.150.p0.1.csv")
e4=read.csv("e4.150.p0.1.csv");e5=read.csv("e5.150.p0.1.csv");e6=read.csv("e6.150.p0.1.csv")
e7=read.csv("e7.150.p0.1.csv");e8=read.csv("e8.150.p0.1.csv");e9=read.csv("e9.150.p0.1.csv")
e10=read.csv("e10.150.p0.1.csv")

data.150.p0.1=rbind(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)
write.csv(data.150.p0.1,"data.150.p0.1.csv")

#+++++++++++++++++++++ 200 NODES +++++++++++++++++++++++++++++++++++++++++++++++++++++++
### 200 nodes with p=runif(200,0.1,0.35)
e1.200.p0.1=networks.hybrid(R=rep(0.15,100),p=runif(100,0.1,0.35),n=200)
write.csv(e1.200.p0.1,"e1.200.p0.1.csv")
e2.200.p0.1=networks.hybrid(R=rep(0.2,100),p=runif(100,0.1,0.35),n=200)
write.csv(e2.200.p0.1,"e2.200.p0.1.csv")
e3.200.p0.1=networks.hybrid(R=rep(0.3,100),p=runif(100,0.1,0.35),n=200)
write.csv(e3.200.p0.1,"e3.200.p0.1.csv")
e4.200.p0.1=networks.hybrid(R=rep(0.4,100),p=runif(100,0.1,0.35),n=200)
write.csv(e4.200.p0.1,"e4.200.p0.1.csv")
e5.200.p0.1=networks.hybrid(R=rep(0.5,100),p=runif(100,0.1,0.35),n=200)
write.csv(e5.200.p0.1,"e5.200.p0.1.csv")
e6.200.p0.1=networks.hybrid(R=rep(0.6,100),p=runif(100,0.1,0.35),n=200)
write.csv(e6.200.p0.1,"e6.200.p0.1.csv")
e7.200.p0.1=networks.hybrid(R=rep(0.7,100),p=runif(100,0.1,0.35),n=200)
write.csv(e7.200.p0.1,"e7.200.p0.1.csv")
e8.200.p0.1=networks.hybrid(R=rep(0.8,100),p=runif(100,0.1,0.35),n=200)
write.csv(e8.200.p0.1,"e8.200.p0.1.csv")
e9.200.p0.1=networks.hybrid(R=rep(0.9,100),p=runif(100,0.1,0.35),n=200)
write.csv(e9.200.p0.1,"e9.200.p0.1.csv")
e10.200.p0.1=networks.hybrid(R=rep(0.99,100),p=runif(100,0.1,0.35),n=200)
write.csv(e10.200.p0.1,"e10.200.p0.1.csv")

######## Saving all 200 nodes p:0.1 data to csv
e1=read.csv("e1.200.p0.1.csv");e2=read.csv("e2.200.p0.1.csv");e3=read.csv("e3.200.p0.1.csv")
e4=read.csv("e4.200.p0.1.csv");e5=read.csv("e5.200.p0.1.csv");e6=read.csv("e6.200.p0.1.csv")
e7=read.csv("e7.200.p0.1.csv");e8=read.csv("e8.200.p0.1.csv");e9=read.csv("e9.200.p0.1.csv")
e10=read.csv("e10.200.p0.1.csv")

data.200.p0.1=rbind(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)
write.csv(data.200.p0.1,"data.200.p0.1.csv")

#+++++++++++++++++++++ 300 NODES +++++++++++++++++++++++++++++++++++++++++++++++++++++++
### 300 nodes with p=runif(300,0.1,0.35)
e1.300.p0.1=networks.hybrid(R=rep(0.15,100),p=runif(100,0.1,0.35),n=300)
write.csv(e1.300.p0.1,"e1.300.p0.1.csv")
e2.300.p0.1=networks.hybrid(R=rep(0.2,100),p=runif(100,0.1,0.35),n=300)
write.csv(e2.300.p0.1,"e2.300.p0.1.csv")
e3.300.p0.1=networks.hybrid(R=rep(0.3,100),p=runif(100,0.1,0.35),n=300)
write.csv(e3.300.p0.1,"e3.300.p0.1.csv")
e4.300.p0.1=networks.hybrid(R=rep(0.4,100),p=runif(100,0.1,0.35),n=300)
write.csv(e4.300.p0.1,"e4.300.p0.1.csv")
e5.300.p0.1=networks.hybrid(R=rep(0.5,100),p=runif(100,0.1,0.35),n=300)
write.csv(e5.300.p0.1,"e5.300.p0.1.csv")
e6.300.p0.1=networks.hybrid(R=rep(0.6,100),p=runif(100,0.1,0.35),n=300)
write.csv(e6.300.p0.1,"e6.300.p0.1.csv")
e7.300.p0.1=networks.hybrid(R=rep(0.7,100),p=runif(100,0.1,0.35),n=300)
write.csv(e7.300.p0.1,"e7.300.p0.1.csv")
e8.300.p0.1=networks.hybrid(R=rep(0.8,100),p=runif(100,0.1,0.35),n=300)
write.csv(e8.300.p0.1,"e8.300.p0.1.csv")
e9.300.p0.1=networks.hybrid(R=rep(0.9,100),p=runif(100,0.1,0.35),n=300)
write.csv(e9.300.p0.1,"e9.300.p0.1.csv")
e10.300.p0.1=networks.hybrid(R=rep(0.99,100),p=runif(100,0.1,0.35),n=300)
write.csv(e10.300.p0.1,"e10.300.p0.1.csv")

######## Saving all 300 nodes p:0.1 data to csv
e1=read.csv("e1.300.p0.1.csv");e2=read.csv("e2.300.p0.1.csv");e3=read.csv("e3.300.p0.1.csv")
e4=read.csv("e4.300.p0.1.csv");e5=read.csv("e5.300.p0.1.csv");e6=read.csv("e6.300.p0.1.csv")
e7=read.csv("e7.300.p0.1.csv");e8=read.csv("e8.300.p0.1.csv");e9=read.csv("e9.300.p0.1.csv")
e10=read.csv("e10.300.p0.1.csv")

data.300.p0.1=rbind(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)
write.csv(data.300.p0.1,"data.300.p0.1.csv")

#+++++++++++++++++++++ 500 NODES +++++++++++++++++++++++++++++++++++++++++++++++++++++++
### 500 nodes with p=runif(500,0.1,0.35)
e1.500.p0.1=networks.hybrid(R=rep(0.15,100),p=runif(100,0.1,0.35),n=500)
write.csv(e1.500.p0.1,"e1.500.p0.1.csv")
e2.500.p0.1=networks.hybrid(R=rep(0.2,100),p=runif(100,0.1,0.35),n=500)
write.csv(e2.500.p0.1,"e2.500.p0.1.csv")
e3.500.p0.1=networks.hybrid(R=rep(0.3,100),p=runif(100,0.1,0.35),n=500)
write.csv(e3.500.p0.1,"e3.500.p0.1.csv")
e4.500.p0.1=networks.hybrid(R=rep(0.4,100),p=runif(100,0.1,0.35),n=500)
write.csv(e4.500.p0.1,"e4.500.p0.1.csv")
e5.500.p0.1=networks.hybrid(R=rep(0.5,100),p=runif(100,0.1,0.35),n=500)
write.csv(e5.500.p0.1,"e5.500.p0.1.csv")
e6.500.p0.1=networks.hybrid(R=rep(0.6,100),p=runif(100,0.1,0.35),n=500)
write.csv(e6.500.p0.1,"e6.500.p0.1.csv")
e7.500.p0.1=networks.hybrid(R=rep(0.7,100),p=runif(100,0.1,0.35),n=500)
write.csv(e7.500.p0.1,"e7.500.p0.1.csv")
e8.500.p0.1=networks.hybrid(R=rep(0.8,100),p=runif(100,0.1,0.35),n=500)
write.csv(e8.500.p0.1,"e8.500.p0.1.csv")
e9.500.p0.1=networks.hybrid(R=rep(0.9,100),p=runif(100,0.1,0.35),n=500)
write.csv(e9.500.p0.1,"e9.500.p0.1.csv")
e10.500.p0.1=networks.hybrid(R=rep(0.99,100),p=runif(100,0.1,0.35),n=500)
write.csv(e10.500.p0.1,"e10.500.p0.1.csv")

######## Saving all 500 nodes p:0.1 data to csv
e1=read.csv("e1.500.p0.1.csv");e2=read.csv("e2.500.p0.1.csv");e3=read.csv("e3.500.p0.1.csv")
e4=read.csv("e4.500.p0.1.csv");e5=read.csv("e5.500.p0.1.csv");e6=read.csv("e6.500.p0.1.csv")
e7=read.csv("e7.500.p0.1.csv");e8=read.csv("e8.500.p0.1.csv");e9=read.csv("e9.500.p0.1.csv")
e10=read.csv("e10.500.p0.1.csv")

data.500.p0.1=rbind(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)
write.csv(data.500.p0.1,"data.500.p0.1.csv")







