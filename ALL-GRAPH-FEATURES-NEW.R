#+++++++++++++++++++++++++++++++++++++++++++++++++++++
# Network properties for multiple graphs
#++++++++++++++++++++++++++++++++++++++++++++++++++

#setwd("C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi-New")
#projectfolder=paste(getwd(), "/", sep = '')
#FolderName=setwd("C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi-New/Animal-Social-Networks")

library(igraph)

##--Deleting isolated vertices
# del_vertices=function(G){
#     Isolated = which(igraph::degree(G)==0)
#     Net=igraph::delete.vertices(G, Isolated)
#     
#     return(Net)
#   }

##--Finding largest component 
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
readtable=function(edgelist="aves-barn-swallow-non-physical.edges"){
x=read.table(edgelist, 
             fill = TRUE , header = FALSE)
return(x)
}

##----Getting graph features on a list of graph
Network.Summary <- function(FolderName){
  
  
     #FolderName="C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi-New/ASN/"  
     Net.lists = list.files(FolderName)
     Net.Feat.Summary = data.frame()
     #filenames <- list.files(FolderName, pattern="*.edges", full.names=TRUE)
     #Net.lists=as.list(filenames) 
     Net.lists=as.list(Net.lists)
     #g = lapply(Net.lists[1:10],read.table)
     g = lapply(Net.lists,read.table,fill = TRUE , header = FALSE)
     GraphNames = cbind(lapply(Net.lists,function(x) gsub(".edges","",x)))
     
     G = lapply(g, largest_comp)
     # 
     g.features=calcGraphFeatures(G)
     all_graphs=cbind(GraphNames, g.features)
  return(all_graphs)
}
#setwd("C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi-New/Animal-Social-Networks")
# file_path <- "Animal-Social-Networks"
# FolderName <- file(file_path, "rt")
#FolderName = c("../Animal-Social-Networks")
#FolderName = c('networks')

#FolderName="C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi-New/ASN/"
FolderName=setwd("C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi-New/ASN/")
z=Network.Summary(FolderName )

z
data.animals=as_tibble(z)

data.animals=data.animals%>%mutate_if(is.character,factor)
data.animals<- apply(data.animals,2,as.character)

write.csv(data.animals,"GraphFeatOnAllAnimNets.csv")






##----Downloadinga numerous network at once from a repository
setwd("C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi-New/ANIMAL-SOCIAL-NETWORS-REPOSITORY")
library(tidyverse)
library(RSelenium)
library(wdman)
library(netstat)

#binman::list_versions("chromedriver") 
#selenium.obj=selenium(retcommand = T,check = F) 


# connecting to selenium server
rs.driver <- rsDriver(
  browser = 'chrome',
  chromever = "112.0.5615.49",
  verbose = F,
  port = free_port()
)

# close the server
#rs.driver$server$stop()

# access the client object
remDr <- rs.driver$client

# open a web browser
remDr$open()

# navigate to the network repository website
remDr$navigate("https://networkrepository.com/asn.php")


# find the 'a' tags within the specified class name using the xpath method
#data_files <- remDr$findElements(using = 'xpath', "//tr[@class='success hrefRow tooltips']")
data_files <- remDr$findElements(using = 'xpath', "//td[@class='tdcell']/a")

file_name=remDr$findElements(using = 'xpath', "//tr[@class='success hrefRow tooltips']")

# Net.lists = data_files 
# Names = cbind(lapply(Net.lists,function(x) gsub(".edges","",x)))

# return the names of the files
data_file_names <- lapply(file_name, function(x) {
  x$getElementText() %>% unlist()
}) %>% flatten_chr() %>% 
  str_remove_all("[:]")

# return the links to the files
data_file_links <- lapply(data_files, function(x) {
  x$getElementAttribute('href') %>% unlist()
}) %>% flatten_chr()

# the loop to download all the files
for (i in 1:length(data_file_names)) {
  download.file(
    url = data_file_links[i],
    destfile = paste0(data_file_names[i], gsub(x = data_file_links[i], pattern = ".*[.]", replacement = "."))
  )
}


### extracting elements of all zip files at once
# 1. type cmd in place of the path in your working directory
# 2. type<- "C:\Program Files\7-Zip\7z.exe" e *.zip and hit enter
# 3. select a for multiple files and done



#calcGraphFeatures(Graphs=m)


# x=makeSpatialGraphs(496,0.06)#m=4,f=0.017,s=9, e=984
# y=calcGraphFeatures(x)
# y
# 
# 
# r1=runif(1000,0.74,0.78);n1=35
# r2=runif(1000,0.44,0.48);n2=17
# r3=runif(1000,0.145,0.15);n3=117
# r4=runif(1000,0.24,0.25);n4=20
# r5=runif(1000,0.17,0.18);n5=62
# r6=runif(1000,0.055,0.06);n6=496
# 
# R=r6
# node.num=n6
#  
#  
#  net=NULL
#  for (i in 1:length(R)){
#    # net=fastSpatialNetwork(n=node.num,r=R[[i]],makeConnected=TRUE, keepCellsSeparate=FALSE)
#    net[i]=makeSpatialGraphs(node.size =node.num,Radius=R[i])
#  }
# 
#  
#  Data=RunSimOnGraphFeatures(net,nreps = 1)
#  
#  DF.DATA=cbind(R,Data)
#  x=write.csv(DF.DATA,"node496-r1000sample.csv")
#  