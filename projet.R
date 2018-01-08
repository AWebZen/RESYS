#Sys.setenv(https_proxy="https://proxy:3128")

##METHODES
#Miic
#Pc -> pcalg
#Aracne -> minet
#Bayesian -> bnlearn

##Installation minet
#source("https:/bioconductor.org/biocLite.R")
#biocLite("minet")

#install.packages("igraph")
library(igraph)

graph1 <- graph(edges = c(4,7,7,2,2,8,2,6,6,1,1,3,1,10,10,5), n=10, directed=F )
#gr1_adj <- as_adj(graph1)#matrice d'adjacence, ligne = entree, colonne = sortie
gr1_adj <- as.matrix(as_adj(graph1),10,10) #matrice d'adjacence, ligne = entree, colonne = sortie
print(rowSums(gr1_adj))

#Verification
degrees = degree(graph1) 