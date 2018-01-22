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

graph1 <- graph(edges = c(4,7,4,8,8,7,7,2,2,8,2,6,6,1,1,3,1,10,10,5), n=10, directed=F )
#gr1_adj <- as_adj(graph1)#matrice d'adjacence, ligne = entree, colonne = sortie
gr1_adj <- as.matrix(as_adj(graph1),10,10) #matrice d'adjacence, ligne = entree, colonne = sortie
print(rowSums(gr1_adj))
print(colSums(gr1_adj))

#Verification
degrees = degree(graph1) 
print(degrees)

g2 <- graph( edges=c(1,2, 2,3, 3,1, 1,4, 4,2,2,4, 3,2, 3,4,4,3,5,7,5,6,6,7,5,8), n=10 )
gr2_adj <- as.matrix(as_adj(g2)) #matrice d'adjacence, ligne = entree, colonne = sortie
print(rowSums(gr2_adj))
print(colSums(gr2_adj))

print(degree(g2, mode = c("out")))
print(degree(g2, mode = c("in")))
print(degree(g2))

#Ajouter distribution noeuds


#Verifier si graphe oriente ou pas
symetric = isSymmetric.matrix(gr1_adj)

#Takes a graph and a boolean as input and gives the degree, and if True, the degree in/out.
degrees = function(gr, symetric, gr_adj = "")
{
  if (! is.matrix(gr_adj))
  {
    gr_adj <- as.matrix(as_adj(gr))
  }
  if (! symetric)
  {
    return(rbind(colSums(gr_adj) + rowSums(gr_adj), colSums(gr_adj), rowSums(gr_adj))) #global, in, out
  }
  else
  {
    return(colSums(gr_adj))
  }
}

#Prend une matrice d'adjacence et la rend symetrique
undirect_graph = function(gr_adj)
{
  undir = gr_adj + t(gr_adj)
  undir[which(undir > 1)] = 1
  return(undir)
}

#local clustering coefficient - GLOBAL A FAIRE
clustering_coeff = function(gr, symetric)
{
  gr_adj <- as.matrix(as_adj(gr))
  diag(gr_adj) <- 0 #Pas d'aretes sur lui meme
  cl = c()
  V = 1:dim(gr_adj)[1]
  if (! symetric)
  {
    gr_adj = undirect_graph(gr_adj)
  }
  
  deg = degrees(gr, TRUE, gr_adj)
  for (node in V)
  {
    print(paste("node", node))
    node.voisins = which(gr_adj[node,] == 1)
    print (paste("voisins", node.voisins))
    aretes.voisins = sum(gr_adj[node.voisins, which(gr_adj[node,] == 1)])/2 # matrice d'adjacence est symetrique donc il faut diviser par 2 le resultat
    print (gr_adj[node.voisins, which(gr_adj[node,] == 1)])
    print(paste("aretes voisins", aretes.voisins))
    cl = c(cl, aretes.voisins * 2/(deg[node] * (deg[node]-1))) 
  }

  # else
  # {
  #   deg = degrees(gr, symetric)[1,]
  #   for (node in V)
  #   {
  #     print(node)
  #       node.voisins = which(gr_adj[node,] == 1| gr_adj[,node] == 1)
  #       print (paste("voisins", node.voisins))
  #       aretes.voisins = sum(gr_adj[node.voisins, which(gr_adj[node,] == 1 | gr_adj[,node] == 1)]) # matrice d'adjacence est non symetrique
  #       print (gr_adj[node.voisins, which(gr_adj[node,] == 1 | gr_adj[,node] == 1)])
  #       print(paste("aretes voisins", aretes.voisins))
  #       cl = c(cl, aretes.voisins/(deg[node] * (deg[node]-1)))
  #   }
  # }
  retour <- cl#c(cl, mean(cl, na.rm = TRUE))
  # names(retour) <- c(V, "global")
  return (retour)
}

#Verification
print(clustering_coeff(g2, FALSE))
print (transitivity(as.undirected(g2, "collapse"), "local"))
