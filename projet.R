#Command line : Rscript --vanilla projet.R -f <network-file> -e <format> [-t]

#!/usr/bin/env Rscript
#install.packages("optparse")
library("optparse")
#install.packages("igraph")
library(igraph)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Network file, required.", metavar="character"),
  make_option(c("-e", "--extension"), type="character", default=NULL, 
              help='Required. File format, choose one of the following : "edgelist", "pajek", "ncol", "lgl", "graphml",
  "dimacs", "graphdb", "gml", "dl".', metavar="character"),
  # make_option(c("-n", "--vertices"), default=0, type = "integer", 
  #             help="Number of network vertices. Ignored if not higher than the largest integer in file."),
  # make_option(c("-d", "--directed"), action="store_true", default=FALSE,
  #             help="When option is present, network is directed. /!\ if not present, network will be considered undirected!"),
  make_option(c("-t", "--traceback"), action="store_true", default=FALSE,
              help="When option is present, outputs the BFS and Floyd-Warshall tracebacks : one of the shortest paths between every 2 vertices.")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file) & is.null(opt$extension) )
{
  print_help(opt_parser)
  stop("At least one input file and one file format must be supplied", call.=FALSE)
}

formats = c("edgelist", "pajek", "ncol", "lgl", "graphml", "dimacs", "graphdb", "gml", "dl")

if ( ! opt$extension %in% formats)
{
  print_help(opt_parser)
  stop('File format must be among the following : "edgelist", "pajek", "ncol", "lgl", "graphml", "dimacs", "graphdb", "gml", "dl".', call.=FALSE)
}




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

#Degree distribution
degree_distrib = function(degree_gr)
{
  png("degree_distribution.png") #Caution! Replaces image if done multiple times
  plot(table(degree_gr)/length(degree_gr), type = "p", ylab = "freq", xlab = "degree") 
  dev.off()
}

#Prend une matrice d'adjacence et la rend symetrique
undirect_graph = function(gr_adj)
{
  undir = gr_adj + t(gr_adj)
  undir[which(undir > 1)] = 1
  return(undir)
}

#Local clustering coefficient
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
    # print(paste("node", node))
    node.voisins = which(gr_adj[node,] != 0)
    # print (paste("voisins", node.voisins))
    aretes.voisins = sum(gr_adj[node.voisins, which(gr_adj[node,] != 0)])/2 # matrice d'adjacence est symetrique donc il faut diviser par 2 le resultat
    # print (gr_adj[node.voisins, which(gr_adj[node,] != 0)])
    # print(paste("aretes voisins", aretes.voisins))
    cl = c(cl, aretes.voisins * 2/(deg[node] * (deg[node]-1))) 
  }

  # else
  # {
  #   deg = degrees(gr, symetric)[1,]
  #   for (node in V)
  #   {
  #     print(node)
  #       node.voisins = which(gr_adj[node,] != 0| gr_adj[,node] != 0)
  #       print (paste("voisins", node.voisins))
  #       aretes.voisins = sum(gr_adj[node.voisins, which(gr_adj[node,] != 0 | gr_adj[,node] != 0)]) # matrice d'adjacence est non symetrique
  #       print (gr_adj[node.voisins, which(gr_adj[node,] != 0 | gr_adj[,node] != 0)])
  #       print(paste("aretes voisins", aretes.voisins))
  #       cl = c(cl, aretes.voisins/(deg[node] * (deg[node]-1)))
  #   }
  # }
  retour <- cl #c(cl, mean(cl, na.rm = TRUE))
  # names(retour) <- c(V, "global")
  return (retour)
}

#Pour oriente ou non, on regarde les degrees out. Pour graphes non ponderes (poids de 1 si arc).
BFS = function(gr_adj, s)
{
  V = dim(gr_adj)[1]
  dist=(rep(Inf, V))
  dist[s]=0
  prev=(rep(-1,V))
  Q=c(s)
  while(length(Q)!=0)
  {
    u=Q[length(Q)]
    Q=Q[-length(Q)]
    for (v in which(gr_adj[u,]!=0))
    {
      if(dist[v]==Inf)
      {
        Q =c(Q, v)
        dist[v]=dist[u]+1
        prev[v]=u
      }
    }
  }
  return(rbind(dist, prev))
}

#Donne un des plus courts chemins entre 2 noeuds
traceback_BFS=function(dist,prev)
{
  chemins = list()
  s=which(dist==0)
  for (u in c(1:length(dist))[-s])
  {
    local = c(u)
    while(prev[u] != -1)
    {
      local = c(local, prev[u])
      u = prev[u]
    }
    local = rev(local)
    chemins[[length(chemins) +1]] = local
  }
  return(chemins)
}


all_BFS=function(gr_adj)
{
  dist_all = c()
  prev_all = c()
  tr_all = list()
  for(s in 1:dim(gr_adj)[1]) #Source
  {
    S = BFS(gr_adj, s)
    tr = traceback_BFS(S[1,], S[2,])
    dist_all = rbind(dist_all, S[1,])
    prev_all = rbind(prev_all, S[2,])
    tr_all[[length(tr_all) +1]] = tr
  }
  return(list(dist_all, prev_all, tr_all))
}

#distance_all = 1er element de la liste retournee par all_BFS
longest_shortest_path_BFS = function(distance_all)
{
  distance_all[distance_all == Inf] = - Inf
  distance_all[distance_all == 0] = - Inf #Non pondere donc si distance = 0 => pas de chemin, ne nous interesse pas
  return (max(distance_all))
}


#Pour oriente ou non, pondere ou non
floyd_warshall = function(gr_adj)
{
  gr2_adj = gr_adj
  V = dim(gr_adj)[1]
  gr_adj[gr_adj == 0] = Inf
  diag(gr_adj) = 0
  Next = t(matrix(1:V, V, V))
  Next[gr_adj == 0] = NA
  Next[gr_adj == Inf] = NA
  # Compt = matrix(0,V,V)
  for(k in 1:V) #noeud obligatoire du chemin
  {
    for(i in 1:V) #source
    {
      for(j in 1:V) #arrivee
      {
        if (gr_adj[i,j] > gr_adj[i,k]+gr_adj[k,j])
        {
          gr_adj[i,j] = gr_adj[i,k]+gr_adj[k,j]
          Next[i,j] = Next[i,k] #update si chemin passant par k entre i et j est plus court que celui qu'on a deja
          # if (i != j)
          # {
          #   Compt[i,j] = 1
          # }
        }
        # else if (k!=i && k !=j && i != j && gr_adj[i,j] != Inf && gr_adj[i,j] == gr_adj[i,k]+gr_adj[k,j])
        # {
        #   Compt[i,j] = Compt[i,j] + 1
        # }
      }
    }
  }
  # Compt[which(gr_adj == gr2_adj)] = Compt[which(gr_adj == gr2_adj)] + 1
  # diag(Compt) = 0
  # return(list(gr_adj,Next, Compt))
  return(list(gr_adj,Next))
}


# Next : 2e element de la liste rendue par floyd_warshall()
traceback_fw = function(Next)
{
  chemins = list()
  V = dim(Next)[1]
  for (u in 1:V)
  {
    for (v in 1:V)
    {
      if(anyNA(Next[u,v]))
        chemins[[length(chemins) +1]] = list()
      else
      {
        local = c(u)
        w = u
        while(w!=v)
        {
          local=c(local,Next[w,v])
          w = Next[w,v]
        }
        # print(local)
        chemins[[length(chemins) +1]] = local
      }
    }
  }
  return (chemins)
}

# traceback_fw = function(Next, compt)
# {
#   chemins = list()
#   V = dim(Next)[1]
#   for (u in 1:V)
#   {
#     for (v in 1:V)
#     {
#       if(anyNA(Next[u,v]))
#         chemins[[length(chemins) +1]] = list()
#       else
#       {
#         local = c(u)
#         w = u
#         while(w!=v)
#         {
#           if (compt[w,v] > compt[u,v])
#           {
#             compt[u,v] = compt[w,v]
#           }
#           local=c(local,Next[w,v])
#           w = Next[w,v]
#         }
#         print(local)
#         chemins[[length(chemins) +1]] = local
#       }
#     }
#   }
#   return (list(chemins, compt))
# }

# Next : 1er element de la liste rendue par floyd_warshall()
longest_shortest_path_fw = function(distance)
{
  distance[distance == Inf] = - Inf
  return (max(distance))
}


#Tests varies des fonctions, verifications par igraph
tests_varies = function()
{
  graph1 <- graph(edges = c(4,7,4,8,8,7,7,2,2,8,2,6,6,1,1,3,1,10,10,5,1,11,11,2), n=11, directed=F )
  #gr1_adj <- as_adj(graph1)#matrice d'adjacence, ligne = entree, colonne = sortie
  gr1_adj <- as.matrix(as_adj(graph1),11,11) #matrice d'adjacence, ligne = entree, colonne = sortie
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
  
  
  #Verifier si graphe oriente ou pas
  symetric = isSymmetric.matrix(gr1_adj)
  
  degree_gr1 = degrees(graph1, TRUE)
  #Verification
  print(clustering_coeff(g2, FALSE))
  print (transitivity(as.undirected(g2, "collapse"), "local"))
  
  S = BFS(gr1_adj,8)
  traceback_BFS(S[1,], S[2,])
  Res_all = all_BFS(gr1_adj)
}

#Output du traceback de BFS
output_traceback_bfs = function(trace, fichier)
{
  write(paste("#BFS Traceback : one of the shortest paths for each pair of vertexes"), append = TRUE, file = fichier)
  for (i in 1:(length(trace)[1]))
  {
    for (j in 1:(length(trace)-1))
    {
      if (length(trace[[i]][[j]]) > 1)
      {
        write(trace[[i]][[j]], sep = "\t", append = TRUE, file =fichier, ncolumns = length(trace[[i]][[j]]))
      }
    }
  }
}

#Output du traceback de F-W
output_traceback_fw = function(trace, fichier)
{
  write(paste("#Floyd-Warshall : one of the shortest paths for each pair of vertexes"), append = TRUE, file = fichier)
  for (i in 1:length(trace))
  {
    if (length(trace[[i]]) > 1)
    {
      write(trace[[i]], sep = "\t",append = TRUE, file =fichier, ncolumns = length(trace[[i]]))
    }
  }
}

# file : network file
main = function(file, tr, extension)
{
  #Input graph
  graphe = read_graph(file, format = extension)
  
  gr_Adj = as.matrix(as_adj(graphe)) #matrice d'adjacence
  
  #Verifier si graphe oriente ou pas
  symetric = isSymmetric.matrix(gr_Adj)
  #Verifier si graphe pondere ou pas
  is_weighted = (length(table(gr_Adj)) > 2)
  
  #Degrees
  graph_degree = degrees(graphe, symetric, gr_adj = gr_Adj)
  degree_distrib(graph_degree) #Outputs png picture
  
  #Local clustering coefficient
  cl_coef = clustering_coeff(graphe, symetric)
  
  #Shortest path via BFS
  if (! is_weighted)
  {
    sh_path_bfs = all_BFS(gr_Adj)
    longest_sh_path_bfs = longest_shortest_path_BFS(sh_path_bfs[[1]])
    print(sh_path_bfs)
  }
  
  #Shortest path via Floyd-Warshall
  sh_path_fw = floyd_warshall(gr_Adj)
  tr_fw = traceback_fw(sh_path_fw[[2]])
  longest_sh_path_fw = longest_shortest_path_fw(sh_path_fw[[1]])
  print(sh_path_fw)
  print(tr_fw)
  
  #Output
  df = as.data.frame(t(rbind(graph_degree, cl_coef)), row.names = paste("Node", 1:dim(gr_Adj)[1]))
  colnames(df) =c("Degree", "Local_transitivity")
  write.table(df, file = "output.txt", quote = FALSE, sep = "\t")
  if (! is_weighted)
  {
    write(c("#BFS Shortest distance matrix"), file ="output.txt", append = TRUE)
    matr = sh_path_bfs[[1]]
    colnames(matr) = paste("Node", 1:dim(gr_Adj)[1])
    rownames(matr) = paste("Node", 1:dim(gr_Adj)[1])
    write.table(matr,file ="output.txt", append = TRUE, quote = FALSE, sep = "\t")
    write(paste("#BFS Longest shortest path", longest_sh_path_bfs), sep = "\t", file ="output.txt", append = TRUE)
  }
  write(c("#Floyd-Warshall Shortest distance matrix"), file ="output.txt", append = TRUE)
  matr = sh_path_fw[[1]]
  colnames(matr) = paste("Node", 1:dim(gr_Adj)[1])
  rownames(matr) = paste("Node", 1:dim(gr_Adj)[1])
  write.table(matr,file ="output.txt", append = TRUE, quote = FALSE, sep = "\t")
  write(paste("#Floyd-Warshall Longest shortest path", longest_sh_path_fw), sep = "\t", file ="output.txt", append = TRUE)
  if (tr)
  {
    if (! is_weighted)
    {
      output_traceback_bfs(sh_path_bfs[[3]], "output_traceback.txt")
    }
    output_traceback_fw(tr_fw, "output_traceback.txt")
  }
}

# MAIN
main(opt$file, opt$traceback, opt$extension)

