library(igraph)
graph1 <- graph(edges = c(4,7,4,8,8,7,7,2,2,8,2,6,6,1,1,3,1,10,10,5,1,11,11,2), n=11, directed=F )
gr1_adj <- as.matrix(as_adj(graph1),11,11)
g1 <- graph(edges = c(1,2,1,3,2,4,3,4,4,5), n=5, directed=F )
g1_adj <- as.matrix(as_adj(g1),5,5)
#plot(graph1)
PlCC=function(gr_adj)
{
	dist=gr_adj
	dist[which(dist == 0)] = Inf
	diag(dist)=0
	d=dim(dist)[1]
	k=1
	PCC=list()
	pcc=list()
	for (x in 1:d)
	{
		for(y in 1:d)
		{
			if (dist[x,y]==k )
			{
				pcc=rbind(pcc,c(x,y))
			}
		}
	}
	PCC[[length(PCC)+1]]=pcc
 	while(length(which(dist == Inf))!=0 && length(pcc)!=0)
 	{
 		pcc=list()
 		for (u in 1:d)
 		{
 			for (v in 1:d)
 			{
 				if (dist[u,v]==k)
 				{
 					for (w in c(1:d))
 					{
 						if (dist[v,w]==1 && w!=u && w!=v)
 						{
 						  for (i in 1:dim(PCC[[k]])[1])
 						  {
 						    if ((PCC[[k]][i,1])==u  && (PCC[[k]][i,length(PCC[[k]][i,])])==v && w%in%PCC[[k]][i,]==FALSE)
 						    {
 						      if (dist[u,w]==k+1|dist[u,w]==Inf)
 						      {
 						        pcc=rbind(pcc,c(PCC[[k]][i,],w))
 						        dist[u,w]=k+1
                  }
 						    }
 						  }
 						}
 					}
 				}
 			}
 		}
 		if (length(pcc)!=0)
 		{
 		  PCC[[length(PCC)+1]]=pcc
 		}
 		
 		k=k+1
	}
  return(list(PCC,dist))
}

A=PlCC(gr1_adj)
A1=PlCC(g1_adj)

# Liste <- list(m1 = matrix(1:4, ncol = 2), m2 = matrix(5:8, ncol = 2), m3 = matrix(9:12, ncol = 2))
Betweenness=function(List_chemins,n)
{
  B=list()
  for (x in 1:n)
  {
  B[length(B)+1]=list(matrix(0,n,n))
  }
  for (i in 2:length(List_chemins[[1]]))
  {
    for (j in 1:dim(List_chemins[[1]][[i]])[1])
    {
      print (List_chemins[[1]][[i]])

      for (n in 3:length(List_chemins[[1]][[i]][j,])-1)
      {
        print (paste('n',n))
        print(List_chemins[[1]][[i]][j,n])
        v=as.numeric(List_chemins[[1]][[i]][j,n])
        s=as.numeric(List_chemins[[1]][[i]][j,1])
        t=as.numeric(List_chemins[[1]][[i]][j,length(List_chemins[[1]][[i]][j,])])
        print (B[[v]])
        B[[v]][s,t]==B[[v]][s,t]+1
      }  
    }
    
  }
print (B)
}

Betweenness(A1[1],5)
print (betweenness(g1))

# print (5%in%List_chemins[[1]][[i]][j,])
# print (List_chemins[[1]][[i]][j,])
# tot=0
# s1=List_chemins[[1]][[i]][j,][1]
# t1=List_chemins[[1]][[i]][j,][length(List_chemins[[1]][[i]][j,])]
# s2=s1
# t2=t1
# b=(rep(0,n))
# if (s1==s2 && t1==t2)
# {
#   for (n in 2:(dim(List_chemins[[1]][[i]][j,])-1))
#   {
#     tot=tot+1
#     b[n]=b[n]+1
#   }
# }
# else
# {
#   tot=0
# }
# # b=(rep(0,n))
# # while (s1==s2 && t1==t2)
# # {
# #   tot=tot+1
# #   b[List_chemins[[1]][[i]][j,][n]]
# #   n+=1
# # }
