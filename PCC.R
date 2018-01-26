library(igraph)
graph1 <- graph(edges = c(4,7,4,8,8,7,7,2,2,8,2,6,6,1,1,3,1,10,10,5,1,11,11,2), n=11, directed=F )
gr1_adj <- as.matrix(as_adj(graph1),11,11)
g1 <- graph(edges = c(1,2,1,3,2,4,3,4,4,5), n=5, directed=F )
g1_adj <- as.matrix(as_adj(g1),5,5)
print (g1_adj)
plot(g1)
PlCC=function(gr_adj)
{
	dist=gr_adj
	dist[which(dist == 0)] = Inf
	diag(dist)=0
	d=dim(dist)[1]
	print (dim(dist)[1])
	k=1
	PCC=c()
	pcc=c()
	for (x in 1:d)
	{
		for(y in 1:d)
		{
			if (dist[x,y]==k )
			{
				#pcc=rbind(pcc,c(x,y))
				#print(pcc)
			}
		}
	PCC=list(PCC,pcc)
	}
	print ('PCC')
	print (PCC)
	
	 
 	
 	while(length(which(dist == Inf))!=0 && k<4)
 	{
 		pcc=list()
 		for (u in 1:d)
 		{
 			for (v in 1:d)
 			{
 				if (dist[u,v]==k)
 				{
 					for (w in c(1:d)[-u][-v])
 					{
 						if (dist[v,w]==k )
 						{
 							dist[u,w]=k+1
 							for (C in PCC)
 							{
 							#print(PCC[1,])
 							}
 							#pcc=rbind(pcc,c(,y))
 						#print (dist[v,w])
 						#print (which(dist[v,w]==k))
 						#print(rbind(u,v,w))	
 						}
 						
 					}
 				}
 			}
 			#print(1 in (c(PCC[1,])))
 		}
 		k=k+1
	}
print(dist)
}

PlCC(g1_adj)