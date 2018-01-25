library(igraph)
graph1 <- graph(edges = c(4,7,4,8,8,7,7,2,2,8,2,6,6,1,1,3,1,10,10,5,1,11,11,2), n=11, directed=F )
gr1_adj <- as.matrix(as_adj(graph1),11,11)
g1 <- graph(edges = c(1,2,1,3,2,4,3,4,4,5), n=5, directed=F )
g1_adj <- as.matrix(as_adj(g1),5,5)

gr_adj=matrix(c(0,1,1,0,0,1,0,0,1,0,1,0,0,1,0,0,1,1,0,1,0,0,0,1,0),5)


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




Betweenness=function(List_chemins,dist,n)
{
    B=list()
    for (x in 1:n)
    {
    B[length(B)+1]=list(matrix(0,n,n))
    }
    Tot=matrix(0,n,n)
    for (i in 2:length(List_chemins[[1]]))
    {
        for (j in 1:dim(List_chemins[[1]][[i]])[1])
        {
            for (x in 3:length(List_chemins[[1]][[i]][j,])-1)
            {
                v=as.numeric(List_chemins[[1]][[i]][j,x])
                s=as.numeric(List_chemins[[1]][[i]][j,1])
                t=as.numeric(List_chemins[[1]][[i]][j,length(List_chemins[[1]][[i]][j,])])
                B[[v]][s,t]=B[[v]][s,t]+1
            }
        }
    }
    for (i in 1:length(B))
    {
        Tot=Tot+(B[[i]])
    }
    for (i in 1:n)
    {
        for (j in 1:n)
        {
            Tot[i,j]=Tot[i,j]/(dist[[1]][i,j]-1)
        }
    }
    b=list()
    for (k in 1:n)
    {
        B[[k]]=B[[k]]/Tot
        for (i in 1:n)
        {
            B[[k]][i,][which(is.na(B[[k]][i,]))]=0
        }
        b=cbind(b,sum(B[[k]])/2)
    }
    return(b)
}
 
#A=PlCC(gr_adj)
#print(Betweenness(A[1],A[2],5))
#print (betweenness(g1))
#B=PlCC(gr1_adj)
#print(Betweenness(B[1],B[2],11))
#graph1 <- graph(edges = c(4,7,4,8,8,7,7,2,2,8,2,6,6,1,1,3,1,10,10,5,1,11,11,2), n=11, directed=F )
#print(betweenness(graph1))
#plot(graph1)

#http://www.dil.univ-mrs.fr/~tichit/rb/tp1/igraph_tutorial.html
g1 <- read.graph("http://www.dil.univ-mrs.fr/~tichit/rb/tp1/depts.txt", format="edgelist")
g2 <- read.graph("http://cneurocvs.rmki.kfki.hu/igraph/karate.net", format="pajek")
g3 <- read.graph("http://www.dil.univ-mrs.fr/~tichit/rb/tp1/depts.lgl", format="ncol")
plot (g3)
#tab<-read.csv("da11.csv", header=TRUE,sep=";",na.strings="NA",encoding="UTF-8")

#print(head(tab))
#print(dim(tab))
