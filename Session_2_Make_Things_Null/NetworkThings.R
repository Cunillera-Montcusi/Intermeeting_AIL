

# Hey! Welcome to the J-AIL intermeeting course!
# I am David Cunillera Montcusí and I am probably talking to you right now. SO focus and look at me!* 
#*for those who already saw this joke before, even though it is not that funny, the same applies... SO focus and look at me!
library(tidyverse);library(viridis)

# We charge the dataset -- Forget about rivers! THE TIME OF PONDS HAS COME! 
xy.Ponds <- read.csv2("Ponds.csv") %>% mutate(UTM_x=as.numeric(UTM_x),UTM_y=as.numeric(UTM_y))
xy.Ponds <- xy.Ponds[1:200,]

### Distance matrix
Dist_matrix <- as.matrix(dist(xy.Ponds[,3:4]))
# For lat and long geosphere::distm()

# Package sna is one of the two packages that we can use to calculate network metrics and create them
library(sna);library(ggnetwork)
# Percolation distance is your friend for local & regional levels
Dist_percol <- ifelse(as.matrix(Dist_matrix)>2800,0,1)
diag(Dist_percol) <- 0

n<- network(Dist_percol, directed=F, diag=F)
ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges( color = "grey60", linewidth=0.1, alpha=0.4) +
  geom_nodes(shape=21,alpha=.75,colour="blue",fill="darkblue")+
  labs(x="",y="")+theme_classic()+
  theme(axis.text = element_blank(),axis.ticks = element_blank(),
        legend.position = "none",panel.background=element_blank())


# sna
# works directly with the matrix
Dist_percol
sna::degree(Dist_percol)

# igraph
g <- igraph::graph.adjacency(Dist_percol)
igraph::degree(g)

sna::degree(Dist_percol)-igraph::degree(g)

## Degree: number of links of each site
# UnWeigheted
sna::degree(Dist_percol,gmode="graph",cmode="freeman",ignore.eval=TRUE)
n %v% "Degree_Unw" <- sna::degree(Dist_percol,gmode="graph",cmode="freeman",ignore.eval=TRUE)

# Weighted degree
sna::degree((Dist_percol*Dist_matrix),gmode="graph", cmode="freeman",ignore.eval=FALSE)
n %v% "Degree_Weig" <- sna::degree((Dist_percol*Dist_matrix),gmode="graph", cmode="freeman",ignore.eval=FALSE)

gridExtra::grid.arrange(
  ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges( color = "grey60", size=0.1, alpha=0.4) +
    geom_nodes(aes(fill=Degree_Unw),shape=21, alpha=.75,size=5)+
    scale_fill_viridis(discrete = F,alpha = 1,begin = 1,end = 0)+
    labs(x="",y="")+theme_classic()+
    theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.background=element_blank()),
  
  ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges( color = "grey60", size=0.1, alpha=0.4) +
    geom_nodes(aes(fill=Degree_Weig),shape=21, alpha=.75,size=5)+
    scale_fill_viridis(discrete = F,alpha = 1,begin = 1,end = 0)+
    labs(x="",y="")+theme_classic()+
    theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.background=element_blank()),
  ncol=2)

# Betweenness
# Unweighted betweenness
sna::betweenness(Dist_percol,gmode="graph",cmode="undirected",ignore.eval=TRUE)
n %v% "Betwee_Unw" <- sna::betweenness(Dist_percol,gmode="graph",cmode="undirected",ignore.eval=TRUE)

# Unweighted betweenness
sna::betweenness((Dist_percol*Dist_matrix),gmode="graph",cmode="undirected",ignore.eval=FALSE)
n %v% "Betwee_Weig" <- sna::betweenness((Dist_percol*Dist_matrix),gmode="graph",cmode="undirected",ignore.eval=FALSE)

gridExtra::grid.arrange(
  ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges( color = "grey60", size=0.1, alpha=0.4) +
    geom_nodes(aes(fill=Betwee_Unw),shape=21, alpha=.75,size=5)+
    scale_fill_viridis(discrete = F,alpha = 1,begin = 1,end = 0)+
    labs(x="",y="")+theme_classic()+
    theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.background=element_blank()),
  
  ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges( color = "grey60", size=0.1, alpha=0.4) +
    geom_nodes(aes(fill=Betwee_Weig),shape=21, alpha=.75,size=5)+
    scale_fill_viridis(discrete = F,alpha = 1,begin = 1,end = 0)+
    labs(x="",y="")+theme_classic()+
    theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.background=element_blank()),
  ncol=2)


##CLOSENESS
# Unweighted closenness
sna::closeness(Dist_percol,gmode="graph",ignore.eval=TRUE,cmode="undirected")
n %v% "Clos_Unw" <- sna::closeness(Dist_percol,gmode="graph",ignore.eval=TRUE,cmode="undirected")

# Weighted closenness
sna::closeness((Dist_percol*Dist_matrix),gmode="graph",ignore.eval=FALSE,cmode="undirected")
n %v% "Clos_Weig" <- sna::closeness((Dist_percol*Dist_matrix),gmode="graph",ignore.eval=FALSE,cmode="undirected")->clos_w

gridExtra::grid.arrange(
  ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges( color = "grey60", size=0.1, alpha=0.4) +
    geom_nodes(aes(fill=Clos_Unw),shape=21, alpha=.75,size=5)+
    scale_fill_viridis(discrete = F,alpha = 1,begin = 1,end = 0)+
    labs(x="",y="")+theme_classic()+
    theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.background=element_blank()),
  
  ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges( color = "grey60", size=0.1, alpha=0.4) +
    geom_nodes(aes(fill=Clos_Weig),shape=21, alpha=.75,size=5)+
    scale_fill_viridis(discrete = F,alpha = 1,begin = 1,end = 0)+
    labs(x="",y="")+theme_classic()+
    theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.background=element_blank()),
  ncol=2)


##EIGENVECTOR centrality
# Unweighted closenness
sna::evcent(Dist_percol,gmode="graph",ignore.eval=TRUE)
n %v% "EgV_Unw" <- sna::evcent(Dist_percol,gmode="graph",ignore.eval=TRUE)

# Weighted closenness
sna::evcent(Dist_percol*Dist_matrix,gmode="graph",ignore.eval=TRUE)
n %v% "EgV_Weig" <- sna::evcent(Dist_percol*Dist_matrix,gmode="graph",ignore.eval=TRUE)

gridExtra::grid.arrange(
  ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges( color = "grey60", size=0.1, alpha=0.4) +
    geom_nodes(aes(fill=EgV_Unw),shape=21, alpha=.75,size=5)+
    scale_fill_viridis(discrete = F,alpha = 1,begin = 1,end = 0)+
    labs(x="",y="")+theme_classic()+
    theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.background=element_blank()),
  
  ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges( color = "grey60", size=0.1, alpha=0.4) +
    geom_nodes(aes(fill=EgV_Weig),shape=21, alpha=.75,size=5)+
    scale_fill_viridis(discrete = F,alpha = 1,begin = 1,end = 0)+
    labs(x="",y="")+theme_classic()+
    theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.background=element_blank()),
  ncol=2)

