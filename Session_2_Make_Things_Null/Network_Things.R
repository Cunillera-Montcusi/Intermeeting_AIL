

# Hey!
# Welcome to the summer school course on river spatiotemporal connectivity!
# I am David Cunillera Montcusí and I am probably talking to you right now. SO focus and look at me! 

library(tidyverse);library(viridis)

# We charge the dataset -- Forget about rivers! THE TIME OF PONDS HAS COME! 
xy.Ponds <- read.csv2("data/Ponds.csv") %>% mutate(UTM_x=as.numeric(UTM_x),UTM_y=as.numeric(UTM_y))
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
  geom_edges( color = "grey60", size=0.1, alpha=0.4) +
  geom_nodes(shape=21,alpha=.75,colour="blue",fill="darkblue")+
  labs(x="",y="")+theme_classic()+
  theme(axis.text = element_blank(),axis.ticks = element_blank(),
        legend.position = "none",panel.background=element_blank())


### MST network - Minimum spanning tree
library(ape)
MST.Dist <- mst(Dist_matrix)
n<- network(MST.Dist, directed=F, diag=F)
ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges( color = "grey60", size=0.1, alpha=0.4) +
  geom_nodes(shape=21,alpha=.75,colour="blue")+
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

##### Diameter #####
# For MST 
g <- igraph::graph.adjacency(MST.Dist)

igraph::diameter(g)
igraph::get.diameter(g)->diam.1
as.vector(diam.1)->id.1

Diam_ID <- rep("NO",nrow(xy.Ponds))
Diam_ID[id.1] <- "YES"

n<- network(MST.Dist, directed=F, diag=F)
n %v% "Diameter" <- Diam_ID

ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges( color = "grey60", size=0.1, alpha=0.4) +
  geom_nodes(aes(fill=Diameter),shape=21, alpha=.75)+
  scale_fill_viridis(discrete = T,alpha = 1,begin = 1,end = 0)+
  labs(x="",y="")+theme_classic()+
  theme(axis.text = element_blank(),axis.ticks = element_blank(),
        legend.position = "none",panel.background=element_blank())

# For percol 
g <- igraph::graph.adjacency(Dist_percol)

igraph::diameter(g)
igraph::get.diameter(g)->diam.2
as.vector(diam.2)->id.2

Diam_ID <- rep("NO",nrow(xy.Ponds))
Diam_ID[id.2] <- "YES"

n<- network(Dist_percol, directed=F, diag=F)
n %v% "Diameter" <- Diam_ID

ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges( color = "grey60", size=0.1, alpha=0.4) +
  geom_nodes(aes(fill=Diameter),shape=21, alpha=.75)+
  scale_fill_viridis(discrete = T,alpha = 1,begin = 1,end = 0)+
  labs(x="",y="")+theme_classic()+
  theme(axis.text = element_blank(),axis.ticks = element_blank(),
        legend.position = "none",panel.background=element_blank())


##### Average path length #####
igraph::average.path.length(g)

# Connectivity
C_obs <- sum(Dist_percol)/2 ## Symetric newtork that is why we divide by two
N<-ncol(Dist_percol)
C_esp <- (N*(N-1))/2 # We calculate all the possible connections (N-1 to remove diagonal and /2 for diag=0)
C_obs/C_esp
# East way
igraph::edge_density(g)


# Modularity 
modul <- igraph::spinglass.community(g)
modul$membership->memb

Memb_ID <- LETTERS[memb]

n<- network(Dist_percol, directed=F, diag=F)
n %v% "Memb_ID" <- Memb_ID

ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges( color = "grey60", size=0.1, alpha=0.4) +
  geom_nodes(aes(fill=Memb_ID),shape=21, alpha=.75, size=4)+
  scale_fill_viridis(discrete = T,alpha = 1,begin = 1,end = 0)+
  labs(x="",y="")+theme_classic()+
  theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.background=element_blank())


# Membership IDs

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


Centr_Values_Pond <- igraph::vertex_attr(intergraph::asIgraph(n))
Matr_Centr_Values_Pond <- matrix(unlist(Centr_Values_Pond),ncol=length(Centr_Values_Pond))
colnames(Matr_Centr_Values_Pond) <- names(Centr_Values_Pond)

Matr_Centr_Values_Pond <- as.data.frame(Matr_Centr_Values_Pond[,1:8])%>% 
                          mutate_if(is.character, as.numeric) %>% 
                          mutate(Modul=Matr_Centr_Values_Pond[,9])

corrmorant::corrmorant(Matr_Centr_Values_Pond)
pairs(Matr_Centr_Values_Pond[,1:8])
dev.off()

Matr_Centr_Values_Pond %>% pivot_longer(cols = 1:8) %>% 
  ggplot()+
  geom_boxplot(aes(x=Modul,y=value))+
  geom_violin(aes(x=Modul,y=value, fill=Modul),alpha=0.5)+
  scale_fill_viridis(discrete = T, direction = -1)+
  theme_classic()+facet_grid(name~.,scales = "free_y")



# Modularity and topological roles
Mod.A <- bipartite::computeModules(Dist_percol)
ListMod <- bipartite::listModuleInformation(Mod.A)
bipartite::printoutModuleInformation(Mod.A)

roles_1 <- bipartite::czvalues(Mod.A, weighted=FALSE, level="lower")
plot(roles_1$z~roles_1$c)
abline(v=0.62) # threshold of Olesen et al. 2007
abline(h=2.5)   # threshold of Olesen et al. 2007

Memb_ID <- rep("NO", nrow(xy.Ponds))
Memb_ID[which(roles_1$c>0.62)] <- "Connector"
Memb_ID[max(roles_1$z)] <- "ModuleHub"

n<- network(Dist_percol, directed=F, diag=F)
n %v% "Memb_ID" <- Memb_ID

ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges( color = "grey60", size=0.1, alpha=0.4) +
  geom_nodes(aes(fill=Memb_ID),shape=21, alpha=.75, size=4)+
  scale_fill_viridis(discrete = T,alpha = 1,begin = 1,end = 0)+
  labs(x="",y="")+theme_classic()+
  theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.background=element_blank())


# STOP about ponds! Lets move to rivers!

# Network building and playing!
#1 
Dendritic <- igraph::graph.tree(n=101,children=2,mode="in")
Regular <- igraph::make_lattice(length = 10, dim = 2)
Random <- igraph::erdos.renyi.game(100,p=0.1,type="gnp")

#2#Visaulizacion
par(mfrow=c(1,3),mar=c(2,2,2,2),bty="l")
igraph::plot.igraph(Dendritic,vertex.size=5,vertex.color="red"
            ,vertex.label=NA,edge.color="black",edge.arrow.size=0.1)
igraph::plot.igraph(Regular,vertex.size=4,vertex.color="red",
            vertex.label=NA,edge.color="black",layout=igraph::layout.grid)
igraph::plot.igraph(Random,vertex.size=4,vertex.color="red",
            vertex.label=NA,edge.color="black",layout=igraph::layout.random)

dev.off()

#3 Centralities 
### Dendritic
degree_Dend <- igraph::degree(Dendritic,mode = "all")
degree_Dend_in <- igraph::degree(Dendritic,mode = "in")
degree_Dend_out <- igraph::degree(Dendritic,mode = "out")
closs_Dend <- igraph::closeness(Dendritic,mode = "all")
betw_Dend <- igraph::betweenness(Dendritic,directed = TRUE)

### Regular
degree_Reg <- igraph::degree(Regular,mode = "all")
degree_Reg_in <- igraph::degree(Regular,mode = "in")
degree_Reg_out <- igraph::degree(Regular,mode = "out")
closs_Reg <- igraph::closeness(Regular,mode = "all")
betw_Reg <- igraph::betweenness(Regular,directed = TRUE)

### Random
degree_Ran <- igraph::degree(Random,mode = "all")
degree_Ran_in <- igraph::degree(Random,mode = "in")
degree_Ran_out <- igraph::degree(Random,mode = "out")
closs_Ran <- igraph::closeness(Random,mode = "all")
betw_Ran <- igraph::betweenness(Random,directed = TRUE)

par(mfrow=c(3,5))
hist(degree_Dend);hist(degree_Dend_in);hist(degree_Dend_out);hist(log(closs_Dend+1));hist(log(betw_Dend+1))
hist(degree_Reg);hist(degree_Reg_in);hist(degree_Reg_out);hist(log(closs_Reg+1));hist(log(betw_Reg+1))
hist(degree_Ran);hist(degree_Ran_in);hist(degree_Ran_out);hist(log(closs_Ran+1));hist(log(betw_Ran+1))
dev.off()
