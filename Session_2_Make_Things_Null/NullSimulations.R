
library(tidyverse);library(viridis)
# START ####
load("CESTES_Null.RData")
databae <- LSmatred[[1]]
databae$coord
databae$envir
databae$comm

storage_outs <- data.frame()
source("function_to_null.R")

# EXERCICE 1 ####


### WE BEGUN TO built all the data that we need to run the simulation
## Distance matrix 
# It is equivalent to the Riv_STconmat 
Dist_Matrix <- dist(databae$coord)

## Community size
J.freshwater <- rep(400,nrow(databae$coord))

# Species tolerances defined in the 2. script
Species_Filter <- matrix(nrow=200,ncol = nrow(databae$coord), data = 0.99)

# Filter of species per site. We will use for tolerance
# Other parameters of the model
id_NOmodule <- rep(1,nrow(databae$coord)) # Modules if we want some sites to belong to the same module. 
pool_200 <- rep(1,nrow(Species_Filter)) # Distribution of the species pool #rlnorm(n = 200,5,1) 

# Distances are related to the distance matrix. 
# Therefore we must see which values correspond to our connections to set the "D50". Dispersal distance is
# considered here as the distance at which probability of dispersal is 0.5. So "higher" or "lower" do not exclude
# other dispersal abilities. They push "overall connectivity" towards higher or lower connections.
summary(as.vector(Dist_Matrix)) 
# We set the dispersal abilities that we want
dispersal <- mean(as.vector(Dist_Matrix))

a <- NULL # We create an output object for each iteration
b <- list()# We create an output object for each iteration
for (it in 1:10) { # We repeat 10 times the same process
output <- Coalescent.exp.Kernel.J(Meta.pool=pool_200,m.pool=0.001, 
                                  Js=J.freshwater, 
                                  filter.env=Species_Filter,
                                  M.dist=as.matrix(Dist_Matrix), 
                                  D50=42406,m.max=1,
                                  id.obs=1:nrow(databae$coord))
a <- rbind(a,output[[1]])
b[[it]] <- output[[2]]
}
Final_out <- resume.out(a)

Final_out$Median
Final_out$`Standard Deviation`

Sim_S <- Final_out$Median[6:(nrow(databae$coord)+5)]
Sim_B <- Final_out$Median[(nrow(databae$coord)+6):((nrow(databae$coord)*2)+5)]
Sim_G <- Final_out$Median[((nrow(databae$coord)*2)+6):length(Final_out$Median)]

True_S <- apply(ifelse(t(databae$comm)>0,1,0),2,sum)
BB<-as.matrix(vegdist(databae$comm, method = "jaccard"))
True_B <- ifelse(is.na(apply(BB,2,mean))==T,0,apply(BB,2,mean))
True_G <- c(length(which(apply(ifelse(t(databae$comm)>0,1,0),1,sum)>1)),
diversity(apply(t(databae$comm),1,sum),"simpson"),
diversity(apply(t(databae$comm),1,sum)))

plot(Sim_S,True_S)
plot(Sim_B,True_B)
plot(Sim_G,True_G)

storage_outs <- storage_outs %>% bind_rows(data.frame("Test"="Ex1",Sim_S,Sim_B,True_S,True_B))

# EXERCICE 2 ####

### WE BEGUN TO built all the data that we need to run the simulation
## Distance matrix 
# It is equivalent to the Riv_STconmat 
Dist_Matrix <- dist(databae$coord)

## Community size
apply(t(databae$comm),2,sum)
J.freshwater <- apply(t(databae$comm),2,sum)

# Species tolerances defined in the 2. script
Species_Filter <- matrix(nrow=200,ncol = nrow(databae$coord), data = 0.99)

# Filter of species per site. We will use for tolerance
# Other parameters of the model
id_NOmodule <- rep(1,nrow(databae$coord)) # Modules if we want some sites to belong to the same module. 
pool_200 <- rep(1,nrow(Species_Filter)) # Distribution of the species pool #rlnorm(n = 200,5,1) 

# Distances are related to the distance matrix. 
# Therefore we must see which values correspond to our connections to set the "D50". Dispersal distance is
# considered here as the distance at which probability of dispersal is 0.5. So "higher" or "lower" do not exclude
# other dispersal abilities. They push "overall connectivity" towards higher or lower connections.
summary(as.vector(Dist_Matrix)) 
# We set the dispersal abilities that we want
dispersal <- mean(as.vector(Dist_Matrix))

a <- NULL # We create an output object for each iteration
b <- list()# We create an output object for each iteration
for (it in 1:10) { # We repeat 10 times the same process
  output <- Coalescent.exp.Kernel.J(Meta.pool=pool_200,m.pool=0.001, 
                                    Js=J.freshwater, 
                                    filter.env=Species_Filter,
                                    M.dist=as.matrix(Dist_Matrix), 
                                    D50=42406,m.max=1,
                                    id.obs=1:nrow(databae$coord))
  a <- rbind(a,output[[1]])
  b[[it]] <- output[[2]]
}
Final_out <- resume.out(a)

Final_out$Median
Final_out$`Standard Deviation`

Sim_S <- Final_out$Median[6:(nrow(databae$coord)+5)]
Sim_B <- Final_out$Median[(nrow(databae$coord)+6):((nrow(databae$coord)*2)+5)]
Sim_G <- Final_out$Median[((nrow(databae$coord)*2)+6):length(Final_out$Median)]

True_S <- apply(ifelse(t(databae$comm)>0,1,0),2,sum)
BB<-as.matrix(vegdist(databae$comm, method = "jaccard"))
True_B <- ifelse(is.na(apply(BB,2,mean))==T,0,apply(BB,2,mean))
True_G <- c(length(which(apply(ifelse(t(databae$comm)>0,1,0),1,sum)>1)),
            diversity(apply(t(databae$comm),1,sum),"simpson"),
            diversity(apply(t(databae$comm),1,sum)))

plot(Sim_S,True_S)
plot(Sim_B,True_B)
plot(Sim_G,True_G)

storage_outs <- storage_outs %>% bind_rows(data.frame("Test"="Ex2",Sim_S,Sim_B,True_S,True_B))

# EXERCICE 3 ####

### WE BEGUN TO built all the data that we need to run the simulation
## Distance matrix 
# It is equivalent to the Riv_STconmat 
Dist_Matrix <- dist(databae$coord)

## Community size
apply(t(databae$comm),2,sum)
J.freshwater <- apply(t(databae$comm),2,sum)

# Species tolerances defined in the 2. script
Species_Filter <- matrix(nrow=nrow(t(databae$comm)),ncol = nrow(databae$coord), data = 0.99)

# Filter of species per site. We will use for tolerance
# Other parameters of the model
id_NOmodule <- rep(1,nrow(databae$coord)) # Modules if we want some sites to belong to the same module. 
pool_200 <- rep(1,nrow(Species_Filter)) # Distribution of the species pool #rlnorm(n = 200,5,1) 

# Distances are related to the distance matrix. 
# Therefore we must see which values correspond to our connections to set the "D50". Dispersal distance is
# considered here as the distance at which probability of dispersal is 0.5. So "higher" or "lower" do not exclude
# other dispersal abilities. They push "overall connectivity" towards higher or lower connections.
summary(as.vector(Dist_Matrix)) 
# We set the dispersal abilities that we want
dispersal <- mean(as.vector(Dist_Matrix))

a <- NULL # We create an output object for each iteration
b <- list()# We create an output object for each iteration
for (it in 1:10) { # We repeat 10 times the same process
  output <- Coalescent.exp.Kernel.J(Meta.pool=pool_200,m.pool=0.001, 
                                    Js=J.freshwater, 
                                    filter.env=Species_Filter,
                                    M.dist=as.matrix(Dist_Matrix), 
                                    D50=42406,m.max=1,
                                    id.obs=1:nrow(databae$coord))
  a <- rbind(a,output[[1]])
  b[[it]] <- output[[2]]
}
Final_out <- resume.out(a)

Final_out$Median
Final_out$`Standard Deviation`

Sim_S <- Final_out$Median[6:(nrow(databae$coord)+5)]
Sim_B <- Final_out$Median[(nrow(databae$coord)+6):((nrow(databae$coord)*2)+5)]
Sim_G <- Final_out$Median[((nrow(databae$coord)*2)+6):length(Final_out$Median)]

True_S <- apply(ifelse(t(databae$comm)>0,1,0),2,sum)
BB<-as.matrix(vegdist(databae$comm, method = "jaccard"))
True_B <- ifelse(is.na(apply(BB,2,mean))==T,0,apply(BB,2,mean))
True_G <- c(length(which(apply(ifelse(t(databae$comm)>0,1,0),1,sum)>1)),
            diversity(apply(t(databae$comm),1,sum),"simpson"),
            diversity(apply(t(databae$comm),1,sum)))

plot(Sim_S,True_S)
plot(Sim_B,True_B)
plot(Sim_G,True_G)
abline(a = 0,b = 1)

storage_outs <- storage_outs %>% bind_rows(data.frame("Test"="Ex3",Sim_S,Sim_B,True_S,True_B))


# EXERCICE 4 ####

### WE BEGUN TO built all the data that we need to run the simulation
## Distance matrix 
# It is equivalent to the Riv_STconmat 
Dist_Matrix <- dist(databae$coord)

## Community size
apply(t(databae$comm),2,sum)
J.freshwater <- apply(t(databae$comm),2,sum)

# Species tolerances defined in the 2. script
env.rda <- vegan::rda(databae$envir, scale = TRUE)
summary(env.rda)

plot(env.rda,scaling = 1,
     display = c("sp", "lc", "cn"),
     main = "PCA")

spe.rda <- vegan::rda(decostand(databae$comm,"hellinger"),databae$envir, scale = TRUE)
summary(spe.rda)

plot(spe.rda,scaling = 1,display = c("sp", "lc", "cn"),
     main = "Triplot RDA spe.hel - scaling 1 - lc scores")

Species_Filter <- matrix(nrow=nrow(t(databae$comm)),ncol = nrow(databae$coord), data = NA)
Pred_Filt <- predict(spe.rda)
for(site_filt in 1:nrow(databae$coord)){
Species_Filter[,site_filt] <- scales::rescale(t(Pred_Filt)[,site_filt],to=c(0.0001,0.99))
}

# Filter of species per site. We will use for tolerance
# Other parameters of the model
pool_200 <- rep(1,nrow(Species_Filter)) # Distribution of the species pool #rlnorm(n = 200,5,1) 

# Distances are related to the distance matrix. 
# Therefore we must see which values correspond to our connections to set the "D50". Dispersal distance is
# considered here as the distance at which probability of dispersal is 0.5. So "higher" or "lower" do not exclude
# other dispersal abilities. They push "overall connectivity" towards higher or lower connections.
summary(as.vector(Dist_Matrix)) 
# We set the dispersal abilities that we want
dispersal <- mean(as.vector(Dist_Matrix))

a <- NULL # We create an output object for each iteration
b <- list()# We create an output object for each iteration
for (it in 1:10) { # We repeat 10 times the same process
  output <- Coalescent.exp.Kernel.J(Meta.pool=pool_200,m.pool=0.001, 
                                    Js=J.freshwater, 
                                    filter.env=Species_Filter,
                                    M.dist=as.matrix(Dist_Matrix), 
                                    D50=42406,m.max=1,
                                    id.obs=1:nrow(databae$coord))
  a <- rbind(a,output[[1]])
  b[[it]] <- output[[2]]
}
Final_out <- resume.out(a)

Final_out$Median
Final_out$`Standard Deviation`

Sim_S <- Final_out$Median[6:(nrow(databae$coord)+5)]
Sim_B <- Final_out$Median[(nrow(databae$coord)+6):((nrow(databae$coord)*2)+5)]
Sim_G <- Final_out$Median[((nrow(databae$coord)*2)+6):length(Final_out$Median)]

True_S <- apply(ifelse(t(databae$comm)>0,1,0),2,sum)
BB<-as.matrix(vegdist(databae$comm, method = "jaccard"))
True_B <- ifelse(is.na(apply(BB,2,mean))==T,0,apply(BB,2,mean))
True_G <- c(length(which(apply(ifelse(t(databae$comm)>0,1,0),1,sum)>1)),
            diversity(apply(t(databae$comm),1,sum),"simpson"),
            diversity(apply(t(databae$comm),1,sum)))

plot(Sim_S,True_S)
plot(Sim_B,True_B)
plot(Sim_G,True_G)

storage_outs <- storage_outs %>% bind_rows(data.frame("Test"="Ex4",Sim_S,Sim_B,True_S,True_B))

# EXERCICE 5 ####
a <- storage_outs %>% 
  ggplot()+geom_point(aes(x=Sim_S,y=True_S,colour=Test))+
  geom_smooth(aes(x=Sim_S,y=True_S),method="lm")+
  facet_grid(.~Test,scales="free_x")
a

b <- storage_outs %>% 
  ggplot()+geom_point(aes(x=Sim_B,y=True_B,colour=Test))+
  geom_smooth(aes(x=Sim_B,y=True_B),method="lm")+
  facet_grid(.~Test,scales="free_x")
b

gridExtra::grid.arrange(a,b,nrow=2)
