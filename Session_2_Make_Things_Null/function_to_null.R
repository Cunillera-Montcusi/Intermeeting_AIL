
###################################################################################################
Coalescent.exp.Kernel.J<-function(Meta.pool, m.pool, Js,filter.env,M.dist, D50, m.max,id.obs){
  
  library(vegan)
  Meta.pool<-Meta.pool/sum(Meta.pool)
  Meta<-NULL
  
  # Function to update: M.migra from Graph to a function based in distance matrix
  M.migra<-migration.kernel.all(M.dist=M.dist, m.pool=m.pool, D50=D50, m.max=m.max)# Funciton defined above. It estimates Migration matrix

  for(i in 1:ncol(M.migra)){Meta<-cbind(Meta, rmultinom(1,1,Meta.pool))}
  
  # If J=0 there is no spp. However, Meta add the possibility of having species there leading to at least 1 species to be found. 
  col_NO_J <- which(Js==0)
  if(length(col_NO_J)==0){col_NO_J <- -(1:length(Js))}
  Meta[,col_NO_J] <- 0
  
  Meta_sml <- Meta[,-col_NO_J]
  M.migra_sml <- M.migra[-col_NO_J,-col_NO_J]
  filter.env_sml <- filter.env[,-col_NO_J]

  for (ii in 2:max(Js)){
    id.j<-which(Js[-col_NO_J]>=ii)
    cat("coalescent construction in J: ", ii," de" ,max(Js),"\n")  
    
    Pool.neighbor<-(Meta_sml%*%M.migra_sml)     # estimates potential reclutants including immigrants for all communities weighted by local abundances
    Pool.neighbor<-Pool.neighbor*filter.env_sml # IMPORTANT: element by element adjustment of species abundances to local filters
    Pool.neighbor <- (Pool.neighbor/max(Pool.neighbor))
    
    if(length(id.j)>1){
      new<-apply(Pool.neighbor[,id.j],2,born, M.pool = Meta.pool, m.pool = m.pool)   # random selection of new individuals from reclutants pool 
      Meta_sml[,id.j]<-Meta_sml[,id.j]+new} else {
      Meta_sml[,id.j]<-Meta_sml[,id.j]+born(n = Pool.neighbor[,id.j], M.pool = Meta.pool, m.pool = m.pool) 
      }                          # upadate communities 
  }

  Meta[,-col_NO_J] <- Meta_sml
  
  MetaCom_t <- Meta
  BB<-as.matrix(vegdist(t(Meta[,id.obs]), method = "jaccard"))
  Bett.all<-apply(BB,2,mean)

  Meta<-list("out"=c("m.pool"=m.pool, "Js.max"=max(Js),"Js.min"=min(Js), "D50"=D50, "m.max"=m.max,
                     "S.loc"=apply(ifelse(Meta[,id.obs]>0,1,0),2,sum),
                     "B.loc.all"=ifelse(is.na(Bett.all)==T,0,Bett.all),
                     "G"=length(which(apply(ifelse(Meta[,id.obs]>0,1,0),1,sum)>1)),
                     "simp"=diversity(apply(Meta[,id.obs],1,sum),"simpson"),
                     "inv.simp"=diversity(apply(Meta[,id.obs],1,sum), "invsimpson")),
  "MetaCom"=MetaCom_t) 
  Meta
}

####


###########################################
# Cells at contact in their edges are assumed the zero distance for migration
# m.max: migration between connected cells
# D50 is the distance at which migration decay yo half m.max

migration.kernel.all<-function(M.dist, m.pool, D50,m.max){
  diag(M.dist)=NA
  M.dist = M.dist-min(M.dist, na.rm=T) # min distance is the distance between neighbour cells. 
  # connected cells has distance zero and migration m.max
  b = -log(0.5)/D50         # b is estimated | m(D50)=m.max*0.5
  M.migra = m.max*exp(-b*M.dist) 
  
  diag(M.migra)<-1                                             # selfrecruitment is considered as 1=m.intra community
  M.migra<-apply(M.migra,2,m_to_1,m.pool)                      # standirize migrations to 1: (m.intra+m.pool+m.neigh=1)
  M.migra
}

############
m_to_1<-function(m, m.pool) (1-m.pool)*m/sum(m) # standarize a vector of migration to add 1 
# also considering migration from an external pool

born<-function(n, M.pool, m.pool)rmultinom(1,1,((1-m.pool)*(n/sum(n))+m.pool*M.pool))
change<-function(n,change)rmultinom(1,change,n)

######################################################################################################################
# function to resume output of simulation
resume.out<-function(out){
  out2<-list()
  out2[[1]]<-apply(out,2,quantile, 0.5, na.rm=T)    # NA originates if only a single module is involved
  out2[[2]]<-apply(out,2,sd, na.rm=T)
  out2[[3]]<-apply(out,2,quantile, 0.975, na.rm=T)
  out2[[4]]<-apply(out,2,quantile, 0.025, na.rm=T)
  names(out2)<-c("Median", "Standard Deviation", "out.IC.up","out.IC.inf")
  out2
}
