
library(FD);library(tidyverse)

load("CESTES_Null.RData")
databae <- LSmatred[[1]]

databae$comm
databae$traits

# Distance-based functional diversity indices
?dbFD
res <-dbFD(
    databae$traits,
    as.matrix(databae$comm),
    w.abun=TRUE,
    stand.FRic = TRUE,
    corr = "cailliez",
    clust.type = "ward.D",
    CWM.type = "all",
    calc.FGR = TRUE)
# g # cut the dendrogram using the number of groups as criterion
# 10 # choose the number of functional groups

# FEVe: Could not be calculated for communities with <3 functionally singular species. 
# FDis: Equals 0 in communities with only one functionally singular species. 
# FRic: To respect s > t, FRic could not be calculated for communities with <3 functionally singular species. 
# FRic: Dimensionality reduction was required. The last 13 PCoA axes (out of 15 in total) were removed. 
# FRic: Quality of the reduced-space representation = 0.3939701 
# FDiv: Could not be calculated for communities with <3 functionally singular species. 


# Using a distance-based approach, Villéger et al. (2008) proposed to distinguish three
# independent components in the functional diversity of a given community: functional
# richness (FRic) represents the amount of functional space filled by a community,
# i.e. the volume of the minimum convex hull that includes all species present in
# the multidimensional space of functional traits; functional evenness (FEve) measures
# the regularity of the abundance distribution of the species along the minimum
# spanning tree that links the species points in multidimensional functional space;
# functional divergence (FDiv) relates to how species abundances are distributed
# within the functional trait space. FEve and FDiv are constrained between 0 and

res

bind_cols("FRic"=res$FRic,"FEve"=res$FEve,"FDiv"=res$FDiv,"FDis"=res$FDis)

Out_FD %>% 
  pivot_longer(cols=4:7) %>% 
  ggplot(aes(x=Etacio,y=value))+geom_boxplot()+theme_classic()+facet_grid(.~name)


