# Packages charging 
library(vegan);library(hillR);library(adespatial);library(ade4)


# START ####
load("CESTES_Null.RData")
databae <- LSmatred[[1]]
databae$coord
databae$envir
databae$comm

# Checking 
str(databae$comm)
summary(databae$comm)

#_____________________________________________________      
# PCA environment ####
#_____________________________________________________      
env.pca <- rda(databae$envir)
summary(env.pca, scaling=2)
summary(env.pca, scaling=1)

screeplot(env.pca,bstick = T,npcs = length(env.pca$CA$eig))

par(mfrow = c(1, 2))
biplot(env.pca,scaling=1, main = "PCA - scaling 1")
biplot(env.pca, main = "PCA - scaling 2")
dev.off()

# Hellinger pre-transformation of the species data
spe.h <- decostand(databae$comm, "hellinger")
spe.h.pca <- rda(spe.h)
spe.h.pca
# Scree plot and broken stick model
screeplot(
  spe.h.pca,
  bstick = TRUE,
  npcs = length(spe.h.pca$CA$eig)
)

par(mfrow = c(1, 2))
biplot(spe.h.pca,scaling=1, main = "PCA - scaling 1")
biplot(spe.h.pca,scaling=2, main = "PCA - scaling 2")
dev.off()

#_____________________________________________________      
# RDA species environment####
#_____________________________________________________      
# Lines to create an RDA species vs environment

spe.h.rda <- rda(spe.h,databae$envir)
spe.h.rda
summary(spe.h.rda)
summary(spe.h.rda,scaling = 1)

# Unadjusted R^2 retrieved from the rda object
RsquareAdj(spe.h.rda)$r.squared
# Adjusted R^2 retrieved from the rda object
RsquareAdj(spe.h.rda)$adj.r.squared

par(mfrow = c(1, 2))
plot(spe.h.rda,
     scaling = 1,
     display = c("sp", "lc", "cn"),
     main = "Triplot RDA spe.hel - scaling 1 - lc scores")
plot(spe.h.rda,
     display = c("sp", "lc", "cn"),type=
     main = "Triplot RDA spe.hel - scaling 2 - lc scores")
dev.off()

## Global test of the RDA result
anova(spe.h.rda, permutations = how(nperm = 999))

R2a.all <- RsquareAdj(spe.h.rda)$adj.r.squared
# Forward selection using forward.sel() {adespatial}
adespatial::forward.sel(spe.h,databae$envir, adjR2thresh = R2a.all)

spe.h.rda_fw <- rda(spe.h~MAP+PercCurrGrass,data=databae$envir)
anova(spe.h.rda_fw, permutations = how(nperm = 999))
RsquareAdj(spe.h.rda_fw)$adj.r.squared

vif.cca(spe.h.rda)
vif.cca(spe.h.rda_fw)

plot(spe.h.rda_fw,
     scaling = 1,
     display = c("sp", "lc", "cn"),
     main = "Triplot RDA spe.hel - scaling 1 - lc scores")




# dbMEM part 
## Step 1. Construct the matrix of dbMEM variables
site.dbmem <- adespatial::dbmem(databae$coord, silent = FALSE)

# Display and count the eigenvalues
attributes(site.dbmem)$values
length(attributes(site.dbmem)$values)

spe.h.rda_mem <- rda(spe.h,site.dbmem)

## Global test of the RDA result
anova(spe.h.rda_mem, permutations = how(nperm = 999))
R2a.all <- RsquareAdj(spe.h.rda_mem)$adj.r.squared
# Forward selection using forward.sel() {adespatial}
adespatial::forward.sel(spe.h,site.dbmem, adjR2thresh = R2a.all)

#_____________________________________________________      
# VARIANCE PARTITIONING ####
#_____________________________________________________      
var.part.spe.h <- varpart(spe.h, databae$envir, site.dbmem)
plot(var.part.spe.h, digits = 2, bg = c("red", "blue"))

env.sel <- databae$envir %>% select(MAP,PercCurrGrass)

var.part.spe.h.envsel <- varpart(spe.h, env.sel, site.dbmem)
plot(var.part.spe.h.envsel, digits = 2, bg = c("red", "blue"))

anova.cca(rda(spe.h, env.sel), perm.max = 999) ## test pure env
anova.cca(rda(spe.h, site.dbmem), perm.max = 999) ## test pure spatial   
anova.cca(rda(spe.h, env.sel, site.dbmem), perm.max = 999) ## test pure spatial  

# New ways to varpart
library(rdacca.hp)
# Lai, J. et al. 2022. Generalizing hierarchical and variation partitioning in multiple regression and 
# canonical analyses using the rdacca.hp R package. - Methods Ecol. Evol. 13: 782–788.
# New approximation for carrying varparts without limitation of variables and specially accouning for the "pure" values

New.RDA <- rdacca.hp(spe.h, cbind(env.sel,site.dbmem), method="RDA", type = "adjR2",var.part = TRUE)
# WARNING: It basically testing different ways of building the RDA based (like permutations) which may take loong times
New.RDA

New.RDA$Hier.part 
New.RDA$Total_explained_variation

# We can also carry permutations to disentangle the significance of the specific variables on the composition of our data
permu.hp(spe.h, cbind(env.sel,site.dbmem), method="RDA", type = "adjR2")


#_____________________________________________________      
# Diversity/ies ####
#_____________________________________________________      

databae$comm

# Alfa diversity - Local scale _______________________

# Compute alpha diversity indices of the fish communities
N0 <- rowSums(databae$comm > 0) # Species richness
N0 <- specnumber(databae$comm) # Species richness (alternate)
H <- diversity(databae$comm) # Shannon entropy (base e)
Hb2 <- diversity(databae$comm, base = 2) # Shannon entropy (base 2)
N1 <- exp(H) # Shannon diversity (base e)
# (number of abundant species)
N1b2 <- 2^Hb2 # Shannon diversity (base 2)
N2 <- diversity(databae$comm, "inv") # Simpson diversity
# (number of dominant species)
J <- H / log(N0) # Pielou evenness
E10 <- N1 / N0 # Shannon evenness (Hill's ratio)
E20 <- N2 / N0 # Simpson evenness (Hill's ratio)
div <- data.frame(N0, H, Hb2, N1, N1b2, N2, E10, E20, J)

pairs(div)
corrmorant::corrmorant(div)# This one will not work... 

# Alfa - Beta - Gamma diversity 
## You can check further information at 
hillR::hill_taxa_parti(databae$comm, q = 0)
hillR::hill_taxa_parti(databae$comm, q = 1)
hillR::hill_taxa_parti(databae$comm, q = 2)

# Beta diversity - Regional scale _______________________

# Jaccard-based Podani indices (presence-absence data)
Jacc_Beta <- beta.div.comp(databae$comm, coef = "J", quant = FALSE)
Ruz_Beta <- beta.div.comp(databae$comm, coef = "J", quant = TRUE)

Sor_Beta <- beta.div.comp(databae$comm, coef = "S", quant = FALSE)
BC_Beta <- beta.div.comp(databae$comm, coef = "S", quant = TRUE)

Jacc_Beta$repl
Jacc_Beta$rich
Jacc_Beta$D
Jacc_Beta$part

tri.plot_Beta <- cbind((1-Jacc_Beta$D),Jacc_Beta$repl,Jacc_Beta$rich)
colnames(tri.plot_Beta) <- c("Similarity", "Repl", "Rich/Ab Diff")
par(mfrow = c(2, 2))
ade4::triangle.plot(as.data.frame(tri.plot_Beta[, c(3, 1, 2)]),min3 = c(0,0,0),max3 = c(1,1,1),
              show = FALSE,labeltriangle = FALSE,addmean = TRUE)
text(-0.45, 0.5, "RichDiff", cex = 1.5);text(0.4, 0.5, "Repl", cex = 1.5);text(0, -0.6, "Jaccard similarity", cex = 1.5)

tri.plot_Beta <- cbind((1-Ruz_Beta$D),Ruz_Beta$repl,Ruz_Beta$rich)
ade4::triangle.plot(as.data.frame(tri.plot_Beta[, c(3, 1, 2)]),min3 = c(0,0,0),max3 = c(1,1,1),
                    show = FALSE,labeltriangle = FALSE,addmean = TRUE)
text(-0.45, 0.5, "RichDiff", cex = 1.5);text(0.4, 0.5, "Repl", cex = 1.5);text(0, -0.6, "Jaccard similarity", cex = 1.5)

tri.plot_Beta <- cbind((1-Sor_Beta$D),Sor_Beta$repl,Sor_Beta$rich)
ade4::triangle.plot(as.data.frame(tri.plot_Beta[, c(3, 1, 2)]),min3 = c(0,0,0),max3 = c(1,1,1),
                    show = FALSE,labeltriangle = FALSE,addmean = TRUE)
text(-0.45, 0.5, "RichDiff", cex = 1.5);text(0.4, 0.5, "Repl", cex = 1.5);text(0, -0.6, "Jaccard similarity", cex = 1.5)

tri.plot_Beta <- cbind((1-BC_Beta$D),BC_Beta$repl,BC_Beta$rich)
ade4::triangle.plot(as.data.frame(tri.plot_Beta[, c(3, 1, 2)]),min3 = c(0,0,0),max3 = c(1,1,1),
                    show = FALSE,labeltriangle = FALSE,addmean = TRUE)
text(-0.45, 0.5, "RichDiff", cex = 1.5);text(0.4, 0.5, "Repl", cex = 1.5);text(0, -0.6, "Jaccard similarity", cex = 1.5)

dev.off()


# New ways to varpart
library(rdacca.hp)
dbRDA.repl <- rdacca.hp(Jacc_Beta$repl, databae$envir, method="dbRDA", type = "adjR2",var.part = TRUE)
dbRDA.repl

dbRDA.rich <- rdacca.hp(Jacc_Beta$rich, databae$envir, method="dbRDA", type = "adjR2",var.part = TRUE)
dbRDA.rich

dbRDA.D <- rdacca.hp(Jacc_Beta$D, databae$envir, method="dbRDA", type = "adjR2",var.part = TRUE)
dbRDA.D

dbRDA.repl$Hier.part
dbRDA.rich$Hier.part
dbRDA.D$Hier.part

dbRDA.repl$Total_explained_variation
dbRDA.rich$Total_explained_variation
dbRDA.D$Total_explained_variation

# We can also carry permutations to disentangle the significance of the specific variables 
#on the composition of our data
permu.hp(Jacc_Beta$repl, databae$envir, method="dbRDA", type = "adjR2")
permu.hp(Jacc_Beta$rich, databae$envir, method="dbRDA", type = "adjR2",permutations = 20)
permu.hp(Jacc_Beta$D, databae$envir, method="dbRDA", type = "adjR2",permutations = 20)

par(mfrow = c(1, 3))
plot(dbrda(Jacc_Beta$repl ~ ., data = databae$envir , add = "cailliez"),scaling = 1,display = c("sp", "lc", "cn"),main = "Repl")
plot( dbrda(Jacc_Beta$rich ~ ., data = databae$envir , add = "cailliez"),scaling = 1,display = c("sp", "lc", "cn"),main = "Richn diff")
plot( dbrda(Jacc_Beta$D ~ ., data = databae$envir , add = "cailliez"),scaling = 1,display = c("sp", "lc", "cn"),main = "D")
dev.off()


