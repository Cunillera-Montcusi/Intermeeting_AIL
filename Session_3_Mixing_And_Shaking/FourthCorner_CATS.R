
library(mvabund);library(lattice)

data(antTraits)

# Single site approximation and trait regression 
y=c(t(antTraits$abund[1,])) #coercing the first row of antTraits$abund into a vector
singleSite = data.frame(y, antTraits$traits,check.names=F) #make a data frame from site abundances and traits
head(singleSite) #here is what the dataset looks like

ft=manyglm(y~.,data=singleSite, family="negative.binomial") #fit a negative binomial regression
summary(ft,nBoot=1) #find which traits are important
plot(ft) #fit looks OK

y=c(t(antTraits$abund[2,])) #coercing the first row of antTraits$abund into a vector
anotherSite = data.frame(y, antTraits$traits,check.names=F) #make a data frame from site abundances and traits
ft=manyglm(y~.,data=anotherSite, family="negative.binomial") #fit a negative binomial regression
summary(ft,nBoot=999) #find which traits are important

### ________________________________________
# Fourth corner and interaction env-trait
### ________________________________________
ft=traitglm(antTraits$abund,antTraits$env,antTraits$traits)
ft$fourth #print fourth corner terms

abSum = apply(antTraits$abund>0,2,mean) # stores the proportion of values that are non-zero for each species.
ab = antTraits$abund[,abSum>0.2] # just keep columns of $abund that are more than 20% presence
tr = antTraits$traits[abSum>0.2,] # and make sure that change is also made to rows of the traits data frame.

# now fit the fourth corner model, only as a function of a couple of traits and env variables:
ftSmall=traitglm(ab,antTraits$env[,1:3],data.frame(tr$Weber,tr$Femur))
anova(ftSmall)

# LASSO penalty incorporated
ft1=traitglm(antTraits$abund,antTraits$env,antTraits$traits,method="glm1path")
ft1$fourth #notice LASSO penalty has shrunk many interactions to zero

library(lattice)
a        = max( abs(ft1$fourth.corner) )
colort   = colorRampPalette(c("blue","white","red")) 
plot.spp = levelplot(t(as.matrix(ft1$fourth.corner)), xlab="Environmental Variables",
                     ylab="Species traits", col.regions=colort(100), at=seq(-a, a, length=100),
                     scales = list( x= list(rot = 45)))
print(plot.spp)


ftspp1=traitglm(antTraits$abund,antTraits$env,method="glm1path")

## No traits matrix entered, so will fit SDMs with different env response for each spp
ftspp1$fourth #notice LASSO penalty has shrunk many interactions to zero
a        = max( abs(ft1$fourth.corner) )
colort   = colorRampPalette(c("blue","white","red")) 
plot.4th = levelplot(t(as.matrix(ft1$fourth.corner)), xlab="Environmental Variables",
                     ylab="Species traits", col.regions=colort(100), at=seq(-a, a, length=100),
                     scales = list( x= list(rot = 45)))
print(plot.4th)


load("CESTES_Null.RData")
databae <- LSmatred[[1]]

databae$comm
databae$traits
databae$envir

# LASSO penalty incorporated
ft1=traitglm(databae$comm,databae$envir,databae$traits,method="glm1path")
ft1$fourth #notice LASSO penalty has shrunk many interactions to zero

library(lattice)
a        = max( abs(ft1$fourth.corner) )
colort   = colorRampPalette(c("blue","white","red")) 
plot.spp = levelplot(t(as.matrix(ft1$fourth.corner)), xlab="Environmental Variables",
                     ylab="Species traits", col.regions=colort(100), at=seq(-a, a, length=100),
                     scales = list( x= list(rot = 45)))
print(plot.spp)

