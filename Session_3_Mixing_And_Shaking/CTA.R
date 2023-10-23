library(ecotraj);library(tidyverse);library(vegan);library(adespatial)

?ecotraj

MyData <- read.csv2("Folder with the data/MyData.csv",dec = ".", sep = ";",
                    col.names = c("Site","Wild","Taxa","AbundL","Time"))

str(MyData)
head(MyData)

MyData_spread <- MyData %>%
  mutate(Time = case_when(
    Time == "1" ~ "1.Autum",
    Time == "2" ~ "2.Winter",
    Time == "3" ~ "3.Spring",
    Time == "4" ~ "4.Summer",
    TRUE ~ "pues otro nombre")) %>% 
  group_by(Site,Time, Wild,Taxa)%>%
  summarise(TotAbb=sum(AbundL)) %>% 
  spread(Taxa, TotAbb, fill=0) 
  
# beta diversity
MyData_spread_tobeta <- MyData_spread%>% 
                        ungroup() %>% 
                        select(-c(Time, Site, Wild))
Dist <- sqrt(vegan::vegdist(MyData_spread_tobeta, "bray"))

# Main trajctory plot
trajectoryPCoA(Dist, MyData_spread$Site, MyData_spread$Time, traj.colors = c("black","red"), lwd = 2,
               survey.labels = T)
legend("topleft", col=c("black","red"), 
       legend=c("Trajectory 1", "Trajectory 2"), bty="n", lty=1, lwd = 2)

trajectoryLengths(Dist, MyData_spread$Site, MyData_spread$Time)
trajectoryAngles(Dist, MyData_spread$Site, MyData_spread$Time)
trajectoryAngles(Dist, MyData_spread$Site, MyData_spread$Time, all=TRUE)
trajectoryDirectionality(Dist, MyData_spread$Site, MyData_spread$Time)

trajectoryConvergence(Dist, MyData_spread$Site, MyData_spread$Time, symmetric = TRUE)


Ds = segmentDistances(Dist, MyData_spread$Site, MyData_spread$Time)$Dseg
mMDS = smacof::mds(Ds)
xret = mMDS$conf
plot(xret, xlab="axis 1", ylab = "axis 2", asp=1, pch=21,
     bg=c(rep("black",3), rep("red",3), rep("blue",3)), 
     xlim=c(-1.5,1), ylim=c(-1,1.5))
text(xret, labels=rep(paste0("s",1:3),3), pos=1)
legend("topleft", pt.bg=c("black","red"), pch=21, bty="n", legend=c("Trajectory 1", "Trajectory 2"))


trajectoryLengths(Dist, MyData_spread$Site, MyData_spread$Time,  relativeToInitial = T)
trajectoryAngles(Dist, MyData_spread$Site, MyData_spread$Time,  relativeToInitial = T)







