
library(tidyverse)
library(babynames)
data(babynames)

# FILTER
filter(babynames, name== "Alice", sex=="F")
# SELECT
select(babynames, year)

# ARRANGE

arrange(babynames, prop)

arrange(babynames, desc(prop)) # asc() is the default

arrange(babynames, desc(name), desc(n))

# MUTATE

mutate(babynames, total=n/prop)

mutate(babynames, prop = prop*100)

mutate(babynames, year_diff=2018-year, months_diff=year_diff-12)

mutate(
  babynames,
  dummy=if_else(sex=="COLEO", "Femenino",
                if_else("HETERO","Masculino",NA))
)
mutate(var2 = case_when(str_detect(var1, '^123') ~ 'ok',
                        TRUE ~ 'not ok'))


mutate(
  babynames,
  grup=if_else(Ordre==c("Coleoptera","Heterp`tera",""),"OCH", 
               if_else(Ordre==c()))
  
  # if_else IS BINNARY
)



mutate(
  babynames,
  dummy = case_when(str_detect(name,"Alic") ~ "Alicia",
                    str_detect(name,"Emm") ~ "Paco",
    TRUE ~ "pues otro nombre"
  )
)

# SUMMARISE 

summarise(
  babynames, 
  mean= mean(year),
  sd= sd(year)
)

# GROUP_BY

by_year <- group_by(babynames, year)

summarise(by_year, mean= mean(n))

plot(summarise(by_year, mean= mean(n)))

# Piping!
babynames%>%
  group_by(year)%>%
  summarise(mean= mean(n))

# Pipes connect each call that you do in "dplyr" is like saying 
#     And then --- 
#     And then --- 

# Not using pipes
summarise(
  group_by(
    filter(babynames, year == 1880), sex
  ),
  max = max(n),prop = max(prop))

# Using them 
year_1880 <- babynames %>%
  filter(year == 1880) %>%
  group_by(sex) %>%
  summarise(max = max(n),prop = max(prop))

babynames%>%
  group_by(year, sex)%>%
  arrange(desc(prop))%>%
  summarise(max= max(prop),
            name=first(name)) %>%
  select(name, sex, year, max)

### PLOT IN THE SAME ENVIRONMENT! 
library(ggplot2)

years <- babynames %>%
  filter(year == c(1880,1910,1950,1960)) %>%
  group_by(year, sex) %>%
  summarise(mean = mean(n))

seq(1880, 2017, 10)
rep(1,231)

years <- babynames %>%
  filter(year == seq(1880, 2017, 10))%>%
  group_by(year, sex) %>%
  summarise(mean = mean(n))

ggplot(data=years,
       aes(x = year, y = mean))+
  geom_point()

# The final value!!!!
babynames%>%
  filter(year==seq(1880, 2017, 10))%>%
  group_by(year, sex)%>%
  summarise(mean = mean(n))%>%
  ggplot(
    aes(x = year, y = mean, color=sex))+
  geom_point(size=5)+
  geom_line()+
  scale_color_manual(values=c("pink","blue"))+
  theme_classic()

#####__________________________________________________________________________________
## TIME TO GO WILD! 

MyData <- read.csv2("Folder with the data/MyData.csv",dec = ".")
colnames(MyData) <- c("ID", "Wildfire","Taxa","AbundL","Time")
str(MyData)
head(MyData)

Matriu_Ab_Spp <- MyData%>%
                    group_by(Wildfire,Time,Taxa)%>%
                    summarise(TotAbb=sum(AbundL))

MyData%>%
  group_by(Time, Wildfire,Taxa)%>%
  summarise(TotAbb=sum(AbundL))%>%
  mutate(PA=1)%>%
  group_by(Time,Wildfire)%>%
  summarise(Rich=sum(PA))%>%
  ggplot(aes(x=Time, y=Rich, fill=Wildfire))+
        geom_line(aes(col=Wildfire),size=2)+
        geom_point(size=5, shape=21)+
  scale_x_discrete(breaks=c("Spring","Summer","Autum","Winter"),
                    limits=c("Spring","Summer","Autum","Winter"))+
  theme_classic()

install.packages("vegan")
library(vegan)
install.packages("adespatial")
library(adespatial)

?specnumber

Matriu_SitxSpp <- Matriu_Ab_Spp%>%
                  spread(key=Taxa, value = TotAbb, fill = 0)

#Richness
richne <-specnumber(Matriu_SitxSpp,MARGIN = 1)

#Beta-diversity partitions
taxa.q_Jaccard <- beta.div.comp(Matriu_SitxSpp[3:ncol(Matriu_SitxSpp)],
                                coef = "J", quant = TRUE)
taxa.q_Jaccard$D
taxa.q_Jaccard$repl
taxa.q_Jaccard$rich

#Beta-diversity - Uniqueness
LCBD_values <- beta.div(Matriu_SitxSpp[3:ncol(Matriu_SitxSpp)], method = "hellinger", nperm = 9999)
LCBD_values$LCBD

# Multivariance NMDS

NMDS_values <- metaMDS(Matriu_SitxSpp[3:ncol(Matriu_SitxSpp)], distance = "bray")

plot(NMDS_values,display = "sites")

x<- NMDS_values$points[,1]
y<- NMDS_values$points[,2]

dataset <- data.frame(x,y,
                      Wildfire=Matriu_SitxSpp$Wildfire)
dataset%>%
ggplot(aes(x=x, y=y, fill=Wildfire))+
           geom_point(shape=21, size=6)+
  theme_classic()





