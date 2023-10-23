
# Modelos nulos para practicar loops i funciones ÈPICAS!

# my first function with David
saluda <- function(x="hola"){print(x)}
saluda(x="Callate, no quiero saludarte")


# Un model simple i senzill, dels de casa----
abundance <- c(1,2,3,1,1,1,6,7,8,9,7,8,9,9,15)
body.size <- c(9,6,4,3,7,9,0.5,4,12,10,2,9,8,7,20)
plot(abundance~body.size)
cor(abundance,body.size)
cor.test(abundance,body.size)

# We create a null model with the sample function that randomizes the values from a given vector
null <- sample(abundance,
               size =  length(abundance),
               replace = F)
a <- cor(null, body.size)

# We can built a function wich reproduces these as many times as wanted
pepino<- function(x,y, permutation=1000){
  out <- c()  # IMPORTANT !!! We create a vector where the output will be loaded 
  for(i in 1:permutation){
    null <- sample(x,
                   size =  length(x),
                   replace = F)
    out[i] <- cor(null, y)
  }
  hist(out)
  abline(v=cor(x,y), col="red")
  a<-1-sum(out<cor(x,y))/length(out)
  e <- ifelse(a>0.05,"LA CAGASTE","****")
  u <- c(a,e)
  u
} 

pepino(x = abundance, y = body.size, permutation = 1000)

#More examples ----
abundance

evennes <- function(x){
  s <- length(x)
  abTOT<-sum(x)
  pi<- x/abTOT
  shannon <- -sum((pi)*log(pi, base=2))
  shannon/log(s, base = 2)
}

evennes(abundance)

#tip: We need to sample 87 individuals and assigng them to 6 species randomly
sum(abundance)

out <- matrix()
for(i in 1:1000){
  a <- sample(1:15, 
              size=sum(abundance),
              replace = T)
  out[i] <- evennes(table(a))
  out
}
hist(out,xlim = c(0.85,1))
abline(v=evennes(abundance), col="red")

for(i in 1:1000){
  antonio <- 1-body.size/sum(body.size)
  a <- sample(1:15, 
              size=sum(abundance),
              replace = T,
              prob =antonio)
  out[i] <- evennes(table(a))
}
hist(out,xlim = c(0.85,1))
abline(v=evennes(abundance), col="red")

