## Set parameters
n= 1000
mu= 170
sd= 7

## Simulate 1000 samples for Y normal distribution
y= rnorm(n, mu, sd)

## Simulate 1000 samples for X multinomial distribution 
x= rmultinom(1000, 1, prob=c(0.04,0.64,0.32)) #p^2, q^2, 2pq

## Four types of models
# Dominant model
i<-1
x.dom=1
for (i in 1:1000){
  if(x[1,i]==1 | x[3,i]==1){
    x.dom[i] <-1
  }
  else {x.dom[i]<-0}
}

# Recessive model
i<-1
x.rec=1
for (i in 1:1000){
  if(x[1,i]==1){
    x.rec[i] <-1
  }  else {x.rec[i]<-0}
}

# Additive model #Has Gradient
i<-1
x.add=1
for (i in 1:1000){
  if(x[1,i]==1){
    x.add[i]<-2
  } else if (x[3,i]==1){
    x.add[i]<- 1}
  else {x.add[i]<-0}
}

# Genotypic model
i<-1
x.gen<-1
for (i in 1:1000){
  if(x[1,i]==1){
    x.gen[i] <-2}
  else if (x[3,i]==1){
    x.gen[i]<-1}
  else {x.gen[i]<-3}
}

x.gen<-as.character(x.gen)

data<- data.frame(y, x.dom, x.rec, x.add, x.gen)

m1<- lm(y~x.dom, data= data)
m2<- lm(y~x.rec, data= data)
m3<- lm(y~x.add, data= data)
m4<- lm(y~x.gen, data= data)

p1<- summary(m1)$coefficient[2,4]
p2<- summary(m2)$coefficient[2,4]
p3<- summary(m3)$coefficient[2,4]
p4<- summary(m4)$coefficient[2,4]
p.dat<-c(p1, p2, p3, p4)
p5<-min(p.dat[1], p.dat[2], p.dat[3], p.dat[4])

y.dist= rnorm(1000, 170, 7)
n.rep=10000
p.dat.rep=matrix(0,nrow=n.rep, ncol=5)
y.dist= matrix(0,nrow=n.rep, ncol=1000)

i.rep <- 1
n.rep <- 10000
p.dat.rep <- matrix(0,nrow=n.rep, ncol=5) 
y.dist <- matrix(0,nrow=n.rep, ncol=1000) 

for (i.rep in 1:n.rep){
  p.dat <- p.dat.rep[i.rep,]
  y.dist[i.rep,] <- rnorm (1000, mean=170, sd=7)
  p.dat[1] <-summary(lm(y.dist[i.rep,]~x.dom))$coefficients[2, 4]
  p.dat[2] <-summary(lm(y.dist[i.rep,]~x.rec))$coefficients[2, 4]
  p.dat[3] <-summary(lm(y.dist[i.rep,]~x.add))$coefficients[2, 4]
  p.dat[4] <-summary(m4 <- lm(y.dist[i.rep,]~x.gen))$coefficients[2, 4]
  p.dat[5] <- min(p.dat[1], p.dat[2], p.dat[3], p.dat[4])
  p.dat.rep[i.rep,]<-p.dat
}

## Plot histogram 
hist(p.dat.rep)
# p value= 0.05
abline(v=0.05,col="blue")
# Bonferroni correction, p= 0.0125
abline(v=0.0125,col="green") 
# 5% data point, p= 0.02
abline(v=quantile(p.dat.rep[,5],0.05),col="red") 

# Pearson correlarion
cor.test(p.dat.rep[,1], p.dat.rep[,2], method= "pearson") 
cor.test(p.dat.rep[,1], p.dat.rep[,3], method= "pearson") 
cor.test(p.dat.rep[,1], p.dat.rep[,4], method= "pearson") 
cor.test(p.dat.rep[,2], p.dat.rep[,3], method= "pearson") 
cor.test(p.dat.rep[,2], p.dat.rep[,4], method= "pearson") 
cor.test(p.dat.rep[,3], p.dat.rep[,4], method= "pearson") 
