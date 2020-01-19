# Input:
#      mu: mean of Guassian distribution
#      sd: standard deviation of Guassian distribution
#      n: number of samples
#      p: minor allele frequency
#      q: major allele frequency
#      n_rep: number of simulation replicates

# Output: 
#      n_rep sets of p-values for each model

# Set parameters
n= 1000
mu= 170
sd= 7
p= 0.2
q= 0.8
n_rep=10000

# Simulate 1000 samples for Y normal distribution
y= rnorm(n, mu, sd)

# Simulate 1000 samples for X multinomial distribution 
x= rmultinom(n, 1, prob=c(p^2,q^2,2*p*q))

# Four types of genetic models
# Dominant model
# Initialization
x_dom<-logical(n)

for (i in 1:n){
  if(x[1,i]==1 | x[3,i]==1){
    x_dom[i] <-TRUE
  }
  else {
    x_dom[i]<-FALSE
  }
}

# Recessive model
# Initialization
x_rec<-logical(n)

for (i in 1:n){
  if(x[1,i]==1){
    x_rec[i] <-TRUE
  }  
  else {
    x_rec[i]<-FALSE
  }
}

# Additive model # Has Gradient
# Initialization
x_add<- numeric(n)

for (i in 1:n){
  if(x[1,i]==1){
    x_add[i]<-2
    }
  else if (x[3,i]==1){
    x_add[i]<- 1
    }
  else{
    x_add[i]<-0
    }
}

# Genotypic model
# Initialization
x_gen<- character(n)

for (i in 1:n){
  if(x[1,i]==1){
    x_gen[i] <-"2"
    }
  else if (x[3,i]==1){
    x_gen[i]<-"1"
    }
  else{
    x_gen[i]<-"3"
    }
}

mydata<- data.frame(y, x_dom, x_rec, x_add, x_gen)

m1<- lm(y~x_dom, data= mydata)
m2<- lm(y~x_rec, data= mydata)
m3<- lm(y~x_add, data= mydata)
m4<- lm(y~x_gen, data= mydata)


p1<- summary(m1)$coefficient[2,4]
p2<- summary(m2)$coefficient[2,4]
p3<- summary(m3)$coefficient[2,4]
p4<- summary(m4)$coefficient[2,4]
p_dat<-c(p1, p2, p3, p4)
p5<-min(p_dat[1], p_dat[2], p_dat[3], p_dat[4])

y_dist= rnorm(n, mu, sd)

p_dat_rep= matrix(0, nrow=n_rep, ncol=5)
y_dist= matrix(0, nrow=n_rep, ncol=n)

p_dat_rep <- matrix(0, nrow=n_rep, ncol=5) 
y_dist <- matrix(0, nrow=n_rep, ncol=n) 

for (i_rep in 1:n_rep){
  p_dat <- p_dat_rep[i_rep,]
  y_dist[i_rep,] <- rnorm (n, mean= mu, sd= sd)
  p_dat[1] <-summary(lm(y_dist[i_rep,]~x_dom))$coefficients[2, 4]
  p_dat[2] <-summary(lm(y_dist[i_rep,]~x_rec))$coefficients[2, 4]
  p_dat[3] <-summary(lm(y_dist[i_rep,]~x_add))$coefficients[2, 4]
  p_dat[4] <-summary(m4 <- lm(y_dist[i_rep,]~x_gen))$coefficients[2, 4]
  p_dat[5] <- min(p_dat[1], p_dat[2], p_dat[3], p_dat[4])
  p_dat_rep[i_rep,]<-p_dat
}

## Plot histogram 
hist(p_dat_rep)
# p value= 0.05
abline(v= 0.05, col= "blue")
# Bonferroni correction, p= 0.0125
abline(v= 0.0125, col= "green") 
# 5% data point, p= 0.02
abline(v=quantile(p_dat_rep[,5], 0.05), col= "red") 

# Pearson correlarion
cor.test(p_dat_rep[,1], p_dat_rep[,2], method= "pearson") 
cor.test(p_dat_rep[,1], p_dat_rep[,3], method= "pearson") 
cor.test(p_dat_rep[,1], p_dat_rep[,4], method= "pearson") 
cor.test(p_dat_rep[,2], p_dat_rep[,3], method= "pearson") 
cor.test(p_dat_rep[,2], p_dat_rep[,4], method= "pearson") 
cor.test(p_dat_rep[,3], p_dat_rep[,4], method= "pearson") 
