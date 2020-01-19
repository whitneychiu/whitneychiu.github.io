# Input: 
#     nA: number of people with A genotype
#     nB: number of people with B genotype
#     nAB: number of people with AB genotype
#     nO: number of people with O genotype

# Output:
#     p: allele frequency of A
#     q: allele frequency of B

# logL~= nA*(p^2+2*p*(1-p-q)) + nB*(q^2+2*q*(1-p-q)) + nAB*(2*p*q) + nO*(1-p-q)^2


# Input parameters
nA = 9123
nB = 2987
nAB = 1269
nO = 7725
n = nA + nB + nAB + nO

### Using EM algorithm ###
# Initialization
p = 0.3
q = 0.4
pOld = Inf
qOld = Inf

totalTime = 0
totalIter = 0
while((abs(p - pOld) > 1e-12) && (abs(q - qOld) > 1e-12))
{ 
  startTime <- Sys.time()
  
  # save old results
  pOld = p
  qOld = q
  
  # E step
  nAA = nA * (p*p)         / (p*p + 2*p*(1-p-q))
  nAO = nA * (2*p*(1-p-q)) / (p*p + 2*p*(1-p-q))
  nBB = nB * (q*q)         / (q*q + 2*q*(1-p-q))
  nBO = nB * (2*q*(1-p-q)) / (q*q + 2*q*(1-p-q))
  
  # M step
  p = (2*nAA + nAO + nAB) / (2*n)
  q = (2*nBB + nBO + nAB) / (2*n)    
  
  endTime <- Sys.time()
  totalTime = (endTime - startTime) + totalTime
  totalIter = totalIter + 1
  
  cat("p=",p,", ")
  cat("q=",q,"\n")
  Sys.sleep(0.01)
}
cat('avg run time:', totalTime/totalIter, '\n')
cat('total iterations:', totalIter, '\n')

### Using Newton-Raphson algorithm ###
library(numDeriv)

f<-function(x)
{
  nA*log(x[1]^2+2*x[1]*(1-x[1]-x[2]))+
    nB*log(x[2]^2+2*x[2]*(1-x[1]-x[2]))+
    nAB*log(2*x[1]*x[2])+nO*log((1-x[1]-x[2])^2)
}

x0=c(0.3,0.4)

xk=x0
xk_next=Inf
dif=abs(xk_next-xk)
print(dif[1])
dif<-c(1,1)

totalTime = 0
totalIter = 0

while((dif[1]>1e-12 && dif[2]>1e-12))
{
  startTime <- Sys.time()
  
  df<-grad(f,xk)
  ddf<-hessian(f,xk)
  inv_ddf<-solve(ddf)
  xk_next=xk-inv_ddf%*%df
  dif=abs(xk_next-xk)
  xk=xk_next
  print(xk_next)
  
  endTime <- Sys.time()
  totalTime = (endTime - startTime) + totalTime
  totalIter = totalIter + 1
  
  cat("dif1 =", dif[1], ", dif2 =", dif[2], "\n")
  Sys.sleep(0.01)
}
end

cat('avg run time:', totalTime/totalIter, '\n')
cat('total iterations:', totalIter, '\n')


### Using Grid searching ###
library("scatterplot3d")

# initialize loglikelihood function
objFunc<-function(p,q)
{
  nA*log(p^2+2*p*(1-p-q))+
    nB*log(q^2+2*q*(1-p-q))+
    nAB*log(2*p*q)+nO*log((1-p-q)^2)
}

# user picked grid size
gridSize = 0.01

# save p, q, val history
pList = c()
qList = c()
vList = c()

p = 0.0 # initial p
while( p <= 1.0 ){
  
  q = 0.0 # initial q
  while( q <= 1.0 && (p+q) <= 1.0){
    val = objFunc(p,q) # compute val 
    
    # save history of p,q,val
    # cat("p=",p,", q=",q, ", v=",val,"\n")
    pList = c(pList, p)
    qList = c(qList, q)
    vList = c(vList, val)
    
    # update q
    q = q + gridSize
  }
  
  # update p
  p = p + gridSize  
}

pList = pList[is.finite(vList)]
qList = qList[is.finite(vList)]
vList = vList[is.finite(vList)]

# print(pList)
# print(qList)
# print(vList)

scatterplot3d(pList,qList,vList)
