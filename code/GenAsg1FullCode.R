#Using EM algorithm

nA = 9123
nB = 2987
nAB = 1269
nO = 7725

# output: nAA, nAO, nBB, nBO, nAB, nOO.

p = 0.3
q = 0.4
n = nA + nB + nAB + nO

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

#Using Newton-Raphson algorithm

nA = 9123
nB = 2987
nAB = 1269
nO = 7725

library(numDeriv)

f<-function(x)
{
  nA*log(x[1]^2+2*x[1]*(1-x[1]-x[2]))+
    nB*log(x[2]^2+2*x[2]*(1-x[1]-x[2]))+
    nAB*log(2*x[1]*x[2])+nO*log((1-x[1]-x[2])^2)
}

x0=c(0.3,0.4)

xk=x0
xkplusone=Inf
dif=abs(xkplusone-xk)
print(dif[1])
dif<-c(1,1)

totalTime = 0
totalIter = 0

while((dif[1]>1e-12 && dif[2]>1e-12))
{
  startTime <- Sys.time()
  
  fprime<-grad(f,xk)
  fpprime<-hessian(f,xk)
  invfpprime<-solve(fpprime)
  xkplusone=xk-invfpprime%*%fprime
  dif=abs(xkplusone-xk)
  xk=xkplusone
  print(xkplusone)
  
  endTime <- Sys.time()
  totalTime = (endTime - startTime) + totalTime
  totalIter = totalIter + 1
  
  cat("dif1 =", dif[1], ", dif2 =", dif[2], "\n")
  Sys.sleep(0.01)
}
end

cat('avg run time:', totalTime/totalIter, '\n')
cat('total iterations:', totalIter, '\n')

#Using Grid searching

library("scatterplot3d")

nA = 9123
nB = 2987
nAB = 1269
nO = 7725

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
