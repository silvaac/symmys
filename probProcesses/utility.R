# test example... generate all the needed things
GBM<-function(T = 2, sigma = 0.01,mu=0,N=1,step=1){
  
  t = seq(from=0,to=T,by=step/252)
  S = matrix(0,length(t),N)
  for(i in 1:N){
    e = mu+sigma*rnorm(n = length(t))
    e[1] = 0
    S[,i] = cumprod(1+e)
  }
  sDf0 = as.data.frame(S)
  colnames(sDf0)= paste("Run",1:N,sep="")
  sDf = list(period=t,value=sDf0)
  return(sDf)
  
}

BM<-function(T = 2, sigma = 0.1,mu=0,N=1,step=1){
  
  dt = step/252
  t = seq(from=0,to=T,by=dt)
  S = matrix(0,length(t),N)
  for(i in 1:N){
    e = mu*dt+sigma*sqrt(dt)*rnorm(n = length(t))
    e[1] = 0
    S[,i] = cumsum(e)
  }
  sDf0 = as.data.frame(S)
  colnames(sDf0)= paste("Run",1:N,sep="")
  sDf = list(period=t,value=sDf0)
  return(sDf)
  
}

VG<-function(T = 2, sigma = 0.01,mu=0,N=1,step=1,kappa=1){
  if(kappa==0){kappa=0.001}  
  dt = step/252
  t = seq(from=0,to=T,by=dt)
  S = matrix(0,length(t),N)
  for(i in 2:length(t)){
    d_tau = kappa*rgamma(n=N, shape=dt/kappa,scale = 1) 
    e = mu*d_tau+sigma*sqrt(d_tau)*rnorm(n = N)
    S[i,] = e
  }
  S = apply(S,2,cumsum)
  sDf0 = as.data.frame(S)
  colnames(sDf0)= paste("Run",1:N,sep="")
  sDf = list(period=t,value=sDf0)
  return(sDf)
}

mertonJumpDiff <- function(T=2,mu=0,sigma=0.2,arrivalRate=3.45,
                           logJumpMu=0,logJumpSigma=0.2,N=1,step=1)
{
  dt = step/252
  t = seq(from=0,to=T,by=dt)
  S = matrix(0,length(t),N)
  jP = rpois(N,T*arrivalRate)
  for(j in 1:N){
    
  }
  
  for(i in 1:length(t)){
    
  }
  
}

