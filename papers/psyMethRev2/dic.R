library(MCMCpack)
library(msm)

loglike=function(y,mu,s2) sum(dnorm(y,mu,sqrt(s2),log=T))

myNormSampler=function(y){
  mu0=0
  s20=1
  a=.5
  b=.5
  
  mu=1:M
  s2=1:M
  ll=1:M
  ybar=mean(y)
  N=length(y)
  for (m in 2:M){
    v=1/(N/s2[m-1]+1/s20)
    c=N*ybar/s2[m-1]+mu0/s20
    mu[m]=rnorm(1,v*c,sqrt(v))
    SSE=sum((y-mu[m])^2)
    s2[m]=rinvgamma(1,b+N/2,SSE/2+a)
    ll[m]=loglike(y,mu[m],s2[m])
  }
  return(cbind(mu,s2,ll))
}


fitNorm=function(dat){
  chains=myNormSampler(dat)
  pm=apply(chains[2:M,],2,mean)
  L=loglike(dat,pm[1],pm[2])
  P=2*(L-pm[3])
  DIC= -2*(L-P)
  bf=2*mean(chains[2:M,1]>0)
  return(list(chains=chains,DIC=DIC,bf=bf))
}


myTNormSampler=function(y){
  mu0=0
  s20=1
  a=.5
  b=.5
  
  mu=1:M
  s2=1:M
  ll=1:M
  ybar=mean(y)
  N=length(y)
  for (m in 2:M){
    v=1/(N/s2[m-1]+1/s20)
    c=N*ybar/s2[m-1]+mu0/s20
    mu[m]=rtnorm(1,v*c,sqrt(v),0,Inf)
    SSE=sum((y-mu[m])^2)
    s2[m]=rinvgamma(1,b+N/2,SSE/2+a)
    ll[m]=loglike(y,mu[m],s2[m])
  }
  return(cbind(mu,s2,ll))
}


fitTNorm=function(dat){
  chains=myTNormSampler(dat)
  pm=apply(chains[2:M,],2,mean)
  L=loglike(dat,pm[1],pm[2])
  P=2*(L-pm[3])
  DIC= -2*(L-P)
  bf=2*mean(chains[2:M,1]>0)
  return(list(chains=chains,DIC=DIC,bf=bf))
}

