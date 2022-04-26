library(survival)
library(rstanarm)
library(rstan)
library(tidyverse)
library(mvtnorm)
library(rstanarm)
source('ExtraFunctions.R')


source('create_environmentalNLMR.R')
#tau=.2

## 2. Simulate trajectories
deltat=1
N=250
times=seq(1,N,by=deltat) # observed times
start=c(50,50) # initial coordinates
TT=length(times) # number of observed times
nsims= 25# number of observed trajectories

track = matrix(NA, nrow=length(times), ncol=3)
track[,1] = sort(times)
delta = diff(track[,1])
track[1,2:3] = start
salida=list()

alpha0=0
alpha=2
beta=1
for (jj in 1:nsims)
{
  J=1000 # number of simulations
  for (k in 2:length(times))
  {
    sigma=sqrt(1/beta*delta[k-1])
    samples=rmvnorm(J,mean=c(track[k-1,2],track[k-1,3]),sigma = diag(sigma^2,2))
    phi=numeric(J)
    
    spatial.sample=SpatialPoints(samples)
    covs= spatial.sample%over% en.c
    quienes=which(!is.na(covs))
    
    for (j in quienes)
    {
      s.sample=covs[j,1]
      pol=s.sample*alpha+alpha0
      phi[j]=exp(pol)
    }
    
    wh=sample(seq(1,J),1,prob=phi)
    track[k,2:3]=samples[wh,]
  }
  
  colnames(track) <- c("time", "x", "y")
  track <- as.data.frame(track)
  
  ## as coordinates
  cctrack=track
  coordinates(cctrack) <- ~x + y
  salida[[jj]]=cctrack
}


## 3. Estimate preliminary beta 

bets.p=list()

stanc("rayleigh_dist.stan")
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
options(mc.cores = parallel::detectCores())

for (i in 1:nsims)
{
  ss=salida[[i]]
  sl=StepLengths(coordinates(ss))
  dlist = list(N=length(sl),y=sl)
  stan.beta <- stan("rayleigh_dist.stan", data=dlist, iter=4000,chains=2)
  beta.p=as.data.frame(stan.beta)
  bets.p[[i]]=1/(beta.p$beta)^2
}


#hist(bets.p[[1]])
# bets.p=numeric(nsims)
# 
# for (i in 1:nsims)
# {
#   ss=salida[[i]]
#   sl=StepLengths(coordinates(ss))
#   bets.p[i]=1/(1/(2*length(sl))*sum(sl^2))
# }
# 
# ebeta=mean(bets.p)


### sample control data

case=c()
x_=c()
y_=c()
indx=c()
nsim=c()
num.traj=c()
ral=c()
JJ=80

for (i in 1:nsims)
{
  trc=salida[[i]]
  N=length(trc)
  x.y=coordinates(trc)
  
  for(t in 1:(N-1))
  {
    case=c(case,1)
    x_=c(x_,x.y[t,1])
    y_=c(y_,x.y[t,2])
    
    loc=x.y[t,]
    #ebeta=runif(1,min=0,max=10) # check this 
    #ebeta=beta
    ebeta=sample(bets.p[[i]],size=1)
    
    sigma=sqrt(1/ebeta*delta[t])
    samples=rmvnorm(JJ,mean=loc,sigma = diag(sigma^2,2))
    
    x_=c(x_,samples[,1])
    y_=c(y_,samples[,2])
    ral=c(ral,c(e.dist(p=x.y[t,],v=x.y[t+1,]),apply(samples,1,e.dist,p=x.y[t,])))
    
    indx=c(indx,rep(t,JJ+1))
    case=c(case,rep(0,JJ))  
  }
  
  num.traj=c(num.traj,rep(seq((JJ+1)*(i-1)+1,(JJ+1)*i),N-1))
  nsim=c(nsim,rep(i,(JJ+1)*(N-1)))
  print(i)
}

control.data=data.frame(case,indx,x_,y_,nsim,num.traj,ral)

# cov values
coords.control=SpatialPoints(control.data[,c('x_','y_')])
covs= coords.control%over% en.c
control.data$cov=covs
control.data$cov=control.data$cov[,1]

#step length
control.data=control.data %>% group_by(num.traj)%>% mutate(sl=c(StepLengths(cbind(x_,y_)),NA))
# control.data=control.data %>% mutate(x2=-sl^2/2)
# control.data=control.data %>% mutate(lsl=log(sl))
# control.data=control.data %>% mutate(sl2=(sl)^2)
control.data=control.data %>% mutate(ral2=-ral^2/2)
control.data$cov=control.data$cov

### set error parameters and add error to covaritates
# covE=control.data$cov+rnorm(length(control.data$cov),0,sd=tau)
# 
# control.data$covE=covE
