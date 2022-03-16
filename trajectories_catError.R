
##### Generate trajectory and fit model considering error in categorical predictors


library(survival)
library(rstanarm)
library(rstan)
library(tidyverse)
library(mvtnorm)
library(rstanarm)
source('ExtraFunctions.R')

## 1. Create categorical Environment 
pi1=0.8
pi0=1-pi1
pii=c(pi0,pi1)
lpii=log(pii)
source('create_environmentalNLMR_categorical.R')

## 2. Simulate trajectories
deltat=1
N=150
times=seq(1,N,by=deltat) # observed times
start=c(50,50) # initial coordinates
TT=length(times) # number of observed times
nsims= 5 # number of observed trajectories

track = matrix(NA, nrow=length(times), ncol=3)
track[,1] = sort(times)
delta = diff(track[,1])
track[1,2:3] = start
salida=list()

alpha0=0
alpha=2
beta=2
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
bets.p=numeric(nsims)

for (i in 1:nsims)
{
  ss=salida[[i]]
  sl=StepLengths(coordinates(ss))
  bets.p[i]=1/(1/(2*length(sl))*sum(sl^2))
}

ebeta=mean(bets.p)

## 4. Generate control sample 

case=c()
x_=c()
y_=c()
indx=c()
nsim=c()
num.traj=c()
ral=c()
JJ=50

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
    ebeta=beta
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
control.data=control.data %>% mutate(x2=-sl^2/2)
control.data=control.data %>% mutate(lsl=log(sl))
control.data=control.data %>% mutate(sl2=(sl)^2)
control.data=control.data %>% mutate(ral2=-ral^2/2)
control.data$cov=control.data$cov

#control.data=control.data %>% filter(!is.na(cov))

### set error parameters and add error to covaritates
covE=control.data$cov
p00=0.95
p11=0.9
p10=1-p00
p01=1-p11

P=matrix(1,ncol=2,nrow=2)
P[1,1]=p00
P[1,2]=p01
P[2,1]=p10
P[2,2]=p11

quienes0=which(control.data$cov==0)
quienes1=which(control.data$cov==1)

q.cambio0=which(rbinom(length(quienes0),1,prob=1-p00)==1)
covE[quienes0[q.cambio0]]=1

q.cambio1=which(rbinom(length(quienes1),1,prob=1-p11)==1)
covE[quienes1[q.cambio1]]=0

control.data$covE=covE

# calculate Q and qi probas 

qi0=sum(c(p00*pi0,p01*pi1))
qi1=sum(c(p10*pi0,p11*pi1))

qi=c(qi0,qi1)
lqi=log(qi)

Q=matrix(1,ncol=2,nrow=2)

q00=p00*pi0/qi0
q11=p11*pi1/qi1
Q[1,1]=q00
Q[1,2]=1-q11
Q[2,2]=q11
Q[2,1]=1-q00
lQ=log(Q)

# fit stan 
#Q=t(Q)
# counts=matrix(1,ncol=2,nrow=2)
# counts[1,1]=70
# counts[1,2]=30
# counts[2,2]=3
# counts[2,1]=97

# datalist = list(N=nrow(control.data),x=covE,x_cat=covE+1,y=control.data$case,Q=Q,counts=counts)
# stanc("clogitCatError2.stan")
# Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
# options(mc.cores = parallel::detectCores())
# 
# fite <- stan("clogitCatError.stan", data=datalist,chains=1)
# clogit_stane=as.data.frame(fite)
# 
# hist(clogit_stane$`beta[2]`)


