#### simulation and test error measurement in categorical data

library(survival)
library(rstanarm)
library(rstan)
library(tidyverse)
library(mvtnorm)
library(rstanarm)

########## ajuste de datos de trajecctoria

# Simulation of the trajectories 
source('trajectories_catError.R')


dat <- control.data[order(control.data$indx), ] # order by strata
dd=dat%>% filter(nsim%in%c(1))
dat <- dd[order(dd$indx), ] # order by strata


## 1. MLE considering and not considering the error
fit=clogit(case~cov+ral2+strata(indx),data=dd)
fit
fitE=clogit(case~covE+ral2+strata(indx),data=dd)
fitE




## Ajuste sin error y feliz

n_cases=dat %>% group_by(indx) %>% summarise(n.cases=sum(case==1))
datalist1 = list(N=nrow(dat),n_grp=N-1,n_coef=2,
                grp=dd$indx, y=dd$case,x=dd[,c('ral2','cov')], 
                n_group=rep(JJ+1,N-1),n_case=rep(1,N-1)) 
                
clogit_stanFELIZ <- stan("clogit2.stan", data=datalist1, iter=4000,chains=2)

outs1=as.data.frame(clogit_stanFELIZ)
hist(outs1$`b[1]`)
hist(outs1$`b[2]`)
abline(v=beta,col='red')

## Ajuste sin error y triste (varaibales con error)

n_cases=dat %>% group_by(indx) %>% summarise(n.cases=sum(case==1))

datalist2 = list(N=nrow(dat),n_grp=N-1,n_coef=2,
                grp=dd$indx, y=dd$case,x=dd[,c('ral2','covE')], 
                n_group=rep(JJ+1,N-1),n_case=rep(1,N-1)) 

clogit_stanTRISTE<- stan("clogit2.stan", data=datalist2, iter=4000,chains=2)

outs2=as.data.frame(clogit_stanTRISTE)
hist(outs2$`b[1]`)
hist(outs2$`b[2]`)
abline(v=beta,col='red')


## Ajuste con error, corregido y feliz (?)
n_cases=dat %>% group_by(indx) %>% summarise(n.cases=sum(case==1))

datalist3 = list(N=nrow(dat),x=dd$covE,x_cat=dd$covE+1,x2=dd$ral2,
                n_grp=N-1,y= dd$case,qgr=dd$indx,lQ=lQ)

stanc("clogitCatError.stan")
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
options(mc.cores = parallel::detectCores())

clogit_stan3 <- stan("clogitCatError.stan", data=datalist3, iter=4000,chains=2)
outs3=as.data.frame(clogit_stan3)
hist(outs3$`beta[1]`)
hist(outs3$`beta[2]`)
abline(v=beta,col='red')


### Save data

Betac=c(outs1$`b[1]`,outs2$`b[1]`,outs3$`beta[1]`)
Alpha=c(outs1$`b[2]`,outs2$`b[2]`,outs3$`beta[2]`)
Model=c(rep('Obs sin error',length(outs1$`b[1]`)),rep('Obs con error',length(outs1$`b[1]`)),
        rep('Obs con error y corregido',length(outs1$`b[1]`)))
Model <- factor(Model, levels = c("Obs sin error", "Obs con error", "Obs con error y corregido"))
realAlpha=rep(alpha,length(Model))

results=data.frame(Betac,Alpha,Model,realAlpha)

filename=paste(paste('results_q00',round(diag(Q)[1],digits=2),sep='='),
      paste('q11',round(diag(Q)[2],digits=2),sep='='),sep='_')
#write_csv(results,file=paste(filename,'.csv',sep=''))

#### create Plots ggplot

Palpha=ggplot(results)+geom_histogram(aes(Alpha))+facet_grid(~Model)+theme_bw()+ 
  geom_vline(xintercept = unique(realAlpha),color= 'red')+
  labs(y='',x=expression(paste('Distribución posterior para ',alpha)))

#ggsave(filename = paste('Alpha',filename,'.pdf'),plot=Palpha)


################################################################################
################################# PRUEBAS ANTIGUAS ##########################
################################################################################
ngroups=100 # number of groups
alphas=runif(ngroups)
n_pgroup=10 # number of observations per group
which.group=rep(seq(1,ngroups),each=n_pgroup)

indx=rep(seq(1,ngroups),each=n_pgroup) #strata
reps=rep(seq(1,n_pgroup),ngroups)
x1=rnorm(n_pgroup*ngroups)
set.seed(988)
x2=rbinom(n_pgroup*ngroups,size=1,prob=.2)
x=cbind(x1,x2)

beta=5
case=numeric(length(x2))
for(i in 1:length(x2))
{
  ps=exp(x2[i]*beta+alphas[which.group[i]])/(1+exp(alphas[which.group[i]]+x2[i]*beta))
  case[i]=sample(c(1,0),1,prob=c(ps,1-ps))
}

df=data.frame(case,x1,x2,indx)
fit1=clogit(case~x2+strata(indx),data=df)


#### common logit 
pi1=0.8
pi0=1-pi1
x2=rbinom(10000,size=1,prob=pi1)
x2e=x2
p00=0.9
p11=0.9
p10=1-p00
p01=1-p11

P=matrix(1,ncol=2,nrow=2)
P[1,1]=p00
P[1,2]=p01
P[2,1]=p10
P[2,2]=p11

quienes0=which(x2==0)
quienes1=which(x2==1)

q.cambio0=which(rbinom(length(quienes0),1,prob=1-p00)==1)
x2e[quienes0[q.cambio0]]=1

q.cambio1=which(rbinom(length(quienes1),1,prob=1-p11)==1)
x2e[quienes1[q.cambio1]]=0

beta0=-1
beta=2
#ps=1/(1+exp(-x2*beta))
ps=exp(beta0+x2*beta)/(1+exp(beta0+x2*beta))
case=numeric(length(ps))

case <- rbinom(n = length(ps), size = 1, prob = ps)

df=data.frame(case,x2,x2e)


fitML=glm(case~x2,data=df,family = binomial)
fitML

fiteML=glm(case~x2e,data=df,family = binomial)
fiteML


# bayesian stan 
Q=matrix(1,ncol=2,nrow=2)

q00=p00*pi0/sum(c(p00*pi0,p01*pi1))
q11=p11*pi1/sum(c(p10*pi0,p11*pi1))

Q[1,1]=q00
Q[1,2]=1-q11
Q[2,2]=q11
Q[2,1]=1-q00
Q=t(Q)
counts=matrix(1,ncol=2,nrow=2)
counts[1,1]=70
counts[1,2]=30
counts[2,2]=97
counts[2,1]=3


dataliste = list(N=nrow(df),x=x2e,x_cat=x2e+1,y=case,Q=Q,counts=counts)
datalist = list(N=nrow(df),x=x2e,y=case)
datalist.base = list(N=nrow(df),x=x2,y=case)

# stan options 
stanc("clogitCatError.stan")
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
options(mc.cores = parallel::detectCores())

# fit the base model: logit model without error in the predictors
fit.base <- stan("logitstan.stan", data=datalist.base,chains=2)
clogit_stan.base=as.data.frame(fit.base)
hist(clogit_stan.base$`beta[1]`)
hist(clogit_stan.base$`beta[2]`)
abline(v=beta,col='red')


fite <- stan("clogitCatError.stan", data=dataliste,chains=2)
clogit_stane=as.data.frame(fite)
# hist(clogit_stane$`beta[1]`)
# hist(clogit_stane$`beta[2]`)

fitNe <- stan("logitstan.stan", data=datalist,chains=2)
clogit_stanNe=as.data.frame(fitNe)
# hist(clogit_stanNe$`beta[2]`)
# hist(clogit_stanNe$`beta[1]`)

values=c(clogit_stane$`beta[2]`,clogit_stane$`beta[1]`,clogit_stanNe$`beta[2]`,clogit_stanNe$`beta[1]`)
real.val=rep(c(rep(beta,length(clogit_stane$`beta[1]`)),rep(beta0,length(clogit_stane$`beta[1]`))),2)
param=rep(c(rep('beta',length(clogit_stane$`beta[1]`)),rep('beta0',length(clogit_stane$`beta[1]`))),2)
model=c(rep('W_error',length(clogit_stane$`beta[1]`)*2),rep('N_error',length(clogit_stane$`beta[1]`)*2))
dfE=data.frame(values,param,real.val,model)

vline.data <- dfE %>%group_by(param) %>% summarize(realval = unique(real.val))

ggplot(dfE)+geom_histogram(aes(x=values,fill=model),alpha=.8)+facet_wrap(~param,scales= 'free')+
  geom_vline(aes(xintercept = realval),vline.data,color='red')+theme_bw()

## Veamos la matriz de clasificacion 
hist(clogit_stane$"Q[1,1]")
abline(v=q00,col='red')
hist(clogit_stane$"Q[2,2]")
abline(v=q11,col='red')

## Conditional Logistic regression 

ngroups=100 # number og groups
alphas=runif(ngroups)
n_pgroup=30 # number of observations per group
which.group=rep(seq(1,ngroups),each=n_pgroup)

indx=rep(seq(1,ngroups),each=n_pgroup) #strata
reps=rep(seq(1,n_pgroup),ngroups)
x1=rnorm(n_pgroup*ngroups)
set.seed(988)
x2=rbinom(n_pgroup*ngroups,size=1,prob=.2)
x=cbind(x1,x2)

beta=2

case=numeric(length(x2))
for(i in 1:length(x2))
{
  ps=exp(x2[i]*beta+alphas[which.group[i]])/(1+exp(alphas[which.group[i]]+x2[i]*beta))
  case[i]=sample(c(1,0),1,prob=c(ps,1-ps))
}

df=data.frame(case,x1,x2,indx)
fit1=clogit(case~x2+strata(indx),data=df)

### Fit data STAN 

## Ajuste sin error 
dat <- df[order(df$indx), ] # order by strata
n_cases=dat %>% group_by(indx) %>% summarise(n.cases=sum(case==1))
datalist = list(N=nrow(dat),n_grp=ngroups,n_coef=2,x=df[,c("x1",'x2')],
                y=df[,'case'],grp=df[,'indx'])

#n_group=rep(n_pgroup,ngroups),n_case=n_cases$n.cases

stanc("clogitCatError2.stan")
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
options(mc.cores = parallel::detectCores())

clogit_stan1 <- stan("clogitCatError2.stan", data=datalist, seed=1001,chains=2)
clogit_stan1=as.data.frame(clogit_stan1)

hist(clogit_stan1$`b[1]`)
hist(clogit_stan1$`b[2]`)
abline(v=beta,col='red')



## Adding error 
pi1=0.8
pi0=1-pi1

x2e=x2
p00=0.9
p11=0.9
p10=1-p00
p01=1-p11

P=matrix(1,ncol=2,nrow=2)
P[1,1]=p00
P[1,2]=p01
P[2,1]=p10
P[2,2]=p11

quienes0=which(x2==0)
quienes1=which(x2==1)

q.cambio0=which(rbinom(length(quienes0),1,prob=1-p00)==1)
x2e[quienes0[q.cambio0]]=1

q.cambio1=which(rbinom(length(quienes1),1,prob=1-p11)==1)
x2e[quienes1[q.cambio1]]=0




