#### Prueba categorica con metropolis hasting Alg

library(matrixStats)
library(tidyverse)
library(extraDistr)
library(dplyr)


# Simulation of the trajectories 
source('trajectories_catError.R')
source('MCMC_alg.R')
dat <- control.data[order(control.data$indx), ] # order by strata
dd=dat%>% filter(nsim%in%c(3))

## 1. MLE considering and not considering the error
fit=clogit(case~cov+ral2+strata(indx),data=dd)
fit
fitE=clogit(case~covE+ral2+strata(indx),data=dd)
fitE

## 2 Bayesian inference with real data

fitB1=MCMC_Nerror(dd,iters=4000,chains=1,
                  priorA=0,priorSDA=3,priorB=0,priorSDB=1,alph=1,b=0,Xcov=dd$cov)

plot(fitB1$out.alpha,type='l')
abline(h=alpha,col='red')
plot(density(rnorm(100000,0,3)),ylim=c(0,1))
lines(density(fitB1$out.alpha[1500:3000]),col='red')
abline(v=alpha,col='blue')

## 2 Bayesian inference with observed data

fitB2=MCMC_Nerror(dd,iters=4000,chains=1,
                  priorA=0,priorSDA=3,priorB=0,priorSDB=1,alph=1,b=0,Xcov=dd$covE)

plot(fitB2$out.alpha,type='l')
abline(h=alpha,col='red')
plot(density(rnorm(100000,0,3)),ylim=c(0,1))
lines(density(fitB2$out.alpha[1500:3000]),col='red')
abline(v=alpha,col='blue')

## 3 Bayesian inference with observed data

fitB=MCMC_categoricalE(dd=dd,iters=1000,chains=1,
                       priorA=0,priorSDA=3,priorB=0,priorSDB=1,
                      Q=Q,qi,alph=2,b=0,Xcov=dd$cov,JJ=JJ)

fitB$ncambios
# PLot results 
plot(fitB$out.alpha,type='l')
abline(h=alpha,col='red')
plot(density(rnorm(100000,0,3)),ylim=c(0,1))
lines(density(fitB$out.alpha[500:1000]),col='red')
abline(v=alpha,col='blue')

plot(fitB$out.beta,type='l')
hist(fitB$out.beta[1500:3000])

# Not considering the error 

# pruebo rcpp code




