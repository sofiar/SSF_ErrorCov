#### simulation and test error measurement in continous  data for nsims

library(rstan)
library(tidyverse)
library(mvtnorm)


# Simulation of the trajectories 
load("~/proyecto_doctoral/Datos_GPS2019/PointProcess/contTray.RData")

# Set error
tau=0.065
covE=control.data$cov+rnorm(length(control.data$cov),0,sd=sqrt(tau))
control.data$covE=covE

stanc("clogitContError.stan")
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
options(mc.cores = parallel::detectCores())

#Save preliminary beta posterior
# for (i in 1:nsims)
# {
#   filename=paste('BETA_results_Contsim',as.character(i),sep='_')
#   write_csv(as.data.frame(bets.p[[i]]),file=paste(filename,'.csv',sep=''))
# }


rm(bets.p,beta.p,dlist,stan.beta)
rm(coords.control,covs,en.c,cctrack,salida,samples,spatial.sample,ss,track,trc,x,x.y)
rm(case,covE,indx,lpii,lqi,lx,nsim,num.traj,ral,s.sample,sl,x_,y_)
rm(delta,phi)
gc()

for (i in 1:nsims)
{
  
  dat <- control.data[order(control.data$indx), ] # order by strata
  dd=dat%>% filter(nsim%in%c(i))
  dat <- dd[order(dd$indx), ] # order by strata
 
  
  # ###### Ajuste sin error y feliz  #########

  # n_cases=dat %>% group_by(indx) %>% summarise(n.cases=sum(case==1))
  # datalist1 = list(N=nrow(dat),n_grp=N-1,n_coef=2,
  #                  grp=dd$indx, y=dd$case,x=dd[,c('ral2','cov')],
  #                  n_group=rep(JJ+1,N-1),n_case=rep(1,N-1))
  # 
  # clogit_stanFELIZ <- stan("clogit2.stan", data=datalist1, iter=4000,chains=2)
  # 
  # outs1=as.data.frame(clogit_stanFELIZ)
  # rm(clogit_stanFELIZ)
  # 
  # # save data
  # Betac=c(outs1$`b[1]`)
  # Alpha=c(outs1$`b[2]`)
  # Model=c(rep('SE-M1',length(outs1$`b[1]`)))
  # realAlpha=rep(alpha,length(Model))
  # 
  # results=data.frame(Betac,Alpha,Model,realAlpha)
  # 
  # filename=paste('results_sim_Cont_NewPrior',as.character(i),sep='_')
  # write_csv(results,file=paste(filename,'.csv',sep=''))
  # rm(outs1)
  # gc()

  ## Ajuste sin error y triste (varaibales con error)
  
  datalist2 = list(N=nrow(dat),n_grp=N-1,n_coef=2,
                   grp=dd$indx, y=dd$case,x=dd[,c('ral2','covE')],
                   n_group=rep(JJ+1,N-1),n_case=rep(1,N-1))

  clogit_stanTRISTE<- stan("clogit2.stan", data=datalist2, iter=3000,chains=2) ### iter=4mil!

  outs2=as.data.frame(clogit_stanTRISTE)
  rm(clogit_stanTRISTE)
  gc()

  # ## Ajuste con error, corregido y feliz (?)


  datalist3 = list(N=nrow(dd),n_grp=N-1,qgr=dd$indx,n_coef=2,
                   x=cbind(dd$covE,dd$ral2),
                   y=dd$case, tau=tau,sigma_x=1)

  # stan options
  clogit_stan3 <- stan("clogitContError.stan", data=datalist3,chains=2,iter=3000)### iter=4mil!
  outs3=as.data.frame(clogit_stan3)
  rm(clogit_stan3)
  gc()

  ### Save data

  Betac=c(outs2$`b[1]`,outs3$`b[1]`)
  Alpha=c(outs2$`b[2]`,outs3$`b[2]`)
  Model=c(rep('CE-M1',length(outs2$`b[1]`)),
          rep("CE-M2",length(outs3$`b[1]`)))
  Model <- factor(Model, levels = c("CE-M1", "CE-M2"))
  realAlpha=rep(alpha,length(Model))

  results=data.frame(Betac,Alpha,Model,realAlpha)

  filename=paste(paste('results_tau_NewPrior_',tau,sep='='),"sim",as.character(i),sep='_')
  write_csv(results,file=paste(filename,'.csv',sep=''))

  ## rm things
  rm(results,filename,Model,realAlpha,Betac,Alpha,dat,dd,outs2,outs3)
  gc()
  print(i)
  
}
