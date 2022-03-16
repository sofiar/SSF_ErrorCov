############################################################################
############################### MCMC algorithms ############################
############################################################################

# LLfunctions

#as clr                
distY=function(alpha,beta,ng,data)
{
  res=0
  for (g in 1:ng){
    data.g=data %>% dplyr::filter(indx==g)
    res=res+alpha*data.g$Xcov[1]+beta*data.g$ral2[1]-logSumExp(alpha*data.g$Xcov+beta*data.g$ral2)
  }
  return(res)
}

#as categorical
distY2=function(alpha,beta,ng,data)
{
  res=0
  for (g in 1:ng){
    data.g=data %>% dplyr::filter(indx==g)
    probs=exp(alpha*data.g$Xcov+beta*data.g$ral2)/sum(exp(alpha*data.g$Xcov+beta*data.g$ral2))
    res=res+dcat(1,probs,log=TRUE)
  }
  return(res)
}
 

#as categorical without data frame
distY3=function(alpha,beta,ng,xcov,ral,quienes.g)
{
  res=0
  for (g in 1:ng){
    X.g=xcov[quienes.g[[g]]]
    ral.g=ral[quienes.g[[g]]]
    probs=exp(alpha*X.g+beta*ral.g)/sum(exp(alpha*X.g+beta*ral.g))
    res=res+dcat(1,probs,log=TRUE)
  }
  return(res)
} 

# Rcpp functions 
Rcpp::sourceCpp('Lrcpp.cpp')


MCMC_categoricalE=function(dd,iters=4000,burn=1000,chains=2,
                           priorA,priorSDA,priorB,priorSDB,Q,qi,
                           alph,b,Xcov,JJ)
{
  d=dd
  N=length(unique(d$indx))+1
  keep.alpha= rep(0,iters)
  keep.beta= rep(0,iters)
  x2=d$ral2
  Xobs=d$covE

  Quienesg=matrix(0,JJ+1,N-1)
  for (i in 1:N-1)
  {Quienesg[,i]=which(dd$indx==i)}
  
  curll=LL(alph,b, N-1 ,Xcov,x2,Quienesg-1,Npg=JJ+1,LOG=TRUE)
  
  ncambios=0
  
  cc=0
  for(j in 1:length(Xcov))
  {cc=cc+lQ[Xcov[j]+1,d$covE[j]+1]}
  
  curxxobs=cc  
  
  
  #MCMC from here. GO!
  for(i in 1:iters){
    
    # Sample alpha
    prior_curll=dnorm(alph,mean=priorA,sd=priorSDA,log=TRUE)
    canalph = rnorm(1,alph,0.6)
    
    canll=LL(canalph,b, N-1 ,Xcov,x2,Quienesg-1,Npg=JJ+1,LOG=TRUE)
    prior_canll=dnorm(canalph,mean=priorA,sd=priorSDA,log=TRUE)
    
    MH <- canll-curll+prior_canll-prior_curll
    
    if(log(runif(1))<MH){
      alph  = canalph
      curll  = canll
    }
    
    # Sample beta
    prior_curll=dnorm(b,mean=priorB,sd=priorSDB,log=TRUE)
    canbet = rnorm(1,b,0.2)
    
    canll=LL(alph,canbet, N-1 ,Xcov,x2,Quienesg-1,Npg=JJ+1,LOG=TRUE)
    prior_canll=dnorm(alph,mean=priorB,sd=priorSDB,log=TRUE)
    
    MH <- canll-curll+prior_canll-prior_curll
    
    if(log(runif(1))<MH){
      b  = canbet
      curll  = canll
    }
    
    #sample Xcov Gibbs 
    
    Xcov=sampleX(alph,b,N-1,Xcov,x2, Xobs,Quienesg,JJ+1,length(Xcov),Q)
    
    # for (j in 1:length(Xcov))
    # {
    # 
    # Xnew=Xcov
    # 
    # ## Xj=0
    # Xcan=Xcov
    # Xcan[j]=0
    # a0=LL(alph,b, N-1 ,Xcan,x2,Quienesg-1,Npg=JJ+1,LOG=FALSE)*Q[1,Xobs[j]+1]
    # 
    # ## Xj=1
    # Xcan=Xcov
    # Xcan[j]=1
    # a1=LL(alph,b, N-1 ,Xcan,x2,Quienesg-1,Npg=JJ+1,LOG=FALSE)*Q[2,Xobs[j]+1]
    # p1=a1/(a0+a1)
    # 
    # # sample new Xj
    # Xnew[j]=rbinom(1,1,prob=p1)
    # Xcov=Xnew
    # 
    # }
    
    # X.newCov=numeric(length(d$Xcov))
    # # puede ir arriba esto
    # quienes.1=which(d$covE==1)
    # quienes.0=which(d$covE==0)
    # 
    # X.newCov[quienes.0]=rbinom(length(quienes.0),size=1,prob=Q[2,1])
    # X.newCov[quienes.1]=rbinom(length(quienes.1),size=1,prob=Q[2,2])
    # 
    # canCov=X.newCov

    
    # one sample a time 
    # Xcan=Xcov
    # for (j in 1:length(Xcov))
    # {
    # Xcan[j]=0
    # if (Xcov[j]==0){Xcan[j]=1}  
    # #Xcan[j]=rbinom(1,1,0.05)    
    # canll=distY3(alph,b,N-1,Xcan,x2,quienes.g)
    # canXXobs=lQ[Xcan[j]+1,Xobs[j]+1]
    # curXXobs=lQ[Xcov[j]+1,Xobs[j]+1]
    # prior_can=lpii[Xcan[j]+1]
    # prior_cur=lpii[Xcov[j]+1]
    # 
    # MH <- canll+canXXobs-curll-curXXobs#+prior_can-prior_cur
    # if(log(runif(1))<MH){
    #   curll  = canll
    #  Xcov=Xcan 
    #   }
    # }    
    
    
    
    
    # sample Xcov (MH)
    # s=rbinom(length(Xcov),1,0.5)
    # canCov= Xcov
    # w.0=canCov[s==1]==0
    # w.1=canCov[s==1]==1
    # canCov[s==1][w.0]=1
    # canCov[s==1][w.1]=0
    # 
    # cc=0
    # for(j in 1:length(Xcov))
    # {cc=cc+lQ[canCov[j]+1,Xobs[j]+1]}
    # 
    # canxxobs=cc  
    # # d.can=d
    # # d.can$Xcv=canCov
    # ##canll=distY(alph,b,N-1,d.can)  
    # canll=LL(alph,b, N-1 ,canCov,x2,Quienesg-1,Npg=JJ+1,LOG=TRUE)
    # MH <- canll+canxxobs-curll-curxxobs
    # 
    # 
    # if(log(runif(1))<MH){
    #   curll  = canll
    #   curxxobs=canxxobs
    #   Xcov=canCov
    #   ncambios=ncambios+1} 
    
    keep.alpha[i]=alph  
    keep.beta[i]=b  
    print(alph)
    print(i) 
    
  }
  
return(list(out.alpha=keep.alpha,out.beta=keep.beta,Xcov=Xcov,ncambios=ncambios))  
  
}
  
MCMC_Nerror=function(dd,iters=4000,chains=2,
                     priorA,priorSDA,priorB,priorSDB,alph,b,N,Xcov)
{
  d=dd
  x2=dd$ral2
  N=length(unique(d$indx))+1
  Quienesg=matrix(0,JJ+1,N-1)
  for (i in 1:N-1)
  {Quienesg[,i]=which(dd$indx==i)}
  
  keep.alpha= rep(0,iters)
  keep.beta= rep(0,iters)
  
  # set first values
  curll=LL(alph,b, N-1 ,Xcov,x2,Quienesg-1,Npg=JJ+1,LOG=TRUE)
  
  #MCMC from here. GO!
  
  for(i in 1:iters){
    
    # Sample alpha
    prior_curll=dnorm(alph,mean=priorA,sd=priorSDA,log=TRUE)
    canalph = rnorm(1,alph,0.5)
    
    canll=LL(canalph,b, N-1 ,Xcov,x2,Quienesg-1,Npg=JJ+1,LOG=TRUE)
    prior_canll=dnorm(canalph,mean=priorA,sd=priorSDA,log=TRUE)
    
    MH <- canll-curll+prior_canll-prior_curll
    
    if(log(runif(1))<MH){
      alph  = canalph
      curll  = canll
    }
    
    # Sample beta
    prior_curll=dnorm(b,mean=priorB,sd=priorSDB,log=TRUE)
    canbet = rnorm(1,b,0.2)
    
    canll=LL(alph,canbet, N-1 ,Xcov,x2,Quienesg-1,Npg=JJ+1,LOG=TRUE)
    prior_canll=dnorm(alph,mean=priorB,sd=priorSDB,log=TRUE)
    
    MH <- canll-curll+prior_canll-prior_curll
    
    if(log(runif(1))<MH){
      b  = canbet
      curll  = canll
    }
    
    keep.alpha[i]=alph  
    keep.beta[i]=b  

    
  }
  return(list(out.alpha=keep.alpha,out.beta=keep.beta))  
  
  }