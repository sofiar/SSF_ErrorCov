## extra functions

# C coordinate matrix
StepLengths=function(C)
{
  X=C[,1]
  Y=C[,2]
  n = length(X)
  adj = X[2:n]-X[1:n-1] 
  op = Y[2:n]-Y[1:n-1]
  step = (adj^2 + op^2)^0.5
  aca = which(step==0)
  step[aca]=0.0000001
  return(step)
}
  
e.dist=function(p,v)
{
  out=sqrt((v[1]-p[1])^2+(v[2]-p[2])^2)
  return(out)
}

CreateQuad=function(track.data,time.int,tile.dim,en.c) #,rad.x,rad.y
{
  # make quad times
  times <- track.data$time - min(track.data$time)
  qt <- seq(0, ceiling(max(times)/time.int)*time.int, time.int)
  qt <- sort(c(times, qt))
  qt <- qt[!duplicated(qt)]
  
  lengths=c(qt[1],diff(qt))
  lengths[2] <- sum(lengths[1:2])
  
  temp.quad <- data.frame(times=qt[-1], length=lengths[-1], 
                          obs=ifelse(qt%in%times, 1, 0)[-1])
  
  temp.quad$x = temp.quad$y = NA
  temp.quad$x[temp.quad$obs==1]= coordinates(track.data)[-1,1]
  temp.quad$y[temp.quad$obs==1] = coordinates(track.data)[-1,2]
  
  
  
  delta <- diff(times)
  d1 <- spDists(track.data)
  pwd <- d1[col(d1)==(row(d1)+1)] #steplength
  disp.per.time <- pwd/delta #vels
  
  disp.rad.lookup <- max(disp.per.time)
  #disp.rad.lookup <- 50
  
  disp.rad.lookup.X <- disp.rad.lookup/tile.dim[1]
  disp.rad.lookup.Y <- (disp.rad.lookup.X*tile.dim[1])/tile.dim[2]
  
  mu <- track.data@coords[1,]
  time.last <- 0
  df <- NULL
  quad.time <- temp.quad$time
  
  pb <- txtProgressBar(min = 0, max = nrow(temp.quad), style = 3)
  
  for(i in 1:nrow(temp.quad)){
    delta.last <- quad.time[i]-time.last
    
    rad.x <- ceiling(disp.rad.lookup.X*delta.last)
    rad.y <- ceiling(disp.rad.lookup.Y*delta.last)
    
  
    x.grid <- c(seq(mu[1]-rad.x*tile.dim[1], mu[1]-tile.dim[1], tile.dim[1]), mu[1], seq(mu[1]+tile.dim[1], 
                                                                                     mu[1]+rad.x*tile.dim[1], tile.dim[1]))
    y.grid <-  c(seq(mu[2]-rad.y*tile.dim[2], mu[2]-tile.dim[2], tile.dim[2]), mu[2], seq(mu[2]+tile.dim[2], mu[2]+rad.y*tile.dim[2], tile.dim[2]))
    
    #para bicho que salta
    # x.grid[x.grid>10]=x.grid[x.grid>10]-10
    # x.grid[x.grid<0]=x.grid[x.grid<0]+10
    # y.grid[y.grid>10]=y.grid[y.grid>10]-10
    # y.grid[y.grid<0]=y.grid[y.grid<0]+10
    # 
    
    spqt <- as(SpatialPoints(expand.grid(x=x.grid, y=y.grid)), "SpatialPixels")
    #spqt2 <- spqt
    
    
    spqtPoly <- as(spqt, "SpatialPolygons")
    if(temp.quad$obs[i]==1){spqt <- rbind(SpatialPoints(temp.quad[i,c("x","y")]), as(spqt,"SpatialPoints"))
    }else {spqt <- as(spqt, "SpatialPoints")}
    
    out <- as(spqt, "data.frame")
    out$t <- temp.quad$time[i]
    out$delta.last <- delta.last
    out$area <- prod(tile.dim)
    out$length<- temp.quad$length[i]
    out$volume <- out$area*temp.quad$length[i]
    out$obs <- rep(0,length(out$area))
    
    # Chequear esto
    if(temp.quad$obs[i]==1){ 
      #out$area <- out$area/c(2, table(spqt%over%spqtPoly))
      
      cord.mu=SpatialPoints(temp.quad[i,c("x","y")])%over%spqtPoly
      div=rep(1,length(out$area))
      div[cord.mu+1]=2
      div[1]=2
      out$area <- out$area/div
      out$obs[1] <- 1
      
      
    } 
    out$volume <- out$area*temp.quad$length[i]
    out$response <- out$obs/out$volume  
    
    out$dist=spDistsN1(spqt, mu)
    
  
    #out$cov=apply(coordinates(spqt),MARGIN=1,xx)
    cov = spqt %over% en.c
    out  <- cbind(out, cov)
    
    out$bmKern <- -0.5*out$dist^2/out$delta.last 
    
    df <- rbind(df, out) # pegoteo a lo anterior
    if(temp.quad$obs[i]==1){
      mu <- as.numeric(temp.quad[i,c("x","y")])
      time.last <- temp.quad$time[i]
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
return(df)
}

xx1=function(v)
{
  dist=(v[1]-2)^2+(v[2]-2)^2
  out=dist
  return(out)
}

xx2=function(v)
{
  dist=(v[1])^2+(v[2])^2
  out=dist
  return(out)
}

pathelements <- function(X,Y){
  n = length(X)
  adj = X[2:n]-X[1:n-1] 
  op = Y[2:n]-Y[1:n-1]
  step = (adj^2 + op^2)^0.5
  aca = which(step==0)
  step[aca]=0.0000001
  
  si<-sign(adj)
  si[si==0]<-1  #corrects for sign(0) == 0
  ang = si*(op<0)*pi+atan(adj/op)
  adif <- ang[2:length(ang)]-ang[1:(length(ang)-1)]
  
  ## corregimos para que quede entre -pi y pi
  adif[which(adif < -pi)] = adif[which(adif < -pi)] + 2*pi
  adif[which(adif > pi)] = adif[which(adif > pi)] - 2*pi  
  
  turns <- adif
  direction <- ang  
  co = adj/step # cosine
  si = op/step  # sine
  res = list(step,turns,direction,co,si)
  names(res) = c("steps","turns","direction","cosine","sine")
  return(res)
}


sim.traj=function(alpha,beta,nsims,delta,times,en.c,J=1000)
{
  for (jj in 1:nsims)
  {
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
        pol=s.sample*alpha
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
  return(salida)
}

create.control=function(JJ=150,salida,ebeta,beta.random=FALSE)
{
  case=c()
  x_=c()
  y_=c()
  indx=c()
  nsim=c()
  num.traj=c()
  ral=c()
  nsims=length(salida)
  
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
      
      if (beta.random)
      {
        ebeta=runif(1,min=0,max=10) # check this   
      }
      
      else{ebeta=beta}
      
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
    #print(i)
  }
  
  control.data=data.frame(case,indx,x_,y_,nsim,num.traj,ral)
  return(control.data)
}


