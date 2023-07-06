### Create environmental covariate 
library(PtProcess)
library(tidyverse)
library(ggplot2)
library(sp)
library(spdep)
library(mvtnorm)
library(geoR)
library(gstat)
library(lctools)
source('ExtraFunctions.R')


xx=function(v)
{
  #if (abs(v[1]-5)<5 & abs(v[2]-5)<5)
  #{
    dist=(v[1]-5)^2+(v[2]-5)^2
  out=dist
  #}
  #out=cos(dist*22)
  #else 
  #{
  #  out=50
  #}´
  return(out)
}


gridd <- GridTopology(c(-5,-5), c(.5,.5), c(40,40)) # create the grid

datafr <- data.frame(cov=apply(coordinates(gridd),MARGIN = 1,xx))# calculate de xx
datafr$cov=(max(datafr$cov)-datafr$cov)/max(datafr$cov)
range(datafr$cov)
en.c <- SpatialGridDataFrame(gridd, datafr) # create the grid data frame

plot(en.c,axes=TRUE,col=hcl.colors(50))


##### Moran and Geary

# set neighborhood
# spqt <- as(en.c, "SpatialPolygons")
# wr <- poly2nb(spqt, queen=FALSE) 
# wm <- nb2mat(wr, style='W', zero.policy = TRUE)
# wrl <- nb2listw(wr, style="W") 

# Mt=moran.test(vals,wrl)
# Gt=geary.test(vals, wrl)
# 
# plot(en.c)
# v=variogram(vals~1, data=en.c)
# plot(v)




# bws=seq(1,20,by=.2)
# vals=en.c$cov
# MIC=moransI.v(coordinates(en.c),bws,vals,family='fixed')
# MIC=data.frame(MIC)
# 
# do=min(which(abs(MIC$Moran.s.I)<0.01))
# d=MIC$k[do]

