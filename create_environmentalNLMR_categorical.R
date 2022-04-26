library(NLMR)
library(raster)

d=10000
x=nlm_gaussianfield(ncol=100,nrow = 100,autocorr_range =d,resolution = 1,mean=0,
                    mag_var=1,nug=0.1,user_seed = 200)

lx=length(values(x))
values(x)=runif(lx,min=0,max=1)
xval=values(x)

xval[values(x)<pi1]=1
xval[values(x)>pi1]=0

values(x)=xval
plot(x,asp=1)

en.c <- as(x, 'SpatialGridDataFrame')
names(en.c)='cov'
en.c$cov
