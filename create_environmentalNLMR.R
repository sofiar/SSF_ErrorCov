library(NLMR)
library(raster)

d=10000
x=nlm_gaussianfield(ncol=100,nrow = 100,autocorr_range =d,resolution = 1,
                    mag_var=1,nug=0.1,user_seed = 200)


values(x)=rnorm(length(values(x)),0,1)
#values(x)=(values(x)-0.5)#/sd(values(x))

plot(x,asp=1)


en.c <- as(x, 'SpatialGridDataFrame')
names(en.c)='cov'
en.c$cov

# bws=seq(1,30,by=1)
# vals=en.c$cov
# MIC=moransI.v(coordinates(en.c),bws,vals,family='fixed')
# MIC=data.frame(MIC)
# do=min(which(abs(MIC$Moran.s.I)<0.01))
# d=MIC$k[do] # Distancia



