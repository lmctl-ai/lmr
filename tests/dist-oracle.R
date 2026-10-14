library(stats)
values<-c(qt(.75,df=7),qchisq(.75,df=7),qgamma(.75,shape=2,rate=3),qbeta(.75,shape1=2,shape2=5),qbinom(.75,size=10,prob=.4),qpois(.75,lambda=3))
writeLines(sprintf('%.17g',values),commandArgs(TRUE)[[1]])
