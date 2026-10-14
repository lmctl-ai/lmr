args <- commandArgs(TRUE)
library(stats)
out <- args[[1]]
save <- function(name, x) writeLines(sprintf("%.17g", as.numeric(x)), file.path(out, paste0(name, ".txt")))
families <- list(
  norm=list(param=list(),d=function(x)dnorm(x),p=function(x)pnorm(x)),
  t=list(param=list(df=7),d=function(x)dt(x,7),p=function(x)pt(x,7)),
  chisq=list(param=list(df=7),d=function(x)dchisq(x,7),p=function(x)pchisq(x,7)),
  gamma=list(param=list(shape=2,rate=3),d=function(x)dgamma(x,2,rate=3),p=function(x)pgamma(x,2,rate=3)),
  beta=list(param=list(a=2,b=5),d=function(x)dbeta(x,2,5),p=function(x)pbeta(x,2,5)),
  binom=list(param=list(size=10,prob=.4),d=function(x)dbinom(x,10,.4),p=function(x)pbinom(x,10,.4)),
  pois=list(param=list(lambda=3),d=function(x)dpois(x,3),p=function(x)ppois(x,3)))
for(n in names(families)){f<-families[[n]];save(paste0("density-",n),f$d(c(0,1,2)));save(paste0("probability-",n),f$p(c(0,1,2)))}
save("probability-log-tail",pnorm(-40,log.p=TRUE))
RNGkind("Mersenne-Twister","Inversion","Rejection")
draws <- list(norm=function()rnorm(3),unif=function()runif(3),exp=function()rexp(3,2),gamma=function()rgamma(3,2,rate=3),beta=function()rbeta(3,2,5),binom=function()rbinom(3,10,.4),pois=function()rpois(3,3),t=function()rt(3,7),chisq=function()rchisq(3,7))
for(n in names(draws)){set.seed(42);save(paste0("rng-",n),draws[[n]]())}
set.seed(42);save("sample",sample(1:8,5))
set.seed(42);save("sample-weighted",sample(1:4,6,replace=TRUE,prob=c(1,2,3,4)))
x<-c(.5,1,2,3)
math<-list(gamma=gamma(x),lgamma=lgamma(x),beta=beta(x,c(2,3,4,5)),lbeta=lbeta(x,c(2,3,4,5)),digamma=digamma(x),trigamma=trigamma(x),choose=choose(c(5,8,9,10),c(2,3,4,5)),lchoose=lchoose(c(5,8,9,10),c(2,3,4,5)),besselJ=besselJ(x,.5),besselI=besselI(x,.5),besselK=besselK(x,.5),besselY=besselY(x,.5))
for(n in names(math))save(paste0("math-",n),math[[n]])
qx<-c(1,2,4,8,16);pr<-c(0,.1,.25,.5,.75,.9,1)
for(t in 1:9)save(paste0("quantile-",t),quantile(qx,pr,type=t,names=FALSE))
z<-fft(c(1,2,3,4));save("fft-re",Re(z));save("fft-im",Im(z));z<-fft(c(1,2,3,4),inverse=TRUE)/4;save("ifft-re",Re(z));save("ifft-im",Im(z))
