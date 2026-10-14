a <- commandArgs(TRUE); outfile <- a[[1]]
library(stats)
q <- qnorm(0.975)
RNGkind('Mersenne-Twister','Inversion','Rejection'); set.seed(42); draws <- rnorm(2)
d <- data.frame(y=c(1,2,4,5),x=c(1,2,3,4)); s<-summary(lm(y~x,data=d)); c<-s$coefficients
values<-c(q,draws,c[,1],c[,2],c[,3],c[,4],s$sigma,s$r.squared,s$df[2])
writeLines(c(as.character(getRversion()),paste(RNGkind(),collapse='|'),sprintf('%.17g',values)),outfile)
