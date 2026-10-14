args <- commandArgs(TRUE)
out <- args[[1L]]
library(stats)
con<-file(out,open='wt')
write_case <- function(name, values) writeLines(paste(c(name, sprintf('%.17g', as.numeric(values))), collapse=' '), con, sep='\n', useBytes=TRUE)
rosen <- function(p) { x<-p[[1L]];y<-p[[2L]];100*(y-x*x)^2+(1-x)^2 }
for (method in c('BFGS','Nelder-Mead')) {
 z<-optim(c(a=-1.2,b=1),rosen,method=method,control=list(maxit=1000,reltol=1e-8))
 write_case(paste0('rosen_',gsub('-','_',method)),c(z$par,z$value,z$convergence,z$counts))
}
box<-optim(c(a=0.5),(function(p)(p[[1L]]-4)^2),method='L-BFGS-B',lower=c(a=0),upper=c(a=2),control=list(maxit=1000,reltol=1e-8))
write_case('box',c(box$par,box$value,box$convergence,box$counts))
RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(42)
sann<-optim(c(a=-1),(function(p)(p[[1L]]-2)^2),method='SANN',control=list(maxit=80,reltol=1e-8))
write_case('sann',c(sann$par,sann$value,sann$convergence,sann$counts))
o<-optimize(function(x)(x-2)^2,c(-5,5),tol=1e-8);write_case('optimize',c(o$minimum,o$objective))
r<-uniroot(function(x)x^2-2,c(0,2),tol=1e-8,maxiter=100);write_case('root',c(r$root,r$f.root,r$iter,r$estim.prec))
i1<-integrate(function(x)exp(-x^2),-Inf,Inf,rel.tol=1e-8);write_case('integral_inf',c(i1$value,i1$abs.error,i1$subdivisions))
i2<-integrate(function(x)exp(-x),0,1,rel.tol=1e-8);write_case('integral_finite',c(i2$value,i2$abs.error,i2$subdivisions))
x<-seq(0,4,length.out=17);y<-4*x/(2+x)+c(0,.01,-.01,.006,-.004,.003,-.006,.002,.004,-.003,.001,-.002,.005,-.004,.002,-.001,.003)
d<-data.frame(x=x,y=y);fit<-nls(y~a*x/(b+x),data=d,start=c(a=3,b=1),control=nls.control(maxiter=100,tol=1e-8,warnOnly=TRUE));sm<-summary(fit)
write_case('nls',c(coef(fit),sm$coefficients[,2],sm$coefficients[,3],sm$coefficients[,4],sm$sigma,df.residual(fit),fit$convInfo$isConv,fit$convInfo$finIter,fit$convInfo$finTol))
close(con)
