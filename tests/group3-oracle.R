args<-commandArgs(TRUE);library(stats);out<-args[[1]]
save<-function(n,x)writeLines(sprintf("%.17g",as.numeric(x)),file.path(out,paste0(n,".txt")))
save_text<-function(n,x)writeLines(as.character(x),file.path(out,paste0(n,".txt")))
dat<-list(x=1:8,yg=c(1,2,3,5,7,8,9,11),yb=c(0,1,0,1,1,0,1,0),yp=c(1,0,2,1,3,4,2,5),yG=c(1.2,2.1,1.7,3,4.2,2.5,5,4.6))
for(f in c("gaussian","binomial","poisson","Gamma")){
 yname<-switch(f,gaussian="yg",binomial="yb",poisson="yp",Gamma="yG");link<-switch(f,gaussian="identity",binomial="logit",poisson="log",Gamma="inverse");fit<-glm(as.formula(paste(yname,"~ x")),data=as.data.frame(dat),family=switch(f,gaussian=gaussian(link),binomial=binomial(link),poisson=poisson(link),Gamma=Gamma(link)),control=glm.control(maxit=50))
 sm<-summary(fit);key<-paste0("glm-",f,"-");save(paste0(key,"estimate"),coef(fit));save(paste0(key,"se"),sm$coefficients[,2]);save(paste0(key,"stat"),sm$coefficients[,3]);save(paste0(key,"p"),sm$coefficients[,4]);save(paste0(key,"deviance"),fit$deviance);save(paste0(key,"null"),fit$null.deviance);save(paste0(key,"aic"),fit$aic);save(paste0(key,"iter"),fit$iter);save_text(paste0(key,"converged"),fit$converged);save_text(paste0(key,"names"),names(coef(fit)))
}
dat$yb<-c(0,0,0,0,1,1,1,1);sep<-glm(as.formula("yb ~ x"),data=as.data.frame(dat),family=binomial(),control=glm.control(maxit=1));save_text("glm-nonconv",sep$converged);save_text("glm-nonconv-status",if(sep$converged)"OK" else "NON_CONVERGED")
x<-1:12;y<-c(0,1,0,2,1,3,2,4,3,5,4,6);sp<-smooth.spline(x,y,df=4);save("spline-fitted",sp$y);save("spline-new",predict(sp,x=c(2.5,6.5,10.5))$y);save("spline-lambda",sp$lambda);save("spline-df",sp$df)
sp<-smooth.spline(x,y,spar=.5);save("spline-spar",sp$y)
sp<-smooth.spline(x,y,cv=TRUE);save("spline-cv",sp$y)
lo<-loess(y~x,data=data.frame(x=x,y=y),span=1,degree=1,family="gaussian");save("loess-fitted",fitted(lo));save("loess-new",predict(lo,newdata=data.frame(x=c(2.5,6.5,10.5))))
sx<-1:6;sy<-c(0,1,0,-1,0,0);xo<-c(1.5,2.5,3.5,4.5,5.5)
for(m in c("fmm","natural","periodic"))save(paste0("is-",m),spline(sx,sy,xout=xo,method=m)$y)
save("is-monoH.FC",splinefun(sx,sy,method="monoH.FC")(xo))
save("approx-linear",approx(sx,sy,xout=xo,method="linear",rule=1)$y);save("approx-constant",approx(sx,sy,xout=c(.5,1.5,5.5,6.5),method="constant",rule=2)$y)
