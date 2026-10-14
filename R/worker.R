`%||%` <- function(x,y) if(is.null(x)) y else x
library(stats)
library(utils)
parse_json <- function(s) {
 ch<-strsplit(enc2utf8(s), '', fixed=TRUE)[[1]]; i<-1L; n<-length(ch); peek<-function() if(i<=n) ch[[i]] else ''
 ws<-function(){while(i<=n&&peek()%in%c(' ','\n','\r','\t'))i<<-i+1L}
 hex4<-function(){if(i+3L>n)stop('bad unicode escape');v<-paste(ch[i:(i+3L)],collapse='');if(!grepl('^[0-9A-Fa-f]{4}$',v))stop('bad unicode escape');i<<-i+4L;strtoi(v,16L)}
 string<-function(){i<<-i+1L; z<-character(n-i+1L);k<-0L;push<-function(v){k<<-k+1L;z[[k]]<<-v};while(i<=n){c<-peek();i<<-i+1L;if(c=='"')return(if(k)paste(z[seq_len(k)],collapse='')else '');if(c=='\\'){e<-peek();i<<-i+1L;if(e=='u'){cp<-hex4();if(cp>=0xD800&&cp<=0xDBFF){if(i+1L>n||ch[[i]]!='\\'||ch[[i+1L]]!='u')stop('bad surrogate pair');i<<-i+2L;lo<-hex4();if(lo<0xDC00||lo>0xDFFF)stop('bad surrogate pair');cp<-0x10000+(cp-0xD800)*1024+(lo-0xDC00)};push(intToUtf8(cp))}else{v<-switch(e,'"'='"','\\'='\\','/'='/',b='\b',f='\f',n='\n',r='\r',t='\t',NULL);if(is.null(v))stop('bad string escape');push(v)}}else push(c)};stop('unterminated JSON string')}
 value<-function(){ws();c<-peek();if(c=='"')return(string());if(c=='{'){i<<-i+1L;obj<-list();ws();if(peek()=='}'){i<<-i+1L;return(obj)};repeat{ws();k<-string();ws();if(peek()!=':')stop('expected colon');i<<-i+1L;obj[k]<-list(value());ws();c<-peek();i<<-i+1L;if(c=='}')break;if(c!=',')stop('expected comma')};return(obj)};if(c=='['){i<<-i+1L;arr<-list();ws();if(peek()==']'){i<<-i+1L;return(arr)};repeat{arr[length(arr)+1L]<-list(value());ws();c<-peek();i<<-i+1L;if(c==']')break;if(c!=',')stop('expected comma')};return(arr)};for(k in c('true','false','null'))if(i+nchar(k)-1L<=n&&paste(ch[i:(i+nchar(k)-1L)],collapse='')==k){i<<-i+nchar(k);return(switch(k,true=TRUE,false=FALSE,null=NULL))};start<-i;if(peek()=='-')i<<-i+1L;while(i<=n&&grepl('^[0-9]$',peek()))i<<-i+1L;if(peek()=='.'){i<<-i+1L;while(i<=n&&grepl('^[0-9]$',peek()))i<<-i+1L};if(peek()%in%c('e','E')){i<<-i+1L;if(peek()%in%c('+','-'))i<<-i+1L;while(i<=n&&grepl('^[0-9]$',peek()))i<<-i+1L};if(i==start)stop('invalid JSON value');as.numeric(paste(ch[start:(i-1L)],collapse=''))}
 value()
}
numeric_input <- function(x) { if(is.null(x)) return(NA_real_); if(is.character(x)&&length(x)==1L) { if(x=='NA') return(NA_real_); return(as.numeric(x)) }; as.numeric(x) }
numeric_file <- function(path) scan(path,what=double(),quiet=TRUE,na.strings='NA')
args <- commandArgs(TRUE); if(length(args)&&args[[1]]=='--') args<-args[-1]; op <- args[[1]]; req <- parse_json(paste(readLines(args[[2]],warn=FALSE),collapse='')); outpath <- args[[3]]
warnings <- character()
encode <- function(x) {
 if (is.double(x)) { if(length(x)!=1L) return(array_encode(x));v<-x[[1L]];if(is.nan(v))return('NaN');if(is.na(v))return('NA');if(is.infinite(v))return(if(v>0)'Inf' else '-Inf');return(v) }
 if (is.list(x)) return(lapply(x, encode))
 x
}
array_encode <- function(x) structure(as.double(x),class='rcli_array')
json_array <- function(x) {
 if(is.numeric(x)) {
  if(!length(x)) return('[]')
  parts<-sprintf('%.17g',x); nan<-is.nan(x); na<-is.na(x)&!nan; inf<-is.infinite(x)
  parts[nan]<-'"NaN"';parts[na]<-'"NA"';parts[inf&x>0]<-'"Inf"';parts[inf&x<0]<-'"-Inf"'
  return(paste0('[',paste(parts,collapse=','),']'))
 }
 paste0('[',paste(vapply(as.list(x),json,''),collapse=','),']')
}
json <- function(x) {
 if(is.null(x)) return('null')
 if(inherits(x,'rcli_array')) return(json_array(unclass(x)))
 if(is.list(x)) {
  n<-names(x)
  if(!is.null(n) && all(n!='')) return(paste0('{',paste(vapply(seq_along(x),function(i) paste0(json(n[[i]]),':',json(x[[i]])),''),collapse=','),'}'))
  return(paste0('[',paste(vapply(x,json,''),collapse=','),']'))
 }
 if(length(x)==0) return('[]')
 if(length(x)>1) return(json_array(x))
 if(is.character(x)) { s<-gsub('\\\\','\\\\\\\\',x);s<-gsub('"','\\\\"',s,fixed=TRUE);s<-gsub('\n','\\\\n',s,fixed=TRUE);s<-gsub('\r','\\\\r',s,fixed=TRUE);s<-gsub('\t','\\\\t',s,fixed=TRUE);return(paste0('"',s,'"')) }
 if(is.logical(x)) return(if(is.na(x)) 'null' else if(x) 'true' else 'false')
 if(is.numeric(x)) { if(is.nan(x)) return('"NaN"'); if(is.na(x)) return('"NA"'); if(is.infinite(x)) return(if(x>0) '"Inf"' else '"-Inf"'); return(sprintf('%.17g',x)) }
 'null'
}
result <- withCallingHandlers(tryCatch({
 if(op=='r_info') list(r_version=as.character(getRversion()),rng_kind=RNGkind())
 else if(op=='dist_quantile') {
  fn <- switch(req$dist,norm=qnorm,t=qt,chisq=qchisq,gamma=qgamma,beta=qbeta,binom=qbinom,pois=qpois)
  if(is.null(fn)) stop('unknown distribution')
  p <- if(!is.null(req$p_file)) scan(req$p_file,what=double(),quiet=TRUE,na.strings='NA') else vapply(req$p,numeric_input,numeric(1)); args2 <- list(p=p,lower.tail=if(is.null(req$lower_tail)) TRUE else req$lower_tail,log.p=if(is.null(req$log_p)) FALSE else req$log_p)
  if(req$dist=='norm') {args2$mean<-if('mean'%in%names(req)) numeric_input(req$mean) else 0;args2$sd<-if('sd'%in%names(req)) numeric_input(req$sd) else 1}
  if(req$dist%in%c('t','chisq')) args2$df<-numeric_input(req$df)
  if(req$dist=='gamma') {args2$shape<-numeric_input(req$shape);if('rate'%in%names(req))args2$rate<-numeric_input(req$rate)}
  if(req$dist=='beta') {args2$shape1<-numeric_input(req$a);args2$shape2<-numeric_input(req$b)}
  if(req$dist=='binom') {args2$size<-numeric_input(req$size);args2$prob<-numeric_input(req$prob)}
  if(req$dist=='pois') args2$lambda<-numeric_input(req$lambda)
  list(value=array_encode(do.call(fn,args2)),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else if(op%in%c('dist_density','dist_probability')) {
  density<-op=='dist_density';family<-req$dist
  fn<-switch(family,norm=if(density)dnorm else pnorm,t=if(density)dt else pt,chisq=if(density)dchisq else pchisq,gamma=if(density)dgamma else pgamma,beta=if(density)dbeta else pbeta,binom=if(density)dbinom else pbinom,pois=if(density)dpois else ppois)
  arg<-list(if(density)numeric_file(req$x_file) else numeric_file(req$q_file));names(arg)<-if(density)'x' else 'q'
  if(family=='norm'){arg$mean<-numeric_input(req$mean);arg$sd<-numeric_input(req$sd)}
  if(family%in%c('t','chisq'))arg$df<-numeric_input(req$df)
  if(family=='gamma'){arg$shape<-numeric_input(req$shape);if(!is.null(req$rate))arg$rate<-numeric_input(req$rate)}
  if(family=='beta'){arg$shape1<-numeric_input(req$a);arg$shape2<-numeric_input(req$b)}
  if(family=='binom'){arg$size<-numeric_input(req$size);arg$prob<-numeric_input(req$prob)}
  if(family=='pois')arg$lambda<-numeric_input(req$lambda)
  if(density)arg$log<-isTRUE(req$log) else {arg$lower.tail<-isTRUE(req$lower_tail);arg$log.p<-isTRUE(req$log_p)}
  list(value=array_encode(do.call(fn,arg)),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else if(op=='rng_draw') {
  RNGkind('Mersenne-Twister','Inversion','Rejection'); set.seed(req$seed)
  n<-req$n; value<-switch(req$dist,norm=rnorm(n,numeric_input(req$mean),numeric_input(req$sd)),unif=runif(n,numeric_input(req$min),numeric_input(req$max)),exp=rexp(n,numeric_input(req$rate)),gamma=rgamma(n,numeric_input(req$shape),rate=numeric_input(req$rate)),beta=rbeta(n,numeric_input(req$a),numeric_input(req$b)),binom=rbinom(n,numeric_input(req$size),numeric_input(req$prob)),pois=rpois(n,numeric_input(req$lambda)),t=rt(n,numeric_input(req$df)),chisq=rchisq(n,numeric_input(req$df)))
  list(value=array_encode(value),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind(),seed=req$seed))
 } else if(op=='rng_sample') {
  RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(req$seed);x<-numeric_file(req$x_file);w<-if(is.null(req$prob_file))NULL else numeric_file(req$prob_file);list(value=array_encode(sample(x,size=req$size,replace=isTRUE(req$replace),prob=w)),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind(),seed=req$seed))
 } else if(op=='math_special') {
  x<-numeric_file(req$x_file);fn<-switch(req$fn,gamma=gamma,lgamma=lgamma,beta=beta,lbeta=lbeta,digamma=digamma,trigamma=trigamma,choose=choose,lchoose=lchoose,besselJ=besselJ,besselI=besselI,besselK=besselK,besselY=besselY)
  value<-if(req$fn%in%c('beta','lbeta','choose','lchoose'))do.call(fn,list(x,numeric_file(req$y_file))) else if(grepl('^bessel',req$fn))do.call(fn,list(x,numeric_input(req$nu))) else fn(x)
  list(value=array_encode(value),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else if(op=='summary_quantile') {
  v<-quantile(numeric_file(req$x_file),probs=numeric_file(req$probs_file),na.rm=isTRUE(req$na_rm),names=FALSE,type=req$type);list(value=array_encode(as.numeric(v)),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind(),type=req$type))
 } else if(op=='num_fft') {
  x<-numeric_file(req$x_file);z<-fft(x,inverse=isTRUE(req$inverse));if(isTRUE(req$inverse))z<-z/length(x);list(re=array_encode(Re(z)),im=array_encode(Im(z)),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind(),inverse=isTRUE(req$inverse)))
 } else if(op=='test_htest') {
  x<-numeric_file(req$x_file);y<-if(is.null(req$y_file))NULL else numeric_file(req$y_file);alt<-req$alternative;cl<-req$conf_level;tst<-req$test
  h<-switch(tst,
   t=if(is.null(y))t.test(x,mu=numeric_input(req$mu),alternative=alt,conf.level=cl) else t.test(x,y,mu=numeric_input(req$mu),alternative=alt,paired=isTRUE(req$paired),var.equal=isTRUE(req$var_equal),conf.level=cl),
   wilcox=wilcox.test(x,y,mu=numeric_input(req$mu),alternative=alt,paired=isTRUE(req$paired),conf.int=isTRUE(req$conf_int),conf.level=cl),
   chisq={xx<-if(!is.null(req$rows))matrix(x,nrow=req$rows,ncol=req$cols) else x;if(is.null(req$expected_file))chisq.test(xx,correct=isTRUE(req$correct)) else chisq.test(xx,p=numeric_file(req$expected_file),correct=isTRUE(req$correct))},
   shapiro=shapiro.test(x),
   ks={fn<-switch(req$distribution,norm=pnorm,unif=punif,exp=pexp,pois=ppois);args<-switch(req$distribution,norm=list(mean=numeric_input(req$mean),sd=numeric_input(req$sd)),unif=list(min=numeric_input(req$min),max=numeric_input(req$max)),exp=list(rate=numeric_input(req$rate)),pois=list(lambda=numeric_input(req$lambda)));do.call(ks.test,c(list(x=x,y=fn,alternative=alt),args))},
   fisher=fisher.test(matrix(x,nrow=2,ncol=2),alternative=alt,conf.level=cl),
   binom=binom.test(as.integer(x[[1L]]),req$n,p=numeric_input(req$p),alternative=alt,conf.level=cl),
   prop={pp<-if(is.null(req$probabilities_file))numeric_input(req$p) else numeric_file(req$probabilities_file);prop.test(as.integer(x),as.integer(numeric_file(req$trials_file)),p=pp,alternative=alt,correct=isTRUE(req$correct),conf.level=cl)},
   cor=cor.test(x,y,alternative=alt,method=req$method,conf.level=cl),
   var=var.test(x,y,ratio=numeric_input(req$ratio),alternative=alt,conf.level=cl))
  list(statistic=if(length(h$statistic))as.numeric(h$statistic)[[1L]]else NULL,parameter=if(!length(h$parameter))NULL else if(length(h$parameter)==1L)as.numeric(h$parameter)[[1L]]else array_encode(as.numeric(h$parameter)),p_value=as.numeric(h$p.value),conf_int=if(is.null(h$conf.int))NULL else array_encode(as.numeric(h$conf.int)),estimate=if(is.null(h$estimate))NULL else array_encode(as.numeric(h$estimate)),method=as.character(h$method),alternative=if(is.null(h$alternative))NULL else as.character(h$alternative),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else if(op=='rng_arima_sim') {
  RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(req$seed);ord<-unlist(req$order);ar<-if(is.null(req$ar_file))numeric()else numeric_file(req$ar_file);ma<-if(is.null(req$ma_file))numeric()else numeric_file(req$ma_file);model<-list(ar=ar,ma=ma);list(value=array_encode(as.numeric(arima.sim(model=model,n=req$n,sd=req$sd))),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind(),seed=req$seed))
 } else if(op=='ts_arima') {
  x<-numeric_file(req$x_file);season<-if(is.null(req$seasonal_order))list(order=c(0L,0L,0L),period=NA_integer_)else list(order=unlist(req$seasonal_order),period=req$period);fit<-stats::arima(x,order=unlist(req$order),seasonal=season,method=req$method);cf<-coef(fit);se<-sqrt(diag(fit$var.coef));pred<-if(req$n_ahead>0)predict(fit,n.ahead=req$n_ahead)$pred else NULL;list(coefficients=as.list(names(cf)),estimate=array_encode(as.numeric(cf)),std_error=array_encode(as.numeric(se)),sigma2=encode(fit$sigma2),loglik=encode(as.numeric(logLik(fit))),aic=encode(fit$aic),residuals=if(isTRUE(req$residuals))array_encode(as.numeric(residuals(fit)))else NULL,forecast=if(is.null(pred))NULL else array_encode(as.numeric(pred)),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else if(op%in%c('ts_acf','ts_pacf')) {
  x<-numeric_file(req$x_file);lag<-if(is.null(req$lag_max))NULL else req$lag_max;z<-if(op=='ts_acf')acf(x,lag.max=lag,plot=FALSE)else pacf(x,lag.max=lag,plot=FALSE);v<-as.numeric(z$acf);list(lags=array_encode(as.numeric(z$lag)),values=array_encode(v),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else if(op=='ts_stl') {
  x<-numeric_file(req$x_file);fit<-stl(ts(x,frequency=req$period),s.window=if(isTRUE(req$periodic))'periodic'else req$seasonal_window,robust=isTRUE(req$robust));list(seasonal=array_encode(fit$time.series[,'seasonal']),trend=array_encode(fit$time.series[,'trend']),remainder=array_encode(fit$time.series[,'remainder']),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else if(op=='cluster_kmeans') {
  d<-as.data.frame(lapply(req$data_files,numeric_file),check.names=FALSE);RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(req$seed);fit<-kmeans(d,centers=req$k,nstart=req$nstart,algorithm='Hartigan-Wong');list(centers=lapply(seq_len(nrow(fit$centers)),function(i)array_encode(fit$centers[i,])),cluster=array_encode(as.integer(fit$cluster)),tot_withinss=encode(fit$tot.withinss),withinss=array_encode(fit$withinss),size=array_encode(fit$size),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind(),seed=req$seed))
 } else if(op=='cluster_hclust') {
  d<-as.data.frame(lapply(req$data_files,numeric_file),check.names=FALSE);dd<-dist(d,method=req$dist_method);fit<-hclust(dd,method=req$method);list(merge=lapply(seq_len(nrow(fit$merge)),function(i)as.list(as.integer(fit$merge[i,]))),height=array_encode(fit$height),order=array_encode(fit$order),cluster=if(is.null(req$cut_k))NULL else array_encode(cutree(fit,k=req$cut_k)),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else if(op=='fit_glm') {
  d<-as.data.frame(lapply(req$data_files,function(path)numeric_file(path)),check.names=FALSE);form<-as.formula(paste(req$response,'~',paste(unlist(req$predictors),collapse=' + ')),env=baseenv());fam<-switch(req$family,gaussian=gaussian(link=if(req$link=='default')'identity' else req$link),binomial=binomial(link=if(req$link=='default')'logit' else req$link),poisson=poisson(link=if(req$link=='default')'log' else req$link),Gamma=Gamma(link=if(req$link=='default')'inverse' else req$link));fit<-glm(form,data=d,family=fam,control=glm.control(maxit=req$maxit));sm<-summary(fit);cf<-coef(fit);nm<-names(cf);tab<-sm$coefficients;column<-function(j){v<-rep(NA_real_,length(cf));names(v)<-nm;common<-intersect(nm,rownames(tab));v[common]<-tab[common,j];v};list(coefficients=as.list(nm),estimate=array_encode(cf),std_error=array_encode(column(2)),statistic=array_encode(column(3)),p_value=array_encode(column(4)),statistic_type='z',deviance=encode(fit$deviance),null_deviance=encode(fit$null.deviance),aic=encode(fit$aic),iterations=fit$iter,converged=isTRUE(fit$converged),domain_status=if(isTRUE(fit$converged))'OK' else 'NON_CONVERGED',rank=fit$rank,meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else if(op=='smooth_spline') {
  x<-numeric_file(req$x_file);y<-numeric_file(req$y_file);args<-list(x=x,y=y,cv=isTRUE(req$cv));if(!is.null(req$df))args$df<-req$df;if(!is.null(req$spar))args$spar<-req$spar;fit<-do.call(smooth.spline,args);new<-if(is.null(req$new_x_file))NULL else predict(fit,x=numeric_file(req$new_x_file))$y;list(fitted=array_encode(fit$y),predicted=if(is.null(new))NULL else array_encode(new),lambda=encode(fit$lambda),df=encode(fit$df),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else if(op=='smooth_loess') {
  x<-numeric_file(req$x_file);y<-numeric_file(req$y_file);fit<-loess(y~x,data=data.frame(x=x,y=y),span=req$span,degree=req$degree,family=req$family);new<-if(is.null(req$new_x_file))NULL else predict(fit,newdata=data.frame(x=numeric_file(req$new_x_file)));list(fitted=array_encode(fitted(fit)),predicted=if(is.null(new))NULL else array_encode(as.numeric(new)),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else if(op=='interp_spline') {
  x<-numeric_file(req$x_file);y<-numeric_file(req$y_file);xo<-numeric_file(req$xout_file);v<-if(req$method=='monoH.FC')splinefun(x,y,method=req$method)(xo)else spline(x,y,xout=xo,method=req$method)$y;list(value=array_encode(v),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind(),method=req$method))
 } else if(op=='interp_approx') {
  z<-approx(numeric_file(req$x_file),numeric_file(req$y_file),xout=numeric_file(req$xout_file),method=req$method,rule=req$rule);list(value=array_encode(z$y),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind(),method=req$method,rule=req$rule))
 } else if(op=='fit_lm') {
  if(!is.null(req$file)) d<-read.csv(req$file,check.names=FALSE) else if(!is.null(req$data_files)) d<-as.data.frame(lapply(req$data_files,function(path) scan(path,what=double(),quiet=TRUE,na.strings='NA')),check.names=FALSE) else d<-as.data.frame(lapply(req$data,function(col) vapply(col,numeric_input,numeric(1))),check.names=FALSE)
  pieces<-strsplit(req$formula,'~',fixed=TRUE)[[1L]];if(length(pieces)!=2L) stop('bad formula');lhs<-trimws(pieces[1L]);rhs<-trimws(pieces[2L]);if(!lhs%in%names(d))stop('response must be a data column');tokens<-unlist(strsplit(gsub('[+-]',' ',rhs), '[[:space:]]+'));tokens<-tokens[nzchar(tokens)&!tokens%in%c('0','1')];if(any(!tokens%in%names(d))||'.'%in%tokens)stop('formula terms must be data columns');f<-as.formula(req$formula,env=baseenv());environment(f)<-baseenv();nafun<-if(identical(req$na_action,'exclude')) na.exclude else na.omit;fit<-lm(f,data=d,na.action=nafun);sm<-summary(fit);cf<-coef(fit);tab<-sm$coefficients;nm<-names(cf);col<-function(j){out<-rep(NA_real_,length(cf));names(out)<-nm;common<-intersect(nm,rownames(tab));out[common]<-tab[common,j];out};
  list(coefficients=as.list(nm),estimate=array_encode(cf),std_error=array_encode(col(2)),statistic=array_encode(col(3)),p_value=array_encode(col(4)),aliased=as.list(unname(nm[is.na(cf)])),rank=fit$rank,sigma=encode(sm$sigma),r_squared=encode(sm$r.squared),df=encode(sm$df[2]),residuals=if(isTRUE(req$residuals)) array_encode(residuals(fit)) else NULL,meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else stop('unknown operation')
},error=function(e) structure(list(message=conditionMessage(e)),class='rcli_error')),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
if(inherits(result,'rcli_error')) {writeLines(json(list(error=result$message,warnings=as.list(warnings))),outpath);quit(status=0)}
result$warnings<-as.list(warnings)
writeLines(json(result),outpath)
