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
 } else if(op=='fit_lm') {
  if(!is.null(req$file)) d<-read.csv(req$file,check.names=FALSE) else if(!is.null(req$data_files)) d<-as.data.frame(lapply(req$data_files,function(path) scan(path,what=double(),quiet=TRUE,na.strings='NA')),check.names=FALSE) else d<-as.data.frame(lapply(req$data,function(col) vapply(col,numeric_input,numeric(1))),check.names=FALSE)
  pieces<-strsplit(req$formula,'~',fixed=TRUE)[[1L]];if(length(pieces)!=2L) stop('bad formula');lhs<-trimws(pieces[1L]);rhs<-trimws(pieces[2L]);if(!lhs%in%names(d))stop('response must be a data column');tokens<-unlist(strsplit(gsub('[+-]',' ',rhs), '[[:space:]]+'));tokens<-tokens[nzchar(tokens)&!tokens%in%c('0','1')];if(any(!tokens%in%names(d))||'.'%in%tokens)stop('formula terms must be data columns');f<-as.formula(req$formula,env=baseenv());environment(f)<-baseenv();nafun<-if(identical(req$na_action,'exclude')) na.exclude else na.omit;fit<-lm(f,data=d,na.action=nafun);sm<-summary(fit);cf<-coef(fit);tab<-sm$coefficients;nm<-names(cf);col<-function(j){out<-rep(NA_real_,length(cf));names(out)<-nm;common<-intersect(nm,rownames(tab));out[common]<-tab[common,j];out};
  list(coefficients=as.list(nm),estimate=array_encode(cf),std_error=array_encode(col(2)),statistic=array_encode(col(3)),p_value=array_encode(col(4)),aliased=as.list(unname(nm[is.na(cf)])),rank=fit$rank,sigma=encode(sm$sigma),r_squared=encode(sm$r.squared),df=encode(sm$df[2]),residuals=if(isTRUE(req$residuals)) array_encode(residuals(fit)) else NULL,meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else stop('unknown operation')
},error=function(e) structure(list(message=conditionMessage(e)),class='rcli_error')),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
if(inherits(result,'rcli_error')) {writeLines(json(list(error=result$message,warnings=as.list(warnings))),outpath);quit(status=0)}
result$warnings<-as.list(warnings)
writeLines(json(result),outpath)
