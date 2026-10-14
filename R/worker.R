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
 } else if(op=='rng_draw') {
  RNGkind('Mersenne-Twister','Inversion','Rejection'); set.seed(req$seed)
  n<-req$n; value<-switch(req$dist,norm=rnorm(n,if('mean'%in%names(req))numeric_input(req$mean) else 0,if('sd'%in%names(req))numeric_input(req$sd) else 1),unif=runif(n,if('min'%in%names(req))numeric_input(req$min) else 0,if('max'%in%names(req))numeric_input(req$max) else 1))
  list(value=array_encode(value),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind(),seed=req$seed))
 } else if(op=='fit_lm') {
  if(!is.null(req$file)) d<-read.csv(req$file,check.names=FALSE) else if(!is.null(req$data_files)) d<-as.data.frame(lapply(req$data_files,function(path) scan(path,what=double(),quiet=TRUE,na.strings='NA')),check.names=FALSE) else d<-as.data.frame(lapply(req$data,function(col) vapply(col,numeric_input,numeric(1))),check.names=FALSE)
  pieces<-strsplit(req$formula,'~',fixed=TRUE)[[1L]];if(length(pieces)!=2L) stop('bad formula');lhs<-trimws(pieces[1L]);rhs<-trimws(pieces[2L]);if(!lhs%in%names(d))stop('response must be a data column');tokens<-unlist(strsplit(gsub('[+-]',' ',rhs), '[[:space:]]+'));tokens<-tokens[nzchar(tokens)&!tokens%in%c('0','1')];if(any(!tokens%in%names(d))||'.'%in%tokens)stop('formula terms must be data columns');f<-as.formula(req$formula,env=baseenv());environment(f)<-baseenv();nafun<-if(identical(req$na_action,'exclude')) na.exclude else na.omit;fit<-lm(f,data=d,na.action=nafun);sm<-summary(fit);cf<-coef(fit);tab<-sm$coefficients;nm<-names(cf);col<-function(j){out<-rep(NA_real_,length(cf));names(out)<-nm;common<-intersect(nm,rownames(tab));out[common]<-tab[common,j];out};
  list(coefficients=as.list(nm),estimate=array_encode(cf),std_error=array_encode(col(2)),statistic=array_encode(col(3)),p_value=array_encode(col(4)),aliased=as.list(unname(nm[is.na(cf)])),rank=fit$rank,sigma=encode(sm$sigma),r_squared=encode(sm$r.squared),df=encode(sm$df[2]),residuals=if(isTRUE(req$residuals)) array_encode(residuals(fit)) else NULL,meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else stop('unknown operation')
},error=function(e) structure(list(message=conditionMessage(e)),class='rcli_error')),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
if(inherits(result,'rcli_error')) {writeLines(json(list(error=result$message,warnings=as.list(warnings))),outpath);quit(status=0)}
result$warnings<-as.list(warnings)
writeLines(json(result),outpath)
