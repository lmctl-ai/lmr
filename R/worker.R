`%||%` <- function(x,y) if(is.null(x)) y else x
library(stats)
library(utils)
parse_json <- function(s) {
 i<-1L; n<-nchar(s); peek<-function() if(i<=n) substr(s,i,i) else ''
 ws<-function(){while(i<=n&&peek()%in%c(' ','\n','\r','\t'))i<<-i+1L}
 string<-function(){i<<-i+1L; z<-'';while(i<=n){c<-peek();i<<-i+1L;if(c=='"')return(z);if(c=='\\'){e<-peek();i<<-i+1L;z<-paste0(z,switch(e,'"'='"','\\'='\\','/'='/',b='\b',f='\f',n='\n',r='\r',t='\t',e))}else z<-paste0(z,c)};stop('unterminated JSON string')}
 value<-function(){ws();c<-peek();if(c=='"')return(string());if(c=='{'){i<<-i+1L;obj<-list();ws();if(peek()=='}'){i<<-i+1L;return(obj)};repeat{ws();k<-string();ws();i<<-i+1L;obj[k]<-list(value());ws();c<-peek();i<<-i+1L;if(c=='}')break};return(obj)};if(c=='['){i<<-i+1L;arr<-list();ws();if(peek()==']'){i<<-i+1L;return(arr)};repeat{arr[length(arr)+1L]<-list(value());ws();c<-peek();i<<-i+1L;if(c==']')break};return(arr)};tail<-substr(s,i,n);if(startsWith(tail,'true')){i<<-i+4L;return(TRUE)};if(startsWith(tail,'false')){i<<-i+5L;return(FALSE)};if(startsWith(tail,'null')){i<<-i+4L;return(NULL)};m<-regexpr('^-?(0|[1-9][0-9]*)(\\.[0-9]+)?([eE][+-]?[0-9]+)?',tail,perl=TRUE);if(m[1]<0)stop('invalid JSON value');v<-regmatches(tail,m);i<<-i+nchar(v);as.numeric(v)}
 value()
}
args <- commandArgs(TRUE); if(length(args)&&args[[1]]=='--') args<-args[-1]; op <- args[[1]]; req <- parse_json(paste(readLines(args[[2]],warn=FALSE),collapse='')); outpath <- args[[3]]
warnings <- character()
encode <- function(x) {
 if (is.double(x)) { x[is.na(x) & !is.nan(x)] <- NA_real_; return(lapply(seq_along(x), function(i) {v<-x[[i]]; if(is.nan(v)) 'NaN' else if(is.na(v)) 'NA' else if(is.infinite(v)) if(v>0) 'Inf' else '-Inf' else v})) }
 if (is.list(x)) return(lapply(x, encode))
 x
}
json <- function(x) {
 if(is.null(x)) return('null')
 if(is.list(x)) {
  n<-names(x)
  if(!is.null(n) && all(n!='')) return(paste0('{',paste(vapply(seq_along(x),function(i) paste0(json(n[[i]]),':',json(x[[i]])),''),collapse=','),'}'))
  return(paste0('[',paste(vapply(x,json,''),collapse=','),']'))
 }
 if(length(x)==0) return('[]')
 if(length(x)>1) return(paste0('[',paste(vapply(as.list(x),json,''),collapse=','),']'))
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
  p <- unlist(req$p); args2 <- list(p=p,lower.tail=if(is.null(req$lower_tail)) TRUE else req$lower_tail,log.p=FALSE)
  if(req$dist=='norm') {args2$mean<-req$mean;args2$sd<-req$sd}
  if(req$dist%in%c('t','chisq')) args2$df<-req$df
  if(req$dist=='gamma') {args2$shape<-req$shape;if(!is.null(req$rate))args2$rate<-req$rate}
  if(req$dist=='beta') {args2$shape1<-req$a;args2$shape2<-req$b}
  if(req$dist=='binom') {args2$size<-req$size;args2$prob<-req$prob}
  if(req$dist=='pois') args2$lambda<-req$lambda
  list(value=encode(do.call(fn,args2)),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else if(op=='rng_draw') {
  RNGkind('Mersenne-Twister','Inversion','Rejection'); set.seed(req$seed)
  n<-req$n; value<-switch(req$dist,norm=rnorm(n,req$mean%||%0,req$sd%||%1),unif=runif(n,req$min%||%0,req$max%||%1))
  list(value=encode(value),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind(),seed=req$seed))
 } else if(op=='fit_lm') {
  if(!is.null(req$file)) d<-read.csv(req$file,check.names=FALSE) else d<-as.data.frame(lapply(req$data,function(col) vapply(col,function(v) if(is.null(v)) NA_real_ else as.numeric(v),numeric(1))),check.names=FALSE)
  if(!grepl('^[A-Za-z.][A-Za-z0-9._]* *~ *[A-Za-z.][A-Za-z0-9._]*( *\\+ *[A-Za-z.][A-Za-z0-9._]*)*$',req$formula)) stop('formula must contain column names and + only')
  f<-as.formula(req$formula); nafun<-if(identical(req$na_action,'exclude')) na.exclude else na.omit; fit<-lm(f,data=d,na.action=nafun); sm<-summary(fit); co<-sm$coefficients
  list(coefficients=rownames(co),estimate=encode(co[,1]),std_error=encode(co[,2]),statistic=encode(co[,3]),p_value=encode(co[,4]),sigma=encode(sm$sigma),r_squared=encode(sm$r.squared),df=sm$df[2],residuals=if(isTRUE(req$residuals)) encode(residuals(fit)) else NULL,meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else stop('unknown operation')
},error=function(e) structure(list(message=conditionMessage(e)),class='rcli_error')),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
if(inherits(result,'rcli_error')) {writeLines(json(list(error=result$message)),outpath);quit(status=0)}
result$warnings<-as.list(warnings)
writeLines(json(result),outpath)
