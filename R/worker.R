`%||%` <- function(x,y) if(is.null(x)) y else x
library(stats)
library(utils)
@RCLI_EXPR_PRELUDE@
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
args <- commandArgs(TRUE); if(length(args)&&args[[1]]=='--') args<-args[-1]; if(length(args)==1L&&args[[1]]=='--expr-spec-hash'){cat(EXPR_SPEC_SHA256,'\n',sep='');quit(status=0)}
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
json_quote <- function(s) {
 cp<-utf8ToInt(enc2utf8(s));escaped<-vapply(cp,function(i){if(i==34L)'\\"'else if(i==92L)'\\\\'else if(i==8L)'\\b'else if(i==9L)'\\t'else if(i==10L)'\\n'else if(i==12L)'\\f'else if(i==13L)'\\r'else if(i<32L)sprintf('\\u%04x',i)else intToUtf8(i)},'');paste0('"',paste(escaped,collapse=''),'"')
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
 if(is.character(x)) return(json_quote(x))
 if(is.logical(x)) return(if(is.na(x)) '"NA"' else if(x) 'true' else 'false')
 if(is.numeric(x)) { if(is.nan(x)) return('"NaN"'); if(is.na(x)) return('"NA"'); if(is.infinite(x)) return(if(x>0) '"Inf"' else '"-Inf"'); return(sprintf('%.17g',x)) }
 'null'
}
expr_condition <- function(code, message=code) stop(structure(list(message=message,tool_code=code),class=c('rcli_expr_error','error','condition')))
expr_keys <- function(x, expected) is.list(x) && !is.null(base::names(x)) && setequal(base::names(x),expected)
build_expr_ast <- function(root, declared, allow_reducers=FALSE, counters=NULL, state=NULL, depth=0L, in_condition=FALSE) {
 if (is.null(state)) { state <- new.env(parent=emptyenv());state$nodes<-0L;state$calls<-0L }
 state$nodes <- state$nodes+1L
 # C++ enforces source nesting depth; a left-associated flat operator chain may
 # have a deeper AST spine, bounded here by the shared total-node cap.
 if (state$nodes>EXPR_CAPS[['nodes']] || depth>EXPR_CAPS[['nodes']]) expr_condition('EXPR_SPEC_MISMATCH')
 if (!is.list(root) || is.null(root$t) || !is.character(root$t) || length(root$t)!=1L) expr_condition('EXPR_SPEC_MISMATCH')
 tick_fun <- function(fun) {
  if (is.null(counters)) return(fun)
  force(fun);force(counters)
  function(...) { counters$evals<-counters$evals+1L;if(counters$evals>counters$limit)expr_condition('EVAL_LIMIT','expression evaluation limit exceeded');fun(...) }
 }
 kind<-root$t
 if(in_condition&&kind!='op')expr_condition('EXPR_SPEC_MISMATCH')
 if(kind=='num'){if(!expr_keys(root,c('t','v'))||!is.numeric(root$v)||length(root$v)!=1L||!is.finite(root$v))expr_condition('EXPR_SPEC_MISMATCH');return(as.double(root$v))}
 if(kind=='const'){if(!expr_keys(root,c('t','n'))||!root$n%in%base::names(EXPR_CONST))expr_condition('EXPR_SPEC_MISMATCH');return(unname(EXPR_CONST[[root$n]]))}
 if(kind=='var'){if(!expr_keys(root,c('t','n'))||!is.character(root$n)||length(root$n)!=1L||!grepl('^[a-z][a-z0-9_]{0,31}$',root$n)||root$n%in%c(EXPR_FUNCTIONS,base::names(EXPR_CONST),EXPR_RESERVED)||!root$n%in%declared)expr_condition('EXPR_SPEC_MISMATCH');return(as.name(root$n))}
 if(kind=='neg'){if(!expr_keys(root,c('t','a'))||length(root$a)!=1L)expr_condition('EXPR_SPEC_MISMATCH');return(as.call(list(tick_fun(FN[['-']]),build_expr_ast(root$a[[1L]],declared,allow_reducers,counters,state,depth+1L,FALSE))))}
 if(kind=='op'){
  if(!expr_keys(root,c('t','o','a'))||length(root$a)!=2L||!root$o%in%EXPR_OPERATORS)expr_condition('EXPR_SPEC_MISMATCH')
  rel<-root$o%in%c('<','<=','>','>=','==','!=','&','|');if(rel!=in_condition)expr_condition('EXPR_SPEC_MISMATCH')
  child_condition<-root$o%in%c('&','|');return(as.call(list(tick_fun(FN[[root$o]]),build_expr_ast(root$a[[1L]],declared,allow_reducers,counters,state,depth+1L,child_condition),build_expr_ast(root$a[[2L]],declared,allow_reducers,counters,state,depth+1L,child_condition))))
 }
 if(kind=='call'){
  if(!expr_keys(root,c('t','f','a'))||!is.character(root$f)||length(root$f)!=1L||!root$f%in%EXPR_FUNCTIONS||root$f=='ifelse'||length(root$a)!=ARITY[[root$f]])expr_condition('EXPR_SPEC_MISMATCH')
  fn<-root$f;if(fn%in%c('sum','mean','prod','min','max')&&!allow_reducers)expr_condition('EXPR_SPEC_MISMATCH')
  state$calls<-state$calls+1L;if(state$calls>EXPR_CAPS[['calls']])expr_condition('EXPR_SPEC_MISMATCH');args<-lapply(root$a,function(a)build_expr_ast(a,declared,allow_reducers,counters,state,depth+1L,FALSE));state$calls<-state$calls-1L
  return(as.call(c(list(tick_fun(FN[[fn]])),args)))
 }
 if(kind=='ifelse'){
  if(!expr_keys(root,c('t','c','a','b')))expr_condition('EXPR_SPEC_MISMATCH')
  state$calls<-state$calls+1L;if(state$calls>EXPR_CAPS[['calls']])expr_condition('EXPR_SPEC_MISMATCH');args<-list(build_expr_ast(root$c,declared,allow_reducers,counters,state,depth+1L,TRUE),build_expr_ast(root$a,declared,allow_reducers,counters,state,depth+1L,FALSE),build_expr_ast(root$b,declared,allow_reducers,counters,state,depth+1L,FALSE));state$calls<-state$calls-1L
  return(as.call(c(list(tick_fun(FN[['ifelse']])),args)))
 }
 expr_condition('EXPR_SPEC_MISMATCH')
}
expr_compile<-function(ast,declared,allow_reducers=FALSE)build_expr_ast(ast,declared,allow_reducers)
expr_counter<-function(limit){z<-new.env(parent=emptyenv());z$evals<-0L;z$limit<-as.integer(limit);z}
expr_runner<-function(ast,declared,allow_reducers,counter){body<-expr_compile(ast,declared,allow_reducers);function(values){counter$evals<-counter$evals+1L;if(counter$evals>counter$limit)expr_condition('EVAL_LIMIT','expression evaluation limit exceeded');env<-list2env(values,parent=emptyenv());eval(body,envir=env)}}
expr_data<-function(req){if(!is.null(req$data_files))return(lapply(req$data_files,numeric_file));if(is.null(req$data))return(list());lapply(req$data,function(x)vapply(x,numeric_input,numeric(1)))}
op <- args[[1]]; req <- parse_json(paste(readLines(args[[2]],warn=FALSE),collapse='')); outpath <- args[[3]]
result <- withCallingHandlers(tryCatch({
 if(op=='expr_builder_check'){
  tryCatch({build_expr_ast(req$ast,unlist(req$names),isTRUE(req$allow_reducers));list(accepted=TRUE)},error=function(e)list(accepted=FALSE,tool_code=if(is.null(e$tool_code))'EXPR_SPEC_MISMATCH' else e$tool_code))
 }
 else if(op=='expr_builder_batch') list(decisions=lapply(req$cases,function(z)tryCatch({build_expr_ast(z$ast,unlist(z$names),isTRUE(z$allow_reducers));list(accepted=TRUE)},error=function(e)list(accepted=FALSE,tool_code=if(is.null(e$tool_code))'EXPR_SPEC_MISMATCH' else e$tool_code))))
 else if(op=='r_info') list(r_version=as.character(getRversion()),rng_kind=RNGkind(),expr_spec_sha256=EXPR_SPEC_SHA256)
 else if(op=='opt_minimize') {
  setTimeLimit(elapsed=45,transient=TRUE);counter<-expr_counter(req$max_evals);dat<-expr_data(req);pnames<-names(req$start);start<-unlist(req$start,use.names=TRUE);method<-req$method
  runner<-expr_runner(req$objective_ast,base::names(c(as.list(start),dat)),TRUE,counter)
  objective<-function(par){vals<-c(as.list(setNames(as.numeric(par),pnames)),dat);v<-runner(vals);if(length(v)!=1L||!is.numeric(v))expr_condition('EXPR_NOT_SCALAR','objective must return one numeric value');v<-as.double(v);if(!is.finite(v)){counter$nonfinite<-if(is.null(counter$nonfinite))1L else counter$nonfinite+1L;if(req$nonfinite=='inf'&&method%in%c('Nelder-Mead','SANN'))return(Inf);expr_condition('NONFINITE_OBJECTIVE','objective returned a non-finite value')};v}
  initial<-tryCatch(objective(start),error=function(e){if(identical(e$tool_code,'NONFINITE_OBJECTIVE'))return(structure(list(),class='rcli_nonfinite_start'));stop(e)});initial_bad<-inherits(initial,'rcli_nonfinite_start')||!is.finite(initial);lower<-rep(-Inf,length(start));upper<-rep(Inf,length(start));names(lower)<-names(upper)<-pnames;if(!is.null(req$lower))lower[names(req$lower)]<-unlist(req$lower);if(!is.null(req$upper))upper[names(req$upper)]<-unlist(req$upper)
  if(initial_bad)list(par=as.list(setNames(as.double(start),pnames)),value='NaN',counts=list(`function`=0L,gradient=NULL),convergence=1L,message='objective is non-finite at starting parameters',domain_status='NONFINITE_OBJECTIVE',hessian=NULL,evals=counter$evals,nonfinite_evals=counter$nonfinite,elapsed_ms=0L,meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind(),seed=if(method=='SANN')req$seed else NULL))else{
   if(method=='SANN'){RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(req$seed)}
   control<-list(maxit=req$maxit);if(method!='L-BFGS-B')control$reltol<-req$reltol
   began<-proc.time()[['elapsed']];fit<-stats::optim(start,objective,method=method,lower=lower,upper=upper,control=control,hessian=isTRUE(req$hessian));elapsed<-round((proc.time()[['elapsed']]-began)*1000);counts<-fit$counts;if(is.null(counts))counts<-integer();conv<-fit$convergence;status<-if(conv==0)'OK'else if(conv==1)'MAXIT'else'NON_CONVERGED';hess<-if(isTRUE(req$hessian)&&!is.null(fit$hessian))lapply(seq_len(nrow(fit$hessian)),function(i)array_encode(fit$hessian[i,]))else NULL
   list(par=as.list(setNames(as.double(fit$par),names(fit$par))),value=encode(fit$value),counts=list(`function`=if(length(counts))as.integer(counts[[1]])else NULL,gradient=if(length(counts)>1L)as.integer(counts[[2]])else NULL),convergence=as.integer(conv),message=if(is.null(fit$message))NULL else as.character(fit$message),domain_status=status,hessian=hess,evals=counter$evals,nonfinite_evals=if(is.null(counter$nonfinite))0L else counter$nonfinite,elapsed_ms=elapsed,meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind(),seed=if(method=='SANN')req$seed else NULL))
  }
 } else if(op=='opt_optimize') {
  setTimeLimit(elapsed=45,transient=TRUE);counter<-expr_counter(req$max_evals);runner<-expr_runner(req$ast,'x',FALSE,counter);f<-function(x){v<-runner(list(x=as.double(x)));if(length(v)!=1L||!is.numeric(v))expr_condition('EXPR_NOT_SCALAR','objective must return one numeric value');as.double(v)};fit<-stats::optimize(f,c(req$lower,req$upper),maximum=isTRUE(req$maximum),tol=req$tol);list(x=encode(if(isTRUE(req$maximum))fit$maximum else fit$minimum),value=encode(fit$objective),evals=counter$evals,meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else if(op=='num_root') {
  setTimeLimit(elapsed=45,transient=TRUE);counter<-expr_counter(req$max_evals);runner<-expr_runner(req$ast,'x',FALSE,counter);f<-function(x){v<-runner(list(x=as.double(x)));if(length(v)==1L&&length(x)>1L)v<-rep(v,length(x));if(length(v)!=length(x)||!is.numeric(v))expr_condition('EXPR_NOT_SCALAR','root function must return one value per x');as.double(v)};extend<-req$extend_int;fl<-f(req$lower);fu<-f(req$upper);if(extend=='no'&&(length(fl)!=1L||length(fu)!=1L||!is.finite(fl)||!is.finite(fu)||fl*fu>0))expr_condition('BAD_INTERVAL','root interval endpoints must be finite and bracket a root');fit<-tryCatch(stats::uniroot(f,c(req$lower,req$upper),tol=req$tol,maxiter=req$maxiter,extendInt=extend),error=function(e){if(!is.null(e$tool_code))stop(e);expr_condition('BAD_INTERVAL',conditionMessage(e))});list(root=encode(fit$root),f_root=encode(fit$f.root),iter=as.integer(fit$iter),estim_prec=encode(fit$estim.prec),evals=counter$evals,meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else if(op=='num_integrate') {
  setTimeLimit(elapsed=45,transient=TRUE);counter<-expr_counter(req$max_evals);runner<-expr_runner(req$ast,'x',FALSE,counter);f<-function(x){v<-runner(list(x=as.double(x)));if(length(v)==1L&&length(x)>1L)v<-rep(v,length(x));if(length(v)!=length(x)||!is.numeric(v))expr_condition('EXPR_NOT_SCALAR','integrand must return one value per x');as.double(v)};lower<-numeric_input(req$lower);upper<-numeric_input(req$upper);fit<-tryCatch(stats::integrate(f,lower=lower,upper=upper,subdivisions=req$subdivisions,rel.tol=req$rel_tol,abs.tol=req$abs_tol,stop.on.error=FALSE),error=function(e){if(grepl('non-finite function value',conditionMessage(e),fixed=TRUE))return(structure(list(message=conditionMessage(e)),class='rcli_nonfinite_integral'));stop(e)});if(!inherits(fit,'rcli_nonfinite_integral')&&fit$message!='OK')warning(fit$message,call.=FALSE);if(inherits(fit,'rcli_nonfinite_integral'))list(value='NaN',abs_error='NaN',subdivisions=0L,message='non-finite function value',domain_status='NON_FINITE',evals=counter$evals,meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))else list(value=encode(fit$value),abs_error=encode(fit$abs.error),subdivisions=as.integer(fit$subdivisions),message=fit$message,domain_status=if(fit$message=='OK')'OK'else'NON_CONVERGED',evals=counter$evals,meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else if(op=='fit_nls') {
  setTimeLimit(elapsed=45,transient=TRUE);counter<-expr_counter(req$max_evals);dat<-expr_data(req);sn<-names(req$start);starts<-unlist(req$start,use.names=TRUE);declared<-base::names(c(as.list(starts),dat));state<-new.env(parent=emptyenv());state$nodes<-0L;state$calls<-0L;body<-build_expr_ast(req$rhs_ast,declared,FALSE,NULL,state);tick<-local({cc<-counter;function(expr){cc$evals<-cc$evals+1L;if(cc$evals>cc$limit)expr_condition('EVAL_LIMIT','expression evaluation limit exceeded');eval(substitute(expr),envir=parent.frame())}});body<-as.call(list(tick,body));safe_parent<-list2env(list(list=base::list),parent=emptyenv());values<-as.list(starts);values_env<-list2env(values,parent=safe_parent);form<-structure(as.call(list(as.name('~'),as.name(req$response),body)),class='formula',.Environment=values_env);d<-as.data.frame(dat,check.names=FALSE);control<-stats::nls.control(maxiter=req$maxiter,tol=req$tol,warnOnly=TRUE)
  fit<-tryCatch(if(req$algorithm=='port')stats::nls(form,data=d,start=starts,algorithm='port',lower={z<-rep(-Inf,length(starts));names(z)<-sn;if(!is.null(req$lower))z[names(req$lower)]<-unlist(req$lower);z},upper={z<-rep(Inf,length(starts));names(z)<-sn;if(!is.null(req$upper))z[names(req$upper)]<-unlist(req$upper);z},control=control)else stats::nls(form,data=d,start=starts,algorithm='default',control=control),error=function(e){if(grepl('singular gradient|step factor|convergence failure|iteration.*exceed',conditionMessage(e),ignore.case=TRUE))return(structure(list(message=conditionMessage(e),iterations=0L,stopCode=if(grepl('iteration.*exceed',conditionMessage(e),ignore.case=TRUE))3L else 1L),class='rcli_nls_nonconvergence'));stop(e)})
  nonconv<-inherits(fit,'rcli_nls_nonconvergence')||!isTRUE(fit$convInfo$isConv);if(nonconv){msg<-if(inherits(fit,'rcli_nls_nonconvergence'))fit$message else fit$convInfo$stopMessage;iters<-if(inherits(fit,'rcli_nls_nonconvergence'))fit$iterations else as.integer(fit$convInfo$finIter);code<-if(grepl('iteration.*exceed',msg,ignore.case=TRUE)||(!inherits(fit,'rcli_nls_nonconvergence')&&fit$convInfo$stopCode==3L))'MAXIT'else'NON_CONVERGED';list(coefficients=as.list(sn),estimate=NULL,std_error=NULL,statistic=NULL,p_value=NULL,sigma='NaN',df=NULL,converged=FALSE,iterations=iters,achieved_tol='NaN',message=msg,domain_status=code,evals=counter$evals,meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))}else{sm<-summary(fit);tab<-sm$coefficients;cf<-coef(fit);list(coefficients=as.list(names(cf)),estimate=array_encode(as.numeric(cf)),std_error=array_encode(tab[,2]),statistic=array_encode(tab[,3]),p_value=array_encode(tab[,4]),sigma=encode(sm$sigma),df=encode(df.residual(fit)),converged=TRUE,iterations=as.integer(fit$convInfo$finIter),achieved_tol=encode(fit$convInfo$finTol),message=NULL,domain_status='OK',evals=counter$evals,fitted=if(isTRUE(req$fitted))array_encode(as.numeric(fitted(fit)))else NULL,residuals=if(isTRUE(req$residuals))array_encode(as.numeric(residuals(fit)))else NULL,meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))}
 }
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
  RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(req$seed);x<-numeric_file(req$x_file);w<-if(is.null(req$prob_file))NULL else numeric_file(req$prob_file);ix<-sample.int(length(x),size=req$size,replace=isTRUE(req$replace),prob=w);list(value=array_encode(x[ix]),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind(),seed=req$seed))
 } else if(op=='math_special') {
  x<-numeric_file(req$x_file);fn<-switch(req$fn,gamma=gamma,lgamma=lgamma,beta=beta,lbeta=lbeta,digamma=digamma,trigamma=trigamma,choose=choose,lchoose=lchoose,besselJ=besselJ,besselI=besselI,besselK=besselK,besselY=besselY)
  value<-if(req$fn%in%c('beta','lbeta','choose','lchoose'))do.call(fn,list(x,numeric_file(req$y_file))) else if(grepl('^bessel',req$fn))do.call(fn,list(x,numeric_input(req$nu))) else fn(x)
  list(value=array_encode(value),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else if(op=='summary_quantile') {
  v<-quantile(numeric_file(req$x_file),probs=numeric_file(req$probs_file),na.rm=isTRUE(req$na_rm),names=FALSE,type=req$type);list(value=array_encode(as.numeric(v)),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind(),type=req$type))
 } else if(op=='num_fft') {
  x<-numeric_file(req$x_file);z<-fft(x,inverse=isTRUE(req$inverse));if(isTRUE(req$normalize)&&isTRUE(req$inverse))z<-z/length(x);list(re=array_encode(Re(z)),im=array_encode(Im(z)),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind(),inverse=isTRUE(req$inverse)))
 } else if(op=='test_htest') {
  x<-numeric_file(req$x_file);y<-if(is.null(req$y_file))NULL else numeric_file(req$y_file);alt<-req$alternative;cl<-req$conf_level;tst<-req$test
  h<-switch(tst,
   t=if(is.null(y))t.test(x,mu=numeric_input(req$mu),alternative=alt,conf.level=cl) else t.test(x,y,mu=numeric_input(req$mu),alternative=alt,paired=isTRUE(req$paired),var.equal=isTRUE(req$var_equal),conf.level=cl),
   wilcox=wilcox.test(x,y,mu=numeric_input(req$mu),alternative=alt,paired=isTRUE(req$paired),conf.int=isTRUE(req$conf_int),conf.level=cl,exact=req$exact,correct=isTRUE(req$correct)),
   chisq={if(!is.null(req$seed)){RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(req$seed)};xx<-if(!is.null(req$rows))matrix(x,nrow=req$rows,ncol=req$cols,byrow=isTRUE(req$byrow)) else x;args<-list(x=xx,correct=isTRUE(req$correct),simulate.p.value=isTRUE(req$simulate_p_value));if(!is.null(req$expected_file))args$p<-numeric_file(req$expected_file);do.call(chisq.test,args)},
   shapiro=shapiro.test(x),
   ks={if(!is.null(req$y_file))ks.test(x,y,alternative=alt,exact=req$exact) else {fn<-switch(req$distribution,norm=pnorm,unif=punif,exp=pexp,pois=ppois);args<-switch(req$distribution,norm=list(mean=numeric_input(req$mean),sd=numeric_input(req$sd)),unif=list(min=numeric_input(req$min),max=numeric_input(req$max)),exp=list(rate=numeric_input(req$rate)),pois=list(lambda=numeric_input(req$lambda)));do.call(ks.test,c(list(x=x,y=fn,alternative=alt,exact=req$exact),args))}},
   fisher={xx<-matrix(x,nrow=if(is.null(req$rows))2 else req$rows,ncol=if(is.null(req$cols))2 else req$cols,byrow=isTRUE(req$byrow));if(nrow(xx)==2&&ncol(xx)==2)fisher.test(xx,alternative=alt,conf.level=cl,or=if(is.null(req$odds_ratio))1 else req$odds_ratio)else fisher.test(xx,alternative=alt,conf.level=cl)},
   binom=binom.test(as.integer(x[[1L]]),req$n,p=if(is.null(req[['p']]))0.5 else numeric_input(req[['p']]),alternative=alt,conf.level=cl),
   prop={args<-list(x=as.integer(x),n=as.integer(numeric_file(req$trials_file)),alternative=alt,correct=isTRUE(req$correct),conf.level=cl);if(!is.null(req$probabilities_file))args$p<-numeric_file(req$probabilities_file) else if(!is.null(req[['p']]))args$p<-numeric_input(req[['p']]);do.call(prop.test,args)},
   cor=cor.test(x,y,alternative=alt,method=req$method,conf.level=cl),
   var=var.test(x,y,ratio=numeric_input(req$ratio),alternative=alt,conf.level=cl))
  list(statistic=if(length(h$statistic))as.numeric(h$statistic)[[1L]]else NULL,statistic_name=if(length(h$statistic))names(h$statistic)[[1L]]else NULL,parameter=if(!length(h$parameter))NULL else if(length(h$parameter)==1L)as.numeric(h$parameter)[[1L]]else array_encode(as.numeric(h$parameter)),p_value=as.numeric(h$p.value),conf_int=if(is.null(h$conf.int))NULL else array_encode(as.numeric(h$conf.int)),estimate=if(is.null(h$estimate))NULL else array_encode(as.numeric(h$estimate)),estimate_names=if(is.null(h$estimate))list()else as.list(names(h$estimate)),null_value=if(is.null(h$null.value))NULL else encode(h$null.value),method=as.character(h$method),alternative=if(is.null(h$alternative))NULL else as.character(h$alternative),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else if(op=='rng_arima_sim') {
  RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(req$seed);ord<-unlist(req$order);ar<-if(is.null(req$ar_file))numeric()else numeric_file(req$ar_file);ma<-if(is.null(req$ma_file))numeric()else numeric_file(req$ma_file);model<-list(ar=ar,ma=ma);list(value=array_encode(as.numeric(arima.sim(model=model,n=req$n,sd=req$sd))),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind(),seed=req$seed))
 } else if(op=='ts_arima') {
  x<-numeric_file(req$x_file);season<-if(is.null(req$seasonal_order))list(order=c(0L,0L,0L),period=NA_integer_)else list(order=unlist(req$seasonal_order),period=req$period);fit<-stats::arima(x,order=unlist(req$order),seasonal=season,method=req$method,include.mean=isTRUE(req$include_mean));cf<-coef(fit);se<-sqrt(diag(fit$var.coef));pr<-if(req$n_ahead>0)predict(fit,n.ahead=req$n_ahead)else NULL;list(coefficients=as.list(names(cf)),estimate=array_encode(as.numeric(cf)),std_error=array_encode(as.numeric(se)),sigma2=encode(fit$sigma2),loglik=encode(as.numeric(logLik(fit))),aic=encode(fit$aic),converged=is.null(fit$code)||fit$code==0,residuals=if(isTRUE(req$residuals))array_encode(as.numeric(residuals(fit)))else NULL,forecast=if(is.null(pr))NULL else array_encode(as.numeric(pr$pred)),forecast_se=if(is.null(pr))NULL else array_encode(as.numeric(pr$se)),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else if(op%in%c('ts_acf','ts_pacf')) {
  x<-numeric_file(req$x_file);lag<-if(is.null(req$lag_max))NULL else req$lag_max;z<-if(op=='ts_acf')acf(x,lag.max=lag,plot=FALSE)else pacf(x,lag.max=lag,plot=FALSE);v<-as.numeric(z$acf);list(lags=array_encode(as.numeric(z$lag)),value=array_encode(v),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else if(op=='ts_stl') {
  x<-numeric_file(req$x_file);fit<-stl(ts(x,frequency=req$period),s.window=if(isTRUE(req$periodic))'periodic'else req$seasonal_window,robust=isTRUE(req$robust));list(seasonal=array_encode(fit$time.series[,'seasonal']),trend=array_encode(fit$time.series[,'trend']),remainder=array_encode(fit$time.series[,'remainder']),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else if(op=='cluster_kmeans') {
  d<-as.data.frame(lapply(req$data_files,numeric_file),check.names=FALSE);RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(req$seed);fit<-kmeans(d,centers=req$k,nstart=req$nstart,iter.max=req$iter_max,algorithm=req$algorithm);list(centers=lapply(seq_len(nrow(fit$centers)),function(i)array_encode(fit$centers[i,])),cluster=array_encode(as.integer(fit$cluster)),tot_withinss=encode(fit$tot.withinss),withinss=array_encode(fit$withinss),size=array_encode(fit$size),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind(),seed=req$seed))
 } else if(op=='cluster_hclust') {
  d<-as.data.frame(lapply(req$data_files,numeric_file),check.names=FALSE);dd<-dist(d,method=req$dist_method,p=req[['p']]);fit<-hclust(dd,method=req$method);list(merge=lapply(seq_len(nrow(fit$merge)),function(i)as.list(as.integer(fit$merge[i,]))),height=array_encode(fit$height),order=array_encode(fit$order),cluster=if(is.null(req$cut_k))NULL else array_encode(cutree(fit,k=req$cut_k)),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else if(op=='fit_glm') {
  d<-as.data.frame(lapply(req$data_files,function(path)numeric_file(path)),check.names=FALSE);form<-as.formula(paste(req$response,'~',paste(unlist(req$predictors),collapse=' + ')),env=baseenv());fam<-switch(req$family,gaussian=gaussian(link=if(req$link=='default')'identity' else req$link),binomial=binomial(link=if(req$link=='default')'logit' else req$link),poisson=poisson(link=if(req$link=='default')'log' else req$link),Gamma=Gamma(link=if(req$link=='default')'inverse' else req$link));gw<-if(is.null(req$weights_file))rep(1,nrow(d))else numeric_file(req$weights_file);go<-if(is.null(req$offset_file))rep(0,nrow(d))else numeric_file(req$offset_file);d$.rcli_weights<-gw;d$.rcli_offset<-go;fit<-glm(form,data=d,family=fam,weights=.rcli_weights,offset=.rcli_offset,control=glm.control(maxit=req$maxit));sm<-summary(fit);cf<-coef(fit);nm<-names(cf);tab<-sm$coefficients;column<-function(j){v<-rep(NA_real_,length(cf));names(v)<-nm;common<-intersect(nm,rownames(tab));v[common]<-tab[common,j];v};stype<-if(req$family%in%c('gaussian','Gamma'))'t'else'z';list(coefficients=as.list(nm),estimate=array_encode(cf),std_error=array_encode(column(2)),statistic=array_encode(column(3)),p_value=array_encode(column(4)),statistic_type=stype,dispersion=encode(sm$dispersion),aliased=as.list(unname(nm[is.na(cf)])),deviance=encode(fit$deviance),null_deviance=encode(fit$null.deviance),aic=encode(fit$aic),iterations=fit$iter,converged=isTRUE(fit$converged),domain_status=if(isTRUE(fit$converged))'OK' else 'NON_CONVERGED',rank=fit$rank,meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else if(op=='smooth_spline') {
  x<-numeric_file(req$x_file);y<-numeric_file(req$y_file);if(anyDuplicated(x)){u<-sort(unique(x));y<-vapply(u,function(v)mean(y[x==v]),numeric(1));x<-u}else if(is.unsorted(x)){o<-order(x);x<-x[o];y<-y[o]};args<-list(x=x,y=y,cv=isTRUE(req$cv));if(!is.null(req$df))args$df<-req$df;if(!is.null(req$spar))args$spar<-req$spar;fit<-do.call(smooth.spline,args);new<-if(is.null(req$new_x_file))NULL else predict(fit,x=numeric_file(req$new_x_file))$y;list(fitted=array_encode(fit$y),predicted=if(is.null(new))NULL else array_encode(new),lambda=encode(fit$lambda),df=encode(fit$df),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else if(op=='smooth_loess') {
  x<-numeric_file(req$x_file);y<-numeric_file(req$y_file);fit<-loess(y~x,data=data.frame(x=x,y=y),span=req$span,degree=req$degree,family=req$family);new<-if(is.null(req$new_x_file))NULL else predict(fit,newdata=data.frame(x=numeric_file(req$new_x_file)));list(fitted=array_encode(fitted(fit)),predicted=if(is.null(new))NULL else array_encode(as.numeric(new)),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else if(op=='interp_spline') {
  x<-numeric_file(req$x_file);y<-numeric_file(req$y_file);if(anyDuplicated(x)){u<-sort(unique(x));y<-vapply(u,function(v)mean(y[x==v]),numeric(1));x<-u}else if(is.unsorted(x)){o<-order(x);x<-x[o];y<-y[o]};xo<-numeric_file(req$xout_file);v<-if(req$method=='monoH.FC')splinefun(x,y,method=req$method)(xo)else spline(x,y,xout=xo,method=req$method)$y;list(value=array_encode(v),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind(),method=req$method))
 } else if(op=='interp_approx') {
  xx<-numeric_file(req$x_file);yy<-numeric_file(req$y_file);if(anyDuplicated(xx)){u<-sort(unique(xx));yy<-vapply(u,function(v)mean(yy[xx==v]),numeric(1));xx<-u};z<-approx(xx,yy,xout=numeric_file(req$xout_file),method=req$method,rule=req$rule);list(value=array_encode(z$y),meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind(),method=req$method,rule=req$rule))
 } else if(op=='fit_lm') {
  if(!is.null(req$file)) d<-read.csv(req$file,check.names=FALSE) else if(!is.null(req$data_files)) d<-as.data.frame(lapply(req$data_files,function(path) scan(path,what=double(),quiet=TRUE,na.strings='NA')),check.names=FALSE) else d<-as.data.frame(lapply(req$data,function(col) vapply(col,numeric_input,numeric(1))),check.names=FALSE)
  pieces<-strsplit(req$formula,'~',fixed=TRUE)[[1L]];if(length(pieces)!=2L) stop('bad formula');lhs<-trimws(pieces[1L]);rhs<-trimws(pieces[2L]);if(!lhs%in%names(d))stop('response must be a data column');tokens<-unlist(strsplit(gsub('[+-]',' ',rhs), '[[:space:]]+'));tokens<-tokens[nzchar(tokens)&!tokens%in%c('0','1')];if(any(!tokens%in%names(d))||'.'%in%tokens)stop('formula terms must be data columns');f<-as.formula(req$formula,env=baseenv());environment(f)<-baseenv();nafun<-if(identical(req$na_action,'exclude')) na.exclude else na.omit;fit<-lm(f,data=d,na.action=nafun);sm<-summary(fit);cf<-coef(fit);tab<-sm$coefficients;nm<-names(cf);col<-function(j){out<-rep(NA_real_,length(cf));names(out)<-nm;common<-intersect(nm,rownames(tab));out[common]<-tab[common,j];out};
  list(coefficients=as.list(nm),estimate=array_encode(cf),std_error=array_encode(col(2)),statistic=array_encode(col(3)),p_value=array_encode(col(4)),aliased=as.list(unname(nm[is.na(cf)])),rank=fit$rank,sigma=encode(sm$sigma),r_squared=encode(sm$r.squared),df=encode(sm$df[2]),residuals=if(isTRUE(req$residuals)) array_encode(residuals(fit)) else NULL,meta=list(r_version=as.character(getRversion()),rng_kind=RNGkind()))
 } else stop('unknown operation')
},error=function(e) structure(list(message=conditionMessage(e),tool_code=if(!is.null(e$tool_code))e$tool_code else if(grepl('reached elapsed time limit',conditionMessage(e),fixed=TRUE))'R_TIMEOUT'else'R_ERROR'),class='rcli_error')),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
if(inherits(result,'rcli_error')) {writeLines(json(list(error=result$message,tool_code=result$tool_code,warnings=as.list(warnings))),outpath);quit(status=0)}
result$warnings<-as.list(warnings)
writeLines(json(result),outpath)
