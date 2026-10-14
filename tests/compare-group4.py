#!/usr/bin/env python3
import json, subprocess, sys, tempfile, pathlib
bin=sys.argv[1]
def invoke(op,params):
 p=subprocess.run([bin,'--input-json',json.dumps({'op':op,'params':params},separators=(',',':'))],text=True,capture_output=True)
 if p.returncode: raise AssertionError(f'{op} failed: {p.stdout} {p.stderr}')
 return json.loads(p.stdout)
r='''args<-commandArgs(TRUE); out<-args[1]; enc<-function(x)paste0('[',paste(sprintf('%.17g',as.numeric(x)),collapse=','),']'); x<-c(1,2,1,3,2,4,2,5,3,5,4,6,5,7,5,8,6,8,7,9); RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(42);sim<-as.numeric(arima.sim(list(ar=.5,ma=numeric()),n=24));fit<-arima(x,order=c(1,0,0),method='CSS-ML'); ac<-acf(x,lag.max=4,plot=FALSE);pa<-pacf(x,lag.max=4,plot=FALSE); st<-stl(ts(1:24,frequency=4),s.window='periodic'); d<-data.frame(x=c(1,1.1,.9,1.2,10,10.2,9.8,10.1));RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(42);km<-kmeans(d,centers=2,nstart=5,algorithm='Hartigan-Wong');hc<-hclust(dist(d),method='complete'); fields<-c(paste0('sim=',enc(sim)),paste0('arima_coef=',enc(coef(fit))),paste0('arima_se=',enc(sqrt(diag(fit$var.coef)))),paste0('arima_fc=',enc(predict(fit,n.ahead=2)$pred)),paste0('arima_sefc=',enc(predict(fit,n.ahead=2)$se)),paste0('acf=',enc(ac$acf)),paste0('pacf=',enc(pa$acf)),paste0('stl=',enc(st$time.series)),paste0('km_centers=',enc(as.vector(t(km$centers)))),paste0('km_cluster=',enc(km$cluster)),paste0('hc_height=',enc(hc$height)),paste0('hc_order=',enc(hc$order)));writeLines(fields,out)'''
with tempfile.NamedTemporaryFile() as f:
 subprocess.run(['Rscript','--vanilla','-e',r,f.name],check=True)
 oracle={k:json.loads(v) for k,v in (line.split('=',1) for line in pathlib.Path(f.name).read_text().splitlines())}
x=[1,2,1,3,2,4,2,5,3,5,4,6,5,7,5,8,6,8,7,9]
sim=invoke('rng_arima_sim',{'n':24,'seed':42,'order':[1,0,0],'ar':[.5]})['value']; assert sim==oracle['sim']
a=invoke('ts_arima',{'x':x,'order':[1,0,0],'n_ahead':2,'residuals':True});assert a['estimate']==oracle['arima_coef'] and a['std_error']==oracle['arima_se'] and a['forecast']==oracle['arima_fc'] and a['forecast_se']==oracle['arima_sefc'] and a['converged'] is True
assert invoke('ts_acf',{'x':x,'lag_max':4})['value']==oracle['acf']
assert invoke('ts_pacf',{'x':x,'lag_max':4})['value']==oracle['pacf']
st=invoke('ts_stl',{'x':list(range(1,25)),'period':4,'periodic':True});assert st['seasonal']+st['trend']+st['remainder']==oracle['stl']
d={'x':[1,1.1,.9,1.2,10,10.2,9.8,10.1]};km=invoke('cluster_kmeans',{'data':d,'k':2,'seed':42,'nstart':5});assert [v for row in km['centers'] for v in row]==oracle['km_centers'] and km['cluster']==oracle['km_cluster']
hc=invoke('cluster_hclust',{'data':d,'method':'complete','cut_k':2});assert hc['height']==oracle['hc_height'] and hc['order']==oracle['hc_order']
print('group4 R oracle: rng_arima_sim, ts_arima, ts_acf/pacf, ts_stl, cluster_kmeans/hclust')
