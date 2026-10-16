import json,subprocess,sys,pathlib,tempfile,os
b=sys.argv[1]
def call(op,p,no_r=False):
 payload=json.dumps({'op':op,'params':p},separators=(',',':'))
 if len(payload)<100000:
  args=[b,'--input-json',payload]
 else:
  with tempfile.NamedTemporaryFile(mode='w') as f:
   f.write(payload);f.flush();args=[b,'--input-file',f.name]
   env=os.environ.copy();env['RCLI_RSCRIPT']='/rcli-test-no-rscript' if no_r else env.get('RCLI_RSCRIPT','/usr/bin/Rscript')
   r=subprocess.run(args,capture_output=True,text=True,env=env);return r.returncode,json.loads(r.stdout)
 env=os.environ.copy();env['RCLI_RSCRIPT']='/rcli-test-no-rscript' if no_r else env.get('RCLI_RSCRIPT','/usr/bin/Rscript')
 r=subprocess.run(args,capture_output=True,text=True,env=env)
 return r.returncode,json.loads(r.stdout)
# R default prop.test permits per-group probabilities to be omitted.
rc,d=call('test_htest',{'test':'prop','x':[15,25],'trials':[40,50]})
assert rc==0 and isinstance(d['p_value'],float),('prop omitted p',rc,d)
# Wilcoxon one-sample signed rank.
rc,d=call('test_htest',{'test':'wilcox','x':[1,2,3,4,6,7,8,9,10,11],'mu':5,'conf_int':True})
assert rc==0 and d['conf_int'] is not None,('one-sample wilcox',rc,d)
# Sampling a one-element vector samples the vector item, not 1:x[1].
rc,d=call('rng_sample',{'x':[10],'size':1,'seed':4})
assert rc==0 and d['value']==[10],('sample singleton',rc,d)
# ks two sample must use y rather than silently test against a named CDF.
rc,d=call('test_htest',{'test':'ks','x':[.1,.5,.9,1.3],'y':[2,3,4,5]})
assert rc==0 and abs(d['p_value']-.0285714285714286)<1e-14,('two-sample KS',rc,d)
# R fft(inverse=TRUE) is unnormalized unless normalize is requested.
rc,d=call('num_fft',{'x':[1,2,3,4,5],'inverse':True})
assert rc==0 and d['re'][0]==15,('inverse fft',rc,d)
rc,d=call('ts_arima',{'x':[1,3,2,4,3,5,4,2,6,3,1,5,4,7,2,6],'order':[1,0,1],'method':'CSS'})
assert rc==0 and d['aic']=='NA' and d['converged'] is False,('logical NA and convergence',rc,d)
print('review2 principal regressions passed')
# Rejection budgets are enforced by C++ before the worker is invoked.
def reject(op,p,code):
 rc,d=call(op,p,no_r=True); got=d.get('error',{}).get('data',{}).get('tool_code')
 assert rc==2 and got==code,(op,code,rc,d)
reject('ts_arima',{'x':[float(i%7) for i in range(40)],'order':[1,0,0],'seasonal_order':[2,0,0],'period':350},'BAD_ORDER')
reject('ts_arima',{'x':[0.0]*200000,'order':[0,0,0],'seasonal_order':[1,0,0],'period':350},'INPUT_TOO_LARGE')
reject('ts_arima',{'x':[0.0]*8,'order':[0,0,0],'seasonal_order':[0,0,0],'period':351},'BAD_ORDER')
reject('ts_acf',{'x':[float(i%7) for i in range(2001)],'lag_max':1001},'BAD_DATA')
reject('ts_pacf',{'x':[float(i%7) for i in range(2001)],'lag_max':1001},'BAD_DATA')
reject('num_fft',{'x':[0.0]*1009},'BAD_FFT_LENGTH')
reject('test_htest',{'test':'cor','x':list(range(5001)),'y':list(range(5001)),'method':'kendall'},'INPUT_TOO_LARGE')
reject('test_htest',{'test':'fisher','x':[250001]*4},'INPUT_TOO_LARGE')
reject('smooth_loess',{'x':list(range(20001)),'y':list(range(20001))},'BAD_DATA')
reject('smooth_spline',{'x':list(range(100001)),'y':list(range(100001))},'BAD_DATA')
reject('num_fft',{'x':[0.0]*100001},'INPUT_TOO_LARGE')
reject('ts_stl',{'x':list(range(80001)),'period':4,'periodic':True},'BAD_PARAMETER')
reject('cluster_kmeans',{'data':{'x':[float(i%2) for i in range(200000)]},'k':500,'seed':42,'nstart':100},'BAD_PARAMETER')
# One-sample and two-sample tests, plus explicitly named table ordering.
rc,d=call('test_htest',{'test':'ks','x':[.1,.5,.9,1.3],'y':[2,3,4,5],'distribution':'norm'})
assert rc==2 and d['error']['data']['tool_code']=='BAD_PARAMETER',('KS mixed inputs',rc,d)
rc,d=call('test_htest',{'test':'fisher','x':[10,20,30,40,50,60],'rows':2,'cols':3,'byrow':True})
assert rc==0 and d['p_value'] is not None,('fisher rxc',rc,d)
# Quoted CSV header names follow CSV quoting rules.
with tempfile.TemporaryDirectory() as temp:
 path=pathlib.Path(temp)/'quoted.csv';path.write_text('"y","x"\n1,1\n2,2\n4,3\n5,4\n')
 # Fit is constrained to cwd, so create it under the example data directory for this CLI test.
 pth=pathlib.Path('examples/data/review2-quoted.csv');pth.write_text(path.read_text())
 try:
  rc,d=call('fit_lm',{'formula':'y ~ x','file':str(pth)})
  # lm() through LAPACK QR: the coefficient error is bounded by about cond(X) * eps with cond([1, x], x = 1..4) = 10, so 64 eps = 1.4e-14
  # relative is the platform-independent bound (the literal bits are R 4.3.3 on aarch64; amd64 may differ in the last bit).
  want=[-0.5,1.4]
  assert rc==0 and len(d['estimate'])==2 and all(abs(g-w)<=64*2.220446049250313e-16*max(1.0,abs(w)) for g,w in zip(d['estimate'],want)),('quoted CSV header',rc,d)
 finally: pth.unlink(missing_ok=True)
