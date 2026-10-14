import json,math,pathlib,subprocess,sys,tempfile
binary=sys.argv[1];root=pathlib.Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory() as td:
 oracle=pathlib.Path(td)/'oracle.txt';subprocess.run(['Rscript','--vanilla','--default-packages=NULL',str(root/'tests/expr-oracle.R'),str(oracle)],check=True,capture_output=True,text=True)
 expected={}
 for line in oracle.read_text().splitlines():
  parts=line.split();expected[parts[0]]=[math.nan if x=='NA' else float(x) for x in parts[1:]]
def call(op,params,ok=True):
 payload=json.dumps({'op':op,'params':params},separators=(',',':'))
 if len(payload)<100000:
  args=[binary,'--input-json',payload]
  p=subprocess.run(args,capture_output=True,text=True)
 else:
  with tempfile.NamedTemporaryFile(mode='w',suffix='.json') as request:
   request.write(payload);request.flush()
   p=subprocess.run([binary,'--input-file',request.name],capture_output=True,text=True)
 doc=json.loads(p.stdout)
 if ok: assert p.returncode==0,(op,params,p.returncode,doc,p.stderr)
 return doc,p.returncode
def near(a,b,tol=2e-12): assert math.isclose(float(a),float(b),rel_tol=tol,abs_tol=tol),(a,b)
def check_min(method,key,params):
 got,_=call('opt_minimize',dict(objective='100*(b-a^2)^2+(1-a)^2',start={'a':-1.2,'b':1},method=method,maxit=1000,reltol=1e-8,**params))
 want=expected[key]
 near(got['par']['a'],want[0]);near(got['par']['b'],want[1]);near(got['value'],want[2]);assert got['convergence']==int(want[3]);assert got['counts']['function']==int(want[4]);assert got['domain_status']=='OK'
check_min('BFGS','rosen_BFGS',{})
check_min('Nelder-Mead','rosen_Nelder_Mead',{})
box,_=call('opt_minimize',{'objective':'(a-4)^2','start':{'a':.5},'method':'L-BFGS-B','lower':{'a':0},'upper':{'a':2}})
for a,b in zip([box['par']['a'],box['value'],box['convergence'],box['counts']['function']],expected['box'][:4]):near(a,b)
sann,_=call('opt_minimize',{'objective':'(a-2)^2','start':{'a':-1},'method':'SANN','seed':42,'maxit':80,'reltol':1e-8})
for a,b in zip([sann['par']['a'],sann['value'],sann['convergence'],sann['counts']['function']],expected['sann'][:4]):near(a,b)
opt,_=call('opt_optimize',{'objective':'(x-2)^2','lower':-5,'upper':5,'tol':1e-8})
near(opt['x'],expected['optimize'][0]);near(opt['value'],expected['optimize'][1])
root,_=call('num_root',{'f':'x^2-2','lower':0,'upper':2,'tol':1e-8,'maxiter':100})
for a,b in zip([root['root'],root['f_root'],root['iter'],root['estim_prec']],expected['root']):near(a,b)
for upper,ref in [('-Inf','integral_inf'),(1,'integral_finite')]:
 lower=-1 if upper==1 else 'Inf'
 if upper=='-Inf': lower=-1;upper=1
 # Cover both a finite interval and the full-line Gaussian integral below.
 if ref=='integral_inf':lower='-Inf';upper='Inf';expr='exp(-x^2)'
 else:lower=0;upper=1;expr='exp(-x)'
 got,_=call('num_integrate',{'f':expr,'lower':lower,'upper':upper,'rel_tol':1e-8})
 near(got['value'],expected[ref][0]);near(got['abs_error'],expected[ref][1]);assert got['domain_status']=='OK'
x=[i*4/16 for i in range(17)];noise=[0,.01,-.01,.006,-.004,.003,-.006,.002,.004,-.003,.001,-.002,.005,-.004,.002,-.001,.003];y=[4*v/(2+v)+e for v,e in zip(x,noise)]
nls,_=call('fit_nls',{'response':'y','rhs':'a*x/(b+x)','data':{'x':x,'y':y},'start':{'a':3,'b':1}})
w=expected['nls'];assert nls['converged'] and nls['domain_status']=='OK'
for a,b in zip(nls['estimate'],w[:2]):near(a,b)
for a,b in zip(nls['std_error'],w[2:4]):near(a,b)
for a,b in zip(nls['statistic'],w[4:6]):near(a,b)
for a,b in zip(nls['p_value'],w[6:8]):near(a,b)
assert nls['df']==w[9] and nls['iterations']==w[11]
# R max-iteration exhaustion is a domain status, not a worker error.
limited,_=call('opt_minimize',{'objective':'100*(b-a^2)^2+(1-a)^2','start':{'a':-1.2,'b':1},'method':'BFGS','maxit':1})
assert limited['domain_status']=='MAXIT' and limited['convergence']==1
# Seeded non-finite simplex evaluations are replaced with +Inf and counted.
nonfinite,_=call('opt_minimize',{'objective':'ifelse(a>0.98,(a-0.5)^2,sqrt(-1))','start':{'a':1},'method':'Nelder-Mead','maxit':100})
assert nonfinite['nonfinite_evals']>0 and nonfinite['domain_status'] in ('OK','MAXIT')
# Evaluation cap is enforced inside the objective closure.
err,rc=call('opt_minimize',{'objective':'(a-2)^2','start':{'a':0},'method':'BFGS','max_evals':1},False)
assert rc==1 and err['error']['data']['tool_code']=='EVAL_LIMIT',err
# Integrable singularity produces R's warning at the low subdivision cap; warnings stay in the result.
singular,_=call('num_integrate',{'f':'1/sqrt(abs(x-0.1234567))','lower':0,'upper':1,'subdivisions':1})
assert singular['warnings'] or singular['domain_status']=='NON_FINITE',singular
# Deliberately expensive vector expression trips max_evals rapidly, well below process timeout.
started=__import__('time').perf_counter();big=call('opt_minimize',{'objective':'sum(exp(sin(exp(sin(x)))))','start':{'a':1},'data':{'x':[i/1000 for i in range(200000)]},'method':'Nelder-Mead','max_evals':2},False);elapsed=__import__('time').perf_counter()-started
assert big[1]==1 and big[0]['error']['data']['tool_code']=='EVAL_LIMIT' and elapsed<10,(big,elapsed)
print(f'expression R-oracle cases passed; 200k-row nested-vector guard: {elapsed:.3f}s')
