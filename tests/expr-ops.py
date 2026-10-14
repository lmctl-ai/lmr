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
maximum,_=call('opt_optimize',{'objective':'sin(x)','lower':0,'upper':3,'maximum':True})
near(maximum['x'],expected['optimize_maximum'][0]);near(maximum['value'],expected['optimize_maximum'][1])
x=[i*4/16 for i in range(17)];noise=[0,.01,-.01,.006,-.004,.003,-.006,.002,.004,-.003,.001,-.002,.005,-.004,.002,-.001,.003];y=[4*v/(2+v)+e for v,e in zip(x,noise)]
# Omitted controls must reproduce the R defaults in this Rscript environment.
for op,params,key,fields in [
 ('opt_optimize',{'objective':'sin(x)','lower':0,'upper':3},'default_optimize',['x','value']),
 ('num_root',{'f':'cos(x)-x','lower':0,'upper':1},'default_root',['root','f_root','iter','estim_prec']),
 ('num_integrate',{'f':'exp(-x^2)','lower':'-Inf','upper':'Inf'},'default_integrate',['value','abs_error','subdivisions']),
 ('fit_nls',{'response':'y','rhs':'a*x/(b+x)','data':{'x':x,'y':y},'start':{'a':3,'b':1}},'default_nls',['estimate','iterations','achieved_tol'])]:
 got,_=call(op,params);want=expected[key]
 actual=[]
 for field in fields:
  value=got[field]
  actual.extend(value if isinstance(value,list) else [value])
 for a,b in zip(actual,want):near(a,b,1e-10)
for method,key in [('Nelder-Mead','default_optim_Nelder_Mead'),('BFGS','default_optim_BFGS'),('CG','default_optim_CG'),('L-BFGS-B','default_optim_L_BFGS_B'),('SANN','default_optim_SANN')]:
 params={'objective':'(a-2)^2','start':{'a':-1},'method':method}
 if method=='SANN':params['seed']=42
 got,_=call('opt_minimize',params);want=expected[key]
 actual=[got['par']['a'],got['value'],got['convergence'],got['counts']['function']]
 for a,b in zip(actual,want[:4]):near(a,b,1e-10)
 if method=='L-BFGS-B':assert got['warnings']==[],got['warnings']
extended,_=call('num_root',{'f':'x^2-2','lower':2,'upper':3,'extend_int':'yes'})
near(extended['root'],expected['extended_root'][0],1e-4)
root,_=call('num_root',{'f':'x^2-2','lower':0,'upper':2,'tol':1e-8,'maxiter':100})
for a,b in zip([root['root'],root['f_root'],root['iter'],root['estim_prec']],expected['root']):near(a,b)
for upper,ref in [('-Inf','integral_inf'),(1,'integral_finite')]:
 lower=-1 if upper==1 else 'Inf'
 if upper=='-Inf': lower=-1;upper=1
 # Cover both a finite interval and the full-line Gaussian integral below.
 if ref=='integral_inf':lower='-Inf';upper='Inf';expr='exp(-x^2)'
 else:lower=0;upper=1;expr='exp(-x)'
 got,_=call('num_integrate',{'f':expr,'lower':lower,'upper':upper,'rel_tol':1e-8,'abs_tol':0})
 near(got['value'],expected[ref][0]);near(got['abs_error'],expected[ref][1]);assert got['domain_status']=='OK'
nls,_=call('fit_nls',{'response':'y','rhs':'a*x/(b+x)','data':{'x':x,'y':y},'start':{'a':3,'b':1},'maxiter':100,'tol':1e-8})
w=expected['nls'];assert nls['converged'] and nls['domain_status']=='OK'
for a,b in zip(nls['estimate'],w[:2]):near(a,b)
for a,b in zip(nls['std_error'],w[2:4]):near(a,b)
for a,b in zip(nls['statistic'],w[4:6]):near(a,b)
for a,b in zip(nls['p_value'],w[6:8]):near(a,b)
assert nls['df']==w[9] and nls['iterations']==w[11]
# Port starts outside the supplied bounds fail in C++ before starting R.
def run_rejected(op,params,tool_code):
 payload=json.dumps({'op':op,'params':params},separators=(',',':'))
 with tempfile.NamedTemporaryFile(mode='w',suffix='.json') as request:
  request.write(payload);request.flush()
  p=subprocess.run([binary,'--input-file',request.name],env={**__import__('os').environ,'RCLI_RSCRIPT':'/rcli-stage7-must-not-start-r'},capture_output=True,text=True)
 d=json.loads(p.stdout);assert p.returncode==2 and d['error']['data']['tool_code']==tool_code,(p.returncode,d,p.stderr)
run_rejected('fit_nls',{'response':'y','rhs':'a*x','data':{'x':[1,2],'y':[2,4]},'start':{'a':10},'algorithm':'port','lower':{'a':0},'upper':{'a':3}},'EXPR_SYNTAX')
# A one-iteration fit returns R's non-convergence message and no fake estimates.
maxed,_=call('fit_nls',{'response':'y','rhs':'a*x/(b+x)','data':{'x':x,'y':y},'start':{'a':3,'b':1},'maxiter':1})
assert maxed['domain_status']=='MAXIT' and maxed['estimate'] is None and maxed['message'] and maxed['iterations']==1,maxed
# Singular-gradient failure before iteration has zero iterations and null estimates.
ts=list(range(10));ys=[2*math.exp(-.3*t) for t in ts]
singular_nls,_=call('fit_nls',{'response':'y','rhs':'a*exp(-b*t)+c*exp(-b*t)','data':{'t':ts,'y':ys},'start':{'a':1,'b':.3,'c':1}})
assert singular_nls['domain_status']=='NON_CONVERGED' and singular_nls['estimate'] is None and singular_nls['iterations']==0 and singular_nls['message'],singular_nls
# R max-iteration exhaustion is a domain status, not a worker error.
limited,_=call('opt_minimize',{'objective':'100*(b-a^2)^2+(1-a)^2','start':{'a':-1.2,'b':1},'method':'BFGS','maxit':1})
assert limited['domain_status']=='MAXIT' and limited['convergence']==1
# Seeded non-finite simplex evaluations are replaced with +Inf and counted.
nonfinite,_=call('opt_minimize',{'objective':'ifelse(a>0.98,(a-0.5)^2,sqrt(-1))','start':{'a':1},'method':'Nelder-Mead','maxit':100})
assert nonfinite['nonfinite_evals']>0 and nonfinite['domain_status'] in ('OK','MAXIT')
bad_start,_=call('opt_minimize',{'objective':'log(a)','start':{'a':-1},'method':'BFGS'})
assert bad_start['domain_status']=='NONFINITE_OBJECTIVE' and bad_start['value']=='NaN' and bad_start['nonfinite_evals']==1,bad_start
special_data,_=call('opt_minimize',{'objective':'a^2','start':{'a':1},'method':'BFGS','data':{'t':[1,'Inf','NaN','-Inf']}})
assert special_data['domain_status']=='OK',special_data
# Evaluation cap is enforced inside the objective closure.
err,rc=call('opt_minimize',{'objective':'(a-2)^2','start':{'a':0},'method':'BFGS','max_evals':1},False)
assert rc==1 and err['error']['data']['tool_code']=='EVAL_LIMIT',err
# Integrable singularity produces R's warning at the low subdivision cap; warnings stay in the result.
singular,_=call('num_integrate',{'f':'1/sqrt(abs(x-0.1234567))','lower':0,'upper':1,'subdivisions':1})
assert singular['warnings'] or singular['domain_status']=='NON_FINITE',singular
# Constants recycle across integrate's vector x; subnormal finite literals parse successfully.
constant,_=call('num_integrate',{'f':'pi','lower':0,'upper':1})
near(constant['value'],math.pi)
precedence,_=call('num_integrate',{'f':'ifelse(x<1 | x>4 & x<0,1,0)','lower':0,'upper':5,'rel_tol':1e-8,'abs_tol':0})
near(precedence['value'],1,1e-10)
nonfinite_integral,_=call('num_integrate',{'f':'sqrt(-1)','lower':0,'upper':1})
assert nonfinite_integral['domain_status']=='NON_FINITE' and nonfinite_integral['meta']['r_version'],nonfinite_integral
subnormal,_=call('opt_optimize',{'objective':'abs(x-1e-310)','lower':0,'upper':1e-300})
# Deliberately expensive vector expression trips max_evals rapidly, well below process timeout.
started=__import__('time').perf_counter();big=call('opt_minimize',{'objective':'sum(exp(sin(exp(sin(x)))))','start':{'a':1},'data':{'x':[i/1000 for i in range(200000)]},'method':'Nelder-Mead','max_evals':2},False);elapsed=__import__('time').perf_counter()-started
assert big[1]==1 and big[0]['error']['data']['tool_code']=='EVAL_LIMIT' and elapsed<10,(big,elapsed)
# Throughput regression: compilation is outside the SANN callback, leaving one AST eval per callback.
started=__import__('time').perf_counter();perf,_=call('opt_minimize',{'objective':'(a-2)^2','start':{'a':-1},'method':'SANN','seed':42,'maxit':10000});perf_s=__import__('time').perf_counter()-started
assert perf["evals"]==10001 and perf_s<50,(perf['evals'],perf_s)
print(f'expression R-oracle cases passed; 200k-row guard: {elapsed:.3f}s; SANN 10001 evaluations: {perf_s:.3f}s ({perf["evals"]/perf_s:.0f} evals/s)')
