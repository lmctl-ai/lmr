import hashlib,json,os,pathlib,subprocess,sys,tempfile
binary,worker,spec_path=sys.argv[1:]
worker=pathlib.Path(worker); spec_path=pathlib.Path(spec_path); root=pathlib.Path(__file__).resolve().parents[1]
corpus=json.loads((root/'tests/expr-corpus.json').read_text())
def cpp_check(case):
    p=subprocess.run([binary,'--expr-check'],input=json.dumps(case,ensure_ascii=False),text=True,capture_output=True,check=True)
    return json.loads(p.stdout)
def r_batch(cases):
    with tempfile.TemporaryDirectory() as td:
        inp=pathlib.Path(td)/'in.json';out=pathlib.Path(td)/'out.json'
        inp.write_text(json.dumps({'cases':cases}))
        p=subprocess.run(['Rscript','--vanilla','--default-packages=NULL',str(worker),'--','expr_builder_batch',str(inp),str(out)],capture_output=True,text=True)
        assert p.returncode==0,(p.stdout,p.stderr)
        return json.loads(out.read_text())['decisions']
sha=hashlib.sha256(spec_path.read_bytes()).hexdigest()
assert subprocess.check_output([binary,'--expr-spec-hash'],text=True).strip()==sha
assert subprocess.check_output(['Rscript','--vanilla','--default-packages=NULL',str(worker),'--expr-spec-hash'],text=True).strip()==sha
accepted=0;accepted_asts=[]
for case in corpus['accept']:
    got=cpp_check(case);assert got['accepted'],(case,got);assert got['sha256']==sha
    accepted_asts.append({'ast':got['ast'],'names':case.get('names',[]),'allow_reducers':case.get('allow_reducers',False)});accepted+=1
assert all(x['accepted'] for x in r_batch(accepted_asts))
rejected=0;bad={'t':'call','f':'system','a':[]}
for case in corpus['reject']:
    assert case.get('expected_tool_code'),case
    got=cpp_check(case);assert not got['accepted'],(case,got);assert got['tool_code']==case['expected_tool_code'],(case,got);rejected+=1
results=r_batch([{'ast':bad,'names':['system']} for _ in corpus['reject']]);assert len(results)==len(corpus['reject']) and all(not x['accepted'] and x['tool_code']=='EXPR_SPEC_MISMATCH' for x in results)
# Force a missing child executable for each invalid expression to prove rejection precedes exec.
for case in corpus['reject']:
    names=set(case.get('names',[]));start={n:1.0 for n in names if n not in {'pi','e'} and n.isascii() and n.islower() and n.replace('_','a').isalnum() and len(n)<=32}
    if not start:start={'a':1.0}
    request={'op':'num_integrate','params':{'f':case['expression'],'lower':0,'upper':1}}
    with tempfile.NamedTemporaryFile('w',delete=False) as f:json.dump(request,f,ensure_ascii=False);path=f.name
    env=os.environ.copy();env['RCLI_RSCRIPT']='/rcli-stage1b-must-not-start-r'
    try:
        p=subprocess.run([binary,'--input-file',path],env=env,capture_output=True,text=True);out=json.loads(p.stdout)
        assert p.returncode==2 and out.get('error',{}).get('data',{}).get('tool_code')!='R_NOT_FOUND',(case,p.returncode,out,p.stderr)
    finally:pathlib.Path(path).unlink(missing_ok=True)
def reject_before_r(req,code=None):
    env=os.environ.copy();env['RCLI_RSCRIPT']='/rcli-stage1b-must-not-start-r'
    with tempfile.NamedTemporaryFile('w',delete=False) as f:json.dump(req,f);path=f.name
    try:
        p=subprocess.run([binary,'--input-file',path],env=env,capture_output=True,text=True);out=json.loads(p.stdout)
        got=out.get('error',{}).get('data',{}).get('tool_code')
        assert p.returncode==2 and got!='R_NOT_FOUND' and (code is None or got==code),(req,p.returncode,out,p.stderr)
    finally:pathlib.Path(path).unlink(missing_ok=True)
for params,code in [
 ({'objective':'a','start':{'a':1},'method':'SANN'},'EXPR_SYNTAX'),
 ({'objective':'a','start':{'a':1},'method':'BFGS','seed':2},'EXPR_SYNTAX'),
 ({'objective':'exp(x)','start':{'a':1},'data':{'x':[1,2]}},'EXPR_NOT_SCALAR'),
 ({'objective':'a','start':{'exp':1}},'EXPR_NAME_COLLISION'),
 ({'f':'sum(x)','lower':0,'upper':1},'EXPR_NOT_SCALAR'),
 ({'objective':'a','start':{'a':1},'max_evals':2000001},'EXPR_TOO_LONG'),
]:
    reject_before_r({'op':'opt_minimize' if 'objective' in params else 'num_integrate','params':params},code)
assert not r_batch([{'ast':{'t':'call','f':'system','a':[{'t':'num','v':1}]},'names':['system']}])[0]['accepted']
# The hidden parser must keep the error excerpt valid UTF-8 when its byte window
# intersects multibyte characters.
near=cpp_check({'expression':'aaaaaaaaéaaaa€b','names':['a']})
assert not near['accepted'] and len(near['near'].encode('utf-8'))<=16,near
print(f'expression corpus: {accepted} accepted, {rejected} rejected; worker builder matched accepted cases and rejected hostile AST probes')
