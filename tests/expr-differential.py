import json, pathlib, random, subprocess, sys, tempfile

binary = sys.argv[1]
rng = random.Random(20261013)
variables = ["x", "y"]
functions = ["abs", "sqrt", "sin", "log1p", "pmin", "pmax"]

def atom(depth):
    if depth <= 0 or rng.random() < 0.28:
        return rng.choice(variables + ["0.25", "0.5", "1.25", "2"])
    if rng.random() < 0.2:
        f = rng.choice(functions)
        if f in ("pmin", "pmax"):
            return f"{f}({atom(depth-1)},{atom(depth-1)})"
        return f"{f}(abs({atom(depth-1)}))" if f == "sqrt" else f"{f}({atom(depth-1)})"
    return f"({atom(depth-1)} {rng.choice(['+','-','*','/'])} {atom(depth-1)})"

expressions=[]
for i in range(2000):
    if i % 4 == 0:
        # Exercise R's distinct logical precedence in varied relational trees.
        expr=f"ifelse({atom(1)}<1 | {atom(1)}>0.5 & {atom(1)}<2,{atom(1)},{atom(1)})"
    else:
        expr=atom(3)
    expressions.append(expr)

reqs=[{"expression":expr,"names":variables} for expr in expressions]
p=subprocess.run([binary,"--expr-check-batch"],input=json.dumps({"cases":reqs}),text=True,capture_output=True,check=True)
decisions=json.loads(p.stdout)["results"]
assert len(decisions)==len(expressions)
cases=[]
for expr,got in zip(expressions,decisions):
    assert got["accepted"],(expr,got)
    cases.append({"expression":expr,"ast":got["ast"]})

def r_literal(v):
    if isinstance(v,dict):
        return "list("+",".join(k+"="+r_literal(x) for k,x in v.items())+")"
    if isinstance(v,list): return "list("+",".join(r_literal(x) for x in v)+")"
    if isinstance(v,str): return json.dumps(v)
    if isinstance(v,bool): return "TRUE" if v else "FALSE"
    return repr(v)

with tempfile.TemporaryDirectory() as td:
    rfile=pathlib.Path(td)/"differential.R"
    lines=[
      "eval_ast <- function(n, env) {",
      " if(n$t=='num') return(n$v); if(n$t=='const') return(if(n$n=='pi') pi else exp(1)); if(n$t=='var') return(env[[n$n]])",
      " if(n$t=='neg') return(-eval_ast(n$a[[1]],env))",
      " if(n$t=='op') { a<-eval_ast(n$a[[1]],env); b<-eval_ast(n$a[[2]],env); return(switch(n$o,'+'=a+b,'-'=a-b,'*'=a*b,'/'=a/b,'^'=a^b,'<'=a<b,'<='=a<=b,'>'=a>b,'>='=a>=b,'=='=a==b,'!='=a!=b,'&'=a&b,'|'=a|b)) }",
      " if(n$t=='call') { a<-lapply(n$a,eval_ast,env=env); return(do.call(get(n$f,baseenv()),a)) }",
      " if(n$t=='ifelse') return(ifelse(eval_ast(n$c,env),eval_ast(n$a,env),eval_ast(n$b,env)))",
      " stop('unexpected AST node')",
      "}",
      "x<-c(-1.25,0.25,1.5); y<-c(0.5,1.75,-0.5); env<-new.env(parent=baseenv());env$x<-x;env$y<-y",
      "cases<-list("+",\n".join("list(expression="+json.dumps(c['expression'])+",ast="+r_literal(c['ast'])+")" for c in cases)+")",
      "for(i in seq_along(cases)){ z<-cases[[i]]; a<-eval_ast(z$ast,env); b<-eval(parse(text=z$expression)[[1]],env); if(!identical(as.vector(a),as.vector(b))) stop(sprintf('case %d differs: %s',i,z$expression)) }",
      "cat(sprintf('differential expressions: %d accepted and bit-identical to R parse/eval\\n',length(cases)))"
    ]
    rfile.write_text("\n".join(lines))
    p=subprocess.run(["Rscript","--vanilla","--default-packages=NULL",str(rfile)],capture_output=True,text=True)
    assert p.returncode==0,(p.stdout,p.stderr)
    print(p.stdout.strip())
