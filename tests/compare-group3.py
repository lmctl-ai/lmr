import json
import pathlib
import subprocess
import sys

binary, oracle = sys.argv[1:]
oracle=pathlib.Path(oracle)
subprocess.run(["Rscript","--vanilla","--default-packages=NULL","tests/group3-oracle.R",str(oracle)],check=True)

def call(op,params,expected):
    argv=[binary,op]
    for key,value in params.items():
        vals=value if isinstance(value,list) else [value]
        for item in vals:
            if isinstance(item,bool): argv.append(("--" if item else "--no-")+key)
            elif isinstance(item,(dict,list)): argv.extend(["--"+key,json.dumps(item,separators=(",",":"))])
            else: argv.extend(["--"+key,str(item)])
    result=json.loads(subprocess.run(argv,check=True,text=True,capture_output=True).stdout)
    for key,name in expected.items():
        want_file=oracle/(name+".txt")
        raw=want_file.read_text().splitlines()
        if key in ("coefficients","converged"):
            want=raw if key=="coefficients" else raw[0]=="TRUE"
            assert result[key]==want,(op,key,result[key],want)
        elif key in ("domain_status",):
            assert result[key]==raw[0],(op,key,result[key],raw)
        else:
            want=[float(x) for x in raw]
            got=result[key]
            if not isinstance(got,list): got=[float(got)]
            else: got=[float(x) for x in got]
            assert got==want,(op,key,got,want)
    return result

data={"x":list(range(1,9)),"yg":[1,2,3,5,7,8,9,11],"yb":[0,1,0,1,1,0,1,0],"yp":[1,0,2,1,3,4,2,5],"yG":[1.2,2.1,1.7,3,4.2,2.5,5,4.6]}
families={"gaussian":("yg","identity"),"binomial":("yb","logit"),"poisson":("yp","log"),"Gamma":("yG","inverse")}
for family,(response,link) in families.items():
    prefix="glm-"+family+"-"
    call("fit_glm",{"response":response,"predictors":["x"],"data":data,"family":family,"link":link},{k:prefix+suffix for k,suffix in (("coefficients","names"),("estimate","estimate"),("std_error","se"),("statistic","stat"),("p_value","p"),("deviance","deviance"),("null_deviance","null"),("aic","aic"),("iterations","iter"),("converged","converged"))})

x=list(range(1,13));y=[0,1,0,2,1,3,2,4,3,5,4,6]
call("smooth_spline",{"x":x,"y":y,"df":4,"new_x":[2.5,6.5,10.5]},{"fitted":"spline-fitted","predicted":"spline-new","lambda":"spline-lambda","df":"spline-df"})
call("smooth_spline",{"x":x,"y":y,"spar":.5},{"fitted":"spline-spar"})
call("smooth_spline",{"x":x,"y":y,"cv":True},{"fitted":"spline-cv"})
call("smooth_loess",{"x":x,"y":y,"span":1,"degree":1,"new_x":[2.5,6.5,10.5]},{"fitted":"loess-fitted","predicted":"loess-new"})
for method in ("fmm","natural","periodic","monoH.FC"):
    call("interp_spline",{"x":[1,2,3,4,5,6],"y":[0,1,0,-1,0,0],"xout":[1.5,2.5,3.5,4.5,5.5],"method":method},{"value":"is-"+method})
call("interp_approx",{"x":[1,2,3,4,5,6],"y":[0,1,0,-1,0,0],"xout":[1.5,2.5,3.5,4.5,5.5]},{"value":"approx-linear"})
call("interp_approx",{"x":[1,2,3,4,5,6],"y":[0,1,0,-1,0,0],"xout":[.5,1.5,5.5,6.5],"method":"constant","rule":2},{"value":"approx-constant"})

for op,params,code in [
    ("fit_glm",{"response":"yg","predictors":["missing"],"data":data,"family":"gaussian"},"BAD_DATA"),
    ("fit_glm",{"response":"yg","predictors":["x"],"data":data,"family":"poisson","link":"cloglog"},"BAD_LINK"),
    ("smooth_spline",{"x":[1,2,3,4],"y":[1,2,3,4],"df":1},"BAD_PARAMETER"),
    ("smooth_loess",{"x":[1,2,3],"y":[1,2,3],"degree":3},"BAD_PARAMETER"),
]:
    argv=[binary,op]
    for key,value in params.items():
        vals=value if isinstance(value,list) else [value]
        for item in vals: argv.extend(["--"+key,json.dumps(item,separators=(",",":")) if isinstance(item,(dict,list)) else str(item)])
    proc=subprocess.run(argv,text=True,capture_output=True)
    doc=json.loads(proc.stdout)
    assert proc.returncode==2 and doc.get("error",{}).get("data",{}).get("tool_code")==code,(op,params,proc.stdout)

# A deliberately short IRLS cap must yield a successful domain result instead of R_ERROR.
nonconv=data.copy();nonconv["y_sep"]=[0,0,0,0,1,1,1,1]
result=call("fit_glm",{"response":"y_sep","predictors":["x"],"data":nonconv,"family":"binomial","maxit":1},{"converged":"glm-nonconv"})
assert result["domain_status"]==(oracle/"glm-nonconv-status.txt").read_text().strip()
