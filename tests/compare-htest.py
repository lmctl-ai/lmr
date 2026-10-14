import json
import pathlib
import subprocess
import sys

binary, expected_dir = sys.argv[1:]
expected_dir = pathlib.Path(expected_dir)
subprocess.run(["Rscript", "--vanilla", "--default-packages=NULL", "tests/htest-oracle.R", str(expected_dir)], check=True)

def call(name, params):
    argv=[binary,"test_htest"]
    for key,value in params.items():
        for item in (value if isinstance(value,list) else [value]):
            if isinstance(item,bool): argv.append(("--" if item else "--no-")+key)
            else: argv.extend(["--"+key,str(item)])
    return json.loads(subprocess.run(argv,check=True,text=True,capture_output=True).stdout)

cases={
 "t":{"test":"t","x":[1,2,3,4]},
 "t-paired":{"test":"t","x":[1,2,4,6],"y":[0,1,2,3],"paired":True},
 "wilcox":{"test":"wilcox","x":list(range(1,9)),"y":list(range(3,11))},
 "chisq":{"test":"chisq","x":[10,20,30],"expected":[.2,.3,.5]},
 "chisq-table":{"test":"chisq","x":[10,20,30,40],"rows":2,"cols":2},
 "shapiro":{"test":"shapiro","x":[1,2,3,5,7]},
 "ks":{"test":"ks","x":[-1,0,1,2],"distribution":"norm"},
 "fisher":{"test":"fisher","x":[1,9,11,3]},
 "binom":{"test":"binom","x":[7],"n":10},
 "prop":{"test":"prop","x":[5,7],"trials":[10,10],"probabilities":[.5,.5]},
 "cor":{"test":"cor","x":[1,2,3,4,5],"y":[2,1,4,3,7],"method":"spearman"},
 "var":{"test":"var","x":[1,2,3,4,5],"y":[2,4,5,7,8]},
}
for name,params in cases.items():
    actual=call(name,params)
    assert "statistic_name" in actual and "estimate_names" in actual and "null_value" in actual
    expected=(expected_dir/(name+".txt")).read_text().splitlines()
    stat,parameter,pvalue,conf,estimate,method,alternative,null_value=expected
    assert (None if stat=="" else float(actual["statistic"]))==(None if stat=="" else float(stat)),(name,"statistic",actual,expected)
    want_param=None if parameter=="" else [float(x) for x in parameter.split(",")]
    got_param=actual["parameter"]
    if got_param is not None and not isinstance(got_param,list): got_param=[float(got_param)]
    assert got_param==want_param,(name,"parameter",actual,expected)
    assert float(actual["p_value"])==float(pvalue),(name,"p_value",actual,expected)
    for key,text in (("conf_int",conf),("estimate",estimate)):
        want=None if not text else [float(x) for x in text.split(",")]
        got=actual[key]
        got=None if got is None else [float(x) for x in got]
        assert got==want,(name,key,got,want)
    assert actual["method"]==method,(name,"method",actual["method"],method)
    assert actual["alternative"]==(alternative or None),(name,"alternative",actual["alternative"],alternative)
    want_null=None if not null_value else [float(x) for x in null_value.split(",")]
    got_null=actual["null_value"]
    got_null=None if got_null is None else [float(got_null)] if not isinstance(got_null,list) else [float(x) for x in got_null]
    assert got_null==want_null,(name,"null_value",got_null,want_null)

# Added R-supported variants: one-sample signed-rank, explicit exact=false, byrow, and seeded simulated chi-square.
wilcox_one=call("wilcox-one",{"test":"wilcox","x":[1,2,3,4,6,7,8,9,10,11],"mu":5,"conf_int":True,"exact":False})
assert wilcox_one["conf_int"] is not None and wilcox_one["null_value"]==5
sim_chisq=call("chisq-sim",{"test":"chisq","x":[10,20,30,40],"rows":2,"cols":2,"byrow":True,"simulate_p_value":True,"seed":42})
assert sim_chisq["p_value"] is not None
for params,code in [
    ({"test":"shapiro","x":[1,2]},"BAD_SIZE"),
    ({"test":"binom","x":[11],"n":10},"BAD_COUNTS"),
    ({"test":"fisher","x":[1,2,3]},"BAD_DIMENSIONS"),
    ({"test":"prop","x":[2,3],"trials":[4]},"BAD_COUNTS"),
]:
    argv=[binary,"test_htest"]
    for key,value in params.items():
        for item in (value if isinstance(value,list) else [value]): argv.extend(["--"+key,str(item)])
    proc=subprocess.run(argv,text=True,capture_output=True)
    doc=json.loads(proc.stdout)
    actual_code=doc.get("error",{}).get("data",{}).get("tool_code")
    assert proc.returncode==2 and actual_code==code,(params,proc.returncode,proc.stdout)
