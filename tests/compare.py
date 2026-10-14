import json, sys
info=json.load(open(sys.argv[1])); q=json.load(open(sys.argv[2])); rng=json.load(open(sys.argv[3])); fit=json.load(open(sys.argv[4])); expected=open(sys.argv[5]).read().splitlines(); vals=list(map(float,expected[2:]))
assert info['r_version'] == expected[0]
assert info['rng_kind'] == expected[1].split('|')
i=0
assert q['value'][0] == vals[i]; i+=1
assert rng['value'] == vals[i:i+2]; i+=2
for key in ('estimate','std_error','statistic','p_value'):
    assert fit[key] == vals[i:i+2], (key,fit[key],vals[i:i+2]); i+=2
for key in ('sigma','r_squared','df'):
    got=fit[key]
    if isinstance(got,list): got=got[0]
    assert got == vals[i], (key,got,vals[i]); i+=1
assert i == len(vals)
