#!/bin/sh
set -eu
bin=$1
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
tmp=$(mktemp -d)
control_file=
cap_file=
row_file=
cleanup(){ rm -rf "$tmp"; [ -z "$control_file" ] || rm -f "$control_file"; [ -z "$cap_file" ] || rm -f "$cap_file"; [ -z "$row_file" ] || rm -f "$row_file"; }
trap cleanup EXIT
"$bin" --schema > "$tmp/schema"
python3 - "$tmp/schema" <<'PY'
import json,sys
s=json.load(open(sys.argv[1])); ops={o['name']:o for o in s['operations']}
assert set(ops)=={'r_info','dist_quantile','dist_density','dist_probability','rng_draw','rng_sample','math_special','summary_quantile','num_fft','test_htest','fit_glm','smooth_spline','smooth_loess','interp_spline','interp_approx','ts_arima','rng_arima_sim','ts_acf','ts_pacf','ts_stl','cluster_kmeans','cluster_hclust','fit_lm'}
assert all(o['access']=='read' for o in ops.values())
p={x['name']:x for x in ops['dist_quantile']['params']}
assert p['p']['type']=='array' and p['p']['items']=='number' and 'log_p' in p
PY
"$bin" r_info > "$tmp/info"
"$bin" dist_quantile --dist norm --p 0.975 > "$tmp/q"
"$bin" rng_draw --dist norm --n 2 --seed 42 > "$tmp/rng"
"$bin" fit_lm --formula 'y ~ x' --data '{"y":[1,2,4,5],"x":[1,2,3,4]}' > "$tmp/lm"
if "$bin" fit_lm --formula 'y ~ letters' --data '{"y":[1,2],"x":[1,2]}' > "$tmp/formula-env"; then exit 1; fi
if "$bin" fit_lm --formula 'y ~ .' --data '{"y":[1,2],"x":[1,2]}' > "$tmp/formula-dot"; then exit 1; fi
if "$bin" fit_lm --formula 'y ~ x' --data '{"y":[1,2],"x":[1]}' > "$tmp/shape-error"; then exit 1; fi
grep -q 'BAD_DATA' "$tmp/shape-error"
for term in args pi T NULL TRUE; do
  if "$bin" fit_lm --formula "y ~ $term" --data '{"y":[1,2],"x":[1,2]}' > "$tmp/formula-env"; then exit 1; fi
done
"$bin" fit_lm --formula 'y ~ 1' --data '{"y":[1,2,3],"x":[1,2,3]}' > "$tmp/intercept"
"$bin" fit_lm --formula 'y ~ x - 1' --data '{"y":[2,4,6],"x":[1,2,3]}' > "$tmp/no-intercept"
"$bin" fit_lm --formula 'y ~ x + 0' --data '{"y":[2,4,6],"x":[1,2,3]}' > "$tmp/no-intercept-zero"
"$bin" fit_lm --formula 'y ~ x + z' --data '{"y":[1,2,4,5],"x":[1,2,3,4],"z":[1,2,3,4]}' > "$tmp/aliased"
python3 - "$tmp/intercept" "$tmp/no-intercept" "$tmp/aliased" "$tmp/no-intercept-zero" <<'PY'
import json,sys
a,b,c,d=(json.load(open(p)) for p in sys.argv[1:])
assert isinstance(a['estimate'],list) and len(a['estimate'])==1
assert isinstance(b['estimate'],list) and len(b['estimate'])==1
assert d['coefficients']==b['coefficients'] and d['estimate']==b['estimate']
assert c['rank']==2 and c['coefficients']==['(Intercept)','x','z'] and c['aliased']==['z'] and c['estimate'][2]=='NA'
assert isinstance(c['std_error'],list) and len(c['std_error'])==3
PY
"$bin" fit_lm --formula 'y ~ x' --file "$root/examples/data/oracle-lm.csv" > "$tmp/lm-csv"
cmp "$tmp/lm" "$tmp/lm-csv"
control_file=$(python3 -c 'p="examples/data/oracle-\x01.csv";open(p,"w").write("y,x\n1,1\n2,2\n4,3\n5,4\n");print(p)')
"$bin" fit_lm --formula 'y ~ x' --file "$control_file" > "$tmp/control-path"
cmp "$tmp/lm" "$tmp/control-path"
Rscript --vanilla --default-packages=NULL "$root/tests/oracle.R" "$tmp/expected"
python3 "$root/tests/compare.py" "$tmp/info" "$tmp/q" "$tmp/rng" "$tmp/lm" "$tmp/expected"
"$bin" fit_lm --formula 'y ~ x' --data '{"y":[1,3,5,7],"x":[1,2,3,4]}' > "$tmp/warn"
python3 -c 'import json,sys; assert json.load(open(sys.argv[1]))["warnings"] == ["essentially perfect fit: summary may be unreliable"]' "$tmp/warn"
"$bin" fit_lm --formula 'y ~ x' --data '{"y":[1,2],"x":[1,2]}' > "$tmp/nonfinite"
python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); assert d["std_error"] == ["NaN","NaN"] and d["sigma"] == "NaN"' "$tmp/nonfinite"
if "$bin" fit_lm --formula 'system(x)' --data '{}' > /dev/null; then exit 1; fi
if "$bin" dist_quantile --dist unknown --p 0.5 > /dev/null; then exit 1; fi
cat > "$tmp/fake-env" <<'SH'
#!/bin/sh
for last do :; done
if [ "${R_PROFILE_USER:-}" = /dev/null ] && [ "${R_ENVIRON_USER:-}" = /dev/null ] && [ "${R_LIBS_USER+x}" = x ] && [ -z "$R_LIBS_USER" ] && [ "${HOME:-}" = /tmp ]; then printf '{"ok":true}' > "$last"; else printf '{"ok":false}' > "$last"; fi
SH
cat > "$tmp/fake-timeout" <<'SH'
#!/bin/sh
for last do :; done
sleep 3
printf '{}' > "$last"
SH
cat > "$tmp/fake-crash" <<'SH'
#!/bin/sh
echo '/tmp/private-rcli-path' >&2
exit 42
SH
cat > "$tmp/fake-r-error" <<'SH'
#!/bin/sh
for last do :; done
printf '{"error":"synthetic failure","warnings":["synthetic warning"]}' > "$last"
SH
chmod +x "$tmp/fake-env" "$tmp/fake-timeout" "$tmp/fake-crash" "$tmp/fake-r-error"
RCLI_RSCRIPT="$tmp/fake-env" "$bin" r_info > "$tmp/env"
grep -q '"ok":true' "$tmp/env"
if RCLI_RSCRIPT="$tmp/fake-timeout" RCLI_TIMEOUT_S=1 "$bin" r_info > "$tmp/timeout"; then exit 1; fi
grep -q 'TIMEOUT' "$tmp/timeout"
if RCLI_RSCRIPT="$tmp/no-such-rscript" "$bin" r_info > "$tmp/no-rscript"; then exit 1; fi
grep -q 'R_NOT_FOUND' "$tmp/no-rscript"
if RCLI_RSCRIPT="$tmp/fake-crash" "$bin" r_info > "$tmp/crash"; then exit 1; fi
python3 - "$tmp/crash" <<'PY'
import json,sys
d=json.load(open(sys.argv[1]))['error']; assert d['data']['tool_code']=='TOOL_CRASH' and 'stderr' not in d['data'] and 'private-rcli-path' not in open(sys.argv[1]).read()
PY
if RCLI_RSCRIPT="$tmp/fake-r-error" "$bin" r_info > "$tmp/r-error-warning"; then exit 1; fi
python3 - "$tmp/r-error-warning" <<'PY'
import json,sys
d=json.load(open(sys.argv[1]))['error']; assert d['data']['warnings']==['synthetic warning']
PY
mkdir "$tmp/private-tmp"
TMPDIR="$tmp/private-tmp" "$bin" r_info > "$tmp/tmpdir-info"
test -z "$(find "$tmp/private-tmp" -type f -print -quit)"
"$bin" dist_quantile --dist norm --p 0 --p 1 > "$tmp/infinities"
python3 -c 'import json,sys; assert json.load(open(sys.argv[1]))["value"] == ["-Inf","Inf"]' "$tmp/infinities"
"$bin" dist_quantile --dist norm --p NA --p NaN --p Inf --p=-Inf > "$tmp/special-p"
python3 -c 'import json,sys; assert json.load(open(sys.argv[1]))["value"] == ["NA","NaN","NaN","NaN"]' "$tmp/special-p"
"$bin" rng_draw --dist norm --n 1 --seed 42 --mean NA > "$tmp/special-mean"
python3 -c 'import json,sys; assert json.load(open(sys.argv[1]))["value"] == ["NaN"]' "$tmp/special-mean"
python3 - "$tmp/null-p.json" <<'PY'
import json,sys
json.dump({'op':'dist_quantile','params':{'dist':'norm','p':[None]}},open(sys.argv[1],'w'))
PY
"$bin" --input-file "$tmp/null-p.json" > "$tmp/null-p"
python3 -c 'import json,sys; assert json.load(open(sys.argv[1]))["value"] == ["NA"]' "$tmp/null-p"
"$bin" fit_lm --formula 'y ~ x' --na-action exclude --residuals --data '{"y":[1,null,4,5],"x":[1,2,3,4]}' > "$tmp/na"
python3 -c 'import json,sys; assert "NA" in json.load(open(sys.argv[1]))["residuals"]' "$tmp/na"
"$bin" fit_lm --formula 'y ~ x' --na-action exclude --data '{"y":[1,"NA",4,5],"x":[1,2,3,4]}' > "$tmp/na-string"
python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); assert len(d["estimate"])==2 and d["warnings"]==[]' "$tmp/na-string"
if "$bin" fit_lm --formula 'y ~ absent' --data '{"y":[1,2,4,5],"x":[1,2,3,4]}' > "$tmp/r-error"; then exit 1; fi
python3 -c 'import json,sys; d=json.load(open(sys.argv[1]))["error"]; assert d["code"] == "VALIDATION_ERROR" and d["data"]["tool_code"] == "BAD_FORMULA"' "$tmp/r-error"
python3 - "$tmp/large.json" <<'PY'
import json,sys
json.dump({"op":"dist_quantile","params":{"dist":"norm","p":[0.5]*200001}},open(sys.argv[1],"w"))
PY
if RCLI_RSCRIPT="$tmp/fake-env" "$bin" --input-file "$tmp/large.json" > "$tmp/large-error"; then exit 1; fi
grep -q 'INPUT_TOO_LARGE' "$tmp/large-error"
if "$bin" rng_draw --dist norm --n 1 --seed 2147483648 > "$tmp/seed-error"; then exit 1; fi
grep -q 'BAD_SEED' "$tmp/seed-error"
python3 - "$tmp/empty.json" <<'PY'
import json,sys
json.dump({'op':'dist_quantile','params':{'dist':'norm','p':[]}},open(sys.argv[1],'w'))
PY
if "$bin" --input-file "$tmp/empty.json" > "$tmp/empty-error"; then exit 1; fi
grep -q 'EMPTY_PROBABILITY' "$tmp/empty-error"
cap_file=$(python3 -c 'p="examples/data/rcli-over-cap.csv";f=open(p,"w");f.write("y,x\n");f.truncate(16777217);f.close();print(p)')
if "$bin" fit_lm --formula 'y ~ x' --file "$cap_file" > "$tmp/csv-error"; then exit 1; fi
grep -q 'INPUT_TOO_LARGE' "$tmp/csv-error"
row_file=$(python3 -c 'p="examples/data/rcli-over-rows.csv";f=open(p,"w");f.write("y,x\n"+"1,1\n"*200001);f.close();print(p)')
if "$bin" fit_lm --formula 'y ~ x' --file "$row_file" > "$tmp/csv-rows-error"; then exit 1; fi
grep -q 'INPUT_TOO_LARGE' "$tmp/csv-rows-error"
"$bin" dist_quantile --dist t --df 7 --p 0.75 > "$tmp/t"
"$bin" dist_quantile --dist chisq --df 7 --p 0.75 > "$tmp/chisq"
"$bin" dist_quantile --dist gamma --shape 2 --rate 3 --p 0.75 > "$tmp/gamma"
"$bin" dist_quantile --dist beta --a 2 --b 5 --p 0.75 > "$tmp/beta"
"$bin" dist_quantile --dist binom --size 10 --prob 0.4 --p 0.75 > "$tmp/binom"
"$bin" dist_quantile --dist pois --lambda 3 --p 0.75 > "$tmp/pois"
Rscript --vanilla --default-packages=NULL "$root/tests/dist-oracle.R" "$tmp/dist-expected"
python3 "$root/tests/compare-dist.py" "$tmp/dist-expected" "$tmp/t" "$tmp/chisq" "$tmp/gamma" "$tmp/beta" "$tmp/binom" "$tmp/pois"
mkdir "$tmp/stage1a-oracle"
python3 "$root/tests/stage1a.py" "$bin" "$tmp/stage1a-oracle"
python3 "$root/tests/compare-htest.py" "$bin" "$tmp/stage1a-oracle"
python3 "$root/tests/compare-group3.py" "$bin" "$tmp/stage1a-oracle"
python3 "$root/tests/compare-group4.py" "$bin" "$tmp/stage1a-oracle"
"$bin" dist_quantile --dist norm --log-p --p -0.025317807984289897 > "$tmp/log-p"
python3 -c 'import json,sys,math; assert abs(json.load(open(sys.argv[1]))["value"][0]-1.9599639845400536)<1e-14' "$tmp/log-p"
python3 - "$bin" "$tmp" <<'PYTIME'
import json,subprocess,sys,time
for n in (20000,200000):
 request=f'{sys.argv[2]}/large-{n}.json'; output=f'{sys.argv[2]}/large-{n}.out'
 with open(request,'w') as f: json.dump({'op':'dist_quantile','params':{'dist':'norm','p':[0.5]*n}},f,separators=(',',':'))
 t=time.perf_counter(); subprocess.run([sys.argv[1],'--input-file',request],stdout=open(output,'w'),check=True); elapsed=time.perf_counter()-t
 data=json.load(open(output)); assert len(data['value'])==n and all(v==0 for v in data['value'])
 print(f'rcli {n} values elapsed={elapsed:.3f} sec')
 if n==200000: assert elapsed<5.0, f'200000-value regression exceeded 5 seconds: {elapsed:.3f}'
PYTIME
