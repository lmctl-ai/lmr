#!/bin/sh
set -eu
bin=$1
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
"$bin" r_info > "$tmp/info"
"$bin" dist_quantile --dist norm --p '[0.975]' > "$tmp/q"
"$bin" rng_draw --dist norm --n 2 --seed 42 > "$tmp/rng"
"$bin" fit_lm --formula 'y ~ x' --data '{"y":[1,2,4,5],"x":[1,2,3,4]}' > "$tmp/lm"
"$bin" fit_lm --formula 'y ~ x' --file "$root/examples/data/oracle-lm.csv" > "$tmp/lm-csv"
cmp "$tmp/lm" "$tmp/lm-csv"
Rscript --vanilla --default-packages=NULL "$root/tests/oracle.R" "$tmp/expected"
python3 "$root/tests/compare.py" "$tmp/info" "$tmp/q" "$tmp/rng" "$tmp/lm" "$tmp/expected"
"$bin" fit_lm --formula 'y ~ x' --data '{"y":[1,3,5,7],"x":[1,2,3,4]}' > "$tmp/warn"
python3 -c 'import json,sys; assert json.load(open(sys.argv[1]))["warnings"] == ["essentially perfect fit: summary may be unreliable"]' "$tmp/warn"
"$bin" fit_lm --formula 'y ~ x' --data '{"y":[1,2],"x":[1,2]}' > "$tmp/nonfinite"
python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); assert d["std_error"] == ["NaN","NaN"] and d["sigma"] == ["NaN"]' "$tmp/nonfinite"
if "$bin" fit_lm --formula 'system(x)' --data '{}' > /dev/null; then exit 1; fi
if "$bin" dist_quantile --dist unknown --p '[0.5]' > /dev/null; then exit 1; fi
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
chmod +x "$tmp/fake-env" "$tmp/fake-timeout"
RCLI_RSCRIPT="$tmp/fake-env" "$bin" r_info > "$tmp/env"
grep -q '"ok":true' "$tmp/env"
if RCLI_RSCRIPT="$tmp/fake-timeout" RCLI_TIMEOUT_S=1 "$bin" r_info > "$tmp/timeout"; then exit 1; fi
grep -q 'TIMEOUT' "$tmp/timeout"
"$bin" dist_quantile --dist norm --p '[0,1]' > "$tmp/infinities"
python3 -c 'import json,sys; assert json.load(open(sys.argv[1]))["value"] == ["-Inf","Inf"]' "$tmp/infinities"
"$bin" fit_lm --formula 'y ~ x' --na-action exclude --residuals --data '{"y":[1,null,4,5],"x":[1,2,3,4]}' > "$tmp/na"
python3 -c 'import json,sys; assert "NA" in json.load(open(sys.argv[1]))["residuals"]' "$tmp/na"
if "$bin" fit_lm --formula 'y ~ absent' --data '{"y":[1,2,4,5],"x":[1,2,3,4]}' > "$tmp/r-error"; then exit 1; fi
python3 -c 'import json,sys; d=json.load(open(sys.argv[1]))["error"]; assert d["code"] == "COMMAND_FAILED" and d["data"]["tool_code"] == "R_ERROR" and d["data"]["r_message"]' "$tmp/r-error"
python3 - "$tmp/large.json" <<'PY'
import json,sys
json.dump({"op":"dist_quantile","params":{"dist":"norm","p":[0.5]*200001}},open(sys.argv[1],"w"))
PY
if RCLI_RSCRIPT="$tmp/fake-env" "$bin" --input-file "$tmp/large.json" > "$tmp/large-error"; then exit 1; fi
grep -q 'INPUT_TOO_LARGE' "$tmp/large-error"
"$bin" dist_quantile --dist t --df 7 --p '[0.75]' > "$tmp/t"
"$bin" dist_quantile --dist chisq --df 7 --p '[0.75]' > "$tmp/chisq"
"$bin" dist_quantile --dist gamma --shape 2 --rate 3 --p '[0.75]' > "$tmp/gamma"
"$bin" dist_quantile --dist beta --a 2 --b 5 --p '[0.75]' > "$tmp/beta"
"$bin" dist_quantile --dist binom --size 10 --prob 0.4 --p '[0.75]' > "$tmp/binom"
"$bin" dist_quantile --dist pois --lambda 3 --p '[0.75]' > "$tmp/pois"
Rscript --vanilla --default-packages=NULL "$root/tests/dist-oracle.R" "$tmp/dist-expected"
python3 "$root/tests/compare-dist.py" "$tmp/dist-expected" "$tmp/t" "$tmp/chisq" "$tmp/gamma" "$tmp/beta" "$tmp/binom" "$tmp/pois"
