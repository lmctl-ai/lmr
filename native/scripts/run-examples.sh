#!/usr/bin/env bash
# Runs every examples/dsl/*.lms and examples/r-vs-dsl/*.lms against the given lmr (cwd must be the repository root)
# and reports each result.
# Needs a DSL v5 lmscript compiler on PATH or in $LMSCRIPT_BIN; without one the test is skipped (exit 0 with a notice).
set -uo pipefail
BIN="${1:?usage: run-examples.sh path/to/lmr [example.lms ...]}"
shift
LMSCRIPT="${LMSCRIPT_BIN:-$(command -v lmscript || true)}"
if [[ -z "$LMSCRIPT" ]]; then echo "lmscript not found: DSL examples skipped"; exit 0; fi
export LMR_BIN="$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN")"
# First on PATH too, so `_dsl_ cli="lmr"` (install-and-run.lms) binds to this build instead of downloading a release.
export PATH="$(dirname "$LMR_BIN"):$PATH"
files=("$@")
[[ ${#files[@]} -gt 0 ]] || files=(examples/dsl/*.lms examples/r-vs-dsl/*.lms)
pass=0; fail=0; failed=()
for f in "${files[@]}"; do
  start=$SECONDS
  if out="$("$LMSCRIPT" "$f" 2>&1)"; then
    pass=$((pass + 1)); status=ok
  else
    fail=$((fail + 1)); status=FAIL; failed+=("$f")
  fi
  printf '%-4s %4ds  %s\n' "$status" "$((SECONDS - start))" "$f"
  [[ $status == ok ]] || printf '%s\n' "$out" | tail -15 | sed 's/^/     | /'
done
echo "examples: $pass passed, $fail failed"
[[ $fail -eq 0 ]]
