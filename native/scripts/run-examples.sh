#!/usr/bin/env bash
# Runs every examples/dsl/*.lms against the given rcli-dsl (cwd must be the repository root) and reports each result.
# Needs a DSL v5 lmscript compiler on PATH or in $LMSCRIPT_BIN; without one the test is skipped (exit 0 with a notice).
set -uo pipefail
BIN="${1:?usage: run-examples.sh path/to/rcli-dsl [example.lms ...]}"
shift
LMSCRIPT="${LMSCRIPT_BIN:-$(command -v lmscript || true)}"
if [[ -z "$LMSCRIPT" ]]; then echo "lmscript not found: DSL examples skipped"; exit 0; fi
export RCLI_DSL="$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN")"
files=("$@")
[[ ${#files[@]} -gt 0 ]] || files=(examples/dsl/*.lms)
pass=0; fail=0; failed=()
for f in "${files[@]}"; do
  start=$(date +%s.%N)
  if out="$("$LMSCRIPT" "$f" 2>&1)"; then
    pass=$((pass + 1)); status=ok
  else
    fail=$((fail + 1)); status=FAIL; failed+=("$f")
  fi
  printf '%-4s %6.2fs  %s\n' "$status" "$(echo "$(date +%s.%N) - $start" | bc)" "$f"
  [[ $status == ok ]] || printf '%s\n' "$out" | tail -15 | sed 's/^/     | /'
done
echo "examples: $pass passed, $fail failed"
[[ $fail -eq 0 ]]
