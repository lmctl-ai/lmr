#!/usr/bin/env bash
# The committed schema must equal what the registry generates, the embedded bytes must equal the file, and (when a
# DSL v5 compiler is available) the compiler's canonical signatures must equal the host registry's.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${1:?usage: check-schema.sh path/to/rcli-dsl}"
SCHEMA="$HERE/schema/rcli.dsl.json"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/rcli-schema-XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
"$BIN" --generate-schema > "$WORK/generated.json"
cmp "$SCHEMA" "$WORK/generated.json" || { echo "schema/rcli.dsl.json is stale: run scripts/update-schema.sh" >&2; exit 1; }
"$BIN" --schema --json > "$WORK/embedded.json"
cmp "$SCHEMA" "$WORK/embedded.json" || { echo "the binary embeds an older schema: rebuild" >&2; exit 1; }
LMSCRIPT="${LMSCRIPT_BIN:-$(command -v lmscript || true)}"
if [[ -n "$LMSCRIPT" ]]; then
  "$LMSCRIPT" --dsl-signatures "$SCHEMA" > "$WORK/compiler.txt"
  "$BIN" --dsl-signatures > "$WORK/host.txt"
  diff -u "$WORK/compiler.txt" "$WORK/host.txt" || { echo "compiler and host signatures differ" >&2; exit 1; }
  echo "schema ok: $(wc -l < "$WORK/host.txt") signatures agree with $LMSCRIPT"
else
  echo "schema ok (lmscript not found: compiler signature parity skipped)"
fi
