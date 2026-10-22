#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$root"

node_bin="$(command -v node)"
lmscript_bin="$(command -v lmscript)"
[[ -n "$node_bin" && -n "$lmscript_bin" ]] || { echo "error: install Node.js and @lmctl-ai/lmctl" >&2; exit 2; }

tmp="$(mktemp -d "${TMPDIR:-/tmp}/lmr-readme.XXXXXX")"
cleanup() { rm -rf "$tmp"; }
trap cleanup EXIT
if [[ -n "${LMSCRIPT_TOOL_CACHE:-}" ]]; then
  cache="$LMSCRIPT_TOOL_CACHE"
else
  cache="$tmp/tool-cache"
fi
mkdir -p "$cache" "$tmp/bin"
export LMSCRIPT_TOOL_CACHE="$cache"
rm -f docs/readme/hero-light.png docs/readme/hero-dark.png docs/readme/architecture.png \
  docs/readme/normal-quantile.png docs/readme/linear-fit.png docs/readme/rosenbrock-minimum.png

# Hide any globally installed lmr/lmdraw so the first render proves the
# published install= declarations on a fresh tool cache.
clean_path=""
IFS=: read -r -a path_entries <<< "${PATH:-}"
for entry in "${path_entries[@]}"; do
  [[ -n "$entry" ]] || continue
  if [[ -x "$entry/lmr" || -x "$entry/lmr.exe" || -x "$entry/lmdraw" || -x "$entry/lmdraw.exe" ]]; then
    continue
  fi
  clean_path="${clean_path:+$clean_path:}$entry"
done
export PATH="$clean_path"

# Hero light installs both released hosts into the isolated cache. Reuse those
# exact binaries for the other variants while preserving install= in each file.
"$node_bin" "$lmscript_bin" docs/readme/hero.lms light
for tool in lmr lmdraw; do
  binary="$cache/github/lmctl-ai/$tool/latest/$tool"
  [[ -x "$binary" ]] || { echo "error: install did not produce $binary" >&2; exit 1; }
  ln -sfn "$binary" "$tmp/bin/$tool"
done
"$tmp/bin/lmr" --version
"$tmp/bin/lmdraw" --version
export PATH="$tmp/bin:$PATH"

for source in docs/readme/hero.lms docs/readme/architecture.lms docs/readme/gallery.lms; do
  "$node_bin" "$lmscript_bin" --check "$source"
done

"$node_bin" "$lmscript_bin" docs/readme/hero.lms dark
"$node_bin" "$lmscript_bin" docs/readme/architecture.lms
"$node_bin" "$lmscript_bin" docs/readme/gallery.lms

python3 - <<'PY'
from pathlib import Path

images = [
    Path("docs/readme/hero-light.png"),
    Path("docs/readme/hero-dark.png"),
    Path("docs/readme/architecture.png"),
    Path("docs/readme/normal-quantile.png"),
    Path("docs/readme/linear-fit.png"),
    Path("docs/readme/rosenbrock-minimum.png"),
]
for image in images:
    if not image.is_file() or image.stat().st_size == 0:
        raise SystemExit(f"missing generated README image: {image}")
    if image.stat().st_size > 400_000:
        raise SystemExit(f"README image exceeds 400 KB: {image} ({image.stat().st_size} bytes)")
    if image.read_bytes()[:8] != b"\x89PNG\r\n\x1a\n":
        raise SystemExit(f"README image is not a PNG: {image}")
print(f"README images ready: {len(images)} PNG files, each under 400 KB")
PY
