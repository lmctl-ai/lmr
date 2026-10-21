#!/usr/bin/env bash
# Re-runs every R recipe (r/*.R) to regenerate data/*.csv and expected/*.json. Needs Rscript; the committed results
# come from R 4.3.3. The lmr side (*.lms) never needs R: it reads the committed files.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
command -v Rscript >/dev/null || { echo "Rscript not found" >&2; exit 1; }
Rscript --version 2>&1 | grep -q "4\.3\.3" || echo "note: the committed expected/*.json came from R 4.3.3; this is $(Rscript --version 2>&1)" >&2
Rscript r/export-data.R
for f in r/*.R; do
  case "$f" in r/_emit.R|r/export-data.R) continue ;; esac
  Rscript "$f" > /dev/null
  echo "ok  $f"
done
