#!/usr/bin/env bash
# Regenerate schema/lmr.dsl.json from the host registry (after changing a bind_*.cpp), then rebuild so the binary
# embeds the new bytes. Bytecode compiled against the old schema must be recompiled (its SHA-256 binding changes).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${BUILD_DIR:-$HERE/../build-native}"
SCHEMA="$HERE/schema/lmr.dsl.json"
[[ -f "$SCHEMA" ]] || printf '{}\n' > "$SCHEMA"
cmake -S "$HERE" -B "$BUILD" -G Ninja -DCMAKE_BUILD_TYPE=Release >/dev/null
cmake --build "$BUILD" --target lmr
"$BUILD/bin/lmr" --generate-schema > "$SCHEMA.tmp"
mv "$SCHEMA.tmp" "$SCHEMA"
cmake --build "$BUILD" --target lmr
"$HERE/scripts/check-schema.sh" "$BUILD/bin/lmr"
