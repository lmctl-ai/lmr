#!/usr/bin/env bash
# Vendor lmcode-cpp (formerly cpp-dsl-sdk): the single-source runtime header, the built-in std catalog it includes
# (stdlib/*.hpp) and the canonical std schema, pinning the commit and every file's content hash in sdk.lock.
#   scripts/sync-sdk.sh            copy from a checkout (LMCODE_CPP, or ../lmcode-cpp, or ../poc/cpp-dsl-sdk) at origin/main
#   scripts/sync-sdk.sh --verify   fail when a vendored file differs from sdk.lock
set -eu
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOCK="$HERE/sdk.lock"
FILES="include/dslvm.hpp:cpp/dslvm.hpp
stdlib/std.hpp:cpp/stdlib/std.hpp
stdlib/string.hpp:cpp/stdlib/string.hpp
stdlib/math.hpp:cpp/stdlib/math.hpp
stdlib/collections.hpp:cpp/stdlib/collections.hpp
stdlib/file.hpp:cpp/stdlib/file.hpp
std/lmcode-std.schema.json:cpp/std/lmcode-std.schema.json"
hash() { sha256sum "$1" | cut -c1-64; }
if [ "${1:-}" = "--verify" ]; then
  [ -f "$LOCK" ] || { echo "sdk.lock missing" >&2; exit 1; }
  bad=0
  while IFS=: read -r _ dst; do
    want="$(sed -n "s|^sha256 $dst=||p" "$LOCK")"
    have="$(hash "$HERE/$dst" 2>/dev/null || echo missing)"
    [ -n "$want" ] && [ "$want" = "$have" ] || { echo "$dst sha256 $have differs from sdk.lock '${want}': run scripts/sync-sdk.sh" >&2; bad=1; }
  done <<< "$FILES"
  [ "$bad" = 0 ] || exit 1
  echo "vendored lmcode-cpp files match sdk.lock ($(sed -n 's/^commit=//p' "$LOCK"))"
  exit 0
fi
SDK="${LMCODE_CPP:-}"
if [ -z "$SDK" ]; then
  for c in "$HERE/../../lmcode-cpp" "$HERE/../../poc/lmcode-cpp" "$HERE/../../poc/cpp-dsl-sdk"; do
    [ -f "$c/include/dslvm.hpp" ] && { SDK="$c"; break; }
  done
fi
[ -n "$SDK" ] && [ -f "$SDK/include/dslvm.hpp" ] || { echo "lmcode-cpp checkout not found (set LMCODE_CPP)" >&2; exit 1; }
ref="${LMCODE_CPP_REF:-origin/main}"
commit="$(git -C "$SDK" rev-parse "$ref")"
{
  printf 'repo=https://github.com/lmctl-ai/lmcode-cpp\ncommit=%s\n' "$commit"
  while IFS=: read -r src dst; do
    mkdir -p "$(dirname "$HERE/$dst")"
    git -C "$SDK" show "$commit:$src" > "$HERE/$dst"
    printf 'sha256 %s=%s\n' "$dst" "$(hash "$HERE/$dst")"
  done <<< "$FILES"
} > "$LOCK.tmp"
mv "$LOCK.tmp" "$LOCK"
echo "synced lmcode-cpp $commit"
