#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
if ! command -v grep >/dev/null 2>&1; then echo 'check-no-libR requires grep' >&2; exit 2; fi
pattern='Rmath\.h|-lRmath?([^[:alnum:]_]|$)|libR(math)?\.so|target_link_libraries.*(libR|Rmath|-lR)'
scan() {
  if grep -E -n "$pattern" "$@"; then echo 'forbidden R library/header reference found' >&2; return 1; fi
}
scan "$root"/R/* "$root"/src/* "$root"/CMakeLists.txt
# Exercise the exact grep expression against the required planted forbidden include.
if printf '%s\n' '#include <Rmath.h>' | grep -E "$pattern" >/dev/null; then :; else echo 'lint self-test failed to detect planted Rmath.h' >&2; exit 1; fi
if [ "${1:-}" ]; then
  if ! command -v ldd >/dev/null 2>&1; then echo 'ldd is required for binary linkage check' >&2; exit 2; fi
  if ldd "$1" | grep -E 'libR(math)?\.so' >/dev/null; then echo 'binary links libR/libRmath' >&2; exit 1; fi
fi
