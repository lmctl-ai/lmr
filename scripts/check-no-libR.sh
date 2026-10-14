#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
if rg --case-sensitive -n 'Rmath\.h|-lRmath?([^[:alnum:]_]|$)|libR(math)?\.so|target_link_libraries.*(libR|Rmath|-lR)' "$root/src" "$root/CMakeLists.txt"; then
  echo 'forbidden R library/header reference found' >&2
  exit 1
fi
