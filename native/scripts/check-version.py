#!/usr/bin/env python3
"""Require exact CLI version bytes and, on a tag build, equality with the release tag."""
import os
from pathlib import Path
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parents[1]
VERSION = re.search(r"^project\(lmr VERSION ([^ )]+)", (HERE / "CMakeLists.txt").read_text(), re.M).group(1)
binary = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE.parent / "build-native" / "bin" / ("lmr.exe" if os.name == "nt" else "lmr")
if os.name == "nt" and binary.suffix.lower() != ".exe":
    binary = binary.with_suffix(".exe")
result = subprocess.run([str(binary), "--version"], capture_output=True)
expected = f"lmr {VERSION}\n".encode()
if result.returncode != 0 or result.stdout != expected or result.stderr:
    raise SystemExit(f"expected --version stdout {expected!r}, empty stderr and exit 0; got {result.stdout!r}, {result.stderr!r}, exit {result.returncode}")

ref = os.environ.get("GITHUB_REF", "")
if ref.startswith("refs/tags/v"):
    tag_version = ref.removeprefix("refs/tags/v")
    if result.stdout != f"lmr {tag_version}\n".encode():
        raise SystemExit(f"release tag {ref.rsplit('/', 1)[-1]} does not match binary output {result.stdout!r}")
print(f"version guard: {result.stdout.decode().strip()}")
