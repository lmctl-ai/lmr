#!/usr/bin/env python3
"""Run the three install-enabled public book examples from isolated working dirs."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
EXPECTED = {
    "normal-quantile": "Normal 97.5% quantile: 1.9599639845\n",
    "linear-fit": "Line fit: intercept 1.0, slope 2.0\n",
    "rosenbrock-minimum": "BFGS minimum: a 0.999804, b 0.999608, value 0.00000004\n",
}
compiler = os.environ.get("LMSCRIPT_BIN", "lmscript")


def compiler_command(script):
    if os.name != "nt":
        executable = shutil.which(compiler)
        if not executable:
            raise SystemExit(f"lmscript compiler not found: {compiler}")
        entry = Path(executable).resolve()
        # lmr may share its PATH directory with both the npm shim and Node.
        # Resolve both before removing that directory from the child's PATH.
        with entry.open("rb") as source:
            shebang = source.readline(256).strip()
        if shebang == b"#!/usr/bin/env node":
            node = shutil.which("node")
            if not node:
                raise SystemExit("lmscript needs Node.js")
            return [str(Path(node).resolve()), str(entry), str(script)]
        return [str(entry), str(script)]
    # npm installs Windows .cmd shims, which CreateProcess cannot launch directly.
    node = shutil.which("node")
    root = subprocess.run("npm root -g", shell=True, check=True, capture_output=True, text=True).stdout.strip()
    entry = Path(root) / "@lmctl-ai" / "lmctl" / "bin" / "lmscript"
    if not node or not entry.is_file():
        raise SystemExit("Windows book checks need Node.js and the globally installed @lmctl-ai/lmctl entry point")
    return [node, str(entry), str(script)]


def run_examples(candidate=None):
    with tempfile.TemporaryDirectory(prefix="lmr-tools-") as cache:
        for name, expected in EXPECTED.items():
            with tempfile.TemporaryDirectory(prefix=f"lmr-book-{name}-") as work:
                env = {**os.environ, "LMSCRIPT_TOOL_CACHE": cache}
                if candidate:
                    env["PATH"] = str(candidate.parent) + os.pathsep + env.get("PATH", "")
                else:
                    env.pop("LMR_BIN", None)
                    # Freeze relative/empty entries before changing cwd, and
                    # exclude Windows command shims as well as native binaries.
                    entries = [Path(entry).resolve() for entry in env.get("PATH", "").split(os.pathsep)]
                    names = ("lmr", "lmr.exe", "lmr.cmd", "lmr.bat", "lmr.com") if os.name == "nt" else ("lmr",)
                    env["PATH"] = os.pathsep.join(
                        str(entry) for entry in entries
                        if not any((entry / name).is_file() for name in names)
                    )
                result = subprocess.run(
                    compiler_command(ROOT / "examples" / f"{name}.lms"),
                    cwd=work, env=env, capture_output=True, timeout=180,
                )
                if result.returncode:
                    detail = result.stdout.decode(errors="replace") + result.stderr.decode(errors="replace")
                    raise SystemExit(f"{name} failed ({result.returncode}):\n{detail}")
                if result.stderr:
                    raise SystemExit(f"{name} wrote unexpected stderr: {result.stderr!r}")
                if result.stdout != expected.encode():
                    raise SystemExit(f"{name} output mismatch:\nexpected {expected.encode()!r}\nactual   {result.stdout!r}")
                artifact = Path(work) / "build" / "rendered" / f"{name}.txt"
                if not artifact.is_file() or artifact.read_bytes() != expected.encode():
                    raise SystemExit(f"{name} rendered text missing or different from expected output")
                if candidate is None:
                    installed = Path(cache) / "github" / "lmctl-ai" / "lmr" / "latest" / ("lmr.exe" if os.name == "nt" else "lmr")
                    if not installed.is_file():
                        raise SystemExit("install check bypassed: remove lmr from PATH to exercise release installation")
                print(f"PASS {name}: {result.stdout.decode().rstrip()}")


# Check the candidate before the network-dependent release installation.
candidate = Path(os.environ.get("LMR_BIN", ROOT / "build-native" / "bin" / ("lmr.exe" if os.name == "nt" else "lmr"))).resolve()
if not candidate.is_file():
    raise SystemExit(f"candidate lmr binary not found: {candidate}")
print(f"Candidate host: {candidate}")
run_examples(candidate)

# Prove release discovery, checksums, and execution from an empty tool cache.
print("Released @latest install:")
run_examples()
