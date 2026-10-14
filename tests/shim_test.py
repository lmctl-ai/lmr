import json
import os
import pathlib
import subprocess
import sys
import tempfile

shim = pathlib.Path(sys.argv[1]).resolve()
home = pathlib.Path.home().resolve()
with tempfile.TemporaryDirectory(prefix="rcli-shim-", dir=home) as td:
    root = pathlib.Path(td)
    fake = root / "fakebin"
    fake.mkdir()
    created = root / "created.txt"
    args_file = root / "args.txt"
    docker = fake / "docker"
    docker.write_text("""#!/bin/sh
printf '%s\\n' "$@" >"$ARGS_FILE"
case "${FAKE_DOCKER_CASE:-ok}" in
  ok) echo '{"success":true,"result":{}}'; exit 0;;
  daemon) echo 'docker: Cannot connect to the Docker daemon' >&2; exit 125;;
  partial) echo '{'; exit 143;;
  crash) exit "${FAKE_DOCKER_RC:-139}";;
esac
""")
    docker.chmod(0o755)
    env = os.environ.copy()
    env.update(PATH=str(fake) + os.pathsep + env["PATH"], ARGS_FILE=str(args_file),
               FAKE_DOCKER_CASE="ok", RCLI_TIMEOUT_S="15", RCLI_SECRET="allowed-prefix",
               SECRET_SHIM_TEST="must-not-pass")
    colon = root / "work:area"
    colon.mkdir()
    p = subprocess.run([str(shim), "--schema"], cwd=colon, env=env, text=True, capture_output=True)
    assert p.returncode == 0 and json.loads(p.stdout)["success"] is True, (p.returncode, p.stdout, p.stderr)
    argv = args_file.read_text().splitlines()
    assert argv[argv.index("--mount") + 1] == f"type=bind,src={colon},dst={colon}", argv
    assert "--network=none" in argv and "--cpus=2" in argv and "--memory=5g" in argv
    assert "-v" not in argv and "TMPDIR=" + str(colon / ".rcli-tmp") in argv
    assert "RCLI_TIMEOUT_S" in argv and "RCLI_SECRET" not in argv and "SECRET_SHIM_TEST" not in argv
    assert not (colon / ".rcli-tmp").exists()

    for case, rc, code in (("crash", 1, "INTERNAL_ERROR"), ("partial", 1, "INTERNAL_ERROR"), ("daemon", 8, "SERVICE_UNAVAILABLE")):
        env["FAKE_DOCKER_CASE"] = case
        env["FAKE_DOCKER_RC"] = "139"
        p = subprocess.run([str(shim), "--schema"], cwd=colon, env=env, text=True, capture_output=True)
        body = json.loads(p.stdout)
        assert p.returncode == rc and body["error"]["code"] == code, (case, p.returncode, p.stdout, p.stderr)

    for unsafe in (pathlib.Path("/"), home, pathlib.Path(tempfile.gettempdir())):
        p = subprocess.run([str(shim), "--schema"], cwd=unsafe, env=env, text=True, capture_output=True)
        assert p.returncode == 2 and json.loads(p.stdout)["error"]["code"] == "INVALID_REQUEST", (unsafe, p.returncode, p.stdout)
print("rcli shim boundary, docker arguments, env filtering, cleanup and failure mapping passed")
