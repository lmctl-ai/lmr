import json
import pathlib
import subprocess
import sys

binary, expected_dir = sys.argv[1:]
expected_dir = pathlib.Path(expected_dir)
subprocess.run(["Rscript", "--vanilla", "--default-packages=NULL", "tests/stage1a-oracle.R", str(expected_dir)], check=True)

def run(op, params, key="value", expected=None):
    argv = [binary, op]
    for name, value in params.items():
        values = value if isinstance(value, list) else [value]
        for item in values:
            if isinstance(item, bool):
                argv.append(("--" if item else "--no-") + name)
            else:
                argv.extend(["--" + name, str(item)])
    got = json.loads(subprocess.run(argv, check=True, text=True, capture_output=True).stdout)[key]
    if expected is not None:
        want = [float(x) for x in (expected_dir / (expected + ".txt")).read_text().splitlines()]
        actual = [float(x) for x in got]
        assert actual == want, (op, params, actual, want)
    return got

def rejects(op, params, code):
    argv = [binary, op]
    for name, value in params.items():
        values = value if isinstance(value, list) else [value]
        for item in values:
            argv.extend(["--" + name, str(item)])
    proc = subprocess.run(argv, text=True, capture_output=True)
    assert proc.returncode == 2, (argv, proc.returncode, proc.stdout, proc.stderr)
    assert json.loads(proc.stdout)["error"]["data"]["tool_code"] == code

dist_params = {
    "norm": {}, "t": {"df": 7}, "chisq": {"df": 7},
    "gamma": {"shape": 2, "rate": 3}, "beta": {"a": 2, "b": 5},
    "binom": {"size": 10, "prob": .4}, "pois": {"lambda": 3},
}
for family, extra in dist_params.items():
    run("dist_density", {"dist": family, "x": [0, 1, 2], **extra}, expected="density-" + family)
    run("dist_probability", {"dist": family, "q": [0, 1, 2], **extra}, expected="probability-" + family)
run("dist_probability", {"dist": "norm", "q": -40, "log_p": True}, expected="probability-log-tail")

rng = {
    "norm": {}, "unif": {}, "exp": {"rate": 2}, "gamma": {"shape": 2, "rate": 3},
    "beta": {"a": 2, "b": 5}, "binom": {"size": 10, "prob": .4}, "pois": {"lambda": 3},
    "t": {"df": 7}, "chisq": {"df": 7},
}
for family, extra in rng.items():
    run("rng_draw", {"dist": family, "n": 3, "seed": 42, **extra}, expected="rng-" + family)
run("rng_sample", {"x": list(range(1, 9)), "size": 5, "seed": 42}, expected="sample")
run("rng_sample", {"x": [1, 2, 3, 4], "size": 6, "seed": 42, "replace": True, "prob": [1, 2, 3, 4]}, expected="sample-weighted")

math = {
    "gamma": {}, "lgamma": {}, "digamma": {}, "trigamma": {},
    "beta": {"y": [2, 3, 4, 5]}, "lbeta": {"y": [2, 3, 4, 5]},
    "choose": {"x": [5, 8, 9, 10], "y": [2, 3, 4, 5]},
    "lchoose": {"x": [5, 8, 9, 10], "y": [2, 3, 4, 5]},
    "besselJ": {"nu": .5}, "besselI": {"nu": .5}, "besselK": {"nu": .5}, "besselY": {"nu": .5},
}
for fn, extra in math.items():
    run("math_special", {"fn": fn, "x": [.5, 1, 2, 3], **extra}, expected="math-" + fn)
for type_ in range(1, 10):
    run("summary_quantile", {"x": [1, 2, 4, 8, 16], "probs": [0, .1, .25, .5, .75, .9, 1], "type": type_}, expected="quantile-" + str(type_))
for inverse, prefix in [(False, "fft"), (True, "ifft")]:
    out = run("num_fft", {"x": [1, 2, 3, 4], "inverse": inverse}, key="re", expected=prefix + "-re")
    # `run` already checks the real part; check the imaginary array from the same invocation too.
    argv = [binary, "num_fft", "--x", "1", "--x", "2", "--x", "3", "--x", "4"]
    if inverse:
        argv.append("--inverse")
    got = json.loads(subprocess.run(argv, check=True, text=True, capture_output=True).stdout)
    want = [float(x) for x in (expected_dir / (prefix + "-im.txt")).read_text().splitlines()]
    assert [float(x) for x in got["im"]] == want, (prefix, got["im"], want)
    assert len(out) == 4

rejects("rng_draw", {"dist": "gamma", "n": 2, "seed": 42}, "MISSING_PARAMETER")
rejects("rng_sample", {"x": [1, 2, 3], "size": 2, "seed": 42, "prob": [1, 2]}, "BAD_PROBABILITY")
rejects("math_special", {"fn": "gamma", "x": [1, 2], "y": [1, 2]}, "BAD_ARGUMENTS")
rejects("summary_quantile", {"x": [1, 2], "probs": [.5], "type": 10}, "BAD_TYPE")
