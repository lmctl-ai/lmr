#!/usr/bin/env python3
"""Randomised differential test of the native core against R (black box). Optional: needs Rscript.

    python3 native/tests/oracle/fuzz_r.py BUILD_DIR [N] [SEED]

Builds fuzz_native from BUILD_DIR/librcli_num.a, generates N random cases per family, evaluates them in one Rscript and
one native process, and reports every mismatch beyond the tolerance (relative, with an absolute floor for values near 0).
"""
import math, os, random, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
INCLUDE = os.path.join(HERE, "..", "..", "include")

FAMILIES = {  # name: (R parameter names, sampler of parameters, sampler of x for d/p, discrete)
    "norm": (["mean", "sd"], lambda r: [r.uniform(-5, 5), 10 ** r.uniform(-3, 3)], None, False),
    "lnorm": (["meanlog", "sdlog"], lambda r: [r.uniform(-2, 2), 10 ** r.uniform(-2, 1)], None, False),
    "t": (["df"], lambda r: [10 ** r.uniform(-1, 4)], None, False),
    "chisq": (["df"], lambda r: [10 ** r.uniform(-1, 3)], None, False),
    "f": (["df1", "df2"], lambda r: [10 ** r.uniform(-1, 3), 10 ** r.uniform(-1, 3)], None, False),
    "gamma": (["shape", "rate"], lambda r: [10 ** r.uniform(-2, 3), 10 ** r.uniform(-2, 2)], None, False),
    "beta": (["shape1", "shape2"], lambda r: [10 ** r.uniform(-2, 2), 10 ** r.uniform(-2, 2)], None, False),
    "exp": (["rate"], lambda r: [10 ** r.uniform(-3, 3)], None, False),
    "weibull": (["shape", "scale"], lambda r: [10 ** r.uniform(-1, 1), 10 ** r.uniform(-2, 2)], None, False),
    "cauchy": (["location", "scale"], lambda r: [r.uniform(-5, 5), 10 ** r.uniform(-2, 2)], None, False),
    "logis": (["location", "scale"], lambda r: [r.uniform(-5, 5), 10 ** r.uniform(-2, 2)], None, False),
    "binom": (["size", "prob"], lambda r: [float(r.randint(0, 2000)), r.random()], None, True),
    "pois": (["lambda"], lambda r: [10 ** r.uniform(-3, 4)], None, True),
    "geom": (["prob"], lambda r: [r.uniform(0.001, 1)], None, True),
    "nbinom": (["size", "prob"], lambda r: [10 ** r.uniform(-2, 3), r.uniform(0.001, 1)], None, True),
    "hyper": (["m", "n", "k"], lambda r: (lambda m, n: [float(m), float(n), float(r.randint(0, m + n))])(r.randint(0, 300), r.randint(0, 300)), None, True),
}


def fmt(v):
    return repr(float(v))


def gen(rng, n):
    cases = []
    for fam, (names, sample, _, discrete) in FAMILIES.items():
        for _ in range(n):
            par = sample(rng)
            fn = rng.choice("dpq")
            lower, lg = rng.random() < 0.6, rng.random() < 0.25
            if fn == "q":
                x = (math.log(rng.uniform(1e-300, 1)) if rng.random() < 0.3 else math.log(rng.random())) if lg else rng.choice([rng.random(), 10 ** -rng.uniform(1, 15), 1 - 10 ** -rng.uniform(1, 15)])
            elif discrete:
                x = float(rng.randint(-2, 3000))
            else:
                x = rng.choice([rng.uniform(-10, 10), 10 ** rng.uniform(-5, 4), -10 ** rng.uniform(-5, 4)])
            if fn == "d":
                lower = True
            p3 = par + [float("nan")] * (3 - len(par))
            rexpr = "%s%s(%s, %s%s%s)" % (fn, fam, fmt(x), ", ".join("%s=%s" % (k, fmt(v)) for k, v in zip(names, par)),
                                         "" if fn == "d" else ", lower.tail=%s" % ("TRUE" if lower else "FALSE"),
                                         ", log=TRUE" if (fn == "d" and lg) else (", log.p=TRUE" if lg else ""))
            native = "dist %s %s %s %s %s %s %d %d" % (fam, fmt(p3[0]), fmt(p3[1]), fmt(p3[2]), fn, repr(x), lower, lg)
            cases.append((rexpr, native, 1))
    for _ in range(n):
        m = rng.randint(1, 40)
        xs = [float(rng.randint(0, 12)) if rng.random() < 0.5 else round(rng.gauss(0, 3), 3) for _ in range(m)]
        t, p = rng.randint(1, 9), rng.choice([0, 1, rng.random()])
        cases.append(("quantile(c(%s), %s, type=%d, names=FALSE)" % (",".join(map(fmt, xs)), fmt(p), t),
                      "quant %d %s %d %s" % (t, repr(p), m, " ".join(map(repr, xs))), 1))
        k = rng.randint(3, 30)
        a = [float(rng.randint(0, 6)) for _ in range(k)]
        b = [float(rng.randint(0, 6)) + (rng.random() if rng.random() < 0.5 else 0) for _ in range(k)]
        meth = rng.choice(["pearson", "spearman", "kendall"])
        vec = lambda v: "%d %s" % (len(v), " ".join(map(repr, v)))
        cases.append(("cor(c(%s), c(%s), method='%s')" % (",".join(map(fmt, a)), ",".join(map(fmt, b)), meth),
                      "cor %s %s %s" % (meth, vec(a), vec(b)), 1))
        kind = rng.choice(["t1", "welch", "w1", "w2", "kendall", "spearman", "pearson"])
        alt = rng.choice(["two.sided", "less", "greater"])
        x = [round(rng.gauss(0.3, 1), rng.choice([1, 3, 6])) for _ in range(rng.randint(3, 60))]
        y = [round(rng.gauss(0, 2), rng.choice([1, 3, 6])) for _ in range(len(x) if kind in ("kendall", "spearman", "pearson") else rng.randint(3, 60))]
        rx, ry = "c(%s)" % ",".join(map(fmt, x)), "c(%s)" % ",".join(map(fmt, y))
        call = {"t1": "t.test(%s, alternative='%s')" % (rx, alt), "welch": "t.test(%s, %s, alternative='%s')" % (rx, ry, alt),
                "w1": "wilcox.test(%s, alternative='%s')" % (rx, alt), "w2": "wilcox.test(%s, %s, alternative='%s')" % (rx, ry, alt)}.get(
            kind, "cor.test(%s, %s, method='%s', alternative='%s')" % (rx, ry, kind, alt))
        native = "test %s %s %s" % (kind, alt, vec(x)) + ("" if kind in ("t1", "w1") else " " + vec(y))
        cases.append(("{h <- %s; c(h$statistic, h$p.value)}" % call, native, 2))
    return cases


def main():
    build = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    work = tempfile.mkdtemp(prefix="rcli-fuzz-")
    exe = os.path.join(work, "fuzz_native")
    subprocess.run(["g++", "-std=c++17", "-O2", "-ffp-contract=off", "-I" + INCLUDE, os.path.join(HERE, "fuzz_native.cpp"),
                    os.path.join(build, "librcli_num.a"), "-pthread", "-o", exe], check=True)
    cases = gen(rng, n)
    rsrc = os.path.join(work, "cases.R")
    with open(rsrc, "w") as f:
        f.write("options(warn=-1)\n")
        for rexpr, _, k in cases:
            f.write("r <- tryCatch(suppressWarnings(%s), error=function(e) NULL); if (is.null(r)) cat('ERROR\\n') else cat(sprintf('%%.17g', r), '\\n')\n" % rexpr)
    rout = subprocess.run(["Rscript", rsrc], capture_output=True, text=True, check=True).stdout.splitlines()
    nout = subprocess.run([exe], input="\n".join(c[1] for c in cases) + "\n", capture_output=True, text=True).stdout.splitlines()
    bad = 0
    for (rexpr, native, _), rl, nl in zip(cases, rout, nout):
        rv, nv = rl.split(), nl.split()
        if rv == ["ERROR"] or nv == ["ERROR"]:
            if (rv == ["ERROR"]) != (nv == ["ERROR"]) and not (rv != ["ERROR"] and all(v in ("NaN", "NA", "nan") for v in rv)):
                bad += 1
                print("ERROR MISMATCH  R=%s native=%s  %s" % (rl.strip(), nl.strip(), rexpr[:200]))
            continue
        for a, b in zip(rv, nv):
            a = float("nan") if a in ("NA", "NaN") else float(a)
            b = float(b)
            if math.isnan(a) and math.isnan(b) or a == b:
                continue
            tol = 1e-9 * max(abs(a), 1e-300) + (1e-300 if "log" not in rexpr else 0)
            if "quantile(c(" in rexpr or rexpr.startswith("cor(") or rexpr.startswith("{"):
                tol = 1e-9 * max(abs(a), 1e-12)
            if not (abs(a - b) <= tol or (math.isinf(a) and a == b)):
                bad += 1
                print("MISMATCH  R=%r native=%r  rel=%.2g  %s" % (a, b, abs(a - b) / max(abs(a), 1e-300), rexpr[:220]))
                break
    print("fuzz: %d cases, %d mismatches" % (len(cases), bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
