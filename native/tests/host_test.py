#!/usr/bin/env python3
"""Host-level tests for lmr: bindings, reverse-communication protocol, CSV, formatting, chunk state, CLI.

Each case compiles a small `_lmcode_` chunk with lmscript and runs it with lmr, then asserts the exit code, values
from the JSON result line, and substrings of stdout/stderr. Positive behaviour is mostly asserted inside the chunk with
out.check (any failure gives exit 4 and names the check). Usage:

    python3 native/tests/host_test.py build-native/bin/lmr [-k substring]

Needs a DSL v5 lmscript ($LMSCRIPT_BIN or on PATH); exits 0 with a notice when none is available.
"""
import concurrent.futures
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCHEMA = os.path.abspath(os.path.join(HERE, "..", "schema", "lmr.dsl.json"))


class Case:
    def __init__(self, name, body, code=0, values=None, out=(), err=(), files=None, chunks=None, args=None):
        self.name = name
        self.body = body          # one chunk body, or None when `chunks` lists several
        self.code = code          # expected exit code
        self.values = values or {}  # expected result values (exact JSON equality)
        self.out = out            # substrings expected on stdout
        self.err = err            # substrings expected on stderr
        self.files = files or {}  # files to create in the work dir before running
        self.chunks = chunks      # several chunk bodies, run in one process
        self.args = args          # raw CLI args instead of --dsl (no compile)


CASES = []


def case(*a, **k):
    CASES.append(Case(*a, **k))


# ---------------------------------------------------------------- Vec
case("vec construction and access", r"""
let v = seq(1.0, 2.0, 0.25);
out.check(v.len() == 5, "seq length");
out.check(v.at(4) == 2.0, "seq end");
out.check(seq_len(3).str() == "[1, 2, 3]", "seq_len");
let l = linspace(0.0, 1.0, 7);
out.check(l.at(6) == 1.0 && l.at(0) == 0.0, "linspace ends exact");
let p = parse_vec("1, 2.5;NA -Inf\t+3 Inf");
out.check(p.len() == 6 && is_nan(p.at(2)) && p.at(3) == -inf() && p.at(4) == 3.0 && p.at(5) == inf(), "parse_vec tokens");
out.check(parse_vec("  ").len() == 0, "parse_vec empty");
let a = vec().push(1.0).push(2.0);
a.append(a);
out.check(a.str() == "[1, 2, 1, 2]", "append to itself");
a.set(0, 9.0);
out.check(a.at(0) == 9.0, "set");
out.check(a.slice(1, 3).str() == "[2, 1]", "slice");
out.check(a.slice(2, 2).len() == 0, "empty slice");
out.check(a.head(10).len() == 4 && a.tail(1).at(0) == 2.0 && a.head(-3).len() == 0, "head/tail clamp");
out.check(a.take(parse_vec("3 0")).str() == "[2, 9]", "take");
out.check(a.where(parse_vec("1 0 NaN 2")).str() == "[9, 2]", "where skips 0 and NaN");
out.check(rep(0.5, 3).sum() == 1.5 && zeros(2).str() == "[0, 0]", "rep zeros");
out.check(vec().push(2598960.0).push(0.1).str() == "[2598960, 0.1]", "str round-trip text");
out.check(parse_vec("10 100 1e20 0.000001 -2500 123456789012").str() == "[10, 100, 1e+20, 1e-06, -2500, 123456789012]", "str avoids needless exponents");
out.check(parse_vec("1 NaN").str() == "[1, NA]", "str NA");
out.check(parse_vec("1.5 2").format("%.2f", "|") == "1.50|2.00", "format");
""")

case("vec summaries", r"""
let x = parse_vec("3 1 2 2 5");
out.check(x.sum() == 13.0 && x.mean() == 2.6, "sum mean");
out.check_close(x.var(), 2.3, 0.000000000001, "var");
out.check(x.median() == 2.0 && x.min() == 1.0 && x.max() == 5.0, "median min max");
out.check(x.sort().str() == "[1, 2, 2, 3, 5]", "sort");
out.check(x.order().str() == "[1, 2, 3, 0, 4]", "order is stable");
out.check(x.rank().str() == "[4, 1, 2.5, 2.5, 5]", "rank average ties");
out.check(x.which_min() == 1 && x.which_max() == 4, "which");
out.check(x.quantile(0.5) == 2.0 && x.quantile(0.25, 1) == 2.0, "quantile types");
out.check(x.quantiles(parse_vec("0 1"), 7).str() == "[1, 5]", "quantiles");
out.check(x.prod() == 60.0 && x.cumsum().str() == "[3, 4, 6, 8, 13]", "prod cumsum");
out.check(x.diff().str() == "[-2, 1, 0, 3]" && x.diff(2).str() == "[-1, 1, 3]", "diff");
let y = parse_vec("1 NaN 3");
out.check(is_nan(y.max()) && is_nan(y.min()), "min/max propagate NaN");
out.check(y.which_max() == 2 && y.count_nan() == 1 && y.drop_nan().len() == 2, "NaN helpers");
out.check_close(parse_vec("1 2 3").cor(parse_vec("2 4 7")), 0.99339926779878, 0.00000000001, "cor pearson");
out.check(parse_vec("1 2 3").cor(parse_vec("1 4 9"), "spearman") == 1.0, "cor spearman");
out.check(parse_vec("1 2 3").dot(parse_vec("4 5 6")) == 32.0, "dot");
let s = parse_vec("2 4 6").standardize();
out.check(abs(s.mean()) < 0.000000000001 && abs(s.sd() - 1.0) < 0.000000000001, "standardize");
""")

case("vec arithmetic and operators", r"""
let a = parse_vec("1 2 3");
let b = parse_vec("4 5 6");
out.check((a + b).str() == "[5, 7, 9]", "vec + vec");
out.check((b - a).str() == "[3, 3, 3]", "vec - vec");
out.check((a * b).str() == "[4, 10, 18]", "vec * vec");
out.check((a + 1.0).str() == "[2, 3, 4]" && (1.0 + a).str() == "[2, 3, 4]", "scalar + both sides");
out.check((10.0 - a).str() == "[9, 8, 7]" && (a - 1.0).str() == "[0, 1, 2]", "scalar - both sides");
out.check((2.0 * a).str() == "[2, 4, 6]" && (-a).str() == "[-1, -2, -3]", "r* and neg");
out.check(b.div(a).at(1) == 2.5 && a.rdiv(6.0).str() == "[6, 3, 2]" && a.div(2.0).at(0) == 0.5, "div rdiv");
out.check(a.pow(2.0).str() == "[1, 4, 9]" && a.sqrt().at(3 - 2) == sqrt(2.0), "pow sqrt");
out.check(a.log().exp().round(12).str() == "[1, 2, 3]", "log exp round");
out.check(parse_vec("-1.5 2").abs().str() == "[1.5, 2]", "abs");
out.check(a.str() == "[1, 2, 3]", "receiver unchanged by arithmetic");
""")

case("round is half to even like R", r"""
out.check(round(2.5) == 2.0 && round(3.5) == 4.0 && round(-0.5) == 0.0, "round half even");
out.check(round(1234.5678, 2) == 1234.57 && round(1250.0, -2) == 1200.0, "round digits");
out.check(parse_vec("0.5 1.5 2.5").round(0).str() == "[0, 2, 2]", "Vec.round");
out.check(is_nan(nan()) && !is_finite(inf()) && is_finite(1.0), "nan inf");
out.check_close(pi(), 3.141592653589793, 0.0, "pi");
out.check(ceil(1.2) == 2.0 && trunc(-1.7) == -1.0 && atan2(1.0, 1.0) == pi() / 4.0, "misc math");
""")

# ---------------------------------------------------------------- strings and formatting
case("fmt and padding", r"""
out.check(fmt("%d", 42) == "42" && fmt("%05d", 7) == "00007", "fmt int");
out.check(fmt("%.3f", 3.14159) == "3.142" && fmt("%10.2e", 12345.0) == "  1.23e+04", "fmt float");
out.check(fmt("%g", 0.0001) == "0.0001" && fmt("x=%.1f%%", 12.34) == "x=12.3%", "fmt %g and %%");
out.check(fmt("[%-5s]", "ab") == "[ab   ]" && fmt("%s!", "hi") == "hi!", "fmt string");
out.check(fmt("%d", 2.0) == "2", "fmt %d with an integral float");
out.check(lpad("ab", 4) == "  ab" && rpad("ab", 4) == "ab  " && lpad("abcdef", 2) == "abcdef", "pad");
let parts = str_split("a,,b,", ",");
out.check(parts.len() == 4 && parts.at(1) == "" && parts.at(3) == "", "split keeps empty fields");
out.check(str_split("a::b", "::").join("|") == "a|b", "multi-char separator");
out.check(parts.contains("b") && parts.index_of("b") == 2 && parts.index_of("zz") == -1, "contains index");
out.check(str_list().add("x").add("y").join("") == "xy", "std str_list add join");
""")

case("lmcode std is provided by the runtime", r"""
out.check(str_trim("  a b ") == "a b" && str_upper("ab") == "AB" && str_replace("a-b-c", "-", "+") == "a+b+c", "str");
out.check(str_starts_with("lmr", "lm") && str_find("abc", "c") == 2 && str_to_int("42") == 42, "str parse");
out.check(math_sqrt(16.0) == 4.0 && math_max(3, 7) == 7 && math_round(2.5) == 3.0, "math");
let m = str_map().set("k", "v").set("j", "w");
out.check(m.len() == 2 && m.get("k") == "v" && m.get_or("zz", "d") == "d" && m.key_at(1) == "j", "StrMap");
let il = int_list().add(3).add(5);
let fl = float_list().add(0.5).set(0, 1.5);
out.check(il.len() == 2 && il.at(1) == 5 && fl.at(0) == 1.5, "IntList and FloatList");
let names = read_csv("in.csv").names();
out.check(names.index_of("y") >= 0 && names.add("extra").len() == names.len(), "host StrList is a std StrList");
file_write("std.txt", "a,b,");
file_append("std.txt", "c");
out.check(file_exists("std.txt") && file_read("std.txt") == "a,b,c" && file_lines("std.txt").len() == 1, "file");
""", files={"in.csv": "x,y\n1,2\n"})

for spec, why in [("%d %d", "exactly one conversion"), ("%n", "unsupported conversion"), ("plain", "no conversion"),
                  ("%5000f", "width too large"), ("%.300f", "precision too large"), ("%x", "unsupported conversion")]:
    case("fmt rejects " + spec, 'out.print(fmt("%s", 1.0));' % spec.replace('"', '\\"'), code=3, out=(why,))
case("fmt %s needs a string", 'out.print(fmt("%s", 1.0));', code=3, out=("needs a string value",))
case("fmt %d rejects a fraction-free huge value", 'out.print(fmt("%d", inf()));', code=3, out=("finite integral",))

# ---------------------------------------------------------------- Frame and CSV
CSV = ('﻿name,score,"note, with comma",empty\r\n'
       'a,1.5,"say ""hi""",\r\n'
       'b,NA,x,\r\n'
       '\r\n'
       'c,-2e3,"multi\nline",\r\n')
case("read_csv quoting, BOM, CRLF, NA, blank lines", r"""
let f = read_csv("in.csv");
out.check(f.nrow() == 3 && f.ncol() == 4, "shape");
out.check(f.names().join("|") == "name|score|note, with comma|empty", "header with BOM stripped and quoted name");
out.check(f.is_text("name") && !f.is_text("score") && !f.is_text("empty"), "type detection; all-empty column is numeric");
out.check(f.col("score").at(0) == 1.5 && is_nan(f.col("score").at(1)) && f.col("score").at(2) == -2000.0, "numbers and NA");
out.check(f.text("note, with comma").at(0) == "say \"hi\"", "doubled quotes");
out.check(f.text("note, with comma").at(2) == "multi\nline", "newline inside quotes");
out.check(f.col("empty").count_nan() == 3, "empty cells are NA");
out.check(f.text("score").at(1) == "NA", "numeric column as text");
""", files={"in.csv": CSV})

case("write_csv round trip", r"""
let f = frame().set("x", parse_vec("1 NaN 0.1")).set("s", str_split("a,b c|q\"t|", "|"));
f.write_csv("sub/out.csv");
let g = read_csv("sub/out.csv");
out.check(g.nrow() == 3 && is_nan(g.col("x").at(1)) && g.col("x").at(2) == 0.1, "numbers round trip exactly");
out.check(g.text("s").at(0) == "a,b c" && g.text("s").at(1) == "q\"t" && g.text("s").at(2) == "", "text round trip");
""", files={"sub/.keep": ""})

case("frame operations", r"""
let f = frame().set("a", parse_vec("1 2 3 4")).set("b", parse_vec("5 NaN 7 8")).set("t", str_split("w,x,y,z", ","));
out.check(f.select("b, a").names().join(",") == "b,a", "select order");
out.check(f.where(parse_vec("0 1 1 0")).text("t").join("") == "xy", "where");
out.check(f.drop_nan().nrow() == 3, "drop_nan");
out.check(f.head(2).nrow() == 2 && f.head(99).nrow() == 4, "head");
f.set("a", parse_vec("9 9 9 9"));
out.check(f.col("a").sum() == 36.0 && f.ncol() == 3, "set replaces in place");
out.check(f.has("t") && !f.has("q"), "has");
""")

for body, msg in [
    ('let f = read_csv("missing.csv");', "cannot open missing.csv"),
    ('let f = read_csv("bad1.csv");', "unterminated quoted field"),
    ('let f = read_csv("bad2.csv");', "has 3 fields; the header has 2"),
    ('let f = read_csv("bad3.csv");', "duplicate column 'a'"),
    ('let f = read_csv("bad4.csv");', "unexpected quote"),
    ('let f = read_csv("bad5.csv");', "empty name"),
    ('frame().set("a", vec().push(1.0)).write_csv("/tmp/x.csv");', "must be relative"),
    ('frame().set("a", vec().push(1.0)).write_csv("../x.csv");', "must not contain '..'"),
    ('frame().set("a", vec().push(1.0)).write_csv("dir/");', "names a directory"),
    ('let f = frame().set("t", str_split("a", ",")); out.print(str(f.col("t").len()));', "is text (use Frame.text)"),
    ('let f = frame().set("a", parse_vec("1 2")); out.print(str(f.col("b").len()));', "no column 'b' (columns: a)"),
    ('let f = frame().set("a", parse_vec("1 2")).set("b", parse_vec("1"));', "has 1 rows; the frame has 2"),
    ('let f = frame().set("a", parse_vec("1 2")).select("a,,b");', "empty column name"),
    ('let f = frame().set("a", parse_vec("1 2")).where(parse_vec("1"));', "mask length 1 differs from nrow 2"),
]:
    case("frame error: " + msg, body, code=3, out=(msg,), files={
        "bad1.csv": 'a,b\n1,"x\n', "bad2.csv": "a,b\n1,2,3\n", "bad3.csv": "a,a\n1,2\n",
        "bad4.csv": 'a,b\n1,x"y\n', "bad5.csv": "a,\n1,2\n"})

# ---------------------------------------------------------------- Vec errors
for body, msg in [
    ("out.print(str(parse_vec(\"1 2\").at(2)));", "Vec.at: index 2 out of range 0..1"),
    ("out.print(str(vec().at(0)));", "Vec.at: index 0 out of range 0..-1"),
    ("parse_vec(\"1 2\").set(-1, 0.0);", "Vec.set: index -1"),
    ("out.print(parse_vec(\"1 2\").add(parse_vec(\"1\")).str());", "Vec.add: lengths differ (2 and 1)"),
    ("out.print(parse_vec(\"1 x\").str());", "parse_vec: not a number: 'x'"),
    ("out.print(str(vec().min()));", "Vec.min: empty vector"),
    ("out.print(parse_vec(\"1 NaN\").sort().str());", "NaN values"),
    ("out.print(parse_vec(\"1 2\").slice(1, 3).str());", "Vec.slice"),
    ("out.print(parse_vec(\"1 2\").take(parse_vec(\"0.5\")).str());", "Vec.take: index"),
    ("out.print(seq(1.0, 2.0, 0.0).str());", "seq: wrong sign or zero 'by'"),
    ("out.print(seq(1.0, 2.0, -1.0).str());", "seq: wrong sign or zero 'by'"),
    ("out.print(zeros(-1).str());", "zeros: n must be between 0 and 1e8"),
    ("out.print(linspace(0.0, 1.0, 1).str());", "linspace: n must be between 2"),
    ("out.print(parse_vec(\"1 2\").diff(0).str());", "lag must be at least 1"),
    ("out.print(str(vec().which_min()));", "no non-NaN values"),
    ("out.print(str_split(\"a\", \"\").join(\"\"));", "str_split: separator must not be empty"),
    ("out.print(lpad(\"a\", 5000));", "lpad: width too large"),
]:
    case("vec error: " + msg, body, code=3, out=(msg,))

# ---------------------------------------------------------------- Out, result line, exit codes
case("out.set types and ordering", r"""
out.set("b", 1.5);
out.set("a", 2);
out.set("s", "q\"uote\n");
out.set("t", true);
out.set("v", parse_vec("1 NaN Inf"));
out.set("l", str_split("x,y", ","));
out.set("b", 3.25);
""", values={"b": 3.25, "a": 2, "s": "q\"uote\n", "t": True, "v": [1, None, None], "l": ["x", "y"]})
case("out.print goes to stderr only", 'out.print("hello report");', err=("hello report",))
case("failed check exits 4 and names the first failure", r"""
out.check(true, "fine");
out.check(1 == 2, "first broken");
out.check(false, "second broken");
out.print("still running");
""", code=4, out=('"code": "APP_FAIL"', "2 of 3 checks failed; first: first broken"),
     err=("CHECK FAILED: first broken", "still running"))
case("check_close reports actual and expected", 'out.check_close(1.0, 2.0, 0.5, "near");', code=4,
     err=("CHECK FAILED: near: actual 1, expected 2, tolerance 0.5",))
case("check_close rejects NaN", 'out.check_close(nan(), 1.0, 1.0, "nan actual");', code=4, err=("nan actual",))
case("out.fail stops with exit 4", 'out.fail("stop here"); out.print("not reached");', code=4,
     out=('"phase": "run"', '"code": "APP_FAIL"', "stop here"))
# The wrapper's `p = _lmcode_` line is line 1 of the compiled file, so the body's second line is source line 3.
case("host errors carry the dsl line", 'let v = vec();\nout.print(str(v.at(3)));', code=3,
     out=('"code": "RUN_HOST"', '"line": 3', '"column": 17', "dsl line 3:17: Vec.at"))
case("failures() counts", 'out.check(false, "x"); out.set("n", out.failures());', code=4)

# ---------------------------------------------------------------- keep / kept across chunks
case("keep and kept across chunks", None, chunks=[
    'keep("v", parse_vec("1 2 3")); keep("f", frame().set("a", parse_vec("4 5")));',
    'out.check(kept_vec("v").sum() == 6.0, "kept vec"); out.check(kept_frame("f").col("a").at(1) == 5.0, "kept frame");'
    ' out.check(is_kept("v") && !is_kept("zz"), "is_kept"); kept_vec("v").push(4.0);',
    'out.check(kept_vec("v").len() == 4, "kept object is shared, not copied"); out.set("ok", true);'],
     values={"ok": True})
case("kept errors", None, chunks=['keep("v", vec());', 'out.print(str(kept_frame("v").nrow()));'], code=3,
     out=("'v' holds a Vec, not a Frame",))
case("nothing kept", 'out.print(str(kept_vec("nope").len()));', code=3, out=("nothing kept under 'nope'",))
case("second chunk error names its file", None, chunks=['out.set("first", 1);', 'out.print(str(vec().at(1)));'],
     code=3, out=("chunk1.dslb: dsl line",))

# ---------------------------------------------------------------- reverse communication
case("minimize, optimize, uniroot, integrate", r"""
let o = minimize(parse_vec("-1.2 1"), "Nelder-Mead");
while (o.next()) {
  o.tell(100.0 * pow(o.x(1) - o.x(0) * o.x(0), 2.0) + pow(1.0 - o.x(0), 2.0));
}
out.check(abs(o.par(0) - 1.0) < 0.01 && o.convergence() == 0 && o.fn_count() == 195, "Nelder-Mead as R (195 evals)");
out.check(o.evaluations() == o.fn_count(), "every evaluation came through tell");
out.check(o.pars().len() == 2 && o.message() == "", "pars message");
let b = minimize(parse_vec("0 0"), "L-BFGS-B").lower(parse_vec("-1 -1")).upper(parse_vec("0.5 2")).hessian(true);
while (b.next()) {
  b.tell(pow(b.x(0) - 1.0, 2.0) + pow(b.x(1) - 1.0, 2.0));
}
out.check(abs(b.par(0) - 0.5) < 0.000001 && abs(b.par(1) - 1.0) < 0.00001, "bound active");
out.check(abs(b.hessian(0, 0) - 2.0) < 0.0001 && abs(b.hessian(0, 1)) < 0.0001, "hessian");
let m = optimize(0.0, 3.0, true).tol(0.000000001);
while (m.next()) {
  m.tell(sin(m.x()));
}
out.check(abs(m.argmin() - pi() / 2.0) < 0.00001 && abs(m.objective() - 1.0) < 0.000000001, "optimize maximum");
let r = uniroot(0.0, 1.0).extend("yes").tol(0.0000000001);
while (r.next()) {
  r.tell(r.x() - 5.0);
}
out.check(abs(r.root() - 5.0) < 0.000001 && r.init_it() > 0 && r.converged(), "uniroot extends the interval");
let q = integrate(-inf(), inf());
while (q.next()) {
  q.tell(exp(-q.x() * q.x()));
}
out.check(abs(q.value() - sqrt(pi())) < 0.000000001 && q.ok() && q.message() == "OK", "Gaussian integral");
let z = integrate(1.0, 1.0);
while (z.next()) {
  z.tell(1.0);
}
out.check(z.value() == 0.0, "empty range");
""")
case("integrate reports a non-finite integrand", r"""
let q = integrate(-1.0, 1.0);
while (q.next()) {
  q.tell(log(q.x()));
}
out.check(!q.ok() && q.ier() == 6, "ier 6");
out.print(q.message());
""", err=("non-finite",))
case("abs_tol defaults to rel_tol", r"""
let g = integrate(0.0, inf()).rel_tol(0.000000000001);
while (g.next()) {
  g.tell(g.x() * g.x() * exp(-g.x()));
}
out.check(abs(g.value() - 2.0) < 0.00000000001, "tight relative tolerance honoured");
""")
case("nested solvers run on separate workers", r"""
# Minimise over a of the integral of (x - a)^2 on [0, 1]: the minimum is at a = 0.5.
let o = optimize(-2.0, 2.0).tol(0.00000001);
while (o.next()) {
  let a = o.x();
  let q = integrate(0.0, 1.0);
  while (q.next()) {
    q.tell(pow(q.x() - a, 2.0));
  }
  o.tell(q.value());
}
out.check(abs(o.argmin() - 0.5) < 0.000001, "outer optimum through inner integrals");
""")
case("abandoned solvers are cancelled cleanly", r"""
let o = minimize(parse_vec("1 1"), "BFGS");
let k = 0;
while (o.next()) {
  o.tell(o.x(0) * o.x(0) + o.x(1) * o.x(1));
  k = k + 1;
  if (k == 3) {
    break;
  }
}
let i = integrate(0.0, 1.0);
i.next();
let r = uniroot(0.0, 1.0);
r.next();
r.tell(-1.0);
r.next();
out.set("k", k);
""", values={"k": 3})
case("many solvers in a loop", r"""
let total = 0.0;
for (j in range(0, 200)) {
  let q = integrate(0.0, 1.0);
  while (q.next()) {
    q.tell(q.x() * float(j));
  }
  total = total + q.value();
}
out.check(abs(total - 9950.0) < 0.000001, "200 integrals");
""")
case("SANN needs a seed", 'let o = minimize(parse_vec("1"), "SANN");\nwhile (o.next()) { o.tell(o.x(0) * o.x(0)); }',
     code=3, out=("SANN needs .seed(n)",))
case("SANN with a seed is reproducible", r"""
let a = minimize(parse_vec("2 2"), "SANN").seed(5).maxit(500);
while (a.next()) {
  a.tell(a.x(0) * a.x(0) + a.x(1) * a.x(1));
}
let b = minimize(parse_vec("2 2"), "SANN").seed(5).maxit(500);
while (b.next()) {
  b.tell(b.x(0) * b.x(0) + b.x(1) * b.x(1));
}
out.check(a.par(0) == b.par(0) && a.value() == b.value(), "same seed, same path");
out.check(a.value() < 0.5, "SANN improves");
""")
for body, msg in [
    ('let o = minimize(parse_vec("1"), "BFGS");\no.tell(1.0);', "no point is waiting (call next() first)"),
    ('let o = minimize(parse_vec("1"), "BFGS");\no.next();\no.next();', "tell(value) was not called for the current point"),
    ('let o = minimize(parse_vec("1"), "BFGS");\no.next();\no.tell(1.0);\no.tell(2.0);', "already has a value"),
    ('let o = minimize(parse_vec("1"), "BFGS");\nout.print(str(o.par(0)));', "results are available after next() returns false"),
    ('let o = minimize(parse_vec("1"), "BFGS");\no.next();\no.maxit(5);', "options must be set before the first next()"),
    ('let o = minimize(parse_vec("1"), "BFGS");\nout.print(str(o.x(0)));', "no current point"),
    ('let o = minimize(parse_vec("1"), "Newton");', "method must be Nelder-Mead, BFGS, CG, L-BFGS-B or SANN"),
    ('let o = minimize(vec(), "BFGS");', "start must have at least one value"),
    ('let o = minimize(parse_vec("1"), "BFGS");\nwhile (o.next()) { o.tell(o.x(0)); }\nout.print(str(o.hessian(0, 0)));', "call .hessian(true)"),
    ('let q = optimize(1.0, 0.0);', "need finite lower < upper"),
    ('let q = uniroot(0.0, inf());', "need finite lower < upper"),
    ('let r = uniroot(0.0, 1.0);\nwhile (r.next()) { r.tell(1.0 + r.x()); }', "opposite sign"),
    ('let r = uniroot(0.0, 1.0).extend("sideways");', "must be \\\"no\\\", \\\"yes\\\""),
    ('let q = integrate(nan(), 1.0);', "bounds must not be NaN"),
    ('let q = integrate(0.0, 1.0).subdivisions(0);', "between 1 and 100000"),
]:
    case("protocol error: " + msg[:60], body, code=3, out=(msg,))

# ---------------------------------------------------------------- nls
case("nls fits Michaelis-Menten", r"""
let x = seq(0.25, 4.0, 0.25);
let y = vec();
for (i in range(0, x.len())) {
  y.push(4.0 * x.at(i) / (2.0 + x.at(i)) + 0.01 * sin(float(i) * 7.0));
}
let m = nls(y).param("a", 3.0).param("b", 1.0);
while (m.next()) {
  for (i in range(0, m.n())) {
    m.tell(i, m.par("a") * x.at(i) / (m.par("b") + x.at(i)));
  }
}
out.check(m.converged() && abs(m.estimate("a") - 4.0) < 0.05 && abs(m.estimate(1) - 2.0) < 0.05, "estimates");
out.check(m.df_residual() == 14 && m.names().join(",") == "a,b", "df and names");
out.check(abs(m.fitted().at(0) + m.residuals().at(0) - y.at(0)) < 0.000000001, "fitted + residual = y");
out.check(m.conf_low("a", 0.95) < m.estimate("a") && m.estimate("a") < m.conf_high("a", 0.95), "interval");
let p = nls(y).param("a", 3.0).param("b", 1.0).algorithm("port").upper("a", 3.5);
while (p.next()) {
  for (i in range(0, p.n())) {
    p.tell(i, p.par(0) * x.at(i) / (p.par(1) + x.at(i)));
  }
}
out.check(abs(p.estimate("a") - 3.5) < 0.000001, "port keeps the bound");
""")
for body, msg in [
    ('let m = nls(parse_vec("1 2")).param("a", 1.0);\nm.next();\nm.tell(0, 1.0);\nm.next();', "tell(i, value) was called for 1 of 2 observations"),
    ('let m = nls(parse_vec("1 2")).param("a", 1.0);\nm.next();\nm.tell(0, 1.0);\nm.tell(0, 1.0);', "observation 0 already has a value"),
    ('let m = nls(parse_vec("1 2"));\nm.next();', "declare parameters with .param"),
    ('let m = nls(parse_vec("1 2")).param("a", 1.0).param("a", 2.0);', "'a' is declared twice"),
    ('let m = nls(parse_vec("1 2")).param("a", 1.0).lower("a", 0.0);\nm.next();', "bounds need .algorithm(\\\"port\\\")"),
    ('let m = nls(parse_vec("1 2")).lower("a", 0.0);', "declare 'a' with .param first"),
    ('let m = nls(parse_vec("1 NaN"));', "non-finite values"),
    ('let m = nls(parse_vec("1 2")).param("a", 1.0);\nm.next();\nout.print(str(m.par("b")));', "no parameter 'b'"),
]:
    case("nls error: " + msg[:60], body, code=3, out=(msg,))

# ---------------------------------------------------------------- models through Frames
case("lm with NA rows, weights and predict", r"""
let f = frame().set("y", parse_vec("1 3 NaN 7 9 12")).set("x", parse_vec("0 1 2 3 4 5")).set("t", str_split("a,b,c,d,e,f", ","));
let m = lm(f, "y ~ x");
out.check(m.n() == 5 && m.df_residual() == 3, "NA row dropped (na.omit)");
out.check(m.has_coef("x") && !m.has_coef("t"), "text column not used");
let nf = frame().set("x", parse_vec("10"));
out.check(abs(m.predict(nf).at(0) - (m.estimate(0) + 10.0 * m.estimate(1))) < 0.000000001, "predict");
out.check(m.predict_se(nf).at(0) > 0.0, "predict_se");
let w = lm(f, "y ~ x", parse_vec("1 1 1 1 1 0"));
out.check(w.n() == 4, "zero weight is not an observation");
let all = lm(frame().set("y", parse_vec("1 2 4 3")).set("a", parse_vec("1 2 3 5")).set("b", parse_vec("2 1 2 1")), "y ~ .");
out.check(all.n_coef() == 3, "dot expands to the other numeric columns");
out.check(m.coef_table() != "", "coef_table");
""")
case("glm poisson with offset", r"""
let f = frame().set("y", parse_vec("2 3 6 7 8 9 10 12 15 20")).set("x", seq_len(10)).set("e", parse_vec("10 12 15 14 13 16 18 17 20 25"));
let m = glm(f, "y ~ x", "poisson", "log", vec(), f.col("e").log());
out.check(abs(m.estimate(0) + 1.3117093525715928) < 0.0000000001 && abs(m.estimate(1) - 0.11445952529770297) < 0.0000000001, "R coefficients");
out.check(abs(m.deviance() - 1.5350120625971155) < 0.000000001, "R deviance");
out.check(m.family() == "poisson" && m.link() == "log" && m.converged(), "meta");
out.check(m.residuals("pearson").len() == 10, "pearson residuals");
""")
for body, msg in [
    ('let m = lm(frame().set("y", parse_vec("1 2")), "y ~ z");', "z"),
    ('let m = lm(frame().set("y", parse_vec("1 2 3")).set("t", str_split("a,b,c", ",")), "t ~ y");', "the response 't' is a factor (text)"),
    ('let m = lm(frame().set("y", parse_vec("1 2")).set("t", str_split("a,a", ",")), "y ~ t");', "contrasts can be applied only to factors with 2 or more levels"),
    ('let m = lm(frame().set("y", parse_vec("1 2")), "y ~ x", parse_vec("1"));', "weights have 1 values"),
    ('let m = glm(frame().set("y", parse_vec("1 2")).set("x", parse_vec("1 2")), "y ~ x", "nonsense");', "family"),
    ('let m = glm(frame().set("y", parse_vec("1 2")).set("x", parse_vec("1 2")), "y ~ x", "poisson", "log", vec(), parse_vec("1"));', "offset has 1 values"),
    ('let m = glm(frame().set("y", parse_vec("1 2")).set("x", parse_vec("1 2")), "y ~ x", "poisson");\nout.print(m.residuals("raw").str());', "deviance, pearson or working"),
    ('let m = lm(frame().set("y", parse_vec("1 2 3")).set("x", parse_vec("1 2 4")), "y ~ x");\nout.print(str(m.estimate("q")));', "no coefficient 'q' (coefficients: (Intercept), x)"),
    ('let m = lm(frame().set("y", parse_vec("1 2 3")).set("x", parse_vec("1 2 4")), "y ~ x");\nout.print(str(m.conf_low(0, 1.5)));', "level must be between 0 and 1"),
]:
    case("model error: " + msg[:60], body, code=3, out=(msg,))

# ---------------------------------------------------------------- tests, dist, rng, fft
case("hypothesis test builder", r"""
let t = test("t").x(parse_vec("5.1 4.9 5.6 5.8 6.0 5.2")).mu(5.0).alt("greater").conf(0.9).run();
out.check(t.statistic_name() == "t" && t.df() == 5.0 && t.alternative() == "greater", "t meta");
out.check(t.has_ci() && t.conf_high() == inf() && t.conf_level() == 0.9, "one-sided interval");
out.check(t.n_params() == 1 && t.param_name(0) == "df" && t.param("df") == 5.0, "params");
let c = test("chisq").table(parse_vec("10 20 30 40"), 2, false).correct(false).run();
out.check(c.expected().len() == 4 && abs(c.expected().at(0) - 12.0) < 0.000000001, "expected counts, filled by column");
let c2 = test("chisq").table(parse_vec("10 30 20 40"), 2, true).correct(false).run();
out.check(abs(c2.statistic() - c.statistic()) < 0.000000000001, "byrow fills the same table");
let g = test("chisq").counts(parse_vec("10 10 20")).run();
out.check(g.df() == 2.0, "goodness of fit, equal probabilities by default");
let b = test("binom").successes(60.0).trials(100.0).run();
out.check(abs(b.p_value() - 0.05688793) < 0.00000001 && b.estimate() == 0.6, "exact binomial");
let pr = test("prop").successes(parse_vec("15 25")).trials(parse_vec("50 50")).run();
out.check(pr.n_estimates() == 2 && pr.estimate_name(0) == "prop 1", "prop estimates");
let k = test("ks").x(parse_vec("0.1 0.4 0.7 0.2 0.9")).dist(unif(0.0, 1.0)).run();
out.check(k.p_value() > 0.5, "ks against a Dist");
out.check(t.summary() != "" && t.method() == "One Sample t-test", "summary and method");
""")
for body, msg in [
    ('let t = test("z");', "unknown kind 'z' (kinds: t, wilcox"),
    ('let t = test("t").run();', "call .x(...) before run()"),
    ('let t = test("var").x(parse_vec("1 2 3")).run();', "call .y(...) before run()"),
    ('let t = test("t").alt("bigger");', "alternative must be"),
    ('let t = test("t").conf(1.0);', "confidence level must be between 0 and 1"),
    ('let t = test("chisq").table(parse_vec("1 2 3"), 2, true);', "3 values do not fill 2 rows"),
    ('let t = test("fisher").run();', "call .table(values, nrow, byrow)"),
    ('let t = test("binom").successes(parse_vec("1 2")).trials(10.0).run();', "take one number each"),
    ('let t = test("shapiro").x(parse_vec("1 2 3")).run();\nout.print(str(t.conf_low()));', "has no confidence interval"),
    ('let t = test("shapiro").x(parse_vec("1 2 3")).run();\nout.print(str(t.df()));', "no parameter 'df'"),
    ('let t = test("shapiro").x(parse_vec("1 2 3")).run();\nout.print(str(t.estimate()));', "has no estimate"),
]:
    case("test error: " + msg[:60], body, code=3, out=(msg,))

case("dist and rng bindings", r"""
let n = norm(1.0, 2.0);
out.check(n.name() == "norm" && n.mean() == 1.0 && n.sd() == 2.0 && !n.is_discrete(), "meta");
out.check(abs(n.q(n.p(0.3)) - 0.3) < 0.000000000001, "round trip");
out.check(n.p_vec(parse_vec("1 1")).str() == "[0.5, 0.5]", "p_vec");
out.check(abs(n.log_p_upper(1.0) - log(0.5)) < 0.000000000001 && n.q_upper(0.5) == 1.0, "tails");
out.check(dist("t", 5.0).name() == "t" && dist("hyper", 5.0, 3.0, 2.0).is_discrete(), "dist by name");
out.check(binom(10, 0.5).q(0.5) == 5.0 && pois(0.0).p(0.0) == 1.0, "discrete");
let r = rng(9);
let a = r.draw(norm(), 1000);
out.check(abs(a.mean()) < 0.15 && abs(a.sd() - 1.0) < 0.1, "normal draws");
out.check(rng(9).draw(norm(), 1000).str() == a.str(), "same seed same stream");
let s = rng(3).sample(seq_len(10), 10, false);
out.check(s.sort().str() == seq_len(10).str(), "permutation sample");
out.check(rng(3).permutation(5).sort().str() == "[0, 1, 2, 3, 4]", "permutation");
let w = rng(4).sample(parse_vec("1 2"), 2000, true, parse_vec("0 1"));
out.check(w.min() == 2.0, "zero weight never drawn");
out.check(r.below(1) == 0, "below 1");
out.check(abs(gamma(5.0) - 24.0) < 0.000000001 && choose(5.0, 2.0) == 10.0 && abs(beta_inc(2.0, 3.0, 1.0) - 1.0) < 0.000000000001, "special");
""")
for body, msg in [
    ("let d = norm(0.0, -1.0);", "sd"),
    ("let d = binom(-1, 0.5);", "binom"),
    ("let d = dist(\"nope\", 1.0);", "nope"),
    ("let r = rng(1);\nout.print(str(r.below(0)));", "n must be at least 1"),
    ("let r = rng(1);\nout.print(r.sample(parse_vec(\"1 2\"), 3, false).str());", ""),
    ("let r = rng(1);\nout.print(r.sample(parse_vec(\"1 2\"), 1, true, parse_vec(\"1\")).str());", "prob must have the length of x"),
    ("let r = rng(1);\nout.print(r.draw(norm(), -1).str());", "n must be between 0 and 1e8"),
]:
    case("dist error: " + (msg or "sample larger than population"), body, code=3, out=(msg,))

case("fft bindings", r"""
let z = fft(parse_vec("1 2 3 4"));
out.check(z.len() == 4 && z.re(0) == 10.0 && z.re(2) == -2.0 && z.im(1) == 2.0, "forward");
let back = z.fft(true).scale(0.25);
out.check(back.re_vec().round(12).str() == "[1, 2, 3, 4]" && back.im_vec().abs().max() < 0.000000000001, "inverse / n");
out.check(abs(z.mod(1) - sqrt(8.0)) < 0.000000000001 && abs(z.arg(1) - 3.0 * pi() / 4.0) < 0.000000000001, "mod arg");
out.check(cvec(parse_vec("1 0"), parse_vec("0 1")).fft(false).len() == 2, "cvec");
out.check(fft(vec().push(5.0)).re(0) == 5.0, "length 1");
""")
case("cvec length mismatch", 'let c = cvec(parse_vec("1 2"), parse_vec("1"));', code=3, out=("differ in length",))

# ---------------------------------------------------------------- smoothing, ts, cluster
case("smoothing bindings", r"""
let x = seq(0.0, 3.0, 0.5);
let y = x.pow(3.0);
let s = spline_fun(x, y);
out.check(abs(s.at(1.25) - 1.953125) < 0.000000001 && abs(s.deriv(1.25, 1) - 4.6875) < 0.000000001, "fmm exact on a cubic");
out.check(s.eval(x).str() == y.str() && s.method() == "fmm", "interpolates the knots");
let ap = approx_fun(parse_vec("0 1"), parse_vec("0 10"));
out.check(ap.at(0.25) == 2.5 && is_nan(ap.at(2.0)), "approx rule 1");
out.check(approx_fun(parse_vec("0 1"), parse_vec("0 10"), "constant", 2).at(5.0) == 10.0, "approx constant rule 2");
let xx = seq(0.0, 10.0, 0.25);
let yy = xx.sqrt();
let ss = smooth_spline_df(xx, yy, 5.0);
out.check(abs(ss.df() - 5.0) < 0.001 && ss.leverages().len() == xx.len(), "df fit");
let lo = loess(xx, yy, 0.5, 1);
out.check(lo.fitted().len() == xx.len() && lo.enp() > 2.0 && is_nan(lo.at(11.0)), "loess");
""")
case("ts bindings", r"""
let r = rng(11);
let x = arima_sim(r, 300, parse_vec("0.6"), vec(), 1.0);
out.check(x.len() == 300 && x.acf(5).len() == 6 && x.acf(5).at(0) == 1.0 && x.pacf(4).len() == 4, "acf pacf shapes");
let m = arima(1, 0, 0).method("ML").fit(x);
out.check(abs(m.estimate("ar1") - 0.6) < 4.0 * m.se("ar1") && m.has_coef("intercept") && m.nobs() == 300, "ar1");
let f = m.forecast(3);
out.check(f.len() == 3 && f.lower(0, 0.95) < f.mean_at(0) && f.mean_at(0) < f.upper(0, 0.95), "forecast interval");
out.check(abs(f.se_at(0) - sqrt(m.sigma2())) < 0.000001, "one-step se");
let i = arima_sim(rng(2), 50, vec(), vec(), 1.0, 1);
out.check(i.len() == 51 && i.at(0) == 0.0, "integrated simulation starts at 0 and has n + d values");
let d = arima(0, 1, 1).fit(i);
out.check(!d.has_coef("intercept"), "no mean for differenced models");
let season = vec();
for (t in range(0, 48)) {
  season.push(sin(2.0 * pi() * float(t) / 12.0) + 0.01 * float(t));
}
let s = stl_periodic(season, 12);
out.check(abs(s.seasonal().at(3) - s.seasonal().at(15)) < 0.000000001, "periodic seasonal repeats");
out.check((s.seasonal() + s.trend() + s.remainder() - season).abs().max() < 0.000000001, "components add up");
""")
for body, msg in [
    ('let a = arima(-1, 0, 0);', "need 0 <= p <= 50"),
    ('let a = arima(1, 0, 0).method("OLS");', "must be CSS-ML, ML or CSS"),
    ('let a = arima(1, 0, 0).seasonal(1, 0, 0, 1);', "period >= 2"),
    ('let x = arima_sim(rng(1), 10, parse_vec("1.5"), vec(), 1.0);', "not stationary"),
    ('let s = stl(parse_vec("1 2 3 4"), 1, 7);', "period must be at least 2"),
    ('let m = arima(1, 0, 0).fit(arima_sim(rng(1), 50, parse_vec("0.5"), vec(), 1.0));\nout.print(str(m.forecast(0).len()));', "horizon must be between 1"),
]:
    case("ts error: " + msg[:60], body, code=3, out=(msg,))

case("cluster bindings", r"""
let f = frame().set("x", parse_vec("0 0.1 0.2 5 5.1 5.2")).set("y", parse_vec("0 0.1 0 5 5 5.1"));
let k = kmeans(f, 2, rng(1), 5);
out.check(k.k() == 2 && k.n() == 6 && k.size(1) + k.size(2) == 6, "sizes");
out.check(k.cluster(0) == k.cluster(2) && k.cluster(0) != k.cluster(3), "partition");
out.check(abs(k.center(k.cluster(3), "x") - 5.1) < 0.000000001, "centre by column name");
out.check(abs(k.totss() - k.tot_withinss() - k.betweenss()) < 0.000000001, "SS decomposition");
let h = hclust(f, "average");
out.check(h.n() == 6 && h.heights().len() == 5 && h.order().len() == 6, "hclust shapes");
out.check(h.cut(2).str() == "[1, 1, 1, 2, 2, 2]" && h.cut(6).max() == 6.0, "cut");
out.check(h.height(4) >= h.height(3), "complete heights monotone for average");
let hm = hclust(f, "single", "manhattan");
out.check(abs(hm.height(0) - 0.1) < 0.000000001, "manhattan distance");
""")
for body, msg in [
    ('let k = kmeans(frame().set("t", str_split("a,b", ",")), 1, rng(1));', "is text; select the numeric columns"),
    ('let k = kmeans(frame().set("x", parse_vec("1 NaN")), 1, rng(1));', "missing values"),
    ('let k = kmeans(frame(), 1, rng(1));', "the frame is empty"),
    ('let k = kmeans(frame().set("x", parse_vec("1 2")), 0, rng(1));', "k must be at least 1"),
    ('let h = hclust(frame().set("x", parse_vec("1")), "average");', "need at least 2 rows"),
    ('let h = hclust(frame().set("x", parse_vec("1 2 3")), "average");\nout.print(h.cut(4).str());', "k must be between 1 and n"),
    ('let k = kmeans(frame().set("x", parse_vec("1 2 9")), 2, rng(1));\nout.print(str(k.size(3)));', "cluster label 3 is outside 1..2"),
]:
    case("cluster error: " + msg[:60], body, code=3, out=(msg,))

# ---------------------------------------------------------------- matrices, PCA, MDS
case("Mat construction, access and algebra", r"""
let m = matrix(parse_vec("1 2 3 4"), 2, false);
out.check(m.at(0, 1) == 3.0 && m.at(1, 0) == 2.0, "matrix fills by column like R");
let r = matrix(parse_vec("1 2 3 4"), 2, true);
out.check(r.at(0, 1) == 2.0 && r.t().at(1, 0) == 2.0, "byrow and transpose");
out.check(m.det() == -2.0, "determinant");
let inv = m.inverse();
let id = m.mul(inv);
out.check(abs(id.at(0, 0) - 1.0) < 0.000000000001 && abs(id.at(0, 1)) < 0.000000000001, "inverse");
let x = m.solve(parse_vec("5 6"));
out.check(abs(x.at(0) + 1.0) < 0.000000000001 && abs(x.at(1) - 2.0) < 0.000000000001, "solve");
out.check(m.values().str() == "[1, 2, 3, 4]" && m.row(1).str() == "[2, 4]" && m.col(1).str() == "[3, 4]", "values row col");
out.check(identity(2).mul_vec(parse_vec("7 8")).str() == "[7, 8]" && diag(parse_vec("2 3")).diag().str() == "[2, 3]", "identity diag");
out.check(m.add(m).sub(m).scale(2.0).at(1, 1) == 8.0, "add sub scale");
out.set("frame_cols", m.to_frame().ncol());
""", values={"frame_cols": 2})
case("eigen, svd, chol", r"""
let a = matrix(parse_vec("4 1 1 3"), 2, false);
let e = a.eigen();
out.check(abs(e.values().at(0) - (3.5 + sqrt(1.25))) < 0.000000000001, "largest eigenvalue");
let v = e.vectors().col(0);
let av = a.mul_vec(v);
out.check(abs(av.at(0) - e.values().at(0) * v.at(0)) < 0.000000000001, "A v = lambda v");
let s = matrix(parse_vec("1 2 3 4 5 6"), 3, false).svd();
out.check(s.d().len() == 2 && s.u().nrow() == 3 && s.v().nrow() == 2, "svd shapes");
let r = a.chol();
out.check(abs(r.t().mul(r).at(0, 1) - 1.0) < 0.000000000001 && r.at(1, 0) == 0.0, "chol is upper, R'R = A");
""")
case("prcomp, cmdscale, mahalanobis", r"""
let f = frame().set("a", parse_vec("1 2 3 4 5")).set("b", parse_vec("2 1 4 3 6"));
let p = prcomp(f, true);
out.check(p.sdev().len() == 2 && abs(p.cumulative().at(1) - 1.0) < 0.000000001, "two components, all variance");
out.check(p.summary() != "", "summary text");
out.check(abs(abs(p.loading("a", 0)) - sqrt(0.5)) < 0.000000000001, "scaled 2-variable PCA loads equally");
let d = as_matrix(f).dist("euclidean");
let m = cmdscale(d, 1);
out.check(m.points().ncol() == 1 && m.eig().len() == 5, "cmdscale");
let x = as_matrix(f);
let md = mahalanobis(x, parse_vec("3 3.2"), x.cov());
out.check(md.len() == 5 && md.min() >= 0.0, "mahalanobis");
""")
for body, msg in [
    ('let m = matrix(parse_vec("1 2 3"), 2, false);', "is not a positive multiple of nrow 2"),
    ('let m = matrix(parse_vec("1 2 2 4"), 2, false).inverse();', "singular"),
    ('let m = matrix(parse_vec("1 2 3 4 5 6"), 2, false).det();', "must be square"),
    ('let m = matrix(parse_vec("1 2 3 4"), 2, false).mul(matrix(parse_vec("1 2 3"), 3, false));', "non-conformable"),
    ('let m = matrix(parse_vec("1 2 3 4"), 2, false).at(2, 0);', "out of range"),
    ('let p = prcomp(frame().set("a", parse_vec("1 1 1")).set("b", parse_vec("1 2 3")), true);', "cannot rescale a constant/zero column"),
    ('let m = as_matrix(frame().set("t", str_split("a,b", ",")));', "is text; select the numeric columns"),
    ('let m = matrix(parse_vec("1 2 2 1"), 2, false).chol();', "not positive"),
]:
    case("matrix error: " + msg, body, code=3, out=(msg,))

# ---------------------------------------------------------------- factors and anova
case("factor formulas and anova", r"""
let f = frame().set("y", parse_vec("1 2 3 5 4 6 8 9")).set("g", str_split("a,a,b,b,c,c,a,b", ","))
  .set("d", parse_vec("1 2 1 2 1 2 1 2"));
let m = lm(f, "y ~ g + factor(d)");
out.check(m.n_coef() == 4 && m.name(1) == "gb" && m.name(3) == "factor(d)2", "R's column names");
let at = m.anova();
out.check(at.n_rows() == 3 && at.row_name(2) == "Residuals" && at.columns().at(4) == "Pr(>F)", "anova shape");
let ss = at.col("Sum Sq").sum();
let y = f.col("y");
out.check(abs(ss - y.var() * 7.0) < 0.000000001, "sequential SS add up to the total");
let small = lm(f, "y ~ g");
let n = anova(small, m);
out.check(n.n_rows() == 2 && n.at(1, "Df") == 1.0, "nested comparison");
let g = glm(f, "y ~ g", "poisson");
out.check(g.anova("Chisq").at("g", "Df") == 2.0, "analysis of deviance");
out.check(m.predict(frame().set("g", str_split("c", ",")).set("d", parse_vec("2"))).len() == 1, "predict with factor levels");
""")
for body, msg in [
    ('let m = lm(frame().set("y", parse_vec("1 2 3")).set("g", str_split("a,b,a", ",")), "y ~ g"); let p = m.predict(frame().set("g", str_split("z", ",")));', "has new level z"),
    ('let m = lm(frame().set("y", parse_vec("1 2 3")).set("x", parse_vec("1 2 3")), "y ~ x"); let a = m.anova().at("nope", "Df");', "nope"),
]:
    case("factor error: " + msg, body, code=3, out=(msg,))

# ---------------------------------------------------------------- k-sample and block tests
case("k-sample and block tests", r"""
let x = parse_vec("1 2 3 4 5 6 7 8 9 10 11 12");
let g = parse_vec("1 1 1 1 2 2 2 2 3 3 3 3");
let k = test("kruskal").x(x).groups(g).run();
out.check(k.method() == "Kruskal-Wallis rank sum test" && k.df() == 2.0 && k.p_value() < 0.02, "kruskal");
let s = test("kruskal").x(x).groups(str_split("b,b,b,b,a,a,a,a,c,c,c,c", ",")).run();
out.check(s.statistic() == k.statistic(), "text groups");
out.check(test("bartlett").x(x.pow(2.0)).groups(g).run().p_value() < 1.0, "bartlett");
out.check(test("fligner").x(x).groups(g).run().statistic_name() == "Fligner-Killeen:med chi-squared", "fligner");
out.check(test("oneway").x(x).groups(g).var_equal(true).run().param("num df") == 2.0, "oneway");
let f = test("friedman").table(parse_vec("1 2 3 2 3 1 3 1 2 1 3 2"), 4, true).run();
out.check(f.df() == 2.0, "friedman from a blocks x groups table");
let f2 = test("friedman").x(parse_vec("1 2 3 2 3 1 3 1 2 1 3 2")).groups(parse_vec("0 1 2 0 1 2 0 1 2 0 1 2")).blocks(parse_vec("0 0 0 1 1 1 2 2 2 3 3 3")).run();
out.check(f2.statistic() == f.statistic(), "friedman from vectors");
out.check(test("quade").table(parse_vec("1 2 3 2 3 1 3 1 2 1 3 2"), 4, true).run().method() == "Quade test", "quade");
out.check(test("mcnemar").table(parse_vec("10 5 2 8"), 2, false).run().df() == 1.0, "mcnemar");
out.check(test("ansari").x(parse_vec("1 2 3 4 5")).y(parse_vec("2 4 6 8 10 12")).run().statistic_name() == "AB", "ansari");
""")
for body, msg in [
    ('let t = test("kruskal").x(parse_vec("1 2")).run();', "call .groups(...) before run()"),
    ('let t = test("friedman").x(parse_vec("1 2")).groups(parse_vec("0 1")).run();', "call .blocks"),
    ('let t = test("mcnemar").x(parse_vec("1 2")).run();', "call .table(values, nrow, byrow) before run()"),
    ('let t = test("kruskal").x(parse_vec("1 2 3")).groups(parse_vec("1 1 1")).run();', "all observations are in the same group"),
]:
    case("test error: " + msg, body, code=3, out=(msg,))

# ---------------------------------------------------------------- density, ecdf, hist, noncentral and tukey
case("density, bandwidths, ecdf, hist, pretty", r"""
let x = parse_vec("1 2 2 3 3 3 4 4 5 7");
let d = density(x);
out.check(d.x().len() == 512 && d.y().len() == 512 && d.bw() == bw(x, "nrd0"), "grid and default bandwidth");
out.check(abs(d.at(3.0) - d.eval(parse_vec("3")).at(0)) < 0.000000000000001 && d.kernel() == "gaussian", "at and eval");
let e = density(x, 0.5, "epanechnikov");
out.check(e.bw() == 0.5 && e.at(100.0) == 0.0, "explicit bandwidth, compact kernel");
out.check(density(x, "SJ", "biweight", 1.5).bw() > 0.0, "rule, kernel and adjust");
let F = ecdf(x);
out.check(F.at(3.0) == 0.6 && F.at(0.0) == 0.0 && F.knots().len() == 6, "ecdf");
let h = hist(x);
out.check(h.counts().sum() == 10.0 && h.breaks().len() == h.counts().len() + 1, "hist");
out.check(hist(x, parse_vec("0 4 8")).counts().str() == "[8, 2]", "hist with breaks (right-closed)");
out.check(pretty(1.0, 15.0, 5).str() == "[0, 2, 4, 6, 8, 10, 12, 14, 16]", "pretty");
""")
case("noncentral and studentized range distributions", r"""
out.check(abs(t_dist(10.0, 0.0).p(1.0) - t_dist(10.0).p(1.0)) < 0.000000000000001, "ncp 0 is central t");
out.check(abs(chisq(3.0, 2.0).mean() - 5.0) < 0.000000000001, "noncentral chi-square mean df + ncp");
let q = f_dist(2.0, 30.0, 3.0).q(0.9);
out.check(abs(f_dist(2.0, 30.0, 3.0).p(q) - 0.9) < 0.000000001, "F quantile round trip");
out.check(abs(tukey(2.0, 10.0).p(3.0) - (2.0 * t_dist(10.0).p(3.0 / sqrt(2.0)) - 1.0)) < 0.00000000001, "two means: ptukey from pt");
out.check(dist("tukey", 3.0, 20.0).name() == "tukey", "by name");
""")
for body, msg in [
    ('let d = tukey(3.0, 20.0).d(1.0);', "dtukey is not available"),
    ('let d = tukey(1.0, 20.0);', "tukey"),
    ('let d = density(parse_vec("1 2 3"), "nonsense");', "nonsense"),
    ('let d = density(parse_vec("1 2 3"), 0.5, "box");', "box"),
    ('let d = density(parse_vec("1 2 3"), -1.0, "gaussian");', "bw must be positive"),
    ('let q = ecdf(parse_vec("1 2")).quantile(2.0);', "p must be in [0, 1]"),
    ('let r = rng(1).draw(tukey(3.0, 10.0));', "studentized range"),
]:
    case("density/dist error: " + msg, body, code=3, out=(msg,))

# ---------------------------------------------------------------- adversarial
case("formula transforms through a frame", r"""
let x = seq(1.0, 20.0, 1.0);
let f = frame().set("x", x).set("my.z", x.sqrt()).set("y", x.pow(2.0).mul(3.0).add(1.0).log());
let m = lm(f, "exp(y) ~ x + I(x^2)");
out.check(abs(m.estimate("I(x^2)") - 3.0) < 0.000001 && abs(m.estimate("(Intercept)") - 1.0) < 0.00001, "I() and exp() response");
let m2 = lm(f, "y ~ log(x) + my.z");
out.check(m2.has_coef("log(x)") && m2.has_coef("my.z"), "dotted column names and log()");
let m3 = lm(f, "y ~ x * my.z");
out.check(m3.has_coef("x:my.z") && m3.n_coef() == 4, "a*b expands");
let m4 = lm(f, "y ~ x - 1");
out.check(m4.n_coef() == 1 && !m4.has_coef("(Intercept)"), "no intercept");
let m5 = lm(f, "y ~ 1");
out.check(m5.n_coef() == 1 && abs(m5.estimate(0) - f.col("y").mean()) < 0.000000001, "intercept only");
""")
case("binomial glm with proportions and weights", r"""
let f = frame().set("p", parse_vec("0.1 0.3 0.5 0.8 0.9")).set("d", parse_vec("1 2 3 4 5")).set("n", parse_vec("10 10 20 10 10"));
let m = glm(f, "p ~ d", "binomial", "logit", f.col("n"));
out.check(abs(m.estimate("d") - 1.1005669234828539) < 0.000000001 && abs(m.estimate(0) + 3.2068016632558991) < 0.000000001, "R: glm(p ~ d, binomial, weights = n)");
out.check(m.df_residual() == 3 && m.dispersion() == 1.0, "df and dispersion");
let pr = m.predict(frame().set("d", parse_vec("3")), "response");
out.check(pr.at(0) > 0.0 && pr.at(0) < 1.0, "response scale");
""")
case("runtime error while a solver waits does not hang", r"""
let o = minimize(parse_vec("1"), "BFGS");
let zero = 0;
while (o.next()) {
  o.tell(o.x(0) * o.x(0));
  let bad = 1 / zero;
}
""", code=3, out=("division by zero",))
case("solver error propagates (BFGS with non-finite objective)", r"""
let o = minimize(parse_vec("1"), "BFGS");
while (o.next()) {
  o.tell(nan());
}
""", code=3, out=('"phase": "run"',))
case("Nelder-Mead tolerates a NaN objective", r"""
let o = minimize(parse_vec("1 1"), "Nelder-Mead").maxit(200);
while (o.next()) {
  let v = o.x(0) * o.x(0) + o.x(1) * o.x(1);
  if (o.x(0) > 1.5) {
    v = nan();
  }
  o.tell(v);
}
out.check(o.value() < 0.001, "NaN treated as +Inf");
""")
case("optimize with a constant function terminates", r"""
let o = optimize(0.0, 1.0);
while (o.next()) {
  o.tell(1.0);
}
out.check(o.evaluations() < 100, "bounded evaluations");
""")
case("uniroot root at an endpoint", r"""
let r = uniroot(0.0, 1.0);
while (r.next()) {
  r.tell(r.x());
}
out.check(r.root() == 0.0 && r.f_root() == 0.0, "exact endpoint root");
""")
case("integrate a function that is NaN only mid-way", r"""
let q = integrate(0.0, 2.0);
while (q.next()) {
  let v = 1.0;
  if (q.x() > 1.0) {
    v = nan();
  }
  q.tell(v);
}
out.check(q.ier() == 6, "flagged");
""")
case("unicode survives CSV and output", r"""
let f = read_csv("u.csv");
out.check(f.text("name").at(0) == "雪 snow ❄", "utf-8 text");
out.set("name", f.text("name").at(0));
f.write_csv("u2.csv");
out.check(read_csv("u2.csv").text("name").at(0) == "雪 snow ❄", "round trip");
""", files={"u.csv": "name,v\n雪 snow ❄,1\n"}, values={"name": "雪 snow ❄"})
case("large data stays fast", r"""
let r = rng(1);
let n = 200000;
let x = r.draw(norm(), n);
let y = x * 2.0 + r.draw(norm(), n);
let f = frame().set("x", x).set("y", y);
f.write_csv("big.csv");
let g = read_csv("big.csv");
let m = lm(g, "y ~ x");
out.check(abs(m.estimate("x") - 2.0) < 0.02 && g.nrow() == n, "200k rows through CSV and lm");
out.check(abs(x.quantile(0.5)) < 0.02 && x.sort().len() == n, "quantile and sort");
let s = test("t").x(x).y(y).run();
out.check(s.p_value() >= 0.0, "t-test on 200k");
""")
case("frame with no rows", r"""
let f = read_csv("h.csv");
out.check(f.nrow() == 0 && f.ncol() == 2 && !f.is_text("a"), "header only");
""", files={"h.csv": "a,b\n"})
case("model on an empty frame errors", 'let m = lm(read_csv("h.csv"), "a ~ b");', code=3, out=("lm",), files={"h.csv": "a,b\n"})
case("t-test drops NA like R", 'let t = test("t").x(parse_vec("1 NaN 3")).run();\n'
     'out.check(t.statistic() == 2.0 && abs(t.p_value() - 0.29516723530086642) < 0.000000000001, "R: t.test(c(1, NA, 3))");')
case("out.set of a large vector", 'out.set("v", seq_len(100000)); out.set("n", 1);', values={"n": 1})
case("deep chain of builder calls", r"""
let t = test("t").x(parse_vec("1 2 3 4")).y(parse_vec("2 3 4 5")).alt("less").alt("two.sided").conf(0.8).conf(0.9).var_equal(true).run();
out.check(t.conf_level() == 0.9 && t.alternative() == "two.sided" && t.df() == 6.0, "last setting wins");
""")

# ---------------------------------------------------------------- CLI
case("cli: no arguments", None, args=[], code=2, err=("usage: lmr",))
case("cli: unknown flag", None, args=["--frobnicate"], code=2, err=("usage",))
case("cli: --dsl without files", None, args=["--dsl"], code=2, err=("usage",))
case("cli: flag among files", None, args=["--dsl", "a.dslb", "--x"], code=2, err=("usage",))
case("cli: version", None, args=["--version"], out=("lmr",))
case("cli: schema summary", None, args=["--dsl-schema"], out=("context out: Out", "Vec:", "fn lm(Frame, string) -> Lm"))
case("cli: --schema alias", None, args=["--schema"], out=("context out: Out",))
case("cli: dsl-schema json", None, args=["--dsl-schema", "--json"], out=('"dsl.schema": 2', '"name": "lmr"'))
case("cli: dsl-schema bad option", None, args=["--dsl-schema", "--yaml"], code=2, err=("usage: lmr",))
case("cli: missing bytecode", None, args=["--dsl", "nope.dslb"], code=3, out=('"phase": "load"', "nope.dslb"))
case("cli: garbage bytecode", None, args=["--dsl", "junk.dslb"], code=3, out=('"phase": "load"',),
     files={"junk.dslb": "not bytecode at all"})


# ---------------------------------------------------------------- runner
def compile_chunk(lmscript, body, path):
    src = path + ".lms"
    with open(src, "w", encoding="utf-8") as f:
        # Forward slashes: a Windows path's backslashes would be escapes inside the lmscript string literals.
        f.write('p = _lmcode_ schema="%s";\n%s\nend _lmcode_\nsavedsl(p, "%s");\n' % (SCHEMA.replace("\\", "/"), body, path.replace("\\", "/")))
    r = subprocess.run([lmscript, src], capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=120)
    if r.returncode != 0:
        raise RuntimeError("compile failed:\n" + r.stdout + r.stderr)


def run_case(binary, lmscript, c, root):
    work = tempfile.mkdtemp(prefix="case-", dir=root)
    try:
        for name, text in c.files.items():
            full = os.path.join(work, name)
            os.makedirs(os.path.dirname(full), exist_ok=True)
            with open(full, "w", newline="", encoding="utf-8") as f:
                f.write(text)
        if c.args is not None:
            argv = [binary] + c.args
        else:
            bodies = c.chunks if c.chunks is not None else [c.body]
            paths = []
            for i, b in enumerate(bodies):
                p = os.path.join(work, "chunk%d.dslb" % i)
                compile_chunk(lmscript, b, p)
                paths.append(os.path.basename(p))
            argv = [binary, "--dsl"] + paths
        r = subprocess.run(argv, capture_output=True, text=True, encoding="utf-8", errors="replace", cwd=work, timeout=120)
        problems = []
        if r.returncode != c.code:
            problems.append("exit %d, expected %d" % (r.returncode, c.code))
        for s in c.out:
            if s not in r.stdout:
                problems.append("stdout lacks %r" % s)
        for s in c.err:
            if s not in r.stderr:
                problems.append("stderr lacks %r" % s)
        if c.args is None:
            lines = [ln for ln in r.stdout.splitlines() if ln.strip()]
            if len(lines) != 1:
                problems.append("expected exactly one stdout line, got %d" % len(lines))
            elif c.values:
                try:
                    got = json.loads(lines[0])["result"]["values"]
                    for k, v in c.values.items():
                        if got.get(k) != v:
                            problems.append("value %s = %r, expected %r" % (k, got.get(k), v))
                except (KeyError, TypeError, ValueError) as e:
                    problems.append("bad result line: %s" % e)
        if problems:
            return c.name, problems + ["stdout: " + r.stdout.strip()[:600], "stderr: " + r.stderr.strip()[:600]]
        return c.name, None
    except Exception as e:  # noqa: BLE001 - report every harness failure as a test failure
        return c.name, [str(e)]
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    binary = os.path.abspath(sys.argv[1])
    only = sys.argv[3] if len(sys.argv) > 3 and sys.argv[2] == "-k" else None
    lmscript = os.environ.get("LMSCRIPT_BIN") or shutil.which("lmscript")
    if not lmscript:
        print("lmscript not found: host tests skipped")
        return 0
    cmake = Path(HERE).parent / "CMakeLists.txt"
    project_version = re.search(r"^project\(lmr VERSION ([^ )]+)", cmake.read_text(), re.M).group(1)
    version_result = subprocess.run([binary, "--version"], capture_output=True)
    if version_result.returncode or version_result.stdout != f"lmr {project_version}\n".encode() or version_result.stderr:
        print("FAIL cli: version bytes", repr(version_result.stdout), repr(version_result.stderr), version_result.returncode)
        return 1
    cases = [c for c in CASES if not only or only in c.name]
    root = tempfile.mkdtemp(prefix="rcli-host-test-")
    failed = 0
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=os.cpu_count() or 2) as pool:
            for name, problems in pool.map(lambda c: run_case(binary, lmscript, c, root), cases):
                if problems:
                    failed += 1
                    print("FAIL", name)
                    for p in problems:
                        print("     ", p.replace("\n", "\n      "))
    finally:
        shutil.rmtree(root, ignore_errors=True)
    print("host tests: %d passed, %d failed" % (len(cases) - failed, failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
