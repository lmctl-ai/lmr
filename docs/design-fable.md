# rcli: R's statistics as lmapi native ops (design)

Decision: build `rcli` as an lmcpp C++17 native tool whose handlers run an **embedded R op catalogue through `Rscript --vanilla` as a child process**. The C++ side owns schema, validation, policy, sandboxing, NA mapping and the JSON envelope; R owns the mathematics. Do **not** link libR or libRmath in v1. Primary value: bit-exact R reference results (quantiles, RNG streams, fits) usable as a differential oracle for lmstats/lmda, plus algorithms those tools lack (nls, smooth.spline, GAM, optim/uniroot/integrate on safe expressions). This is a design, not a build report.

## 1. Which R layer, and why (verified on this host)

Facts checked on this aarch64 Ubuntu host (R 4.3.3-2build2, `/usr/lib/R/lib/libR.so`, `/usr/share/R/include`):
- `stats.so` exports only Fortran entry points (`dqrls_`, `lminfl_`, `hclust_`, loess `ehg*`, PORT `d*7*`); the `.Call` routines (`C_Cdqrls`, `C_optim`, `C_ARIMA_Like`, `C_kmeans_Hartigan`, …) are registered, not exported. Reaching them means vendoring R's GPL-2 C sources, and the algorithms' glue (glm IRLS loop, `arima` CSS-ML driver, `loess` surface logic, `kmeans` restarts, `hclust` merge ordering) is R code. A "C layer" wrapper would reimplement R and then need R as its oracle anyway. **(d) rejected.**
- `libR.so` exports `Rf_initEmbeddedR`, `R_tryEval`, `R_ParseVector`, `vmmin/nmmin/lbfgsb`, `dqrls_`, `Rdqags`. Embedding works but costs SEXP/longjmp handling, `R_HOME` setup, single-thread constraints, and makes the binary a GPL-2 combined work. It saves only the process spawn: measured `Rscript --vanilla -e 'invisible(1)'` is 0.19 s with default packages and **0.06 s with `R_DEFAULT_PACKAGES=NULL`** (0.11 s including `library(stats)` + `lm` on 1e4 rows, 20 MB RSS). **(b) rejected.**
- `r-mathlib` (apt candidate 4.3.3-2build2, not installed) ships `libRmath.so.1.0.0` (264 KB) and `libRmath.a` (456 KB), 155 symbols. `Rmath.h` is LGPL-2.1 but the nmath sources are GPL-2 (R's COPYRIGHTS lists only `R.h Rdefines.h Rgraphics.h Rinternals.h Rmath.h S.h R_ext/*.h` as LGPL). Its standalone `unif_rand` is Marsaglia-multicarry with `set_seed(u1,u2)`, **not** R's Mersenne-Twister, so it cannot reproduce `set.seed(42); rnorm()` streams. Good for d/p/q/special functions only. **(a) is the optional tier 2, not v1.**
- BLAS/LAPACK/LINPACK (c): generic LAPACK adds nothing R-specific; lmda already has faer/nalgebra. Only relevant as an oracle for decompositions; defer to stage 2 as one op.
- Rscript as child (e): is R itself, so every result is bit-exact R by construction; MIT-clean; 60 ms startup is below lmscript's own per-call overhead. **Chosen.**

Tier 2 (optional, after stage 2, only if a measured lmscript loop needs sub-millisecond `dist_*` calls): a separate `rmath` build target that statically links `libRmath.a` for `dist_d/p/q` and `math_special` in-process. That binary is a GPL-2 combined work and must live in its own GPL-2 repo or stay private. Cut from v1.

### Overlap with lmstats (45 ops) and lmda (199 ops)
Duplicated as *features* and therefore included only as **oracle ops**: lm (lmda `model_regress`, lmstats `regression`), glm (lmstats `glm`, lmda `model_logit/poisreg`), arima (lmstats `arima`), loess (lmda `model_loess`), classical htests (lmda `test_*`), fft (lmda `series_spectrum`), quantiles (lmda `summary_describe`). Neither tool has: `nls`, `smooth.spline`, `mgcv::gam`, `optim/nlm/uniroot/integrate` over a user expression, R's 9 quantile types, `ptukey/qtukey/pwilcox/psignrank/pnbeta/pnt`, Bessel/`digamma`/`lchoose`, `set.seed` streams and `sample()` (Rejection sampling), `arima.sim`, `approx/spline/splinefun` interpolation. v1 is 21 ops: 9 oracle-grade reference ops and 12 new-algorithm ops.

| Op | R core | Why it earns a slot |
|---|---|---|
| `dist_d`, `dist_p`, `dist_q` | `d/p/q{norm,t,chisq,f,beta,gamma,binom,pois,nbinom,hyper,unif,exp,weibull,lnorm,logis,cauchy,geom,signrank,wilcox,tukey}` + `ncp` | bit-exact quantiles incl. distributions missing from lmda/lmstats; `lower_tail`, `log_p`, vectorized `x` |
| `rng_draw` | `RNGkind("Mersenne-Twister","Inversion","Rejection"); set.seed(seed); r<dist>(n, …)` | the R RNG stream oracle; `rmultinom` too |
| `rng_sample` | `sample(x, size, replace, prob)` | seeded permutations/bootstraps reproducible from any R session |
| `rng_arima_sim` | `arima.sim` | seeded ARIMA generators for lmstats differential tests |
| `math_special` | `gamma lgamma digamma trigamma beta lbeta choose lchoose besselI/J/K/Y log1p expm1` | special functions as oracle for lmstats' hand-written ones |
| `summary_quantile` | `quantile(type=1..9)`, `fivenum`, `mad`, `IQR` | all nine quantile definitions, exact |
| `fit_lm` | `lm` (formula, weights, offset, contrasts) + `summary`, `anova`, `influence` | oracle for `model_regress`/`regression`; QR rank, aliased terms, hat/Cook's |
| `fit_glm` | `glm` families gaussian/binomial/poisson/Gamma/inverse.gaussian/quasi*, links, offsets | oracle for lmstats `glm`; IRLS trace, deviance, AIC, `converged` as domain status |
| `fit_nls` | `nls` (Gauss-Newton/port/plinear), `start`, bounds | new |
| `fit_arima` | `arima` (CSS-ML, Kalman exact), `predict(n.ahead)` | oracle for lmstats `arima` loglik/coef/se |
| `fit_gam` | `mgcv::gam` with `s()/te()`, REML/GCV, `predict` | new; mgcv ships with r-base recommended |
| `smooth_spline` | `smooth.spline` (GCV, `spar`, `df`), `predict(deriv)` | new |
| `smooth_loess` | `loess` (`span`, `degree`, `family`, surface exact/interpolate) | oracle for lmda `model_loess` |
| `smooth_interp` | `approx`, `spline`/`splinefun` (fmm/natural/periodic/monoH.FC/hyman) | new |
| `opt_minimize` | `optim` (Nelder-Mead/BFGS/CG/L-BFGS-B/SANN), `nlm` | new; objective is a sandboxed expression (§2) |
| `num_root` | `uniroot`, `polyroot` | new |
| `num_integrate` | `integrate` (QUADPACK dqags/dqagi, infinite bounds) | new |
| `num_fft` | `fft`, `spec.pgram` | oracle for lmda `series_spectrum` |
| `test_htest` | `t/wilcox(exact,conf.int)/ks/shapiro/fisher(r×c)/chisq/cor(kendall exact)/var/bartlett/kruskal/mantelhaen/mcnemar/prop/binom/poisson.test` | one op returning R's uniform `htest` shape; oracle for lmda `test_*` |

Stage 2 candidates: `cluster_kmeans` (Hartigan-Wong, seeded `nstart`), `cluster_hclust` + `cutree`, `linalg_decomp` (qr/svd/eigen/chol/solve/det via R's LAPACK), `fit_coxph` (survival), `fit_lme` (nlme), `multi_prcomp`.

## 2. Interface on lmcpp

**Process model.** `rcli <op> …` → lmcpp validates → handler serializes params to `Context::temp_file_path()` (`args.json`), writes the embedded R program (one C++ raw-string literal per op plus `prelude.R`) to a second temp file, execs `Rscript --vanilla --default-packages=NULL <prog.R> <args.json> <out.json>` with a scrubbed environment (`LC_ALL=C.UTF-8`, `OMP_NUM_THREADS=1`, `R_LIBS_USER=` empty, `HOME` set to a temp dir, no network), waits with a timeout (default 50 s, `RCLI_TIMEOUT_S`), then `read_result_file(out.json,"json")` and prints. No shell strings ever; argv arrays only.

**Naming.** Flat snake_case `group_name` as in the table; `tool:"rcli"`, `inputs:["flags","json"]`, `returns:"object"`, `product` omitted (native tool). All 21 v1 ops are `access:"read"`, none dangerous.

**Data handoff.** Inputs: numeric vectors/matrices inline as JSON (`Type::Array` of numbers for vectors; `Type::Json` for row-major matrices and `data` records) **or** `file` (CSV path, header row) + `formula` string referencing column names. CSV is the cross-tool table format: lmda reads/writes it, lmcreator `render_chart` reads it, R reads it natively with `read.csv`. Arrow needs the R `arrow` package (large, not in Debian main); safetensors has no R reader and is tensor-shaped, not table-shaped. JSON inline is fine for the sizes v1 targets (lmstats already does this). Outputs: inline JSON. Tabular results use one envelope everywhere: `{"columns":["term","estimate","std_error","t","p"],"types":["string","number",…],"rows":[[…],…],"n_rows":5}`; vectors are plain arrays; matrices are row-major arrays of arrays with `dim`. Results are capped at lmcpp's 16 MiB file bound; handlers additionally enforce `max_values` (default 200 000 numbers) and fail with `COMMAND_FAILED`/`tool_code:"OUTPUT_TOO_LARGE"` before serialization. Large fitted/residual/prediction outputs are opt-in via `return:["coefficients","fitted","residuals"]`. Writing result files is deferred to stage 2 as explicit `write` ops (`table_export`) copied from torchcli's atomic-install pattern; keeping all v1 ops `read` keeps `--readonly` useful.

**Numbers.** R emits doubles with `sprintf("%.17g")` (round-trip exact) via an in-repo 60-line emitter; input parsing uses `jsonlite::fromJSON(simplifyVector=TRUE)` (apt `r-cran-jsonlite`, MIT). Input `null` inside numeric arrays means `NA`; strings are rejected in numeric positions. Output: `NA`, `NaN`, `Inf`, `-Inf` all become JSON `null` (same rule as torchcli) with a sibling `nonfinite:{"na":k,"nan":k,"inf":k,"ninf":k}` per vector/column so nothing is silently lost and columns stay numerically typed. Integers from R stay integers; logical stays boolean.

**Seed and reproducibility.** Every stochastic op has required integer `seed` and optional `sample_kind` (`Rejection` default, `Rounding` for pre-3.6 compatibility). The prelude always sets `RNGkind` explicitly before `set.seed`; it never inherits `.Random.seed`. Every result carries `meta:{"r_version":"4.3.3","platform":"aarch64-unknown-linux-gnu","blas":"<from sessionInfo()>","rng":["Mersenne-Twister","Inversion","Rejection"],"seed":42,"elapsed_ms":n}`. Reproducibility promise: bit-exact for the pinned (R version, platform, BLAS) triple recorded in `meta`; last-ulp differences across BLAS builds are documented, not promised away.

**Errors and warnings.** The prelude wraps the op in `withCallingHandlers(tryCatch(...))`. R `stop()` → `COMMAND_FAILED` exit 1 with `data.tool_code:"R_ERROR"`, `data.r_message`, `data.r_call`. Warnings are collected into `result.warnings:[string]`, never dropped to stderr. Non-convergence (`glm$converged`, `optim$convergence`, `nls` step failures, `integrate` message) is a domain status inside the result per protocol §6.1, not an exit code. Domain errors we can detect before R (empty vectors, `df<=0`, `p` outside [0,1] when `log_p` false, non-square matrix) are lmcpp `VALIDATION_ERROR` with precise `tool_code`. Rscript missing → `SERVICE_UNAVAILABLE` exit 8 `tool_code:"R_NOT_FOUND"`; timeout → `TIMEOUT` exit 7; signal/crash → `INTERNAL_ERROR` `tool_code:"TOOL_CRASH"`.

**Expression sandbox** (`opt_minimize`, `num_root`, `num_integrate`). The objective is a string in R expression syntax over named parameters, e.g. `"(x[1]-1)^2 + (x[2]-2)^2"`. The prelude parses it with `parse(text=, keep.source=FALSE)`, walks `all.names()` and rejects anything outside an allowlist (arithmetic, comparison, `[`, `(`, `sum prod mean length abs sqrt exp log log1p expm1 sin cos tan atan tanh pmin pmax ifelse dnorm pnorm …`), rejects assignment, `function`, `::`, backticks, strings and calls with more than 10 000 tokens, then evaluates in `new.env(parent=emptyenv())` populated only with allowlisted functions and the `data` vectors. C++ rejects the string if it contains any byte outside printable ASCII. This is the single security-sensitive surface in rcli and gets its own negative test corpus.

**lmcpp gaps to add (lmcpp repo, kept generic):**
1. `lmcpp::run_child(argv, env, stdin_bytes, timeout, limits)` — fork/exec without shell, env allowlist, bounded stdout/stderr capture, timeout kill, exit/signal → `Error` mapping. Listed as deferred "subprocess helpers" in lmcpp parity.md; this is its first consumer.
2. `Type::Array` must accept `null` items for numeric arrays (NA); currently item typing is strict. Add `nullable_items` on `Param`.
3. A second RAII temp path per `Context` (currently one `temp_file_path()`); make it `temp_file_path(name)` or a temp directory.
4. Table envelope helper `lmcpp::table(columns, types, rows)` and a documented convention in lmcpp docs; no schema change.
5. Not needed: GraphQL, CSV result adapter, CLI11 bridge, large-result files (stage 2 only).

## 3. Build, packaging, testing, CI

**Binary.** `rcli` links only `lmcpp::core` and nlohmann-json; it is a plain host binary (native build on this host in under a minute). R scripts are embedded at compile time (CMake `file(READ)` → generated `ops_embedded.cpp`), so one executable carries its catalogue and `--schema` needs no R.

**Docker (canonical, pinned).** Base `ubuntu:24.04` arm64, digest-pinned, so the container carries exactly the host's `r-base-core 4.3.3-2build2` and reference `libblas3/liblapack3`; add `r-cran-jsonlite`, `r-cran-mgcv`, `r-cran-nlme`, `r-cran-survival` (all GPL, unchanged Debian binaries). Builder stage adds `g++ cmake ninja-build nlohmann-json3-dev`, builds lmcpp from the sibling checkout (same `~/repos` context + `.dockerignore` allowlist as torchcli), then rcli. Test stage runs ctest; runtime stage is ubuntu + R + `/usr/local/bin/rcli` as entrypoint. Budget: runtime image 400–500 MB unpacked (R + ICU + Tcl/Tk pulled by r-base-core), build 4–6 min cold on this host, dominated by apt; no compiler heavier than lmcpp. Record real numbers at stage 1.

**Shim.** Copy `bin/torchcli` verbatim as `bin/rcli` (`--network=none --cpus=2 --memory=2g`, cwd mount below `$HOME`, uid/gid, `RCLI_READONLY/RCLI_DRYRUN` passthrough, Docker-failure → `SERVICE_UNAVAILABLE`, crash → `TOOL_CRASH`). Because the host has the same R build, developers may also put `build/bin/rcli` on PATH directly; CI and recipes use the shim. Docker adds about 0.3–0.5 s per call; the recipe budget (60 s per call) is unaffected.

**Testing.** (1) *Glue oracle*: for every op, a CTest Rscript fixture computes the same thing directly in R and compares to rcli's JSON after `%.17g` round-trip: **exact equality** for `dist_*`, `rng_*`, `math_special`, `summary_quantile`, `num_fft`, `smooth_interp`, and for deterministic fits on the same image; `all.equal(tolerance=1e-12)` only where R's own iteration order can differ (`optim` SANN, `fit_gam` with REML). (2) *Envelope tests*: NA/NaN/Inf mapping, warnings capture, `R_ERROR` mapping, timeout, missing Rscript, `OUTPUT_TOO_LARGE`, readonly/dryrun floors, schema offline. (3) *Sandbox corpus*: 40+ rejected expressions (`system("id")`, `eval`, `base::`, `<-`, `function(x)`, `.Internal`, `\u` escapes, lambda `\(x)`, deep nesting) must fail with `tool_code:"EXPR_REJECTED"` and never spawn. (4) *Differential vs lmstats/lmda* (stage 2): JSON report, not a gate, with policy coefficients/statistics 1e-8 relative, p-values 1e-6 absolute, iterative fits (glm/arima/nls) 1e-4 relative on log-likelihood; each documented divergence gets an issue in the Rust repo. (5) lmapi conformance inherited from lmcpp.

**CI.** GitHub Actions on `ubuntu-24.04-arm` hosted runners (public repos): build test + runtime stages, run ctest, run `lmscript --check` on every recipe, upload `cli.schema.dev.json` and the differential report as artifacts. A second job builds and runs natively (apt R) to prove the no-Docker path.

## 4. Stages, recipes, risks, cuts

**Stage 0 — vertical slice (3 ops + recipe, native host, no Docker).** Ops: `dist_q`, `rng_draw`, `fit_lm`. Deliver `run_child` in lmcpp, prelude + emitter, embedded scripts, NA mapping, warnings, `R_ERROR`. Recipe `examples/oracle-lm.lms`: load `examples/data/cars.csv`, call `lmda.model_regress` and `rcli.fit_lm` on the same formula, print max |Δcoef| and exit 1 if above 1e-8; handle `status:"simulated"` like torchcli recipes.
```sh
cd ~/repos/lmctl-ai/lmcpp && cmake -S . -B build -G Ninja -DLMAPI_DIR=../lmapi -DLMCPP_BUILD_TESTS=ON && cmake --build build --parallel 2 && ctest --test-dir build --output-on-failure
cd ~/repos/cppwrapper/rcli && cmake -S . -B build -G Ninja -Dlmcpp_DIR=../../lmctl-ai/lmcpp/stage/lib/cmake/lmcpp && cmake --build build --parallel 2
build/bin/rcli --schema > cli.schema.dev.json
build/bin/rcli dist_q --dist t --p 0.975 --df 7            # {"value":[2.364624251592784],...}
build/bin/rcli rng_draw --dist norm --n 3 --seed 42         # [1.3709584471466685,-0.5646981713960888,0.36312841133733914]
build/bin/rcli fit_lm --file examples/data/cars.csv --formula 'dist ~ speed'
ctest --test-dir build --output-on-failure
PATH="$PWD/build/bin:$HOME/repos/lmctlhq/lmda/target/release:$PATH" lmscript --check examples/oracle-lm.lms && lmscript examples/oracle-lm.lms
```
Gate: the three outputs match the literal values above (host R 4.3.3); recipe prints a Δ below 1e-8; `--schema` works with `Rscript` removed from PATH; `RCLI_READONLY=1` changes nothing (all read).

**Stage 1 — all 21 ops, Docker, shim, oracle tests, 3 recipes.**
```sh
docker build --platform linux/arm64 -f docker/Dockerfile --target test -t rcli:test ~/repos
docker run --rm --network=none rcli:test ctest --test-dir /build --output-on-failure
docker build --platform linux/arm64 -f docker/Dockerfile --target runtime -t rcli:dev ~/repos
export PATH="$PWD/bin:$PATH" TMPDIR="$PWD/.tmp"; mkdir -p "$TMPDIR"
rcli --schema | diff - cli.schema.dev.json
for r in oracle-lm decay-nls spline-physics; do lmscript --check examples/$r.lms && lmscript examples/$r.lms artifacts/$r; done
```
Gate: 21 oracle fixtures exact/tolerance-pass; sandbox corpus rejected; native and Docker outputs byte-identical for the whole fixture set; each call under 2 s in Docker.

**Stage 2 — differential suite, write ops, extra ops, CI.** `tests/differential/` report against lmstats and lmda (gate: report generated, divergences filed); `table_export` write op with atomic install; `cluster_*`, `linalg_decomp`, `fit_coxph`; GitHub Actions. Tier-2 libRmath only after a measured need.

**Recipes.**
1. `oracle-lm.lms` (stage 0): lmda regression and t-test vs `rcli.fit_lm` and `rcli.test_htest` on one CSV; prints an agreement table; nonzero exit on disagreement. Story: R certifies the Rust tool.
2. `decay-nls.lms`: `lmsim.gen_decay` (seeded SSA event stream) → bin counts in lmscript → `rcli.fit_nls` (`N0*exp(-lambda*t)`) and `rcli.fit_glm` (Poisson log-link) → compare recovered λ with the simulation input → `lmcreator.render_chart` of counts, both fitted curves and residuals. Story: simulate → fit → present.
3. `spline-physics.lms`: `physim.run_simple_pendulum` trajectory → `rcli.rng_draw` Gaussian noise (seeded) added in lmscript → `rcli.smooth_spline` with GCV → `rcli.num_root` on the spline's first derivative for the turning time and `rcli.num_integrate` for the swept area → `lmstats.arima` vs `rcli.fit_arima` on the residual series (loglik Δ) → `lmcreator.render_chart`. Story: physics → smoothing → calculus → time-series oracle → chart.

**Risks (ordered).** (1) *Expression sandbox escape* in the three expression ops: mitigated by allowlist on parsed AST, empty-parent environment, Docker `--network=none`, cwd-only mount, scrubbed env, and the negative corpus; if any doubt remains at stage 1 review, ship those three ops Docker-only. (2) *Licensing of the embedded R scripts*: the C++ is clean (R is a separately invoked program; GPL FAQ "mere aggregation"/pipes), but FSF's interpreter-bindings reading could treat R scripts that call GPL packages (mgcv, nlme, survival) as linked; community practice (thousands of MIT CRAN packages) accepts this, and we ship no modified R code. Keep the R scripts clearly separable under `r/` with an explicit note; drop `fit_gam/fit_coxph/fit_lme` to a private repo if counsel objects. (3) *Reproducibility is only as strong as the pinned triple*: Ubuntu 24.04 `r-base-core` will eventually move; pin the `.deb` versions and digest, record `meta`, and treat a version bump as a fixture-regeneration event. (4) Per-call overhead in tight lmscript loops (≈60 ms native, ≈400 ms Docker): mitigated by vectorized ops (`x` arrays, `n` draws); tier 2 if ever measured as blocking. (5) 16 MiB inline cap for large fits: `return` selection now, `table_export` in stage 2. (6) R's locale/encoding quirks in messages: forced `C.UTF-8`, messages passed through as strings only.

**Cuts for v1.** No libR embedding, no libRmath link, no vendored R sources, no arbitrary R code execution, no plotting/graphics, no CRAN installs beyond Debian packages, no write ops, no Arrow/safetensors I/O, no daemon/REPL session, no formula features beyond what `model.frame` parses from CSV columns, no Windows/macOS.

## 5. Naming, placement, licence

Repo and tool name **`rcli`** (binary `rcli`, schema `tool:"rcli"`), public at **`github.com/lmctl-ai/rcli`**, MIT, holder Lattice Systems LLC, commits as Mike Ma <mike.ma@lmctl.com>. Rationale: the shipped source is our C++ plus our R scripts; R is consumed as an installed program over argv and files, exactly like `git` behind `lmgit`; the Docker image redistributes unmodified Debian R packages, which the GPL permits with upstream source references in `NOTICE.md`. Add a `LICENSE-THIRD-PARTY.md` stating R (GPL-2), mgcv/nlme/survival (GPL-2/3), jsonlite (MIT) and that no R source is vendored. The optional tier-2 `libRmath` build, if ever made, goes to a separate repo `rmath-cli` under GPL-2 (or private `lmctlhq`) because it is a GPL combined work. `lmr` is rejected as a name: it collides with existing R packages and reads as "lmctl R", which overstates the scope; `rcli` matches `torchcli`/`physim` naming.
