# rcli: R as an lmapi native tool (design, Opus)

Decision: `rcli` is an **MIT C++17 lmcpp front end that runs a fixed catalogue of R calls in a separate
`Rscript` worker process** (layer e). It wraps no R C code, links no R library and evaluates no caller-supplied R
code. R 4.5.0 from Debian trixie runs in a pinned arm64 image, called through a torchcli-style shim. Its two jobs:
(1) a **bit-reproducible R oracle** for lmstats/lmda (and later Rust ports), and (2) R algorithms our Rust tools
lack. Tier 2 adds ops from R's recommended packages (survival, mgcv, nlme, MASS, rpart), using the same worker.

## 1. Verified facts (this host, 2026-10-13)

- Host R 4.3.3 (Ubuntu noble, aarch64), `/usr/lib/R/lib/libR.so`, headers in `/usr/share/R/include`
  (`Rmath.h Rinternals.h Rembedded.h Rinterface.h R_ext/`). `r-mathlib` (libRmath) is **not installed**; candidate
  4.3.3. R uses Debian **reference** BLAS/LAPACK 3.12 (`La_library()`), RNG Mersenne-Twister/Inversion/Rejection.
- Recommended packages are already installed: MASS, survival, mgcv, nlme, cluster, rpart, KernSmooth, boot, Matrix,
  nnet, spatial, class. jsonlite is not installed (`r-cran-jsonlite` 1.8.8 on noble, 1.9.1 on trixie).
- Debian trixie-slim offers `r-base-core`/`r-mathlib` **4.5.0-3** (r-base-core Installed-Size 46.6 MB). Bookworm
  offers only 4.2.2.
- Licences: `Rmath.h` is LGPL-2.1+, but nmath sources (`qnorm.c`, `pgamma.c`, `snorm.c`) and stats sources
  (`kmeans.c`, `optim.c`) are **GPL-2+** (Debian: "GPL-2, Portions LGPL"): libRmath is effectively GPL-2+.
- Stats kernels in `stats.so` are mostly `.Call(SEXP)` entry points plus Fortran (`lowesb_`, `stl_`, `hclust_`,
  `rbart_`). They link libR, and **much of each algorithm lives in R code** (glm.fit IRLS, the optim driver in
  arima, nls, smooth.spline setup, loess parameter setup).
- Timings: `Rscript --vanilla` takes ~0.19 s per start; `docker run --rm` of an existing image takes ~0.69 s.
  `write.csv` keeps only 15 significant digits (`0.333333333333333`), while `sprintf("%.17g")` round-trips.

## 2. Layer decision

| Layer | Verdict | Reason |
|---|---|---|
| (a) libRmath | Not primary; maybe a later GPL sidecar | GPL-2+ code. Covers only d/p/q/r and special functions, which the worker can vectorise in one call |
| (b) embed libR | Reject | Same results as (e), but our binary becomes GPL, and longjmp across C++ frames is hazardous. Saves only ~0.2 s/call |
| (c) LAPACK/BLAS | Reject | Not R-specific; reference LAPACK is MIT/BSD-like and Rust has faer/nalgebra |
| (d) stats C/Fortran kernels | Reject | Kernels without R's glue do not reproduce R results, which defeats the oracle; they also need SEXP/libR and are GPL |
| (e) Rscript worker | **Primary** | Exact R semantics, the whole R ecosystem, a clean licence boundary (separate process), and low effort |

R's value is *being R*: bit-exact outputs need R's own glue, defaults, RNG and BLAS, so the process boundary
costs no fidelity.

**Tier 2 (stage 3):** same worker, Debian `r-cran-*` only (pinned with the image, no CRAN installs): `fit_coxph`,
`fit_survreg`, `fit_gam` (mgcv), `fit_lme` (nlme), `fit_rpart`, `fit_glm_nb` (MASS).

## 3. Overlap analysis and v1 catalogue (21 ops)

lmstats (45 ops) is a statsmodels port: time series, GLM/discrete, survival phreg, GARCH, Prophet. lmda (199 ops)
covers tables, tests, posthoc, multivariate, clustering, SPC and power. Neither exposes raw **d/p/q/r distribution
functions, R's RNG streams, sample-quantile types, nls, smoothing splines or R-exact htest outputs**. Ops marked
"oracle" duplicate a lmda/lmstats feature on purpose: they are reference implementations for differential tests,
not features for end users.

| # | op | R call | Role (what it checks or adds) |
|---|---|---|---|
| 1 | `r_info` | R.version, RNGkind, La_library, packageVersion | provenance stamp for every oracle run |
| 2 | `dist_density` | d<family>(x, ..., log) | new; oracle for lmda dist_* / lmstats likelihoods |
| 3 | `dist_cdf` | p<family>(q, ..., lower.tail, log.p) | new; tail accuracy (pnorm(-40,log=T) = -804.6) |
| 4 | `dist_quantile` | q<family>(p, ...) | new; oracle for every critical value in lmda |
| 5 | `dist_random` | set.seed(seed, kind…); r<family>(n, ...) | new; R-identical streams to replay R papers/tests |
| 6 | `rng_sample` | sample(n, size, replace, prob) | oracle for lmda resample_*, permutation tests |
| 7 | `special_fn` | gamma lgamma beta lbeta digamma trigamma psigamma choose lchoose besselI/J/K/Y | new |
| 8 | `sample_quantile` | quantile(x, probs, type = 1..9) | oracle for summary_describe/hdi/ci |
| 9 | `kde` | density(x, bw = nrd0/nrd/ucv/bcv/SJ, kernel, n) | oracle for lmda summary_kde, lmstats kde |
| 10 | `htest` | test ∈ t, wilcox (exact, conf.int), ks (exact), fisher (r×c), shapiro, ansari, binom, poisson, prop, chisq (simulate.p.value+seed), mantelhaen, mcnemar, bartlett, fligner, kruskal, friedman, mood, var, cor | one envelope that checks ~25 lmda test_* ops |
| 11 | `p_adjust` | p.adjust(p, method) | oracle for test_adjust |
| 12 | `fit_lm` | lm.fit through a programmatic formula + summary | oracle for model_regress / lmstats regression |
| 13 | `fit_glm` | glm(family, link, weights, offset) | oracle for lmstats glm/logit/poisson |
| 14 | `fit_nls` | nls with a selfStart catalogue: SSlogis SSmicmen SSasymp SSgompertz SSfpl SSbiexp | **new** (no nls in Rust tools) |
| 15 | `smooth_spline` | smooth.spline(x, y, df / spar / GCV) + predict | **new** |
| 16 | `smooth_loess` | loess(span, degree, family) + predict | oracle for lmda model_loess (R is the reference loess) |
| 17 | `ts_arima` | arima(order, seasonal, method = ML/CSS-ML/CSS) + predict | oracle for lmstats arima (compare with method=ML) |
| 18 | `ts_stl` | stl(s.window, robust) | oracle for lmstats stl |
| 19 | `fft` | fft / mvfft (inverse) | oracle for lmda series_spectrum, physim FFT op |
| 20 | `cluster_kmeans` | kmeans(algorithm = Hartigan-Wong/Lloyd/MacQueen, nstart, seed) | oracle; adds Hartigan-Wong |
| 21 | `cluster_hclust` | hclust(dist(method), method incl. ward.D2, centroid) + cutree(k/h) | oracle for cluster_hclust |

Distribution families (v1): norm, lnorm, t, chisq, f, gamma, beta, exp, weibull, cauchy, logis, unif, binom, pois,
nbinom, geom, hyper, signrank, wilcox, tukey (p/q only), with `ncp` where R supports it. Cut from v1: optim/nlm
with a generic objective, because an objective would mean evaluating caller-supplied code; spec.pgram,
power.*.test and ppr also wait.

## 4. Interface on lmcpp

- **Tool `rcli`, flat snake_case ops** as above, `returns:"object"`, inputs flags+json, no GraphQL. `--schema`,
  `--list` and `--help` never start R.
- **Data in:** reuse lmda's table triple verbatim: `file` (CSV path, confined to cwd) | `records` (inline array of
  objects) | `csv` (inline string), plus column-name params (`y`, `x: array<string>`, `column`). Vector ops also
  accept `x: array<number>` directly. lmscript can then pipe lmda/lmsim records straight into rcli with no
  conversion. Skip Arrow (the R arrow package is not in Debian, and the image cost is large). Skip safetensors in
  v1 (no R reader; add `x_ref` later if torchcli interop is needed).
- **Data out:** inline JSON only in v1, so **all 21 ops are `access:"read"`**: they write nothing and stay usable
  under `--readonly`. Per-row vectors (residuals, fitted values, draws, clusters, STL components) are columnar
  arrays, capped at 200k numbers (well under lmcpp's 16 MiB result_file cap). Above the cap the op fails with
  `COMMAND_FAILED/OUTPUT_TOO_LARGE`. Stage 3 may add a `write` op `dist_random_csv` with an exclusive `out`, using
  torchcli's no-replace atomic writer.
- **Result envelope** (always present): `{ r: {version, rng_kind}, call: "<canonical R call>", warnings: [...],
  ... }`. Tabular results use a columnar `table: {columns: [..], data: {col: [..]}}`; coefficient tables use
  `coefficients: {term, estimate, std_error, statistic, p_value}`. htest is uniform: `statistic, parameter,
  p_value, conf_int, estimate, null_value, alternative, method`.
- **Numbers:** the worker encodes every double with `sprintf("%.17g")` through its own encoder (jsonlite is used
  only for parsing input), and nlohmann re-emits the shortest round-trip form, so values cross the boundary
  bit-exactly. **NA/NaN/±Inf are encoded as JSON strings `"NA"`, `"NaN"`, `"Inf"`, `"-Inf"`**. These are lossless;
  NA≠NaN matters to an oracle. Inputs accept `null` or `"NA"` as NA. CSV input: empty or `NA` cells are NA.
- **Seeds/RNG:** every stochastic op (`dist_random`, `rng_sample`, `kmeans`, `chisq` simulate) requires `seed`
  (an integer). Optional `rng_kind` / `normal_kind` / `sample_kind` default to R's defaults. The worker calls
  `set.seed(seed, kind, normal.kind, sample.kind)` first, and the result echoes all of them plus `r.version`.
  Reproducibility is guaranteed only within a pinned image. R's RNG has been stable since 3.6, but nmath fixes do
  change tails, so goldens are tied to the image digest.
- **Errors:** front-end validation (types, enums, family/param arity, column existence, n ≤ caps) fails with exit 2
  before R starts. The worker wraps each op in `withCallingHandlers` (collecting warnings into `warnings[]`) inside
  `tryCatch`. An R error is `COMMAND_FAILED` exit 1, `tool_code:"R_ERROR"`, `data:{r_message, r_call}`.
  Non-convergence (nls/arima/glm `converged=FALSE`) is `tool_code:"NOT_CONVERGED"` for nls (where R itself
  errors). For glm/arima it is a success with `converged:false`. Domain issues like qnorm(2) give `"NaN"` plus the
  R warning "NaNs produced". The global `--strict-warnings` turns any warning into `R_WARNING` exit 1. The worker
  timeout (default 50 s, below lmscript's 60 s) kills the process group and returns TIMEOUT 7. A crash or empty
  response is `INTERNAL_ERROR`/`TOOL_CRASH`.
- **No code injection:** no formula, expression or function strings are accepted. Column names must match
  `^[A-Za-z.][A-Za-z0-9._]{0,63}$`, and formulas are built with `reformulate(backquoted names)`. The worker
  dispatches through a closed `switch(op)` in `R/ops.R`; the request is plain data from `jsonlite::fromJSON`. The
  worker runs as `Rscript --vanilla --default-packages=stats,utils` with a cleared env (`LANG=C.UTF-8`,
  `R_LIBS_USER=` empty, `OMP_NUM_THREADS=1`), `--network=none`, and a read-only rootfs.
- **Worker protocol:** the front end writes `req.json` to `Context::temp_file_path()`, then execs (argv array, no
  shell) `Rscript … /opt/rcli/R/main.R <op> req.json resp.json` and reads the response with
  `read_result_file(json)`. Each call runs one R process, so batching comes from vector params, not from loops.

**lmcpp gaps to fill (upstream, MIT):** (1) a `subprocess` helper: argv exec, cleared env, timeout +
process-group kill, bounded stdout/stderr, exit/signal mapping (deferred in lmcpp parity "generic subprocess
helpers"). (2) A `table_input` helper implementing lmda's file/records/csv triple with typed column extraction and
NA. (3) `nonfinite` JSON encode/decode helpers. (4) Lifting torchcli's cwd-confined path canonicaliser and
no-replace atomic writer into lmcpp. (5) A reusable docker shim template (`scripts/docker-shim.sh`) parameterised
by image/tool/env-prefix. (6) Optional: `Type::Integer` items in arrays (`seed`, `order`). GraphQL is not needed.

## 5. Build, packaging, testing

- **Image:** `debian:trixie-slim@sha256:<arm64 digest>`. Pin apt through a `snapshot.debian.org` date (fixing
  torchcli's unpinned-apt gap: oracle values depend on the R build). Builder: `g++ cmake ninja-build
  nlohmann-json3-dev` plus lmcpp from the repos-root context (same allowlisted `.dockerignore` pattern as
  torchcli). Runtime: `r-base-core r-cran-jsonlite --no-install-recommends` (stage 3 adds `r-cran-{survival,mgcv,
  nlme,mass,rpart}`), the rcli binary and `/opt/rcli/R/*.R`. Keep reference BLAS; never install OpenBLAS,
  because threaded BLAS breaks bit-reproducibility.
- **Estimates (to be measured at stage 0):** runtime image 250–350 MB, cold build 4–8 min (apt dominates), C++
  rebuild under 1 min, per call ~0.9–1.1 s (0.7 s docker + 0.2 s R).
- **Shim:** copy torchcli's `bin/torchcli` as `bin/rcli`: same cwd-under-HOME rule, `--mount` of cwd, uid/gid,
  `--network=none --cpus=1 --memory=2g --read-only --tmpfs /tmp`, and forwarding of `RCLI_READONLY/RCLI_DRYRUN`.
  **Dev mode:** `RCLI_RSCRIPT=/usr/bin/Rscript build/rcli …` runs natively on the host (R 4.3.3) for protocol and
  unit work. Golden tests run only in the image.
- **Tests:** (1) lmcpp protocol vectors plus CLI regressions (validation, readonly, timeout, crash, injection
  attempts such as a column named ``a`); system("x")``, oversized output). (2) **Oracle goldens per op:**
  `tests/oracle/<op>.R` calls the R function directly (independently of `ops.R`) on 3–6 fixtures, including edge
  cases (extreme tails, ties, NA, n=1). It writes expected JSON in the same image, and CTest asserts **bit
  equality** of rcli output. This proves the plumbing, not R. (3) A round-trip test of 10⁵ random doubles plus all
  nonfinite values through req→R→resp. (4) **Consumer tolerance policy** (published in `docs/oracle.md` for
  lmstats/lmda differential tests): d/p/q and special functions rel 1e-12 (abs floor 1e-300; log scale abs 1e-10);
  closed-form tests/quantiles/fft rel 1e-10; IRLS/QR fits rel 1e-8; iterative optimisers (arima, nls,
  smooth.spline GCV, loess) rel 1e-5 on estimates and 1e-4 on SEs; cluster labels must give ARI = 1 after
  relabelling; RNG draws exact or the test fails.
- **CI:** `ubuntu-24.04-arm` Actions: build test image, ctest, `lmscript --check examples/*.lms`; nightly golden
  regeneration diff (any diff = reviewed image bump).

## 6. Licensing, naming, placement

- Linking libRmath or libR, or copying, translating or porting R C/Fortran/R source (including nmath → Rust),
  makes the result GPL-2+. **Running R as a separate program** that exchanges JSON files over argv is
  mere aggregation/arm's-length communication. CRAN itself hosts many MIT-licensed R packages that call GPL R
  functions, so our `R/*.R` op scripts and the C++ front end can be **MIT ("Lattice Systems LLC")**.
- CONTRIBUTING rules: never include `Rmath.h` or link `-lR/-lRmath`, never paste R source; Rust ports use rcli only as a black-box oracle.
- **The image contains unmodified Debian GPL binaries.** Building locally is fine. Publishing the image requires a
  source offer, which is satisfied by pointing at the pinned snapshot.debian.org source packages in a NOTICE.
  Do not push images in v1.
- **Name/placement:** repo and tool `rcli` (parallels torchcli). Start **private `lmctlhq/rcli`** (beside lmstats/lmda)
  with MIT LICENSE from commit 1; move to **public `lmctl-ai/rcli`** at the stage-3 gate with no relicensing.

## 7. Staged plan with acceptance

**Stage 0 — thin vertical slice (3 ops + recipe):** `r_info`, `dist_quantile`, `dist_random`, the shim, the
worker protocol, the nonfinite encoding and lmcpp subprocess helper (1).
```sh
cd ~/repos/lmctlhq/rcli
docker build --platform linux/arm64 -f docker/Dockerfile --target test -t rcli:test ~/repos
docker run --rm --network=none rcli:test ctest --test-dir /build --output-on-failure
docker build --platform linux/arm64 -f docker/Dockerfile --target runtime -t rcli:dev ~/repos
export PATH="$PWD/bin:$PATH" TMPDIR="$PWD/.tmp"; mkdir -p .tmp
rcli --schema > cli.schema.dev.json
rcli dist_quantile --family norm --p 0.975 --p 2      # → [1.9599639845400534,"NaN"] + warning "NaNs produced"
rcli dist_random --family norm --n 2 --seed 42         # → [1.3709584471466685,-0.5646981713960888] (host-verified)
rcli --readonly r_info                                  # works (read)
lmscript --check examples/stage0-ci.lms && lmscript examples/stage0-ci.lms
```
Gate: the recipe draws 500 values with `rcli dist_random(seed=7)`, then runs `lmda summary_describe(records)` and
`rcli dist_quantile(t, df=499)`. It prints the 95% CI and checks it against lmda `summary_ci` to within 1e-10. The
schema works offline; the bit-equality goldens pass; a call takes under 1.5 s.

- **Stage 1 — distributions + oracle core (ops 2–11).** Gate: goldens for every family×function incl. extreme
  tails; `lmscript examples/oracle-tests.lms` passes.
- **Stage 2 — models (ops 12–21) + `docs/oracle.md`.** Gate: goldens; `oracle-arima.lms`, `sim-fit-chart.lms` pass;
  timing table recorded.
- **Stage 3 — tier 2 + publish:** recommended-package ops, `dist_random_csv`, lmcpp helpers upstreamed, CI, NOTICE,
  transfer to lmctl-ai.

## 8. Recipes (stories across tools)
1. `examples/oracle-tests.lms` — **R checks lmda.** For a fixed CSV, call lmda `test_ttest`, `test_mannwhitney`,
   `test_ks`, `test_fisher`, `test_adjust` and the matching `rcli htest`/`p_adjust`. Print a table of
   |Δstatistic| and |Δp|, and fail on any value over the tolerance policy. A nightly run of this recipe turns R
   into lmda's regression oracle.
2. `examples/oracle-arima.lms` — **R checks lmstats.** Run `lmsim gen_timeseries(seed)` to make an ARMA series,
   then `lmstats arima(order=[1,0,1], method=statespace)` and `rcli ts_arima(method="ML")`. Compare coefficients
   (rel 1e-5), log-likelihood and 12-step forecasts, then `lmcreator render_chart` with both forecast bands overlaid.
3. `examples/sim-fit-chart.lms` — **simulate → fit → present.** Run `physim run_michaelis_menten_enzyme_kinetics`
   to get substrate/product trajectories and add seeded noise with `rcli dist_random`. Then fit `rcli fit_nls`
   (SSmicmen) for rate vs. substrate and `smooth_spline` on the product curve. Make
   `lmcreator render_chart` (data, fit, residuals) and `export_document` a DOCX with the recovered Vmax/Km
   alongside the physim ground-truth parameters.

## 9. Risks and cuts
- **Per-call latency (~1 s)** makes loop-heavy recipes slow. Mitigate with vectorised params and measure first.
  A `rcli serve` session or native host mode come only after measurement. Never embed libR to save time.
- **Oracle drift/false confidence:** R 4.5.0 (image) ≠ 4.3.3 (host) ≠ what users run, and nmath/arima fixes
  change digits. Pin by digest + snapshot, stamp `r_info` in every result, and treat golden diffs as reviewed bumps.
- **Licence leakage:** someone links Rmath "for speed" or ports nmath into lmstats. Enforce with a CI grep for
  `Rmath.h|-lR` and the CONTRIBUTING rules.
- Smaller: non-bug semantic mismatches (arima ML vs statespace, ward vs ward.D2) — recipes name the matching R method.
- **Cut if needed:** tier 2, `smooth_loess`, `kde`, `fft`, Bessel fns, nightly diff, publication. **Never cut:**
  goldens, seed echo, nonfinite encoding, the process boundary.
