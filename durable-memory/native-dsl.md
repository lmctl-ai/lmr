# Native rcli: `_dsl_` host on a C++17 core (no R, no Docker)

Status (2026-10-20): implemented in `native/`. The owner asked for:
1. a build with no Docker and no R runtime, using pure C/C++ libraries only, with no GUI or display code;
2. a defined API to expose;
3. a `_dsl_` version of every example;
4. multi-step examples.

This supersedes the Rust `lmnum` route in `native-plan.md` for rcli. That plan's op mapping, clean-room rules and
tolerance tiers still apply.

## Shape

- **`native/include/rcli/*.hpp`: the contracts.** These are the numeric core API in namespace `rcli`: special, dist,
  rng, summary, fft, interp, smooth, optim, htest, models, ts, cluster. Each module has
  `native/src/num/<module>.cpp` and a unit test `native/tests/test_<module>.cpp`. The tests use `native/tests/check.hpp`,
  and their expectations come from Rscript run as a black-box oracle, with the R expression in a comment.
- **`native/src/host/`: the DSL host `rcli-dsl`.** It follows the lmcad/lmhtml pattern.
  - `host.hpp` holds the arenas, the `Binder` (type-checks arguments, maps `rcli::Error` to `dslvm::Error`) and the
    `Driver`. `bind_<domain>.cpp` registers the API. `main.cpp` holds the CLI and generates the schema.
- **`native/cpp/dslvm.hpp`: the vendored SDK.** It is pinned by `native/sdk.lock` and refreshed with
  `native/scripts/sync-sdk.sh`. The source is `../poc/cpp-dsl-sdk` at commit ef86184, bytecode v5.
- **`native/schema/rcli.dsl.json`: GENERATED from the registry.** Never edit it by hand. Change a `bind_*.cpp`, then run
  `native/scripts/update-schema.sh`. The ctest `schema_matches_registry` diffs the file, the embedded bytes, and
  `lmscript --dsl-signatures` against `rcli-dsl --dsl-signatures`.
- **Build**: `cmake -S native -B build-native -G Ninja && cmake --build build-native`. The binary is one executable at
  `build-native/bin/rcli-dsl`: libstdc++ and libgcc are linked statically, glibc dynamically (the owner's choice,
  2026-10-20; `ldd` shows only libc, libm and the loader). It needs glibc >= the build host's (2.38 on the Ubuntu 24.04
  dev box: `__isoc23_strtoul`, `fmod`); build on an older distro for wider reach. `-DRCLI_DSL_FULL_STATIC=ON` restores
  a fully static binary. With dynamic libm, last-ulp results of exp/log/pow etc. can follow the target's glibc. The flag `-ffp-contract=off` is required (aarch64 GCC
  otherwise fuses multiply-adds, which changes roundoff and some optimiser counts).
- **Examples**: `examples/dsl/*.lms`, sharing `examples/dsl/lib/rclidsl.lms` (save, exec and check). Run them with
  `native/scripts/run-examples.sh build-native/bin/rcli-dsl` from the repo root.

## Decisions (and why)

- **Reverse communication for callbacks.** The DSL cannot pass functions, so `minimize`, `optimize`, `uniroot`,
  `integrate` and `nls` run the native algorithm on a worker thread. The script answers `next()` / `x()` / `tell()`.
  - Only one side runs at a time, so results are deterministic.
  - `Driver` is declared last in every solver struct so it is destroyed (cancelled and joined) first.
  - Native algorithms must be exception-safe and must never swallow exceptions with `catch (...)`: cancellation is a
    non-`std::exception` type thrown through them.
- **No expression-string grammar.** The old `expr-spec.json` evaluator was R-side. Objectives are now DSL code, with
  `exp`, `log` and the special functions added as host functions because the DSL only has
  sin/cos/sqrt/abs/floor/pow/min/max.
- **Output protocol.** stdout carries exactly one JSON line (`values` from `out.set`, plus checks and failures). Report
  text goes to stderr. A failed `out.check` gives exit 4. lmscript orchestration reads `values`.
- **Integrate does not throw on a NaN integrand.** It returns `ier` 6 with a message, because the DSL has no try/catch.
- **The RNG is not R's stream.** It is xoshiro256** with exact variate methods, which the owner accepted. The examples'
  checks must be statistical (truth within k SE), never "R's first draw". Changing a seed is fine; weakening a check is
  not.
- **Index convention.** Indices are 0-based. Labels keep R's values: k-means/cutree labels 1..k, hclust merge signs, and
  1-based `rank()` and `order()` of hclust leaves.

## Clean-room rule (licensing)

R is GPL. The rules:

- Never read R/CRAN sources or print R function bodies. Do not reproduce nmath or other R code from memory.
- Rscript may be run only as a black-box oracle, and R help pages may be read for documented behaviour.
- Algorithms come from papers (cited in comments).
- Public-domain QUADPACK and netlib STL may be ported. No Numerical Recipes code.
- `special.hpp` exposes `stirlerr`, `bd0`, `dpois_raw` and `dbinom_raw`. These are Loader (2000) saddle-point terms,
  named as in that paper.

## Measured agreement with R 4.3.3 (from the module reports)

| Module | Agreement |
|---|---|
| cluster | exact merge, order and cutree including ties; heights 1e-12; kmeans with given centres identical (cluster, iter, ifault) |
| optim | integrate ~1e-15 with the same ier; optimize bit-identical; Nelder-Mead bit-identical including counts; BFGS same par (1e-12), counts ±1; L-BFGS-B same optimum, counts differ; uniroot roots within tol, some iteration counts differ (Brent's published zeroin vs R's netlib variant) |
| models | lm 1e-14; glm 4e-15 (offset too); nls estimates 1e-9, se 3e-8; nls port is not nl2sol (counts differ) |
| htest | 1e-15 for t/var/binom/chisq/prop/pearson; exact wilcox identical; shapiro W 2e-14, p 2e-12; ks 1e-12; fisher p 1.5e-11. Known: prop.test null_value with k > 1 and given p is not representable |
| ts | acf/pacf 1e-16; ARIMA (ML, CSS, CSS-ML; mean, differenced, seasonal) coef <= 1e-11 and se <= 6e-8 for most models, airline models 1e-6..5e-5; forecasts <= 1e-6; arima_sim identical to arima.sim given the same innovations; STL 1e-13. Two findings by black box: R's ARIMA se use the Hessian in the transformed (Jones) parameters mapped back by a 1e-3 forward-difference Jacobian; the intercept is optimised scaled by 10 sd(x)/sqrt(n). Known: R's robust STL on even-length series uses an off-median robustness scale at some iterations (we use the true median; up to 5e-2 on remainders); near-non-stationary seasonal fits differ (R's Gardner start vs our exact covariance) |
| summary, fft, interp, smooth | summary 1e-15 (all 9 quantile types); fft 4e-15 normalised error; splines/approx 1e-13; loess direct 1e-15 (R's default interpolate surface differs by up to 0.28 on noisy n = 50); smooth.spline same df within 1e-5 of the data range (nx <= 49), GCV df within 0.02. R's monoH.FC third derivative is wrong (6dh instead of 6d); ours is correct |
| dist, special, rng | d/p <= 4.4e-13 over ~15,400 grid points; q <= 4e-13 (qt 2e-11, qf 4e-10 at extremes); special <= 1e-13; Bessel J 1.1e-11, Y 1.3e-12, I/K 1e-13. 58 grid points deviate where R is the less accurate side (underflow to 0, e.g. pnorm(-38) = 2.885e-316 natively; extreme-tail qt/qf; lgamma near 1 and 2, checked at 50 digits). Binomial variates use BTRS. RNG integer streams are bit-identical across platforms; libm-dependent variates may differ in the last ulp |

## Tests (16 ctest targets, all passing on 2026-10-20)

- **Unit tests per module**: `native/tests/test_<module>.cpp`, about 19,000 R-oracle checks.
- **`test_edge.cpp`**: about 200 checks.
  - Hostile callbacks (NaN, ±Inf, constant, huge or alternating values; a throw at evaluation k) for every
    optim/optimize/uniroot/integrate/nls variant: each must terminate within an evaluation budget or raise
    `rcli::Error`.
  - Degenerate data (empty, n = 1, constant, all ties, separation, more columns than rows, extreme parameters).
  - Exact R answers and messages for degenerate cases.
- **`host_test.py`** (ctest `host_bindings`): about 160 cases that compile small `_dsl_` chunks and run `rcli-dsl`. It
  asserts exit codes, the JSON result line, error messages and DSL line/column. It covers every binding family, CSV
  edge cases, `fmt`, keep/kept across chunks, the solver protocol and its misuse errors, nested and abandoned solvers,
  and the CLI. Filter with `-k substring`.
- **Sanitizers**: `cmake -S native -B build-asan -DRCLI_SANITIZE=address,undefined`. The unit tests, `test_edge`,
  `host_test.py` and all examples run clean under ASan and UBSan.
- **`tests/oracle/fuzz_r.py BUILD_DIR [N] [SEED]`**: optional, needs Rscript, and is not part of ctest.
  - It generates random d/p/q cases for all 16 families (both tails, log scale), quantile types, cor with ties and
    t/wilcox/cor tests, then compares one Rscript process with the native driver `fuzz_native.cpp`.
  - Expect about 6-15 mismatches per 2850 cases. All are extreme-tail quantiles: continuous qt/qf, where R does not
    invert its own CDF; and discrete quantiles with p within about 1e-14 of 1, where both sides are at double
    resolution. Investigate anything else.
- **Bugs found and fixed this round**:
  - `json_number`/`Vec.str` printed 10 as `1e+01`.
  - `integrate` did not default `abs_tol` to `rel_tol` (R does).
  - `binom.test` used a different order and wording for its argument errors.
  - Upper-tail discrete quantiles with p near 1 returned the support minimum: the fuzz p·(1 + 8ε) reached 1, so
    `qbinom(1 - 1.3e-15, 533, 0.885, lower.tail = FALSE)` gave 0 instead of 407. Now the fuzz applies to 1 − p.
- **Working note**: background jobs and subagents in this environment are killed when the Claude session restarts. A
  hardening round run as three background agents was lost twice. Run long builds and tests in the foreground in steps
  under 10 minutes.

## Not done / next

- The lmapi JSON ops (`rcli fit_lm {...}`) still run R in Docker. Re-pointing them at `rcli_num` would give a static
  `rcli` with the old interface. Everything they need exists in the core, except the expression-string evaluator.
- The `oracle-*` recipes need `statsrust` and were not ported.
- No CI leg builds `native/` yet.
