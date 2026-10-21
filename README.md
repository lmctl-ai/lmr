# rcli

`rcli` is an lmcpp native tool that runs a fixed catalogue of R computations in a fresh `Rscript --vanilla` child process per call. C++ validates the published schema and policy; the worker accepts data only and dispatches through fixed operations. The worker is embedded in the binary and materialized in a private temporary file for each call. Large numeric inputs travel through private one-number-per-line files and are loaded with R's vectorized `scan()`; small options and strings use the in-repo JSON codec. It never evaluates caller-supplied R code.

R is licensed under GPL. This project invokes the unmodified R program as a separate process and exchanges JSON through temporary files; it does not link `libR` or `libRmath`. See [CONTRIBUTING.md](CONTRIBUTING.md) for the dependency boundary.

## Native `_dsl_` host (no R)

`native/` builds `lmr`, a single executable (no R runtime; libstdc++ linked in, only glibc linked dynamically). It runs lmscript `_dsl_` bytecode against a
clean-room C++17 statistics core that uses only the standard library: no R, no Docker, no Python and no third-party
libraries. It covers the same domains as the operations below (distributions, random draws, special functions,
quantiles, FFT, hypothesis tests, lm/glm/nls, splines, loess, ARIMA, STL, k-means, hclust, optimisation, roots and
integration) as a typed object API. Objectives and models are written in the DSL itself.

```sh
cmake -S native -B build-native -G Ninja && cmake --build build-native && (cd build-native && ctest)
lmscript examples/dsl/tut-lm.lms
```

Prebuilt binaries for Linux (x64, arm64; glibc 2.35 or newer), macOS (arm64, x64) and Windows (x64) are attached to
each [release](https://github.com/lmctl-ai/lmr/releases). Put `lmr` on `PATH` (or set `LMR_BIN`) and install
the compiler with `npm i -g @lmctl-ai/lmctl`. Or let lmscript fetch it: `_dsl_ cli="lmr"
install="github:lmctl-ai/lmr@v0.3.0"` downloads the release binary when `lmr` is not on PATH (see
[examples/dsl/install-and-run.lms](examples/dsl/install-and-run.lms)).

The API is in [docs/dsl-api.md](docs/dsl-api.md), and the examples, including multi-step flows, are in
[examples/dsl/](examples/dsl/README.md). [examples/r-vs-dsl/](examples/r-vs-dsl/README.md) solves 44 common R recipes
(regression, GLMs, factors and ANOVA, classical and rank tests, power, distributions, PCA, MDS, clustering, ARIMA,
STL, optimisation, MLE, density estimation, smoothing, FFT) both in R and in lmr on the same data, and compares 3304
numbers side by side. R is used only as an oracle process for the unit tests' expected values.

## Operations

| Operation | Description |
|---|---|
| `r_info` | R version and RNG kinds. |
| `dist_quantile` | Vectorized quantiles for norm, t, chisq, gamma, beta, binom, and pois. |
| `dist_density` | Vectorized density functions for the supported distribution families. |
| `dist_probability` | Vectorized CDFs for the supported distribution families. |
| `rng_draw` | Seeded draws for norm, unif, exp, gamma, beta, binom, pois, t, and chisq. |
| `rng_sample` | Seeded sample from an explicit population using sampled indices. |
| `math_special` | Vectorized gamma/log-gamma, beta/log-beta, digamma, trigamma, choose, and Bessel functions. |
| `summary_quantile` | R sample quantiles, types 1 through 9. |
| `num_fft` | Real and imaginary parts of forward or inverse FFT. |
| `test_htest` | t, Wilcoxon, chi-square, Shapiro, KS, Fisher, binomial, proportion, correlation, and variance tests. |
| `fit_glm` | Gaussian, binomial, Poisson, and Gamma generalized linear models. |
| `smooth_spline` | Smoothing spline fits and optional predictions. |
| `smooth_loess` | LOESS fits and optional predictions. |
| `interp_spline` | Spline interpolation with selected R methods. |
| `interp_approx` | Piecewise linear or constant interpolation. |
| `ts_arima` | ARIMA estimates, optional residuals, and forecasts. |
| `rng_arima_sim` | Seeded ARMA simulation. |
| `ts_acf` | Autocorrelation summary. |
| `ts_pacf` | Partial autocorrelation summary. |
| `ts_stl` | STL seasonal, trend, and remainder decomposition. |
| `cluster_kmeans` | Seeded Hartigan-Wong k-means. |
| `cluster_hclust` | Hierarchical clustering and optional cutree assignments. |
| `fit_lm` | Linear regression on inline columns or a confined CSV file. |
| `opt_minimize` | Fixed-grammar objective optimization with supported optim methods. |
| `opt_optimize` | One-dimensional bounded optimization. |
| `num_root` | One-dimensional root finding with optional interval extension. |
| `num_integrate` | Numerical integration for validated expressions. |
| `fit_nls` | Nonlinear least squares using validated expressions and named starts. |

Smooth-spline and interpolation inputs sort unsorted x values and combine ties by their mean; smooth-spline fitted values follow sorted unique x order. Loess preserves the input row order, including ties, as R does. Vector inputs are capped at 200,000 items unless an operation documents a lower complexity or output cap. `smooth_loess` and Kendall correlation cap at 20,000 and 5,000 observations, respectively; STL caps at 80,000, FFT caps the largest prime factor at 1,000, ACF/PACF lags cap at 1,000, and seasonal ARIMA effective AR/MA orders cap at 350. Density, probability, quantile, special-function, summary, and FFT results preserve R's double values; scalars remain scalars and vector results remain arrays. `rng_sample` samples from explicit `x` values by drawing indices, so a one-element population remains that element. All these operations accept numeric vectors and never evaluate caller-supplied expressions.

All operations are read-only. Output numbers preserve round-trip doubles, and non-finite values use the strings `NA`, `NaN`, `Inf`, and `-Inf`. Fit statistic scalars remain scalars; estimates and other vector fields remain arrays, including one-element coefficient vectors. Results include R version/RNG metadata and a `warnings` array, including on R errors. Numeric inputs accept null or the strings `NA`, `NaN`, `Inf`, and `-Inf`; R determines each operation-specific result. R errors return `COMMAND_FAILED` with `tool_code: R_ERROR`; input validation happens in C++ before the worker starts. Calls time out after 50 seconds (`RCLI_TIMEOUT_S` can lower the limit for development, never raise it above 50 seconds).

## Run

R 4.3.3 is the runtime policy for development and production (`r-base-core=4.3.3-2build2` on Ubuntu 24.04). The image pins the Ubuntu 24.04 multi-architecture manifest digest and lmcpp commit. Oracle tests still calculate version-sensitive expectations with the R executable in the same environment.

For a local native build, install lmcpp and configure rcli with its CMake package prefix, then invoke `build/bin/rcli` directly. The production container and host shim are built with:

```sh
docker build -f docker/Dockerfile --target runtime -t rcli:runtime .
docker build -f docker/Dockerfile --target test -t rcli:test .
```

Use `build/bin/rcli` for a native host binary, or put the absolute path to `bin` on `PATH` to run the isolated container shim. The shim requires a working directory below `$HOME`; it mounts only that directory read-write, disables network access, and starts a fresh container for each call. Build `rcli:runtime` first; `bin/rcli --schema` runs offline from that local image and can regenerate `cli.schema.dev.json`. lmscript recipes import that schema and work when `bin` is on `PATH`; see [examples/README.md](examples/README.md) for the fourteen tutorial conversions and their synthetic fixtures. The shim passes `RCLI_TIMEOUT_S` when set; its value can lower the 50-second call timeout.

Array flags take one item per occurrence (`--p 0.975`; repeat `--p` for additional values). A JSON array passed as one flag value is not accepted.

```sh
build/bin/rcli --schema
bin/rcli --schema
bin/rcli dist_quantile --dist norm --p 0.975
bin/rcli rng_draw --dist norm --n 2 --seed 42
bin/rcli fit_lm --formula 'y ~ x' --file examples/data/oracle-lm.csv
```

The test target runs CTest, including R-computed oracles, expression conformance, no-libR checks, and host-shim contract tests. To check schema compatibility with lmts, generate `cli.schema.dev.json` from the binary and run the pinned lmts `parsePublishedSchema` check (also run in CI). CI (`.github/workflows/ci.yml`) builds and runs the suite with R 4.3.3 on `ubuntu-24.04-arm` and `ubuntu-24.04`: GitHub Actions run 38934174891 (commit 31ce63b, 2026-10-16) passed every gate on both architectures (7/7 CTest groups, no-libR, schema, lmts parse); the earlier amd64 failure was a one-ulp difference between two R code paths for the GLM null deviance, now compared with a derived tolerance (docs/amd64-notes.md). The R oracles are bit-exact for the same R binary and code path; across BLAS builds and architectures only the documented tolerances are promised. See [docs/parity.md](docs/parity.md) for the operation contribution checklist and container behavior.

Recipes: `examples/oracle-lm.lms`, `examples/oracle-tests.lms`, `examples/fit-and-forecast.lms`, `examples/fit-curve.lms`, and `examples/minimize-rosenbrock.lms`. Each imports the schema, checks tool status, and asserts numeric answers. The optional lmstats cross-tool recipes are listed below.

### Cross-tool oracle recipes

`examples/oracle-regression.lms`, `examples/oracle-acf.lms`, and `examples/oracle-arima.lms`
run the same inputs through rcli and lmstats `statsrust`, print both values and absolute differences,
and exit nonzero when a stated tolerance is exceeded. They import rcli's generated schema and the
pinned lmstats schema snapshot in `examples/schemas/`; that snapshot records its source commit.
The OLS and ACF comparisons use `1e-9` absolute tolerance; iterative AR(1) coefficient and variance
comparisons use `1e-3`.

To build `statsrust` without filling the lmstats checkout with build output, use a scratch target
directory and put the resulting executable on an absolute `PATH` entry:

```sh
cd ~/repos/lmctlhq/lmstats
CARGO_TARGET_DIR=/tmp/oracle-target cargo build --locked --bin statsrust
mkdir -p /tmp/oracle-bin
cp /tmp/oracle-target/debug/statsrust /tmp/oracle-bin/statsrust
export PATH="$(realpath ~/repos/lmctlhq/rcli/bin):/tmp/oracle-bin:$PATH"
cd ~/repos/lmctlhq/rcli
lmscript examples/oracle-regression.lms
lmscript examples/oracle-acf.lms
lmscript examples/oracle-arima.lms
```

Remove `/tmp/oracle-target` after building if disk space matters. These cross-tool recipes require
`statsrust` and fail with a setup hint when it is missing; the existing R-only recipes remain
independent of lmstats.

### Expression operations

- `opt_minimize`: `optim` methods Nelder-Mead, BFGS, CG, L-BFGS-B, and seeded SANN; named starts, optional named box bounds for L-BFGS-B, optional numeric data columns, Hessian, and `nonfinite` policy. R defaults are `maxit` 500 (Nelder-Mead), 100 (BFGS/CG/L-BFGS-B), or 10,000 (SANN), and `reltol = sqrt(.Machine$double.eps)`; L-BFGS-B does not receive `reltol`. Returns named parameters, value, counts, convergence/domain status, evaluation count, non-finite count, elapsed time, and warnings.
- `opt_optimize`: one-dimensional `optimize` over a finite interval in `x`, with optional maximization. Defaults to R's `tol = .Machine$double.eps^0.25`; maximum mode returns its maximizing `x`.
- `num_root`: `uniroot` over a finite interval in `x`, with R defaults `tol = .Machine$double.eps^0.25` and `maxiter = 1000`; `extend_int` enables R's bracket extension. Reports the root, residual, iteration count, and precision.
- `num_integrate`: `integrate` over finite bounds or `Inf`/`-Inf`; defaults `rel_tol` and `abs_tol` to `.Machine$double.eps^0.25`; reports value, error estimate, subdivisions, status, and warnings. Scalar integrands recycle to the requested evaluation-vector length, matching R's numerical integration contract.
- `fit_nls`: `nls` default or port algorithm with named starts, optional bounds, and optional fitted/residual vectors. Defaults to `maxiter=50`, `tol=1e-5`, matching `nls.control`; port starts must be inside bounds. Non-convergence returns null estimates and R's message. The response and data columns are ordinary numeric columns; `plinear` and weights are not supported.

### Safety model

Callers send expression text only to the hand-written C++ recursive-descent parser. It enforces the shared `expr-spec.json` allowlist, declared names, and byte/token/node/depth/call/identifier caps, then transports a canonical AST. R never parses or evaluates caller text: its builder revalidates every AST node and arity, inserts function objects from the generated allowlist table, and evaluates in an environment whose parent is `emptyenv()`. `fit_nls` uses the same rebuilt call tree in an empty-parent formula environment, with the safe `list` primitive bound because `nls`'s model-frame construction requires it; expression names cannot access that binding unless it is declared (and `list` is reserved). No caller gradients, callbacks, indexing, assignments, strings, or arbitrary R code are accepted.

Each expression closure has `max_evals` (default 500,000, maximum 2,000,000); it is compiled once per call and its evaluation counter advances once per objective evaluation. R's transient 45-second `setTimeLimit` is the second guard and lmcpp's 50-second process timeout is the last. SANN requires a seed and sets RNGkind explicitly. Invalid syntax/name/arity/cap errors are returned as positional `EXPR_*` validation errors before R starts. Worker-side AST mismatches are internal errors; `EVAL_LIMIT` and non-finite objective errors retain distinct tool codes/domain status. Non-convergence is returned as a domain status.

### Stage 1a supported parameter limits

Complexity caps are validated in C++ before R starts: `test_htest` Kendall uses at most 5,000 pairs; Fisher tables total at most 1,000,000; `ts_arima` requires p + period·P and q + period·Q ≤ 350 and bounds observations × effective-order² at 2e8; `ts_acf`/`ts_pacf` lag_max ≤ 1,000; `smooth_loess` uses at most 20,000 rows; `ts_stl` returns three vectors and uses at most 80,000 rows (keeping serialized output below roughly 8 MiB) and period ≤ 350; `smooth_spline` caps combined input/prediction rows at 100,000; `num_fft` caps its input at 100,000; `cluster_kmeans` requires n·k·nstart·iter_max ≤ 2e8 and k ≤ 500; FFT rejects a largest prime factor over 1,000. Distribution parameter vector recycling and Bessel `expon.scaled` remain unsupported. Wilcoxon is capped at 20,000 samples per side, binomial tests at n=1,000,000, simulated chi-square at 1,000 cells, Fisher at 25 cells and total count 1,000,000. `test_htest` exposes `exact`, `correct`, `byrow`, Fisher `odds_ratio`, and seeded simulated chi-square. `fit_glm` accepts weights and offset vectors.


### Host cap timings

Host measurements use the host binary and R 4.3.3, including JSON transport and process startup. Each row uses a valid input at the operation cap (or the largest cap combination shown); results remained below 11 seconds in these runs.

| Operation / cap fixture | Wall time |
|---|---:|
| `dist_quantile`, `dist_density`, `dist_probability` (200k values) | 3.05 / 2.76 / 2.46 s |
| `rng_draw`, `rng_sample` (200k draws / population and output) | 1.62 / 2.89 s |
| `math_special`, `summary_quantile` (200k values / values plus 200k probs) | 3.00 / 4.60 s |
| `num_fft` (100k values) | 2.10 s |
| `test_htest`: t (200k × 2), Wilcoxon (20k × 2), Fisher (1m total), Kendall (5k), Shapiro (5k), KS (200k), binom (n=1m), prop (200k groups), simulated chisq (1k cells), Pearson/Spearman/variance (200k × 2) | 3.61 / 0.87 / 1.22 / 1.60 / 0.56 / 2.34 / 0.43 / 10.53 / 5.45 / 3.49 / 3.54 / 3.38 s |
| `fit_lm` (200k rows, residuals) / `fit_glm` (200k rows, 1 predictor) | 4.75 / 3.64 s |
| `smooth_spline` (100k combined input and predictions) / `smooth_loess` (20k input and 20k predictions) | 2.33 / 4.05 s |
| `interp_spline` / `interp_approx` (200k x, y and xout) | 6.55 / 6.04 s |
| `ts_arima` (200k observations from seeded AR(1)) / `rng_arima_sim` (200k output) | 3.29 / 1.80 s |
| `ts_acf` / `ts_pacf` (200k observations, lag_max=1000) / `ts_stl` (80k, period=350, robust) | 2.13 / 2.18 / 2.98 s |
| `cluster_kmeans` (200k × 1, k=100, nstart=10, iter_max=1; 200m full budget) / `cluster_hclust` (5k × 1) | 2.93 / 1.33 s |

The JSON document transport is capped at 16 MiB by lmcpp; CSV inputs are capped at 16 MiB and 200,000 rows. Timings are machine-specific measurements, not runtime guarantees.
