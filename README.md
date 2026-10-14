# rcli

`rcli` is an lmcpp native tool that runs a fixed catalogue of R computations in a fresh `Rscript --vanilla` child process per call. C++ validates the published schema and policy; the worker accepts data only and dispatches through fixed operations. The worker is embedded in the binary and materialized in a private temporary file for each call. Large numeric inputs travel through private one-number-per-line files and are loaded with R's vectorized `scan()`; small options and strings use the in-repo JSON codec. It never evaluates caller-supplied R code.

R is licensed under GPL. This project invokes the unmodified R program as a separate process and exchanges JSON through temporary files; it does not link `libR` or `libRmath`. See [CONTRIBUTING.md](CONTRIBUTING.md) for the dependency boundary.

## Operations

- `r_info`: R version and RNG kinds.
- `dist_quantile`: vectorized quantiles for `norm`, `t`, `chisq`, `gamma`, `beta`, `binom`, and `pois`, with `lower_tail` and `log_p` controls.
- `dist_density`: vectorized `dnorm`/`dt`/`dchisq`/`dgamma`/`dbeta`/`dbinom`/`dpois`, with `log`.
- `dist_probability`: vectorized matching CDFs, with `lower_tail` and `log_p`.
- `rng_draw`: seeded `norm`, `unif`, `exp`, `gamma`, `beta`, `binom`, `pois`, `t`, and `chisq` draws. Each stochastic call sets R's RNG kind to Mersenne-Twister / Inversion / Rejection before `set.seed`.
- `rng_sample`: seeded sample from an explicit numeric population via sampled indices (`sample.int`), optionally with replacement and probability weights.
- `math_special`: vectorized `gamma`, `lgamma`, `beta`, `lbeta`, `digamma`, `trigamma`, `choose`, `lchoose`, and Bessel J/I/K/Y functions.
- `summary_quantile`: R's nine sample-quantile definitions for supplied probabilities.
- `num_fft`: real-valued FFT input with separate real and imaginary arrays; inverse mode matches R’s unnormalized `fft(inverse=TRUE)`; `normalize=true` divides by input length.
- `test_htest`: `t`, `wilcox`, `chisq`, `shapiro`, `ks`, `fisher`, `binom`, `prop`, `cor`, and `var` tests. It uses one enum operation because R htests share the statistic/parameter/p-value/interval/estimate/method shape; named statistics and estimates include companion name fields. Tables default to R column-major filling (`byrow=false`); Fisher supports explicit r×c dimensions. Simulated chi-square requires a seed.
- `fit_glm`: base R Gaussian, binomial, Poisson, and Gamma GLMs with supported links, weights/offsets, coefficients, standard errors, t or z inference, dispersion, aliases, deviance, AIC, iterations, and convergence status.
- `smooth_spline` and `smooth_loess`: base R smoothers with fitted values and optional predictions at `new_x`.
- `interp_spline` and `interp_approx`: R spline and piecewise approximation methods evaluated at `xout`.
- `ts_arima`: base R ARIMA estimation (default method CSS-ML, matching `stats::arima()`), optional seasonal order and mean term, coefficient inference, optional residuals, convergence flag, forecasts and forecast standard errors.
- `rng_arima_sim`: seeded `arima.sim` for ARMA models (d must be zero; AR/MA coefficient counts match p/q).
- `ts_acf` and `ts_pacf`: autocorrelation and partial autocorrelation summaries, with vector output under `value`.
- `ts_stl`: STL seasonal/trend/remainder decomposition.
- `cluster_kmeans`: seeded Hartigan-Wong k-means with centers, memberships, within-cluster sums, and sizes.
- `cluster_hclust`: distance plus hierarchical clustering with merge, height, order, and optional `cutree` assignments. Hclust is capped at 5,000 rows due to quadratic memory use.
- `fit_lm`: formulas using named data columns, plus intercept-only `y ~ 1` and intercept removal with `-1` or `+0`, with inline JSON columns or a cwd-confined CSV path. Dot formulas and expressions are rejected. CSV inputs are capped at 16 MiB and 200,000 rows.

Smooth-spline and interpolation inputs sort unsorted x values and combine ties by their mean; smooth-spline fitted values follow sorted unique x order. Loess preserves the input row order, including ties, as R does. Vector inputs are capped at 200,000 items unless an operation documents a lower complexity or output cap. `smooth_loess` and Kendall correlation cap at 20,000 and 5,000 observations, respectively; STL caps at 80,000, FFT caps the largest prime factor at 1,000, ACF/PACF lags cap at 1,000, and seasonal ARIMA effective AR/MA orders cap at 350. Density, probability, quantile, special-function, summary, and FFT results preserve R's double values; scalars remain scalars and vector results remain arrays. `rng_sample` samples from explicit `x` values by drawing indices, so a one-element population remains that element. All these operations accept numeric vectors and never evaluate caller-supplied expressions.

All operations are read-only. Output numbers preserve round-trip doubles, and non-finite values use the strings `NA`, `NaN`, `Inf`, and `-Inf`. Fit statistic scalars remain scalars; estimates and other vector fields remain arrays, including one-element coefficient vectors. Results include R version/RNG metadata and a `warnings` array, including on R errors. Numeric inputs accept null or the strings `NA`, `NaN`, `Inf`, and `-Inf`; R determines each operation-specific result. R errors return `COMMAND_FAILED` with `tool_code: R_ERROR`; input validation happens in C++ before the worker starts. Calls time out after 50 seconds (`RCLI_TIMEOUT_S` can lower the limit for development, never raise it above 50 seconds).

## Build and run

The host has no CMake/Ninja/nlohmann packages, so the development image extends `lmcpp-dev` with `r-base-core` and builds against the installed lmcpp package:

```sh
docker build -f Dockerfile.dev -t rcli-dev .
# In lmcpp, build and install to a local prefix:
cmake -S . -B build-rcli -G Ninja -DCMAKE_INSTALL_PREFIX="$PWD/stage-rcli" -DLMCPP_BUILD_TESTS=ON
cmake --build build-rcli --parallel 2 && ctest --test-dir build-rcli --output-on-failure
cmake --install build-rcli
# Here, configure rcli with -DCMAKE_PREFIX_PATH=/path/to/lmcpp/stage-rcli
```

The Debian bookworm apt candidate in `rcli-dev` is R 4.2.2; the host runtime is R 4.3.3. Version-sensitive test expectations are computed by the `Rscript` in the same environment. R version is included in every result. No production Docker image or shim is part of this slice.

```sh
build/bin/rcli --schema
build/bin/rcli dist_quantile --dist norm --p 0.975
build/bin/rcli rng_draw --dist norm --n 2 --seed 42
build/bin/rcli fit_lm --formula 'y ~ x' --file examples/data/oracle-lm.csv
```

Build and test in the dev image with `cmake --build build --parallel 2 && ctest --test-dir build --output-on-failure`. R fixture scripts calculate expected values directly with the same environment's R; the lmscript recipes import the schema and assert against R-computed fixtures. The smoke suite covers R-oracle values, schema descriptions, schema-safe formulas, warnings, non-finite values, timeout/signal cleanup, TMPDIR forwarding, complexity-cap rejection, result shapes, and large vector input. Array flags take one item per occurrence (`--p 0.975`; repeat `--p` for more); JSON-array text as one flag value is not accepted. `scripts/check-no-libR.sh` enforces the no-link boundary and checks its detector against a planted `Rmath.h` include.

Alternate RNG kinds and normal-generation kinds are not configurable; calls use the explicit defaults listed above. The `examples/oracle-lm.lms`, `examples/oracle-tests.lms`, and `examples/fit-and-forecast.lms` recipes import the published schema, handle call status, and check R-computed numeric fixtures. Oracle-test fixtures use absolute tolerance 1e-8; the forecast recipe uses the same tolerance.

Stage 1a delivers the seven-family distribution functions, expanded seeded RNG catalogue, sampling, math-special functions, sample quantiles, FFT, common hypothesis tests, GLM, smoothing/interpolation, ARIMA/time-series summaries, STL, and clustering. It leaves the design's larger distribution family list for later. The fixed worker uses a small in-repo JSON codec, so host R does not need `jsonlite` installed. The dev image is based on Debian bookworm because that is the available `lmcpp-dev` image; its R is 4.2.2 rather than the host's 4.3.3, and every result records that version.

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
