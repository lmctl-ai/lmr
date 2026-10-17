# Examples and tutorial conversions

Run from the repository root after building `rcli:runtime` and putting the absolute `bin` directory on `PATH`:

```sh
export PATH="$PWD/bin:$PATH"
lmscript --check examples/tut-lm.lms
lmscript examples/tut-lm.lms
```

Each `tut-*.lms` recipe imports `cli.schema.dev.json`, uses only documented operations, prints its result tables, and
exits 1 when any of its self-checks fails. The checks compare against known truth: the generating parameters of
seeded synthetic data (with tolerances stated in the recipe and derived from the standard errors or the noise level),
or closed-form values for analytic problems. `lmscript --dry-run` previews every recipe (the first call is simulated
and the recipe stops with a notice). Fixtures are original seeded synthetic data; see [data/README.md](data/README.md).
The `oracle-*.lms` recipes compare rcli with lmstats/lmda and need the cross-tool setup in the top-level README.

| Recipe | Tutorial topic | Operations | What is checked |
|---|---|---|---|
| `tut-lm.lms` | Multiple regression; after R `?lm`, `?confint` | `fit_lm` → `dist_quantile`, `test_htest` | coefficient table with t intervals; truth inside every CI; residuals sum to 0, reproduce sigma and R²; Shapiro on residuals; aliasing of a collinear column |
| `tut-glm-binomial.lms` | Logistic regression; after R `?glm` | `rng_draw` → `fit_glm` → `dist_quantile`, `dist_probability` | odds ratios with Wald intervals containing the truth; deviance likelihood-ratio test; fitted P at the true ED50 |
| `tut-glm-poisson.lms` | Count regression; after R `?glm` | `fit_glm` → `dist_quantile`, `dist_probability` | rate ratios with intervals containing the truth; LR test; deviance/df near 1; fitted mean vs Poisson quartiles |
| `tut-tests.lms` | Classical tests; after R `?t.test`, `?wilcox.test`, `?chisq.test`, `?fisher.test`, `?binom.test`, `?cor.test`, `?shapiro.test` | `rng_draw` → `test_htest` (9 tests) | a built-in shift detected with its CI; equal variances kept; balanced table gives p = 1; exact binomial p = 0.05688793; correlation recovered; normality kept/rejected as constructed |
| `tut-nls.lms` | Nonlinear least squares; after R `?nls` | `rng_draw` → `fit_nls` → `dist_quantile` | Michaelis-Menten Vmax, Km inside t intervals; sigma vs noise; fitted curve error; port and default agree |
| `tut-timeseries.lms` | Seasonal decomposition and AR(1) forecasting; after R `?stl`, `?arima` | `rng_arima_sim` → `ts_acf`, `ts_stl`, `ts_arima`, `ts_pacf`, `dist_quantile` | seasonal ACF peak; STL recovers the cycle; ar1 within 3 SE; PACF cut-off; forecast SEs match the AR(1) formula |
| `tut-cluster.lms` | k-means and hierarchical clustering; after R `?kmeans`, `?hclust` | `rng_draw` → `cluster_kmeans`, `cluster_hclust` | purity against the true groups; centre recovery within 3 SE; within-SS vs 2σ²; dendrogram gap; seeded determinism |
| `tut-smoothing.lms` | Smoothing and interpolation; after R `?smooth.spline`, `?loess`, `?spline`, `?approx` | `rng_draw` → `smooth_spline`, `smooth_loess`, `interp_spline`, `interp_approx` | RMSE against the true cubic below half the noise; prediction at new x; fmm spline exact on a cubic; linear interpolation exact |
| `tut-optimization.lms` | Optimisation, roots, quadrature; after R `?optim`, `?optimize`, `?uniroot`, `?integrate` | `opt_minimize`, `opt_optimize`, `num_root`, `num_integrate` | four methods reach Rosenbrock's minimum; active bound solution; Hessian [[4,1],[1,6]]; closed forms for 1-D problems |
| `tut-distributions.lms` | Distribution tables and simulation; after R `?Distributions`, `?quantile` | `dist_quantile`, `dist_probability`, `dist_density`, `num_integrate`, `rng_draw`, `summary_quantile` | reference quantiles to 1e-9; χ²(1) = normal²; p(q(p)) round trips; density integrates to 1; sample mean within 3 SE; seed reproducibility |
| `tut-pacf-order.lms` | AR order identification; after R `?pacf`, `?ar` | `rng_arima_sim` → `ts_acf`, `ts_pacf`, `ts_arima` (orders 1-3) | PACF cut-off after lag 2; AIC selects AR(2); coefficients within 3 SE; ACF(1) vs φ1/(1-φ2) |
| `tut-arima-roundtrip.lms` | ARMA simulation and estimation; after R `?arima.sim`, `?arima` | `rng_arima_sim` → `ts_arima` (ML and CSS-ML) → `ts_acf` | estimates within 3 SE; sigma²; white residuals; forecast SEs follow the psi weights; seed reproducibility |
| `tut-integrate.lms` | Quadrature against closed forms; after R `?integrate` | `num_integrate`, `math_special`, `num_root`, `opt_optimize` | seven integrals inside their error bounds; divergent and non-finite integrands flagged; Γ(3), Γ(½)² = π; π/2 from cos and sin |
| `tut-special-fft.lms` | Special functions and the DFT; after R `?Special`, `?besselJ`, `?fft` | `math_special`, `num_fft` | factorials, half-integer Γ, beta ratios, Pascal's rule, Euler's constant, Bessel Wronskian; impulse/constant/cosine spectra, Parseval, inverse round trip |

Expression-taking operations use rcli's constrained expression grammar; none of these recipes passes arbitrary R code.
Two client-side limits shaped the recipes: `lmscript` types numeric params as numbers, so infinite integration bounds
are written as wide finite ones, and arrays above roughly 100 KB are handed over as a temp file outside the shim's
mount, so samples stay at a few thousand values.
