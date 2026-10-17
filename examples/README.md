# Examples and tutorial conversions

Run from the repository root after building `rcli:runtime` and putting the absolute `bin` directory on `PATH`:

```sh
export PATH="$PWD/bin:$PATH"
lmscript --check examples/tut-lm.lms
lmscript examples/tut-lm.lms
```

Each `tut-*.lms` recipe imports `cli.schema.dev.json`, uses only documented operations, prints a result, and exits nonzero when its stated check fails. Fixtures are original seeded synthetic data; see [data/README.md](data/README.md).

| Recipe | Tutorial topic | Operations |
|---|---|---|
| `tut-lm.lms` | Multiple regression and residual summaries; after R `?lm` | `fit_lm` → `summary_quantile` |
| `tut-glm-binomial.lms` | Logistic dose response; after R `?glm` | `rng_draw` → `fit_glm` |
| `tut-glm-poisson.lms` | Poisson counts; after R `?glm` / Dobson-style examples | `fit_glm` → `dist_quantile` |
| `tut-tests.lms` | Group tests and normality; after R `?t.test` / `?wilcox.test` | `rng_draw` → `test_htest` (t, Wilcoxon, Shapiro) |
| `tut-nls.lms` | Nonlinear saturation curve; after R `?nls` | `rng_draw` → `fit_nls` |
| `tut-timeseries.lms` | Seasonal decomposition, ARIMA, and diagnostics; after R `?arima` | `rng_arima_sim` → `ts_arima`, `ts_acf`, `ts_pacf`, `ts_stl` |
| `tut-cluster.lms` | Partitioning and hierarchical clustering; after scikit-learn clustering guide | `cluster_kmeans` → `cluster_hclust` |
| `tut-smoothing.lms` | Smooth curves; after scipy-lectures statistics chapter | `rng_draw` → `smooth_spline`, `smooth_loess` |
| `tut-optimization.lms` | Scalar optimization, roots, and quadrature; after scipy-lectures optimization chapter | `opt_optimize`, `num_root`, `num_integrate` |
| `tut-distributions.lms` | Distribution tables and simulation; after scipy-lectures statistics chapter | `dist_quantile` → `dist_density`, `rng_draw` → `summary_quantile` |

Expression-taking operations use rcli's constrained expression grammar; none of these recipes passes arbitrary R code. Time-series input is a generated AR(1) noise sequence with an added fixed seasonal cycle rather than an external AirPassengers table.
