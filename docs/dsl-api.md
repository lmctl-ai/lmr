# rcli `_dsl_` API (native, no R)

`lmr` is a single executable (no R runtime; it needs only the system glibc) that runs lmscript `_dsl_` bytecode against a native C++17 statistics core.
It needs no R, Python, Docker or third-party libraries. It replaces the R-backed operations of `rcli` with an
object API. A script describes a whole multi-step analysis in one typed chunk; lmscript compiles that chunk once, and
the host runs it in one process.

```
lmscript examples/dsl/tut-lm.lms             # compiles the chunk to build/dsl/tut-lm.dslb, runs lmr, checks
lmr --dsl build/dsl/a.dslb build/dsl/b.dslb  # several chunks, one process, shared host state
lmr --dsl-schema                             # every signature, grouped by type
lmr --dsl-schema --json > lmr.dsl.json       # the schema the compiler needs (byte-identical to native/schema/)
```

The machine-checked list of all signatures is `native/schema/lmr.dsl.json`. The test `schema_matches_registry` checks
it against the host registry and against `lmscript --dsl-signatures`. This page explains the design and each type;
`lmr --dsl-schema` prints the full list (`--schema` is an older alias).

## Why a DSL rather than one process per call

The schema-module form (`import "cli.schema.dev.json" as rcli; rcli.fit_lm({...})`) costs a process start, a JSON
encode/decode of every array, and (formerly) an R start per call. A tutorial with ten tests made ten round trips. A
`_dsl_` chunk is compiled and type-checked by lmscript before anything runs, then executed in one process. Vectors stay
native, so loops over observations (simulating data, residual checks, bootstraps) run at bytecode speed. Objectives,
integrands and model curves are written in the chunk itself (see "Reverse communication"), so no expression-string
grammar is needed.

## Process contract

| Item | Contract |
|---|---|
| stdout | exactly one JSON line: `{"ok": true, "result": {"values": {...}, "checks": N, "failures": 0}}` or `{"ok": false, "phase": "load"/"run"/"app", "code": ..., "message": ..., "line": L, "column": C}` |
| stderr | the chunk's `out.print` text, then `CHECK FAILED: ...` lines |
| exit | 0 ok; 2 usage; 3 load or run error (bad bytecode, schema mismatch, host error such as a singular fit); 4 application failure (`out.fail`, or any failed `out.check`/`out.check_close`) |
| state | `--dsl A B C` loads and verifies every file first, then runs them in order. Objects stored with `keep(name, obj)` stay reachable from later chunks through `kept_vec`/`kept_frame` |
| files | `read_csv` reads any path; `Frame.write_csv` writes only relative paths below the working directory (no `..`), atomically |
| numbers | IEEE binary64; NaN is R's NA. JSON output writes NaN/Inf as `null` |

The schema hash binds bytecode to the exact schema bytes. After the API changes, recompile the `.lms` files; the loader
refuses stale bytecode with `LOAD_SCHEMA_MISMATCH`.

## Conventions

- **Indices are 0-based** (`v.at(0)`, `fit.name(0)`, `Hclust.height(0)`). **Labels keep R's values**: k-means and
  `cutree` cluster labels are 1..k, `Hclust.merge_a/b` use R's sign convention, and `Vec.rank` gives 1-based ranks.
- **R parameterisations**: `gamma_dist(shape, rate)`, `weibull(shape, scale)`, `nbinom(size, prob)`,
  `hyper(m, n, k)`, `binom(size, prob)`.
- **Builders** return themselves, so calls chain: `test("t").x(a).y(b).var_equal(true).run()`.
- **New objects**: arithmetic and transforms return a new `Vec`. Only `push`, `set`, `append` (Vec) and `set` (Frame)
  mutate their receiver.
- **Errors** are host errors with a message that names the call, e.g. `dsl line 12:9: Vec.at: index 30 out of range 0..29`.
- DSL literals: floats need `digits.digits` (`0.000001`, not `1e-6`); widen ints with `float(i)`.

## Types

### `out` (the context object, type `Out`)

| Method | Meaning |
|---|---|
| `print(string)` | one line of report text (stderr) |
| `check(bool, string) -> bool` | records a self-check; a false one prints `CHECK FAILED: msg` and makes the run exit 4 |
| `check_close(actual, expected, tol, label) -> bool` | `abs(actual - expected) <= tol` (and actual finite) |
| `set(key, float/int/string/bool/Vec/StrList)` | a value in the result JSON for the calling lmscript (multi-step flows) |
| `fail(string)` | stop now with exit 4 |
| `failures() -> int` | failed checks so far |

### `Vec`: numeric vector

- **Construction**: `vec()`, `zeros(n)`, `rep(x, n)`, `seq(from, to, by)`, `seq_len(n)` (1..n), `linspace(a, b, n)`,
  `parse_vec("1 2.5 NA -Inf")`, `Frame.col(name)`, `Rng.draw(dist, n)`.
- **Access**: `len`, `at(i)`, `set(i, x)`, `push(x)`, `append(v)`, `copy`, `slice(start, end)`, `head(k)`, `tail(k)`,
  `take(indices)`, `where(mask)` (keep where mask != 0), `rev`, `sort`, `order` (0-based indices), `rank`.
- **Summaries**: `sum`, `mean`, `var`, `sd`, `median`, `mad`, `min`, `max`, `prod`, `which_min`, `which_max`,
  `count_nan`, `drop_nan`, `quantile(p[, type])`, `quantiles(probs, type)` (Hyndman & Fan types 1-9),
  `cor(v[, "pearson"|"spearman"|"kendall"])`, `cov`, `dot`.
- **Arithmetic**: operators `+ - *` and unary `-` work Vec-Vec (equal lengths) and Vec-float, with floats on either
  side (`2.0 * x`, `1.0 - x`). Methods: `add`, `sub`, `rsub`, `mul`, `div`, `rdiv`, `neg`, `pow(p)`, `exp`, `log`,
  `sqrt`, `abs`, `round(digits)`, `cumsum`, `diff([lag])`, `standardize`.
- **Time series**: `acf([lag_max])` (lags 0..lag_max), `pacf([lag_max])` (lags 1..lag_max).
- **Text**: `str()` gives `[1, 2.5, 2598960]` (shortest round-trip text, NaN as `NA`); `format("%.3f", ", ")` uses an explicit format.

### `StrList`, `Frame`

- **StrList**: `strs()`, `split(text, sep)`, then `len`, `at`, `push`, `join(sep)`, `contains`, `index` (-1 when absent).
- **Frame**: `read_csv(path)` (RFC 4180; a column is numeric when every cell is a number or NA/empty) and `frame()`.
  Methods: `nrow`, `ncol`, `names`, `has`, `is_text`, `col(name) -> Vec`, `text(name) -> StrList`,
  `set(name, Vec|StrList)`, `select("a,b")`, `drop(name)`, `where(mask)`, `head(k)`, `drop_nan`, `write_csv(path)`.
- **Keeping objects**: `keep(name, Vec|Frame)`, `kept_vec(name)`, `kept_frame(name)`, `is_kept(name)` carry objects
  across chunks of one `--dsl` run.

### Math functions

- **Elementary**: `exp`, `log`, `log10`, `log2`, `log1p`, `expm1`, `tan`, `asin`, `acos`, `atan`, `atan2`, `sinh`,
  `cosh`, `tanh`, `ceil`, `trunc`, `round(x[, digits])` (half to even, like R), `is_nan`, `is_finite`, `nan()`,
  `inf()`, `pi()`. The DSL itself already has `sin`, `cos`, `sqrt`, `abs`, `floor`, `pow`, `min`, `max`, `clamp`.
- **Special**: `gamma`, `lgamma`, `digamma`, `trigamma`, `beta`, `lbeta`, `choose`, `lchoose`, `factorial`,
  `lfactorial`, `erf`, `erfc`, `bessel_j(x, nu)`, `bessel_y`, `bessel_i(x, nu[, scaled])`, `bessel_k`,
  `gamma_p(a, x)`, `beta_inc(a, b, x)`.
- **Formatting**: `fmt("%.4f", x)` takes exactly one `%f %e %g %d %s` conversion with flags, width and precision, plus
  any literal text. `lpad(s, w)`, `rpad(s, w)`.

### `Dist`: 21 distributions with R semantics

- **Constructors**: `norm()`, `norm(mean, sd)`, `lnorm`, `t_dist(df)`, `chisq(df)`, `f_dist(df1, df2)`,
  `gamma_dist(shape, rate)`, `beta_dist(a, b)`, `exp_dist(rate)`, `unif(min, max)`, `weibull`, `cauchy`, `logis`,
  `binom(size, prob)`, `pois(lambda)`, `geom(prob)`, `nbinom(size, prob)`, `hyper(m, n, k)`, and
  `dist("name", p1[, p2[, p3]])`.
- **Noncentral and studentized range**: `t_dist(df, ncp)`, `chisq(df, ncp)`, `f_dist(df1, df2, ncp)`, and
  `tukey(nmeans, df[, nranges])` (R's `ptukey`/`qtukey`: `p` and `q` only; `df` may be `inf()`). These power
  `power.t.test`-style calculations and Tukey HSD intervals; see `examples/r-vs-dsl/power-noncentral.lms` and
  `tukey-hsd.lms`.
- **Methods**: `d`, `log_d`, `p`, `p_upper`, `log_p`, `log_p_upper`, `q`, `q_upper`, `q_log`, `mean`, `var`, `sd`,
  `is_discrete`, `name`, and vectorised `d_vec`, `p_vec`, `q_vec`.
- Upper tails and logs are computed directly, as in R, so `norm().p_upper(40.0)` and `t_dist(5.0).log_p(-200.0)` are
  accurate.

### `Rng`: seeded random numbers

- **Construction**: `rng(seed)`.
- **Methods**: `unif()`, `norm()`, `below(n)`, `draw(dist) -> float`, `draw(dist, n) -> Vec`,
  `sample(x, size, replace[, prob])`, `permutation(n)`.
- **Generator**: xoshiro256** seeded through splitmix64, with exact published variate methods. Streams are
  reproducible across platforms but are **not R's streams**, a decision recorded in `durable-memory/status.md`. Checks
  must therefore be statistical, never "R's first draw".

### Hypothesis tests: `test(kind)` → `Test` → `run()` → `Htest`

Kinds: `t`, `wilcox`, `var`, `cor`, `shapiro`, `ks`, `chisq`, `fisher`, `binom`, `prop`, `kruskal`, `bartlett`,
`fligner`, `oneway` (Welch unless `var_equal(true)`), `friedman`, `quade`, `mcnemar`, `ansari`.

| Setter | Used by |
|---|---|
| `x(Vec)`, `y(Vec)` | samples (one or two) |
| `dist(Dist)` | one-sample `ks` |
| `alt("two.sided"\|"less"\|"greater")`, `conf(level)` | all with intervals |
| `mu(float)`, `paired(bool)`, `var_equal(bool)` | t, wilcox |
| `exact(bool)`, `correct(bool)`, `conf_int(bool)` | wilcox, ks, cor; chisq/prop continuity |
| `ratio(float)` | var |
| `method("pearson"\|"kendall"\|"spearman")` | cor |
| `table(values, nrow, byrow)` | chisq (independence), fisher (R fills by column unless byrow) |
| `counts(Vec)`, `probs(Vec)` | chisq goodness of fit |
| `successes(x)`, `trials(n)`, `p(float \| Vec)` | binom (one value each), prop (vectors) |
| `odds_ratio(float)` | fisher null |
| `groups(Vec \| StrList)` | kruskal, bartlett, fligner, oneway; friedman, quade (with `blocks`); text labels sort in byte order |
| `blocks(Vec \| StrList)` | friedman, quade; or give `.table(values, nblocks, byrow)` with blocks as rows, groups as columns |
| `table(values, nrow, byrow)` | mcnemar (k x k) |

`Htest` methods:

- `statistic`, `statistic_name`, `p_value`, and the parameters `df()`, `param("num df")`, `n_params`,
  `param_name(i)`, `param(i)`.
- `has_ci`, `conf_low`, `conf_high`, `conf_level`.
- `estimate()` and `estimate(i)`, `estimate_name(i)`, `n_estimates`.
- `null_value`, `method`, `alternative`, `expected()` (chisq/prop), and `summary()` (R-like print text).

The method strings are byte-identical to R's.

### Models

- **`lm(frame, formula[, weights])` → `Lm`.** Formulas: `y ~ a + b`, `a:b`, `a*b`, `(a+b)^2`, `- 1` / `0 +`, `.`, and
  `log()`, `exp()`, `sqrt()`, `I(x^k)` of a column or of the response.
  - **Factors**: a text column is a factor, and `factor(v)` makes a numeric column one. Levels sort in byte order (as
    `read.csv` + R's C locale), treatment contrasts use the first level as the baseline, and columns are named as in R
    (`suppVC`, `factor(dose)2`, `suppVC:factor(dose)2`). Predicting with a level not seen in the fit is an error.
  - **ANOVA**: `Lm.anova()` (sequential sums of squares, R's `anova.lm`), `Glm.anova(test)` (analysis of deviance;
    test `"none"`, `"Chisq"`/`"LRT"` or `"F"`), and `anova(smaller, bigger)` (nested F test). An `Anova` has
    `n_rows`, `row_name(i)`, `columns`, `col(name)` → Vec, `at(row name | index, column)` and `table()` (printable),
    with R's column names (`"Df"`, `"Sum Sq"`, `"F value"`, `"Pr(>F)"`, `"Deviance"`, `"Resid. Dev"`, `"Pr(>Chi)"`). Rows with NaN in a used column are dropped
  (R's `na.omit`). Aliased (collinear) terms are reported as in R, through `aliased(name)` and `rank()`.
- **`glm(frame, formula, family[, link[, weights[, offset]]])` → `Glm`.** Families: gaussian, binomial, poisson and
  Gamma, with R's allowed links. A binomial response is 0/1, or proportions with the trials as weights. An empty `vec()`
  means no weights. The offset is added to the linear predictor, e.g. `f.col("exposure").log()` for rates.
- **Shared coefficient methods**, for `Lm`, `Glm` and `Nls`:
  - `n_coef`, `name(i)`, `names`, `has_coef`.
  - `estimate`, `se`, `statistic`, `p_value` and `aliased`, each indexed by an int or a term name.
  - `conf_low`/`conf_high(term, level)`: t intervals for lm, nls and gaussian/Gamma glm; Wald normal intervals for
    binomial and poisson.
  - `estimates()` → Vec, and `coef_table()` (printable).
- **`Lm` methods**: `sigma`, `r_squared`, `adj_r_squared`, `f_statistic`, `f_df1`, `f_df2`, `f_p_value`, `rss`,
  `loglik`, `aic`, `bic`, `n`, `rank`, `df_residual`, `residuals`, `fitted`, `predict(frame)`, `predict_se(frame)`.
- **`Glm` methods**: `deviance`, `null_deviance`, `aic`, `dispersion`, `family`, `link`, `n`, `df_residual`,
  `df_null`, `iter`, `converged`, `fitted` (response scale), `linear_predictors`,
  `residuals("deviance"|"pearson"|"working")`, `predict(frame, "link"|"response")`.
- **`nls(y)` → `Nls`.** This uses reverse communication (below). Configure it with `.param(name, start)`,
  `.lower(name, x)` and `.upper(name, x)` (which need `.algorithm("port")`), `.maxiter(n)` and `.tol(x)`. Then
  `while (m.next()) { for (i ...) { m.tell(i, model(m.par("Vmax"), m.par("Km"), x.at(i))); } }`.
  - Results: the coefficient methods plus `sigma`, `rss`, `df_residual`, `iterations`, `converged`, `message`,
    `conv_criterion`, `fitted`, `residuals`, `evaluations`.

### Smoothing and interpolation

- **`spline_fun(x, y[, "fmm"|"natural"|"periodic"|"monoH.FC"])` → `Spline`**: `at`, `deriv(x, order)`, `eval(Vec)`.
- **`approx_fun(x, y[, "linear"|"constant", rule[, f]])` → `ApproxFun`**: `at`, `eval`.
- **Smoothing splines**: `smooth_spline(x, y)` chooses lambda by GCV; `smooth_spline(x, y, true)` uses
  leave-one-out CV; `smooth_spline_df(x, y, df)` and `smooth_spline_lambda(x, y, lambda)` fix the smoothing.
  - `SmoothSpline` methods: `at`, `deriv`, `eval`, `x`, `y`, `leverages`, `df`, `lambda`, `cv_crit`, `pen_crit`.
  - This is a penalised natural cubic spline on all unique x (Reinsch). R's `spar` scale is deliberately not
    reproduced. Results with a given `df` agree with R where R uses all knots (nx ≤ 49).
- **`loess(x, y[, span[, degree[, family]]])` → `Loess`**: `at`, `eval`, `fitted`, `residuals`,
  `robustness_weights`, `enp`, `s`, `trace_hat`, `one_delta`, `two_delta`.
  - This uses exact local fits (R's `surface = "direct"`, which it matches to about 1e-15). R's default
    interpolation surface is an approximation that can differ by several percent of the noise on small noisy samples.
    `enp` is tr(L'L), as R's `$enp`, and `trace_hat` is tr(L). The cost is O(n³): about 1 s at n = 1000.

### Time series

- **ARIMA**: `arima(p, d, q)` → `ArimaSpec`, configured with `.seasonal(P, D, Q, period)`,
  `.method("CSS-ML"|"ML"|"CSS")` and `.include_mean(bool)`, then `.fit(Vec)` → `Arima`.
  - `Arima` methods: `n_coef`, `name`, `estimate`, `se`, `statistic`, `p_value`, `has_coef`, `sigma2`, `loglik`,
    `aic`, `nobs`, `convergence`, `residuals`.
  - `forecast(h)` → `Forecast`, with `mean`, `se`, `mean_at`, `se_at`, `lower(i, level)`, `upper(i, level)`.
- **`arima_sim(rng, n, ar, ma, sd[, d])` → Vec.** This rejects a non-stationary AR part and uses R's burn-in rule.
- **STL**: `stl(x, period, s_window[, robust])` and `stl_periodic(x, period[, robust])` → `Stl`, with `seasonal`,
  `trend`, `remainder`, `weights`.

### Clustering

The data is the numeric columns of a `Frame`; use `select` first if it also has text columns.

- **`kmeans(frame, k, rng[, nstart[, iter_max, algorithm]])` → `Kmeans`.** Algorithms: Hartigan-Wong, Lloyd and
  MacQueen.
  - Methods: `k`, `n`, `cluster(row)` (labels 1..k), `clusters()`, `size(label)`, `withinss(label)`,
    `center(label, col | "name")`, `tot_withinss`, `betweenss`, `totss`, `iter`, `ifault`.
- **`hclust(frame, method[, dist_method[, p]])` → `Hclust`.**
  - Methods: `n`, `height(i)`, `heights`, `merge_a(i)`, `merge_b(i)`, `order` (1-based leaves), `cut(k)` (labels by
    first appearance).
  - Merges, heights and order equal R's, including R's tie order.

### Matrices and multivariate analysis

- **`Mat`** (dense): `matrix(values, nrow, byrow)` (R fills by column unless `byrow`), `as_matrix(frame)` (numeric
  columns), `diag(Vec)`, `identity(n)`.
  - Access: `nrow`, `ncol`, `at(i, j)`, `set(i, j, x)`, `row(i)`, `col(j)`, `values()` (column-major, as R's
    `as.vector`), `diag()`, `to_frame()` (columns V1..Vp), `str()`, `format(spec)`.
  - Algebra: `t`, `mul(Mat)`, `mul_vec(Vec)`, `add`, `sub`, `scale(x)`, `solve(Vec)` and `inverse` (LU with partial
    pivoting), `det`, `chol` (upper R with R'R = A, as R's `chol`), `cov`, `cor`, `dist(method)`.
  - `eigen()` → `Eigen` (`values` decreasing, `vectors` as columns; symmetric matrices) and `svd()` → `Svd` (`d`, `u`,
    `v`). Each vector is unit length with its largest component positive (R's signs are arbitrary).
- **`prcomp(frame | Mat, scale[, …])`** → `Pca`: `prcomp(frame, scale)` centres; `prcomp(frame, center, scale)`.
  Methods `sdev`, `rotation` (Mat), `scores` (Mat), `center`, `scale`, `proportion`, `cumulative`,
  `loading(variable, k)`, `summary()`.
- **`cmdscale(distances, k)`** → `Mds` (`points`, `eig`, `gof`, `warning`), classical multidimensional scaling.
- **`mahalanobis(x, center, cov)`** → squared distances, as R's `mahalanobis`.

### Density estimation and histograms

- **`bw(x, rule)`**: `"nrd0"`, `"nrd"`, `"ucv"`, `"bcv"`, `"SJ"` (`"SJ-ste"`), `"SJ-dpi"`, as R's `bw.*`.
- **`density(x[, rule | bw, kernel[, adjust]])`** → `Density` on R's 512-point grid (from min − 3 bw to max + 3 bw),
  kernels gaussian, epanechnikov, rectangular, triangular, biweight, cosine, optcosine with bw the kernel's sd.
  Methods `x`, `y`, `bw`, `kernel`, `at(x0)`, `eval(Vec)`. Values are exact kernel sums; R bins and uses an FFT, so it
  differs by its binning error (about 1e-3 of the peak).
- **`ecdf(x)`** → `Ecdf` (`at`, `eval`, `knots`, `quantile(p)`); **`hist(x[, rule | breaks])`** → `Hist` (`breaks`,
  `counts`, `density`, `mids`; R's Sturges/Scott/FD rules and right-closed bins); **`pretty(lo, hi, n)`**.

### FFT

`fft(Vec)` returns a `CVec`; `cvec(re, im)` builds one directly. `CVec` methods: `fft(inverse)`, `len`, `re(i)`,
`im(i)`, `mod(i)`, `arg(i)`, `re_vec`, `im_vec`, `mod_vec`, `scale(x)`. Neither direction is normalised, as in R's
`fft()`.

### Reverse communication: `minimize`, `optimize`, `uniroot`, `integrate`, `nls`

The objective is ordinary DSL code. The host runs the native algorithm on a worker thread. Each `next()` returns `true`
when the algorithm needs a value at a new point and `false` when it has finished. Exactly one side runs at a time, so
results are deterministic.

```
let o = minimize(parse_vec("-1.2 1"), "BFGS").maxit(500);
while (o.next()) {
  let a = o.x(0);
  let b = o.x(1);
  o.tell(100.0 * pow(b - a * a, 2.0) + pow(1.0 - a, 2.0));
}
out.print(fmt("%.6f", o.par(0)) + fmt(" %.6f", o.par(1)) + fmt("  f=%.3e", o.value()));
```

| Builder | Options (before the first `next`) | Results (after `next` is false) |
|---|---|---|
| `minimize(start, method)`: Nelder-Mead, BFGS, CG, L-BFGS-B, SANN | `lower(Vec)`, `upper(Vec)`, `maxit`, `reltol`, `ndeps`, `hessian(bool)`, `seed(int)` (SANN) | `par(i)`, `pars`, `value`, `fn_count`, `gr_count`, `convergence`, `message`, `hessian(i, j)`, `evaluations` |
| `optimize(lower, upper[, maximum])` | `tol` | `argmin` (the maximiser when maximum), `objective`, `evaluations` |
| `uniroot(lower, upper)` | `tol`, `maxiter`, `extend("no"\|"yes"\|"upX"\|"downX")` | `root`, `f_root`, `iter`, `init_it`, `estim_prec`, `converged` |
| `integrate(lower, upper)`, bounds may be `inf()`/`-inf()` | `subdivisions`, `rel_tol`, `abs_tol` (defaults to `rel_tol`, as R) | `value`, `abs_error`, `subdivisions_used`, `ier`, `ok`, `message`, `evaluations` |

Protocol errors are run errors with a precise message: a `tell` without a pending point, a second `tell`, `next()`
without a `tell`, reading results early, or setting options after the start. A non-finite integrand stops integration
with `ier` 6 and a message rather than an error, so a chunk can test for it. If a chunk abandons a solver part way
through, the worker is cancelled cleanly when the host exits.

## Mapping from the old R-backed operations

| `rcli` op (R) | `_dsl_` API |
|---|---|
| `dist_quantile` / `dist_density` / `dist_probability` | `Dist.q` / `d` / `p` (+ `_upper`, `log_`, `_vec`) |
| `rng_draw`, `rng_sample` | `Rng.draw(dist, n)`, `Rng.sample` |
| `math_special` | `gamma`, `lgamma`, `beta`, `lbeta`, `digamma`, `trigamma`, `choose`, `lchoose`, `bessel_*` |
| `summary_quantile` | `Vec.quantiles(probs, type)` |
| `num_fft` | `fft(Vec)`, `CVec.fft(inverse)` |
| `test_htest` | `test(kind)...run()` |
| `fit_lm`, `fit_glm`, `fit_nls` | `lm`, `glm`, `nls` (model written in the chunk) |
| `smooth_spline`, `smooth_loess`, `interp_spline`, `interp_approx` | `smooth_spline*`, `loess`, `spline_fun`, `approx_fun` |
| `ts_arima`, `rng_arima_sim`, `ts_acf`, `ts_pacf`, `ts_stl` | `arima(...).fit`, `arima_sim`, `Vec.acf`, `Vec.pacf`, `stl` |
| `cluster_kmeans`, `cluster_hclust` | `kmeans`, `hclust` |
| `opt_minimize`, `opt_optimize`, `num_root`, `num_integrate` | `minimize`, `optimize`, `uniroot`, `integrate` (objective in DSL) |
| `r_info` | `lmr --version` |

## Agreement with R

Each module has a unit test (`native/tests/test_*.cpp`) whose expectations were produced by running R 4.3.3 as a black
box. Measured agreement is recorded in `durable-memory/native-dsl.md`. The deliberate differences are:

- the random streams;
- smoothing-spline lambda selection (R's `spar` and knot subset are not reproduced);
- loess, which uses exact rather than interpolated surfaces;
- iteration counts of some optimiser paths (L-BFGS-B, uniroot, CG Polak-Ribiere);
- robust STL on even-length series (R uses an off-median scale at some iterations);
- near-non-stationary seasonal ARIMA fits. Ordinary ARIMA fits agree with R to about 1e-11.

R itself is never linked or read: it was used only as an oracle process.
