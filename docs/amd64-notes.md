# amd64 notes: the CI failure and the float-comparison audit

rcli is developed on aarch64 (R 4.3.3-2build2). CI runs the same R on `ubuntu-24.04-arm` and `ubuntu-24.04`. Every R oracle in
`tests/` is computed by the same R binary on the same machine as the tool, so values produced by the same R code path are
bit-identical on either architecture; only values produced by two different code paths can differ by platform rounding.

## The amd64 failure (runs 38523688440 and 38840663416)

`rcli_ops` / `tests/compare-group3.py`: `fit_glm` `null_deviance` 1.695819109974191 (tool) vs 1.6958191099741908 (oracle), one ulp,
under exact equality. The worker always passes `weights` (all 1) and an `offset` (all 0) to `glm()`; with an offset present, R's
`glm()` recomputes `null.deviance` by an intercept-only IRLS fit (`fit2 <- method(x = X[, "(Intercept)"], ..., offset = offset)`),
while the oracle's plain `glm()` takes the closed form `sum(dev.resids(y, weighted mean, weights))` from `glm.fit`. Two evaluation
orders of the same number: identical on aarch64, one ulp apart on amd64 (FMA / vectorised sums). Not a tool defect.

Fix: a relative tolerance of 1e-8 (plus 1e-9 absolute) for `null_deviance` only, derived from the IRLS stopping rule that produces the
tool's value (`glm.control(epsilon = 1e-8)`: `|dev - devold| / (|dev| + 0.1) < epsilon`). Every other `fit_glm` field (estimates,
standard errors, statistics, p-values, deviance, AIC, iterations) comes from the same IRLS iterations in both paths (adding a zero offset
and unit weights is exact) and stays under exact equality. Mutation check at the tool boundary (`null_deviance x (1 + 1e-7)`): the test
fails; `x (1 + 1e-9)`: passes, as derived.

## Audit of exact and sub-1e-12 float comparisons

| Where | Comparison | Classification | Action |
|---|---|---|---|
| `tests/compare.py` (`fit_lm`, `dist_quantile`, `rng_draw`) | exact vs `tests/oracle.R` | same R code path (`lm`, `qnorm`, `rnorm` with the same RNG kind), same binary: deterministic | kept exact |
| `tests/compare-dist.py` | exact vs `dist-oracle.R` | same quantile functions | kept exact |
| `tests/compare-htest.py` | exact vs `htest-oracle.R` (statistics, parameters, p-values, intervals) | same `*.test` calls | kept exact |
| `tests/compare-group3.py` | exact vs `group3-oracle.R` | same path for smooth_spline, loess, interp_spline, interp_approx and all `fit_glm` fields except `null_deviance` (above) | `null_deviance` 1e-8 relative; rest exact |
| `tests/compare-group4.py` | exact vs an inline R oracle (arima, acf, pacf, stl, kmeans, hclust) | same path, seeded RNG kind set explicitly | kept exact |
| `tests/stage1a.py` | exact vs `stage1a-oracle.R` (math_special, summary_quantile, num_fft) | same path | kept exact |
| `tests/review2-regressions.py` quoted-CSV `fit_lm` | exact vs literal doubles `[-0.4999999999999992, 1.3999999999999997]` | **platform-dependent**: the literal bits are LAPACK QR on aarch64; amd64 may round the last bit differently | relative 64 eps = 1.4e-14 against the exact `[-0.5, 1.4]` (cond([1, x], x = 1..4) is 10, so cond * eps = 2e-15 is the backward-error bound); mutation `x (1 + 1e-13)` fails, `x (1 + 1e-15)` passes |
| `tests/review2-regressions.py` two-sample KS p-value | absolute 1e-14 vs 0.0285714285714286 (= 1/35, a closed-form rational) | legitimately tight (exact rational reached by a fixed computation) | kept |
| `tests/smoke.sh` `qnorm` log-p | absolute 1e-14 vs 1.9599639845400536 | same R `qnorm`, deterministic | kept |
| `tests/expr-ops.py` | `near` with 2e-12 relative / absolute | iterative optimisers against R references computed by the same R | kept |
| `examples/oracle-lm.lms` | max absolute 1e-12 against the same literal `lm` coefficients | 1e-12 is 500 x the backward-error bound 2e-15: platform-safe | kept (documented) |
| `examples/oracle-*.lms`, `fit-*.lms`, `minimize-rosenbrock.lms` | absolute 1e-8 / 1e-9 / 1e-3 | cross-tool tolerances, stated in each recipe header | kept |

## CI

`.github/workflows/ci.yml` now runs the build, CTest (`--output-on-failure`), the no-libR check, the schema gate (`rcli --schema`
against the checked-in `cli.schema.dev.json`) and the lmts parse as separate gates with `if: !cancelled()`, with a gate table in the
job summary; the job still fails when any gate failed.

## Status

amd64 CI status: see README "Build and test"; updated from the CI run, never from local reasoning.
