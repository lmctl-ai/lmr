# R vs lmr: common R recipes, same problem, same data

Each recipe here is a well-known R example, taken from the R help pages (`?optim`, `?arima`, `?chisq.test`, ...) or
the standard tutorial exercises (`lm(mpg ~ wt + hp, mtcars)`, `kmeans(iris)`, `hclust(dist(scale(USArrests)))`, ...).
It is solved twice:

- **R**: `r/NAME.R` is the recipe as you would write it in R. It saves R's results with full precision to
  `expected/NAME.json`.
- **lmr**: `NAME.lms` solves the same problem on the same data in a `_dsl_` chunk run by the native `lmr` binary, then
  prints R's value and lmr's value side by side and checks each one against a stated tolerance.

R is only needed to regenerate the expected results. Running the comparison needs only lmscript and `lmr`.

```sh
cmake -S native -B build-native -G Ninja && cmake --build build-native   # once, from the repository root
lmscript examples/r-vs-dsl/lm-mtcars.lms                                 # one recipe
native/scripts/run-examples.sh build-native/bin/lmr                      # every example, as ctest does
examples/r-vs-dsl/refresh-r.sh                                           # regenerate R's side (needs Rscript)
```

Sample output (`lm-mtcars`):

```
lm-mtcars: R (R version 4.3.3 (2024-02-29)) vs lmr
value                                      R                 lmr      diff      tol
estimate[3]                      37.22727012         37.22727012   4.6e-16    1e-12  ok
std_error[3]                     1.598787538         1.598787538   4.4e-16    1e-12  ok
p_value[3]                   2.565458512e-20     2.565458512e-20   6.1e-18    1e-12  ok
...
residuals[32]                   -2.572329403        -2.572329403   2.6e-14    1e-12  ok
17 values compared, 0 outside tolerance
lm-mtcars: lmr agrees with R.
```

`diff` is |lmr − R| / max(1, |R|): an absolute difference for values below 1 and a relative one above. For an array
(`[32]`) the row shows the first element and the largest difference over all elements.

## Recipes

Results are from R 4.3.3 and lmr 0.3.0, grouped by topic. "Numbers" counts every compared element; "worst" is the
largest `diff`.

| Recipe | R | Numbers | Worst |
|---|---|---:|---:|
| [lm-mtcars](lm-mtcars.lms) | `lm(mpg ~ wt + hp, mtcars)`: summary, confint, predict with se | 62 | 2.6e-14 |
| [lm-cars-poly](lm-cars-poly.lms) | `lm(dist ~ speed + I(speed^2), cars)` | 68 | 2.7e-14 |
| [anova-plantgrowth](anova-plantgrowth.lms) | `aov(weight ~ group, PlantGrowth)` (dummies built in the chunk) | 17 | 1.3e-15 |
| [glm-logistic-mtcars](glm-logistic-mtcars.lms) | `glm(am ~ wt + hp, binomial, mtcars)` | 50 | 3.0e-15 |
| [glm-poisson-warpbreaks](glm-poisson-warpbreaks.lms) | `glm(breaks ~ wool + tension, poisson, warpbreaks)` | 16 | 4.2e-15 |
| [lm-interaction-weights](lm-interaction-weights.lms) | `lm(mpg ~ wt * hp)`, `lm(log(dist) ~ log(speed))`, weighted `lm` | 19 | 2.6e-15 |
| [lm-prediction-intervals](lm-prediction-intervals.lms) | `predict(interval = "confidence"/"prediction")` | 25 | 3.9e-15 |
| [anova-twoway-toothgrowth](anova-twoway-toothgrowth.lms) | `aov(len ~ supp * dose, ToothGrowth)`, sequential sums of squares | 24 | 3.9e-15 |
| [tukey-hsd](tukey-hsd.lms) | `TukeyHSD(aov(breaks ~ tension))`, `ptukey`, `qtukey` | 21 | 1.6e-8 |
| [pairwise-t-plantgrowth](pairwise-t-plantgrowth.lms) | `pairwise.t.test`: pooled SD, Holm, Bonferroni | 9 | 3.9e-16 |
| [factor-models](factor-models.lms) | `lm(len ~ supp * factor(dose))`, `anova()`, nested `anova(m1, m2)`, `anova(glm, test = "Chisq")` | 46 | 2.2e-14 |
| [glm-dobson-poisson](glm-dobson-poisson.lms) | the Poisson example of `?glm` (Dobson 1990) | 41 | 7.3e-15 |
| [glm-gamma-clotting](glm-gamma-clotting.lms) | the Gamma example of `?glm` (McCullagh & Nelder, clotting) | 13 | 2.1e-14 |
| [glm-binomial-links](glm-binomial-links.lms) | logit/probit/cloglog by AIC; grouped binomial with weights, LD50 | 18 | 3.1e-15 |
| [logistic-odds-ratios](logistic-odds-ratios.lms) | `exp(coef)`, `exp(confint.default())`, likelihood-ratio test | 11 | 1.5e-15 |
| [nls-puromycin](nls-puromycin.lms) | `nls(rate ~ Vm * conc / (K + conc), Puromycin)` (treated) | 20 | 8.3e-9 |
| [t-test-sleep](t-test-sleep.lms) | `t.test` on `sleep`: Welch, pooled, paired | 16 | 6.0e-16 |
| [wilcox-tests](wilcox-tests.lms) | `wilcox.test`: exact paired, ties with correction, Hodges-Lehmann interval | 9 | 5.0e-16 |
| [chisq-tests](chisq-tests.lms) | `chisq.test`: 2 × 3 independence, goodness of fit, Yates 2 × 2 | 13 | 6.9e-18 |
| [fisher-binom-prop](fisher-binom-prop.lms) | `fisher.test(TeaTasting)`, `binom.test(682, 925)`, `prop.test(smokers, patients)` | 20 | 4.4e-16 |
| [k-sample-tests](k-sample-tests.lms) | `kruskal.test`, `bartlett.test`, `fligner.test`, `oneway.test` (equal variances and Welch) | 13 | 1.7e-18 |
| [block-paired-tests](block-paired-tests.lms) | `friedman.test`, `quade.test`, `mcnemar.test` (help-page examples) | 11 | 3.5e-17 |
| [power-noncentral](power-noncentral.lms) | `power.t.test`, `power.prop.test`, noncentral `pt`/`qt`/`dt`/`pchisq`/`qchisq`/`pf`/`qf` | 26 | 9.7e-8 * |
| [cor-test-mtcars](cor-test-mtcars.lms) | `cor.test(mpg, wt)`: Pearson, Spearman, Kendall | 11 | 7.4e-16 |
| [cor-test-exact](cor-test-exact.lms) | exact Kendall and Spearman tests (Hollander & Wolfe, `?cor.test`) | 10 | 1.1e-16 |
| [contingency-mtcars](contingency-mtcars.lms) | `chisq.test(table(cyl, gear))`, `fisher.test(table(am, vs))` | 22 | 2.3e-16 |
| [normality-variance](normality-variance.lms) | `shapiro.test`, `var.test(mpg ~ am)`, `ks.test(mpg, "pnorm", 20, 6)` | 10 | 1.3e-10 * |
| [summary-quantiles](summary-quantiles.lms) | `summary`, `quantile(type = 1..9)`, `sd`, `IQR`, `mad`, `cor` | 54 | 1.1e-16 |
| [distributions](distributions.lms) | `dnorm`, `pt`, `qt`, `qchisq`, `qf`, `dbinom`, `dpois`, `qgamma`, `pbeta`, ... | 47 | 4.3e-16 |
| [pca-usarrests](pca-usarrests.lms) | `prcomp(USArrests, scale. = TRUE)`, `prcomp(iris)`, `eigen(cor(mtcars))` | 67 | 2.5e-6 * |
| [mds-mahalanobis](mds-mahalanobis.lms) | `cmdscale(eurodist)`, `mahalanobis`, `cor` | 116 | 1.4e-13 |
| [kmeans-iris](kmeans-iris.lms) | `kmeans(iris[, 1:4], 3, nstart = 25)` (label-free results) | 12 | 1.9e-15 |
| [hclust-usarrests](hclust-usarrests.lms) | `hclust(dist(scale(USArrests)))`, complete and average, `cutree(4)` | 252 | 1.9e-16 |
| [arima-lh](arima-lh.lms) | `arima(lh, c(1,0,0))`, `c(3,0,0)`, `c(1,0,1)`, `predict(n.ahead = 12)` | 44 | 4.0e-8 |
| [arima-airline](arima-airline.lms) | airline model on `log(AirPassengers)` with a 12-month forecast | 43 | 4.4e-7 |
| [arima-nile-diagnostics](arima-nile-diagnostics.lms) | `arima(Nile, c(1,1,1))` by CSS-ML, ML, CSS; Ljung-Box on residuals | 18 | 2.3e-10 |
| [acf-pacf](acf-pacf.lms) | `acf`, `pacf` of `lh` and `Nile` | 74 | 3.5e-16 |
| [stl-decomposition](stl-decomposition.lms) | `stl(nottem, "periodic")`, `stl(log(co2), s.window = 21)` | 1656 | 5.8e-13 |
| [optim-rosenbrock](optim-rosenbrock.lms) | `optim` on Rosenbrock: Nelder-Mead, BFGS, L-BFGS-B | 11 | see below * |
| [mle-gamma-optim](mle-gamma-optim.lms) | gamma MLE with `optim(BFGS, hessian = TRUE)`, standard errors | 8 | 3.6e-11 |
| [integrate-roots](integrate-roots.lms) | `integrate` (incl. infinite bounds), `uniroot`, `optimize` | 8 | 1.4e-16 |
| [density-ecdf](density-ecdf.lms) | `density()`, `bw.nrd0/nrd/ucv/bcv/SJ`, `ecdf`, `hist` | 42 | 1.3e-14 * |
| [smoothing-cars](smoothing-cars.lms) | `smooth.spline(df = 5)`, `loess`, `splinefun`, `approx` | 118 | 9.9e-6 * |
| [fft-spectrum](fft-spectrum.lms) | `fft(1:8)`, amplitude spectrum, inverse transform | 113 | 1.4e-15 |

In total, 44 recipes compare 3304 numbers.

Several recipes also show how to do in a few DSL lines what R does with a helper lmr does not have: R's treatment
dummies for factors, sequential ANOVA sums of squares, Holm's adjustment, `Box.test`, prediction intervals, Wald odds
ratios, and the contingency table itself.

## Known differences (*)

The comparison is a bug detector, not a byte-identity requirement: tolerances are set where a real bug would show.
Each known difference is explained in the recipe's header and compared with its own tolerance.

- **`density()`** (density-ecdf): R bins the data and convolves with an FFT; lmr sums the kernel exactly. Density
  values differ by R's binning error (about 3e-4 here) and are compared at 2e-3; bandwidths, grid, ECDF and
  histogram agree to 1e-14.
- **`summary(prcomp)`** (pca-usarrests): R rounds the proportions of variance to 5 decimals.
- **Sample sizes from `power.*.test`** (power-noncentral): R's `uniroot` stops at tol = 1.2e-4; lmr solves to 1e-10.

- **`ks.test` with ties** (normality-variance): both sides use the asymptotic Kolmogorov distribution. lmr's p-value,
  0.74573634804, is the limit series summed to double precision; R stops its series early and reports a value
  1.3e-10 higher.
- **`optim` BFGS evaluation count** (optim-rosenbrock): the iterates, the minimum and the gradient count agree to
  1e-13, but R reports 120 function evaluations and lmr 118.
- **`optim` L-BFGS-B** (optim-rosenbrock): a different implementation of the same algorithm. Both stop short of the
  exact minimum (1, 1), R at 0.99980 and lmr at 0.99991.
- **`smooth.spline(df = 5)`** (smoothing-cars): R's search for the smoothing parameter stops at df 5.00055, while lmr
  hits a requested df exactly. Asked for R's achieved df, lmr's curve agrees with R's to about 1e-5.
- **`loess` default surface** (smoothing-cars): lmr always computes the exact local fit (R's `surface = "direct"`, which
  it matches to 3e-15). R's default `surface = "interpolate"` is an approximation, and differs from the exact fit by up
  to 1.6% here; the recipe shows that difference rather than hiding it.
- **Random numbers**: lmr's generator is not R's, so `kmeans` is compared on results that do not depend on the random
  starts or on label numbering. With 25 starts both reach the same optimum.
- **Factors**: lmr's formulas take numeric columns, so the ANOVA, Poisson GLM and Puromycin recipes build R's
  treatment dummies or subsets from the text columns in the chunk.

## Files

- `r/NAME.R`, the R recipe. `r/_emit.R` writes the JSON with 17 significant digits. `r/export-data.R` writes `data/`.
- `data/*.csv`: the R datasets used (mtcars, cars, sleep, PlantGrowth, ToothGrowth, warpbreaks, iris, USArrests,
  faithful, Puromycin, lh, AirPassengers, nottem, co2, Nile), exported from R's `datasets` package by `r/export-data.R`.
- `expected/NAME.json`: R's results (R 4.3.3).
- `NAME.lms`: the lmr version and the comparison. `lib/compare.lms` prints the table and applies the tolerances.
