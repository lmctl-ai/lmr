## 0.5.1

- Omit built-in StrList methods from the generated host schema; the lmcode standard registry provides them by default.

## 0.5.0

- Adopt lmcode-cpp (renamed from cpp-dsl-sdk). Its built-in lmcode std (72 signatures: `str_*`, `math_*`, `file_*`, `StrList`, `IntList`, `FloatList`, `StrMap`) is available in every chunk without lmr registering it; the schema carries the std declarations until the compiler injects them.
- `StrList` is now std's type. Renamed lmr-only names: `strs()` → `str_list()`, `split` → `str_split`, `StrList.push` → `add`, `StrList.index` → `index_of`.
- Blocks are spelled `_lmcode_ ... end _lmcode_` (same bytecode as `_dsl_`); needs lmscript 0.2.420 or newer.

## Unreleased

- Make `lmr --version` print the exact `lmr X.Y.Z` release identity and guard every release tag against the built binary version.
- Add compact install-enabled LMScript examples for quantiles, linear regression, and BFGS optimization, with expected-output checks.
- Match R expression semantics for logical precedence, defaults, interval extension, scalar integration, and NLS non-convergence; compile expression ASTs once per call.
- Add 2,000-case deterministic parser differential coverage, expression guard/throughput regressions, method-specific optimizer defaults, and the R third-party notice.
- Add the shared expression grammar, C++ parser, R-side AST validation, and the optimization, root-finding, integration, and nonlinear least-squares operations.
- Add expression conformance/oracle tests and the noisy curve-fit and Rosenbrock lmscript recipes.
- Added base R ARIMA fitting and forecasting, seeded ARMA simulation, ACF/PACF, STL decomposition, Hartigan-Wong k-means, and hierarchical clustering.
- Added the `fit-and-forecast.lms` seeded simulation-to-forecast recipe.

# Changelog

## Unreleased

- Add vectorized density/CDF operations, expanded seeded draws and sampling, special math functions, quantile types 1–9, and FFT/inverse FFT.
- Add the enum-based `test_htest` operation for common base R tests, with shared typed results and pre-R validation.
- Add GLM fitting with convergence status and base R spline/loess smoothers plus interpolation methods.
- Stream large quantile probabilities and inline model columns through temporary numeric files, load them with vectorized `scan()`, and encode numeric result arrays in one vectorized pass. The 200,000-probability smoke test enforces a five-second ceiling.
