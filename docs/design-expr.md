# rcli stage 1b: caller-supplied expressions (opt_minimize, opt_optimize, num_root, num_integrate, fit_nls)

Decision: callers never send R code. They send an arithmetic **expression string** that a hand-written C++ recursive-descent parser turns into a small JSON AST; R rebuilds that AST into a call object from a fixed table of *function objects* (never names, never `parse(text=)`), and evaluates it in an environment whose parent is `emptyenv()`. The allowlist, arities and caps live in one spec file embedded into both sides at build time.

## 1. Grammar (ASCII only, whitespace = space/tab/newline, max 2048 bytes)
```
expr    := sum
sum     := product (('+' | '-') product)*
product := unary (('*' | '/') unary)*
unary   := '-' unary | power
power   := atom ('^' unary)?                      # right-assoc, as in R
atom    := NUMBER | IDENT | FUNC '(' expr (',' expr)* ')' | '(' expr ')' | 'ifelse' '(' cond ',' expr ',' expr ')'
cond    := rel (('&' | '|') rel)*                 # only inside ifelse
rel     := expr ('<' | '<=' | '>' | '>=' | '==' | '!=') expr
NUMBER  := [0-9]+ ('.' [0-9]*)? ([eE] [+-]? [0-9]+)?   # must be a finite double; no hex, no 1L, no Inf/NA/NaN/TRUE
IDENT   := [a-z][a-z0-9_]{0,31}                        # lowercase; must be a declared name (see §3)
```
`FUNC` is exactly the allowlist below, each with a fixed arity; `pi` and `e` are built-in constants (idents that are not declarable). Nothing else: no `[`, `[[`, `$`, `@`, `::`, strings, backticks, `=`/`<-`, `;`, comments, `%%`/`%/%`/`%in%`, `**`, `&&`/`||`, `!`, `function`/`\`, `c()`, `if`, `NULL`.

Allowlist (name → arity → R binding): `exp log log1p expm1 log2 log10 sqrt abs sin cos tan asin acos atan sinh cosh tanh floor ceiling trunc round sign gamma lgamma digamma pnorm dnorm qnorm plogis` → 1; `atan2 pmin pmax beta lbeta choose lchoose` → 2; `ifelse` → 3; reducers `sum mean prod min max` → 1 (produce a scalar from a vector; allowed only where the op expects a scalar objective, see §4). Operators bind to `base::"+" "-" "*" "/" "^" "<" "<=" ">" ">=" "==" "!=" "&" "|"`.

Caps enforced by the C++ parser before anything is spawned: ≤2048 bytes, ≤512 tokens, ≤512 AST nodes, nesting depth ≤32, call nesting ≤16, ≤64 distinct identifiers, number literals that `strtod` cannot round-trip as finite are rejected. Byte values outside 0x20–0x7E (plus `\t\n`) are rejected before tokenizing (kills Unicode homoglyphs and NUL).

## 2. AST transport and rebuild in R
C++ emits canonical JSON into the request file, e.g. `a*exp(-b*t)` → `{"t":"op","o":"*","a":[{"t":"var","n":"a"},{"t":"call","f":"exp","a":[{"t":"op","o":"*","a":[{"t":"neg","a":[{"t":"var","n":"b"}]},{"t":"var","n":"t"}]}]}]}`. Node kinds: `num{v}`, `var{n}`, `const{n}`, `neg{a}`, `op{o,a[2]}`, `call{f,a[]}`, `ifelse{c,a,b}`.

R rebuilds with a builder that is itself allowlist-driven:
```r
FN <- list(exp=base::exp, log=base::log, pnorm=stats::pnorm, ifelse=base::ifelse, `+`=base::`+`, ...)  # function OBJECTS from the spec
build <- function(n, names) switch(n$t,
  num   = n$v,
  const = c(pi=base::pi, e=base::exp(1))[[n$n]],
  var   = { if(!n$n %in% names) stop("EXPR_SPEC_MISMATCH"); as.name(n$n) },
  neg   = as.call(list(FN[["-"]], build(n$a[[1]], names))),
  op    = { f <- FN[[n$o]]; if(is.null(f)) stop("EXPR_SPEC_MISMATCH"); as.call(list(f, build(n$a[[1]],names), build(n$a[[2]],names))) },
  call  = { f <- FN[[n$f]]; if(is.null(f) || length(n$a)!=ARITY[[n$f]]) stop("EXPR_SPEC_MISMATCH"); as.call(c(list(f), lapply(n$a, build, names))) },
  ifelse= as.call(list(FN$ifelse, build(n$c,names), build(n$a,names), build(n$b,names))),
  stop("EXPR_SPEC_MISMATCH"))
compile <- function(ast, names) { body <- build(ast, names); function(env) eval(body, env) }   # env: list2env(vars, parent=emptyenv())
```
Because every function in the call object is the function value itself, symbol lookup can only ever resolve the declared variables; `eval` with `parent=emptyenv()` means even `pi` is unreachable unless the builder inserted its numeric value. The R side re-checks node kinds, arity and names (a mismatch is `INTERNAL_ERROR` with `tool_code:"EXPR_SPEC_MISMATCH"`, never a validation message, because C++ already validated). No `parse`, `eval(parse())`, `get`, `match.fun`, `do.call` with caller strings, or `as.formula` on caller text anywhere in the worker.

## 3. Names, vectors, data
- **Named scalar parameters only**, no `par[i]`: `start: {"a": 1, "b": 0.5}` declares the parameter names; the R wrapper converts between the named vector `optim` wants and the environment the expression reads. Names must match IDENT, be distinct, not collide with a function, constant or data column, and number ≤32.
- **Data columns** (`data: {"t": [...], "y": [...]}`, numeric-file transport, ≤200 000 rows, ≤32 columns, same name rules) are bound as vectors; the expression is evaluated vectorised over rows exactly as R does. For `num_root`/`num_integrate`/`opt_optimize` the single free variable is fixed to `x` and no data is accepted; the wrapper asserts `length(result) == length(x)` so a non-vectorised expression cannot silently return a scalar to `integrate`.
- **Scalar objective check**: for `opt_minimize` the wrapper asserts `length(value) == 1`; an expression over data columns without a reducer therefore fails with `tool_code:"EXPR_NOT_SCALAR"` (validation, from a dry evaluation at `start` before optimisation).
- `fit_nls` takes `response` (column name) and `rhs` (expression over parameters and columns); R gets `y ~ <rebuilt call>` constructed with `as.call(list(as.name("~"), as.name(response), body))` and `environment(form) <- emptyenv()`-derived env holding only the data; `nls(..., start=, algorithm = "default"|"port", lower/upper for port, control=nls.control(maxiter, tol, warnOnly=TRUE))`. `plinear` and `weights` are out of scope.

## 4. Operations
| Op | Params | Result |
|---|---|---|
| `opt_minimize` | `objective`, `start{}`, `method` ∈ Nelder-Mead/BFGS/CG/L-BFGS-B/SANN, `lower{}`/`upper{}` (L-BFGS-B only), `maxit` ≤10000, `reltol`, `seed` (required iff SANN, rejected otherwise), `hessian` bool, `data{}` optional, `nonfinite` ∈ `inf`(default)/`error` | `par{}`, `value`, `counts{function,gradient}`, `convergence` (R code), `message`, `domain_status` OK/NON_CONVERGED/MAXIT, `hessian` rows, `evals`, `nonfinite_evals`, `elapsed_ms`, `warnings[]` |
| `opt_optimize` | `objective` in `x`, `lower`, `upper`, `maximum` bool, `tol` | `x`, `value`, `evals` |
| `num_root` | `f` in `x`, `lower`, `upper`, `tol`, `maxiter` ≤1000, `extend_int` ∈ no/yes/upX/downX | `root`, `f_root`, `iter`, `estim_prec`, `evals` |
| `num_integrate` | `f` in `x`, `lower`/`upper` numbers or the special strings `"Inf"`/`"-Inf"`, `subdivisions` ≤1000, `rel_tol`, `abs_tol` | `value`, `abs_error`, `subdivisions`, `message`, `evals` |
| `fit_nls` | `response`, `rhs`, `data{}`, `start{}`, `algorithm`, `lower{}`/`upper{}` (port), `maxiter` ≤1000, `tol`, `fitted`/`residuals` bool | coefficient table (`coefficients[]`, `estimate[]`, `std_error[]`, `statistic[]`, `p_value[]`), `sigma`, `df`, `converged`, `iterations`, `achieved_tol`, `domain_status`, `evals` |

## 5. Runaway and non-finite guards
- Every compiled closure increments an evaluation counter; above `max_evals` (default 500 000, cap 2 000 000) it `stop()`s with a sentinel that maps to `COMMAND_FAILED`/`tool_code:"EVAL_LIMIT"`. `setTimeLimit(elapsed = timeout_s - 5, transient = TRUE)` is set in the prelude as a second layer; the lmcpp 50 s process kill remains the last.
- `nonfinite:"inf"` replaces NA/NaN in a minimisation objective by `+Inf` (so Nelder-Mead/SANN can step away) and counts it; L-BFGS-B and BFGS cannot accept Inf, so for those methods a non-finite value raises `tool_code:"NONFINITE_OBJECTIVE"` with the offending `par`. `num_root` validates finite, opposite-sign endpoints up front (`tool_code:"BAD_INTERVAL"`); `num_integrate` maps R's "non-finite function value" to `domain_status:"NON_FINITE"`. `fit_nls` singular gradient / step-factor failures are `domain_status` values with `converged:false`, not failures (protocol §6.1).
- SANN: `RNGkind("Mersenne-Twister","Inversion","Rejection"); set.seed(seed)` as in `rng_draw`; `meta.seed` echoed. All other methods are deterministic and `seed` is rejected (`UNKNOWN_PARAM` semantics via validation) to keep results reproducible by construction.
- `^` with huge exponents, `exp(1e6)` and division by zero are legal and produce Inf/NaN handled by the rules above; `gamma(1e6)` etc. are bounded by R's own overflow handling (returns Inf with a warning captured into `warnings[]`).

## 6. Error mapping
C++ parser: `VALIDATION_ERROR` with `tool_code` ∈ `EXPR_SYNTAX`, `EXPR_BAD_BYTE`, `EXPR_UNKNOWN_IDENT`, `EXPR_UNKNOWN_FUNCTION`, `EXPR_ARITY`, `EXPR_TOO_LONG`, `EXPR_TOO_DEEP`, `EXPR_TOO_MANY_NODES`, `EXPR_NAME_COLLISION`, `EXPR_NOT_SCALAR`, plus `position` (byte offset) and `near` (≤16 chars) in `data`. R rebuild mismatch: `INTERNAL_ERROR`/`EXPR_SPEC_MISMATCH`. Runtime: `COMMAND_FAILED` with `R_ERROR`, `EVAL_LIMIT`, `NONFINITE_OBJECTIVE`; timeouts and crashes as today.

## 7. Single source of truth and tests
`expr-spec.json` in the repo: `{functions:{exp:{arity:1,r:"base::exp"},…}, operators:[…], constants:{pi,e}, caps:{bytes,tokens,nodes,depth,calls,idents}}`. CMake embeds it into C++ (`configure_file` → `expr_spec.hpp`) and into `worker.R` (generated prelude defining `FN`/`ARITY`/caps from the same JSON), and a CTest asserts both generated artifacts carry the spec's SHA-256. A second CTest runs `tests/expr-corpus.json` (accept and reject cases with expected `tool_code`) through `rcli --expr-check` (a hidden read op that parses without spawning R) and through the R builder in isolation, so both sides agree on every corpus entry.

Negative corpus outline (≥45 cases): injection `system("id")`, `eval(parse(text="1"))`, `base::exp(1)`, backtick names, `x <- 1`, `x = 1`, `function(x) x`, `\(x) x`, `.Internal(…)`, `get("system")`, `match.fun("exp")`, `do.call(…)`, `Sys.setenv(…)`, `quote(x)`, `substitute(x)`, `environment()`, `x$y`, `x@y`, `x[[1]]`, `x[1]`, `c(1,2)`, `"a"`, `'a'`, `# comment`, `1;2`, `1\nsystem("id")`, Cyrillic `ехр(1)`, NUL byte, `1e999`, `0x10`, `1L`, `TRUE`, `NA`, `Inf`, `NULL`, `exp()`, `exp(1,2)`, `log(x,2)`, `pi(1)`, `2**3`, `7%%2`, `x %in% y`, `a && b`, `!x`, `--x`, `x ? y`, 33 nested parentheses, 513 nodes, 2049 bytes, undeclared `y`, parameter named `exp`, duplicate parameter, parameter named `x` in `num_root`, `ifelse(x, 1, 2)` (non-relational condition), empty string, whitespace only, `sum(x)` in `num_integrate` (non-vectorised), `exp(x)` in `opt_minimize` with data (non-scalar), SANN without seed, BFGS with seed.

## 8. Not supported (by decision)
User-defined functions or recursion, vector literals and indexing, string or logical literals, matrices, gradients supplied by the caller (optim uses finite differences; nls uses numeric derivatives), `plinear`, weights, constraints other than box bounds, callbacks, more than one free variable for the 1-D ops, Unicode identifiers, ODE solvers (deSolve is not in base R). These keep the parser under ~400 lines of C++ and the R builder under 40, which is the entire surface that has to be audited.
