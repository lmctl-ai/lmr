# rcli status and hand-off (2026-10-20, updated after the testing round)

## State
- `main` (legacy, unchanged): a C++17 lmcpp front end that spawns `Rscript` per call inside a pinned offline docker
  image.
  - It comes with the hardened host shim `bin/rcli`, 28 ops, and the lmapi schema via `--schema` / `--schema --json`.
  - It is now also the parity oracle for the native core.
- `native/` (NEW, committed 2026-10-20): the R-free implementation the owner asked for. Details are in
  [native-dsl.md](native-dsl.md).
  - **Core**: a C++17 core using only the standard library: 12 modules with unit tests whose expectations come from R
    run as a black-box oracle.
  - **Host**: the single-executable `_dsl_` host `lmr` (cpp-dsl-sdk bytecode v5), exposing 546 signatures.
  - **Build**: `cmake -S native -B build-native -G Ninja && cmake --build build-native && (cd build-native && ctest)`.
    All 16 ctest targets pass: 12 module unit tests, `test_edge`, the host suite `host_test.py`, schema parity with
    `lmscript --dsl-signatures`, and all 22 DSL examples. They are also clean under ASan/UBSan. Testing details and
    the bugs fixed are in native-dsl.md ("Tests").
- `examples/`: 14 `tut-*.lms` tutorial conversions for the R-backed tool.
  - They use `tut-vehicles.csv` and `scripts/generate-tutorial-data.py`, and each recipe asserts a known answer.
  - Three `oracle-*` recipes need `statsrust` on PATH and fail without it (expected).
  - Known limitation: the lm "power" parameter is unidentifiable.
  - Run them as `lmscript examples/x.lms` with the shim on PATH (never `lmscript run`).
- `examples/dsl/` (NEW): `_dsl_` ports of all 17 non-oracle recipes, plus 5 multi-step `flow-*` recipes. The shared
  runner is `lib/lmr.lms`; see its README.
- Docs:
  - `docs/dsl-api.md` (NEW): the native API guide and the mapping from the 28 ops.
  - `docs/parity.md` and `docs/oracle-findings.md`.
  - R is GPL. It is only ever spawned as a process; never read or translate R/CRAN source.
- The repo is PRIVATE. Per the CI rule, private repos are verified locally; keep paid macOS/Windows legs gated with
  `github.event.repository.private == false`. There is no CI job for `native/` yet.

## Not done
- **Old op names on the native core.** The lmapi op surface (`rcli fit_lm {...}`) still runs R. Serving those op names
  from `rcli_num` would need a native evaluator for `expr-spec.json` (the grammar exists in `src/expr_parser.hpp`; the
  evaluator was R-side) plus thin adapters. `native-plan.md` (the Rust `lmnum` route) is superseded for rcli by the
  C++ DSL route the owner chose on 2026-10-19.
- **Known numeric gaps** (each documented in native-dsl.md and in the module tests):
  - smooth.spline does not reproduce R's `spar` or knot subset;
  - loess uses exact local fits only, not R's kd-tree interpolation;
  - the nls "port" algorithm is not nl2sol;
  - some uniroot and L-BFGS-B iteration counts differ from R's;
  - robust STL on even-length series and near-non-stationary seasonal ARIMA fits differ (ordinary ARIMA agrees to
    about 1e-11).
- **RNG decision**: R-identical streams are NOT required. The native core uses xoshiro256** with exact variate methods.
  Do not copy GPL R/Rmath source (clean-room; scipy/statsmodels/sklearn are BSD and may be consulted, with their notices
  kept).

## Rules that apply
- Commit and push after each finished piece of work without asking (owner, 2026-10-20, for this early phase).
- Commit identity: Mike Ma <mike.ma@lmctl.com>, one Co-Authored-By per contributing model, explicit `git add`, and a
  push after every commit. Staging and deployment belong to another team.
