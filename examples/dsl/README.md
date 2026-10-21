# `_dsl_` examples (native lmr, no R)

Every recipe here compiles one or more lmscript `_dsl_` chunks against `native/schema/lmr.dsl.json` and runs them with
the single-executable `lmr` host. No R, Docker or Python is involved. The API is in [docs/dsl-api.md](../../docs/dsl-api.md).

```sh
cmake -S native -B build-native -G Ninja && cmake --build build-native   # once
lmscript examples/dsl/tut-lm.lms                                         # cwd = repository root
native/scripts/run-examples.sh build-native/bin/lmr                      # all of them, as ctest does
```

## Zero setup: `cli=` and `install=`

[`install-and-run.lms`](install-and-run.lms) needs only Node.js and lmscript (`npm i -g @lmctl-ai/lmctl`, 0.2.359 or
newer). No checkout, build or schema file:

```
p = _dsl_ cli="lmr" install="github:lmctl-ai/lmr@v0.4.0";
...
end _dsl_
savedsl(p, "job.dslb");
r = system ${dsl_cli(p)} --dsl job.dslb;
```

- `cli="lmr"`: lmscript compiles against the bytes `lmr --dsl-schema --json` prints, so the bytecode is bound to the
  binary that will run it.
- `install="github:lmctl-ai/lmr@v0.4.0"`: only when `lmr` is not on PATH, lmscript downloads that release's
  `lmr-OS-ARCH[.exe]` once into its tool cache and checks it against the SHA-256 digest GitHub reports. Nothing is
  installed without this line.
- `dsl_cli(p)`: the argv prefix of that exact binary.

The other recipes use `schema="../../native/schema/lmr.dsl.json"` and `lib/lmr.lms`, so they test the checkout's own
build.

`lib/lmr.lms` finds the binary: `$LMR_BIN`, then `lmr` on PATH, then `build-native/bin/lmr`. It writes
the bytecode to `build/dsl/NAME.dslb`. You can also run a saved chunk directly with
`build-native/bin/lmr --dsl build/dsl/tut-lm.dslb`. Each recipe exits 0 only when every self-check holds. A failed
`out.check` makes `lmr` exit 4, and the recipe passes that code on.

## Ports of the tutorial recipes

Each `tut-*.lms` and the three small recipes (`fit-and-forecast`, `fit-curve`, `minimize-rosenbrock`) is a `_dsl_`
version of the recipe with the same name in `examples/`. Each keeps the same topic, data and self-checks. The
differences:

- **One chunk instead of a round trip per call.** The analysis is a typed `_dsl_` chunk that runs in one native
  process. The original made one schema-module call (one R process) per step.
- **Functions are DSL code.** Objectives, integrands and nonlinear models are written in the chunk and driven by
  reverse communication. The original passed expression strings.
- **Different random streams.** The random generator is not R's, so seeded data differ from the originals. The checks
  are statistical (truth within k standard errors), and a recipe's header notes any seed that was changed.
- **`fit-and-forecast` checks parity with R** without depending on any random stream. It embeds the series R generated
  and compares the native ARIMA forecasts with R's.
- **No `oracle-*` ports.** Those recipes compare against `statsrust`, which is not part of this build.

## Multi-step flows

| Recipe | Pattern | Steps |
|---|---|---|
| `flow-arima-selection.lms` | values back to lmscript, then `${}` into the next chunk, CSV hand-off | explore (PACF, AIC table) → lmscript picks the order → fit, residual whiteness, 12-step forecast to CSV → lmscript report |
| `flow-ab-test.lms` | the test is chosen from diagnostics | Shapiro and F test → lmscript picks Welch t / pooled t / Wilcoxon → the chosen test plus a 2000-resample bootstrap written as a DSL loop |
| `flow-dose-response.lms` | three chunks in **one** process, sharing state through `keep`/`kept_*` | simulate and clean → compare logit/probit/cloglog GLMs by AIC → ED50/ED90 by `uniroot` on the fitted curve |
| `flow-cluster-k.lms` | a decision rule in lmscript between scans | k-means for k = 1..8 → elbow rule → final k-means vs Ward (Rand index), profiles, labels to CSV |
| `flow-anomaly-stl.lms` | an iterative loop between lmscript and a chunk | robust STL plus a MAD threshold, re-run at a lower threshold until every injected spike is found |

They show where each kind of work belongs. Numerical work, loops over observations and resampling go in the native
chunk. Decisions, logging and file orchestration go in lmscript. The values each chunk publishes with `out.set` are the
interface between the two.
