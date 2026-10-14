# rcli

`rcli` is an lmcpp native tool that runs a fixed catalogue of R computations in a fresh `Rscript --vanilla` child process per call. C++ validates the published schema and policy; the worker accepts data only and dispatches through fixed operations. The worker is embedded in the binary and materialized in a private temporary file for each call. It never evaluates caller-supplied R code.

R is licensed under GPL. This project invokes the unmodified R program as a separate process and exchanges JSON through temporary files; it does not link `libR` or `libRmath`. See [CONTRIBUTING.md](CONTRIBUTING.md) for the dependency boundary.

## Stage 0 operations

- `r_info`: R version and RNG kinds.
- `dist_quantile`: vectorized probabilities for `norm`, `t`, `chisq`, `gamma`, `beta`, `binom`, and `pois`, with `lower_tail` and `log_p` controls.
- `rng_draw`: seeded `norm` and `unif` draws. Each call sets R's RNG kind to Mersenne-Twister / Inversion / Rejection before `set.seed`.
- `fit_lm`: formulas using named data columns, plus intercept-only `y ~ 1` and intercept removal with `-1` or `+0`, with inline JSON columns or a cwd-confined CSV path. Dot formulas and expressions are rejected. CSV inputs are capped at 16 MiB and 200,000 rows.

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

Build and test in the dev image with `cmake --build build --parallel 2 && ctest --test-dir build --output-on-failure`. `tests/oracle.R` calculates expected values directly with the same environment's R. The smoke suite covers exact oracle values, schema-safe formulas, warnings, non-finite values, timeouts, environment scrubbing, result shapes, and 200,000-value input. Array flags take one item per occurrence (`--p 0.975`; repeat `--p` for more); JSON-array text as one flag value is not accepted. `scripts/check-no-libR.sh` enforces the no-link boundary and checks its detector against a planted `Rmath.h` include.

Alternate RNG kinds and normal-generation kinds are not configurable; calls use the explicit defaults listed above. The `examples/oracle-lm.lms` recipe is self-contained because no prebuilt lmda binary was available and its release build was intentionally not started. It checks the rcli fit against the known R result.

Stage 0 differences from the design notes: this slice exposes the four requested operations, with the selected seven quantile families and an `na_action` enum for `fit_lm`; it does not add a Docker runtime/shim or lmda cross-tool assertion. Its worker uses a small in-repo JSON codec, so host R does not need `jsonlite` installed. The dev image is based on Debian bookworm because that is the available `lmcpp-dev` image; its R is 4.2.2 rather than the host's 4.3.3, and every result records that version.
